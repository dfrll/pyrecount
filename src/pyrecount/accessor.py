#! /usr/bin/env python3
import re
import logging
from os import path
import polars as pl
from scipy.io import mmread
from typing import Optional, Union, List, Tuple
from pybiocfilecache import BiocFileCache, models
from pybiocfilecache.cache import BiocFileCache as BiocFileCacheType

from .utils import *
from .cache import QCache
from .api import EndpointConnector
from .models import Dtype, Annotation, Extensions
from .locator import ProjectLocator, MetadataLocator

from pprint import pprint

log = logging.getLogger()

class Project():
    def __init__(
            self,
            metadata: pl.DataFrame,
            dbase: str,
            dtype: Dtype,
            cache_location: str,
            annotation: Optional[Annotation] = None,
            jxn_format: str = None,
            root_url: str = 'http://duffel.rail.bio/recount3/'
    ):

        if dtype == Dtype.JXN and jxn_format is None:
            raise ValueError(f'Parameter `jxn_format` is required when `dtype` is {Dtype.JXN}.')
        if dtype == Dtype.GENE and annotation is None:
            raise ValueError(f'Parameter `annotations` is required when `dtype` is {Dtype.GENE}.')

        self.metadata: pl.DataFrame = metadata
        self.dbase: str = dbase
        self.dtype: Dtype = dtype
        self.cache_location: str = cache_location
        self.annotation: Annotation = annotation
        self.jxn_format: Optional[str] = jxn_format
        self.root_url: str = root_url

        self.organism = self.metadata['organism'].unique().first()
        self.organism_cache_location: str = path.join(self.cache_location, self.organism) if self.cache_location else None
        self.endpoints = EndpointConnector(root_url=self.root_url, organism=self.organism)
        self.project = self.metadata['project'].unique().to_list()
        self.sample = self.metadata['external_id'].unique().to_list()


    def cache(self):

        project = ProjectLocator(
            root_organism_url = self.endpoints.root_organism_url,
            data_sources = self.endpoints.data_sources,
            dbase = self.dbase,
            dtype = self.dtype,
            annotation = self.annotation,
            project = self.project,
            sample = self.sample,
            jxn_format = self.jxn_format
        )

        qcache = QCache(
            fpaths = project.fpaths,
            organism_cache_location = self.organism_cache_location
        )

        match self.dtype:
            case Dtype.METADATA | Dtype.JXN | Dtype.BW:
                qcache.biocache()
            case Dtype.GENE | Dtype.EXON:
                # not spawning threads for so few data sources.
                qcache.biocache_serial()
            case _:
                raise ValueError(f'Invalid dtype: {self.dtype}')


    def load(self) -> Union[pl.DataFrame, Tuple[pl.DataFrame, pl.DataFrame]]:

        cache: BiocFileCacheType = BiocFileCache(self.organism_cache_location)
        cache_resources: List[models.Resource] = cache.list_resources()

        match self.dtype:
            case Dtype.METADATA:
                return self._metadata_load(cache_resources)
            case Dtype.JXN:
                # TODO: sending `cache` not necessary
                return self._jxn_load(cache, cache_resources)
            case Dtype.GENE:
                return self._gene_load(cache_resources)
            case Dtype.EXON:
                return self._exon_load(cache_resources)
            #case Dtype.BW:
                ## TODO: sending cache not necessary
                ## XXX: expose BigWig URLs rather than caching
                #return self._bw_load(cache, cache_resources)
            case _:
                raise ValueError(f'Invalid dtype: {self.dtype}')


    def _get_jxn_ids(self, cache: BiocFileCacheType) -> List:
        # XXX: redundant procedure for reading cache.list_resources()
        cache_items = {item.rname: item.rpath for item in cache.list_resources()}
        ids_path = [rpath for rname, rpath in cache_items.items() if 'ID.gz' in rname][0]
        return pl.read_csv(ids_path, infer_schema=False)['rail_id'].to_list()#.cast(str).to_list()


    def _id_matrix_market(self, mm_dataframe: pl.DataFrame, cache: BiocFileCacheType) -> pl.DataFrame:
        ''' set column names to sample identifiers '''
        ids: List[str] = self._get_jxn_ids(cache)
        return mm_dataframe.rename(dict(zip(mm_dataframe.columns, ids)))


    def _metadata_load(self, cache_resources: List[models.Resource]):

        join_cols = ['rail_id', 'external_id', 'study']
        cache_dataframe = None

        for resource in cache_resources:
            if self.dbase not in resource.rname:
                continue

            # read files associated with self.dtype only.
            ext_vals = getattr(Extensions, self.dtype.name).value
            if not any(ext in resource.rname for ext in ext_vals):
                continue

            try:
                current_dataframe = pl.read_csv(resource.rpath, separator='\t', infer_schema=False)
            except Exception as e:
                logging.error(f'Error reading file {resource.rpath}: {e}')
                return

            try:
                if cache_dataframe is None:
                    cache_dataframe = current_dataframe
                else:
                    # use suffix to avoid name clashes in column names
                    cache_dataframe = cache_dataframe.join(current_dataframe, on=join_cols, suffix=resource.rid)
            except Exception as e:
                logging.error(f'Error joining resource {resource.rid} to metadata: {e}')
                return

        # drop duplicate columns
        cache_dataframe = cache_dataframe.drop([col for col in cache_dataframe.columns if 'BFC' in col])

        return replace_organism(cache_dataframe).unique()


    def _jxn_load(self, cache: BiocFileCacheType, cache_resources: List[models.Resource]) -> Tuple[pl.DataFrame, pl.DataFrame]:

        for resource in cache_resources:
            if self.dbase not in resource.rname:
                continue
            if self.dtype.value not in resource.rname:
                continue

            # ids skipped until call to self._id_matrix_market()
            if 'ID' in resource.rname:
                continue

            if 'MM' in resource.rname:
                # the samples in the MM jxn table are not in the same order as the metadata.
                # this is the reason for calling self._get_jxn_ids().
                mm_array = mmread(resource.rpath).toarray()
                mm_dataframe = self._id_matrix_market(pl.from_numpy(mm_array), cache)
            else:
                try:
                    current_dataframe = pl.read_csv(resource.rpath, separator='\t', infer_schema=False)
                except Exception as e:
                    logging.error(f'Error reading file {resource.rpath}: {e}')
                    return 

                # ignores {ID}_{dbase}.recount_project.MD.gz
                if all(col in current_dataframe.columns for col in self.metadata.columns):
                    continue
                else:
                    cache_dataframe = current_dataframe

        # possible to concatenate instead.
        # return pl.concat([cache_dataframe, mm_dataframe])
        return mm_dataframe, cache_dataframe


    # TODO:
    def _bw_load(self, cache: BiocFileCacheType, cache_resources: List[models.Resource]) -> Tuple[pl.DataFrame, pl.DataFrame]:
        # XXX: expose BigWig URLs rather than caching
        for resource in cache_resources:
            if self.dbase not in resource.rname:
                continue
            #bw: bigWigFile = pyBigWig.open(resource.rpath)
        return pl.DataFrame()


    def _read_gtf(self, rpath: str) -> pl.DataFrame:

        annotation_dataframe = pl.read_csv(
            rpath,
            comment_prefix = '#',
            separator = '\t',
            new_columns = ['seqname', 'source', 'feature', 'start', 'end',
                            'score', 'strand', 'frame', 'attribute']
        )

        fields = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_source',
                    'gene_biotype', 'transcript_name', 'transcript_source',
                    'transcript_biotype', 'protein_id', 'exon_id', 'tag'
        ]

        return annotation_dataframe.with_columns([
                annotation_dataframe['attribute'].map_elements(
                    lambda x: re.findall(rf'{field} "([^"]*)"', x)[0] if rf'{field} "' in x else '',
                    return_dtype=pl.Utf8  
                ).alias(field) for field in fields]
        )


    def _read_counts(self, rname: str):
        # TODO: extract first column (chromosome|start_1base|end_1ba…)
        counts_dataframe = pl.read_csv(
            rname,
            comment_prefix = '#',
            separator = '\t',
        )
        return counts_dataframe 


    def _gene_load(self, cache_resources: List[models.Resource]) -> pl.DataFrame:
        for resource in cache_resources:
            if self.annotation.value in resource.rname:
                if any(resource.rname.endswith(ext) for ext in Extensions.EXON.value):
                    annotation = self._read_gtf(resource.rpath)
                if resource.rname.endswith(f'{self.annotation.value}.gz'):
                    counts = self._read_counts(resource.rpath)

        return annotation, counts


    def _exon_load(self, cache_resources: List[models.Resource]) -> pl.DataFrame:
        for resource in cache_resources:
            if self.annotation.value in resource.rname:
                if any(resource.rname.endswith(ext) for ext in Extensions.EXON.value):
                    annotation = self._read_gtf(resource.rpath)
                if resource.rname.endswith(f'{self.annotation.value}.gz'):
                    counts = self._read_counts(resource.rpath)

        return annotation, counts


class Metadata():
    def __init__(
            self,
            organism: str,
            cache_location: str,
            root_url: str = 'http://duffel.rail.bio/recount3/'
    ):
        self.organism: str = organism
        self.organism_cache_location: str = path.join(cache_location, organism) if cache_location else None
        self.root_url: str = root_url


    def cache(self) -> None:
        try:
            endpoints = EndpointConnector(root_url=self.root_url, organism=self.organism)
            metadata = MetadataLocator(root_organism_url=endpoints.root_organism_url, data_sources=endpoints.data_sources)

            QCache(
                fpaths = metadata.fpaths,
                organism_cache_location = self.organism_cache_location
            # not spawning threads for so few data sources.
            ).biocache_serial()

        except Exception as e:
            log.error(e)


    def load(self) -> pl.DataFrame:
        pattern = '.recount_project.'
        sep = '\t'

        cache: BiocFileCacheType = BiocFileCache(self.organism_cache_location)
        cache_resources: List[models.Resource] = cache.list_resources()

        dataframes = list()
        for resource in cache_resources:
            if pattern not in resource.rname:
                continue

            try:
                logging.info(f'Reading file: {resource.rpath}')
                current_dataframe = pl.read_csv(resource.rpath, separator=sep)
                dataframes.append(current_dataframe)
            except Exception as e:
                logging.error(f'Error reading file {resource.rpath}: {e}')
                return

        if not dataframes:
            logging.error('All file reads failed.')
            return pl.DataFrame()

        meta_dataframe = pl.concat(dataframes, how='vertical')
        return replace_organism(meta_dataframe).unique()
