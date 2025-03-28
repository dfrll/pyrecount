#! /usr/bin/env python3
import re
import logging
from os import path
import polars as pl
from scipy.io import mmread
from typing import Optional, Union, List, Tuple, Dict
from pybiocfilecache import BiocFileCache, models
from pybiocfilecache.cache import BiocFileCache as BiocFileCacheType

from .utils import *
from .cache import Cache
from .api import EndpointConnector
from .models import Dtype, Tags, Annotation, Extensions
from .locator import ProjectLocator, MetadataLocator

from pprint import pprint

log = logging.getLogger()


class Project():
    def __init__(
            self,
            metadata: pl.DataFrame,
            dbase: str,
            organism: str,
            dtype: List[Dtype],
            cache_location: str,
            annotation: Optional[Annotation] = None,
            jxn_format: str = None,
            root_url: str = 'http://duffel.rail.bio/recount3/'
    ):

        if not isinstance(metadata, pl.DataFrame):
            raise TypeError('Parameter value for `metadata` must be a Polars DataFrame.')
        if not isinstance(dtype, List):
            raise TypeError('Parameter value for `dtype` must be a list.')
        if dtype == Dtype.JXN and jxn_format is None:
            raise ValueError(f'Parameter `jxn_format` is required when `dtype` is {Dtype.JXN}.')
        if dtype in (Dtype.GENE, Dtype.EXON) and annotation is None:
            raise ValueError(f'Parameter `annotations` is required when `dtype` is {dtype}.')

        self.metadata: pl.DataFrame = metadata
        self.dbase: str = dbase
        self.organism: str = organism
        self.dtype: List[Dtype] = dtype
        self.cache_location: str = cache_location
        self.annotation: Annotation = annotation
        self.jxn_format: Optional[str] = jxn_format
        self.root_url: str = root_url

        self.project: List[str] = self.metadata['project'].unique().to_list()
        self.sample: List[str] = self.metadata['external_id'].unique().to_list()
        self.endpoints = EndpointConnector(root_url=self.root_url, organism=self.organism)


    def _get_project_fpaths(self, dtype) -> List[str]:
        project = ProjectLocator(
            root_organism_url = self.endpoints.root_organism_url,
            data_sources = self.endpoints.data_sources,
            dbase = self.dbase,
            dtype = dtype,
            annotation = self.annotation,
            project = self.project,
            sample = self.sample,
            jxn_format = self.jxn_format)
        return project.fpaths


    def cache(self):

        for dtype in self.dtype:
            fpaths = self._get_project_fpaths(dtype)
            cache = Cache(fpaths = fpaths, cache_location = self.cache_location)
            match dtype:
                case Dtype.METADATA | Dtype.JXN:
                    cache.biocache_threadpool()
                case Dtype.GENE | Dtype.EXON | Dtype.BW:
                    cache.biocache_serial()
                case _:
                    raise ValueError(f'Invalid dtype: {self.dtype}')


    def load(self) -> Union[pl.DataFrame, Tuple[pl.DataFrame, pl.DataFrame]]:

        cache: BiocFileCacheType = BiocFileCache(self.cache_location)
        self.cache_resources: List[models.Resource] = cache.list_resources()

        for dtype in self.dtype:
            match dtype:
                case Dtype.METADATA:
                    return self._metadata_load()
                case Dtype.JXN:
                    return self._jxn_load()
                case Dtype.GENE:
                    return self._gene_load()
                case Dtype.EXON:
                    return self._exon_load()
                #case Dtype.BW:
                    ## XXX: expose BigWig URLs rather than caching
                    #return self._bw_load()
                case _:
                    raise ValueError(f'Invalid dtype: {self.dtype}')


    # XXX: implement multi-project support before this method
    #def get_read_counts(self, raw_counts: pl.DataFrame, meta_dataframe: pl.DataFrame) -> pl.DataFrame:
        ## use pivot rather than transpose for external_id tracking
        ## XXX: use variable_name and value_name in pivot function calls
        #counts = raw_counts \
            #.unpivot(index='gene_id') \
            #.pivot(on='gene_id', index='variable', values='value') \
            #.sort(by='variable')

        #meta_dataframe = meta_dataframe \
            #.with_columns(
                #pl.col('star.average_mapped_length').cast(pl.Float64)
            #).sort(by='external_id')

        #counts = counts.with_columns(
            #(counts[col] / meta_dataframe['star.average_mapped_length']).alias(col)
            #for col in counts.columns if col != 'variable'
        #).unpivot(index='variable', variable_name='gene_id').\
        #pivot(on='variable', index='gene_id', values='value')

        #return counts


    def _metadata_load(self):

        cache_dataframe = None
        join_cols = ['rail_id', 'external_id', 'study']

        for project_id in self.project:
            project_dataframe = None
            for resource in self.cache_resources:
                if project_id not in resource.rname:
                    continue
                if self.dbase not in resource.rname:
                    continue
                if Dtype.METADATA.value not in resource.rname:
                    continue

                try:
                    current_dataframe = pl.read_csv(resource.rpath, separator='\t', infer_schema=False)
                except Exception as e:
                    logging.error(f'Error reading file {resource.rpath}: {e}')
                    return

                try:
                    if project_dataframe is None:
                        project_dataframe = current_dataframe
                    else:
                        # use suffix to avoid name clashes in column names
                        project_dataframe = project_dataframe.join(current_dataframe, on=join_cols, suffix=resource.rid)
                except Exception as e:
                    logging.error(f'Error joining resource {resource.rid} to metadata: {e}')
                    return

            # drop duplicate columns
            project_dataframe = project_dataframe.drop([col for col in project_dataframe.columns if 'BFC' in col])

            if cache_dataframe is None:
                cache_dataframe = project_dataframe
            else:
                cache_dataframe, project_dataframe = self._add_missing_columns(cache_dataframe, project_dataframe)
                cache_dataframe = pl.concat([cache_dataframe, project_dataframe], how='vertical')

        return replace_organism(cache_dataframe).unique()


    # TODO: refactor: abstract method
    def _add_missing_columns(self, df1: pl.DataFrame, df2: pl.DataFrame) -> Tuple[pl.DataFrame, pl.DataFrame]:

        missing_in_df1 = [col for col in df2.columns if col not in df1.columns]
        missing_in_df2 = [col for col in df1.columns if col not in df2.columns]

        df1_types = {col: df1.schema[col] for col in missing_in_df2}
        df2_types = {col: df2.schema[col] for col in missing_in_df1}

        for col in missing_in_df1:
            df1 = df1.with_columns(pl.lit(None).cast(df2_types[col]).alias(col))
        for col in missing_in_df2:
            df2 = df2.with_columns(pl.lit(None).cast(df1_types[col]).alias(col))
        return df1, df2


    def _id_matrix_market(self, mm_dataframe: pl.DataFrame, cache_resources: BiocFileCacheType) -> pl.DataFrame:
        ''' set column names to sample identifiers '''
        ids: List[str] = list()
        for resource in cache_resources:
            if 'ID' not in resource.rname:
                continue
            else:
                # inaccurate if multiple ID files match in cache_resources per project.
                # this shouldn't happen based on file structure provided by the api.
                ids = pl.read_csv(resource.rpath)['rail_id'].cast(pl.String).to_list()
        return mm_dataframe.rename(dict(zip(mm_dataframe.columns, ids)))


    def _jxn_load(self) -> Tuple[pl.DataFrame, pl.DataFrame]:

        for project_id in self.project:
            for resource in self.cache_resources:
                if project_id not in resource.rname:
                    continue
                if self.dbase not in resource.rname:
                    continue
                if Dtype.JXN.value not in resource.rname:
                    continue

                # ids skipped until call to self._id_matrix_market()
                if 'ID' in resource.rname:
                    continue

                if 'MM' in resource.rname:
                    # the samples in the MM jxn table are not in the same order as the metadata.
                    # this is the reason for calling self._id_matrix_market().
                    mm_array = mmread(resource.rpath).toarray()
                    mm_dataframe = self._id_matrix_market(pl.from_numpy(mm_array), self.cache_resources)
                else:
                    try:
                        current_dataframe = pl.read_csv(resource.rpath, separator='\t', infer_schema=False)
                    except Exception as e:
                        logging.error(f'Error reading file {resource.rpath}: {e}')
                        return

                    # ignores {ID}_{dbase}.recount_project.MD.gz, which seems to be a copy of self.metadata.
                    if all(col in current_dataframe.columns for col in self.metadata.columns):
                        continue
                    else:
                        cache_dataframe = current_dataframe

        # possible to concatenate instead.
        # return pl.concat([cache_dataframe, mm_dataframe])
        return mm_dataframe, cache_dataframe


    def _bw_load(self) -> pl.DataFrame:
        # XXX: expose BigWig URLs rather than caching
        for resource in self.cache_resources:
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

        fields = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 
                  'gene_source', 'gene_biotype', 'transcript_name', 
                  'transcript_source', 'transcript_biotype', 'protein_id', 
                  'exon_id', 'tag']

        return annotation_dataframe.with_columns([
                annotation_dataframe['attribute'].map_elements(
                    lambda x: re.findall(rf'{field} "([^"]*)"', x)[0] if rf'{field} "' in x else '',
                    return_dtype=pl.Utf8  
                ).alias(field) for field in fields]
        )


    # TODO: extract first column (chromosome|start_1base|end_1ba…)
    def _read_counts(self, rname: str):
        counts_dataframe = pl.read_csv(
            rname,
            comment_prefix = '#',
            separator = '\t',
        )
        return counts_dataframe


    # TODO: abstract method GENE, EXON load
    def _gene_load(self) -> pl.DataFrame:
        for resource in self.cache_resources:
            if self.annotation.value in resource.rname:
                if any(resource.rname.endswith(ext) for ext in Extensions.GENE.value):
                    annotation = self._read_gtf(resource.rpath)
                if resource.rname.endswith(f'{self.annotation.value}.gz'):
                    counts = self._read_counts(resource.rpath)
        return annotation, counts


    # TODO: abstract method GENE, EXON load
    def _exon_load(self) -> pl.DataFrame:
        for resource in self.cache_resources:
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
        self.cache_location: str = cache_location
        self.root_url: str = root_url

    def cache(self) -> None:
        try:
            endpoints = EndpointConnector(root_url=self.root_url, organism=self.organism)
            metadata = MetadataLocator(root_organism_url=endpoints.root_organism_url, data_sources=endpoints.data_sources)

            Cache(
                fpaths = metadata.fpaths,
                cache_location = self.cache_location
            # not spawning threads for so few data sources.
            ).biocache_serial()

        except Exception as e:
            log.error(e)


    def load(self) -> pl.DataFrame:
        pattern = '.recount_project.'
        sep = '\t'

        cache: BiocFileCacheType = BiocFileCache(self.cache_location)
        cache_resources: List[models.Resource] = cache.list_resources()

        dataframes = list()
        for resource in cache_resources:
            if pattern not in resource.rname:
                continue
            if self.organism not in resource.rname:
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
