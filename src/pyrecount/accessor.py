#! /usr/bin/env python3
import logging
from os import path
import polars as pl
from glob import glob
from scipy.io import mmread
from typing import Optional, Union, List, Tuple
from pybiocfilecache import BiocFileCache, models
from pybiocfilecache.cache import BiocFileCache as BiocFileCacheType

from .utils import *
from .cache import QCache
from .api import EndpointConnector
from .models import Dtype, Annotation
from .locator import ProjectLocator, MetadataLocator

log = logging.getLogger()

class Project():
    def __init__(
            self,
            metadata: pl.DataFrame,
            dbase: str, 
            # TODO: Union[Dtype, List[Dtype]]: allow multiple dtypes
            dtype: Dtype,
            cache_location: Optional[str],
            # TODO: handle annotation argument
            annotation: Optional[Annotation] = None,
            jxn_format: str = None,
            root_url: str = 'http://duffel.rail.bio/recount3/'
    ):

        if dtype == Dtype.JXN and jxn_format is None:
            raise ValueError(f'Parameter `jxn_format` is required when `dtype` is {Dtype.JXN}.')

        self.metadata: pl.DataFrame = metadata
        self.dbase: str = dbase
        self.dtype: Dtype = dtype
        self.cache_location: str = path.join(cache_location, self.dtype.value)
        self.annotation: Annotation = annotation
        self.jxn_format: Optional[str] = jxn_format
        self.root_url: str = root_url

    def load(self) -> Union[pl.DataFrame, Tuple[pl.DataFrame, pl.DataFrame]]:
        cache: BiocFileCacheType = BiocFileCache(self.cache_location)
        cache_resources: List[models.Resource] = cache.list_resources()

        match self.dtype:
            case Dtype.METADATA:
                return self._metadata(cache_resources)
            case Dtype.JXN:
                return self._jxn(cache, cache_resources)
            case _:
                raise ValueError(f'Invalid dtype: {self.dtype}')

    def cache(self):
        organism = self.metadata['organism'].unique()
        endpoints = EndpointConnector(root_url=self.root_url, organism=organism.first())

        self.project = self.metadata['project'].unique().to_list()
        sample = self.metadata['external_id'].unique().to_list()

        project_metadata = ProjectLocator(
            root_organism_url = endpoints.root_organism_url,
            data_sources = endpoints.data_sources,
            dbase = self.dbase,
            dtype = self.dtype,
            # TODO: handle Annotation type
            annotation = self.annotation,
            project = self.project,
            sample = sample,
            jxn_format = self.jxn_format
        )

        QCache(
            fpaths=project_metadata.fpaths,
            cache_location=self.cache_location
        ).biocache()

    def _get_jxn_ids(self, cache: BiocFileCacheType) -> List:
        cache_items = {item.rname: item.rpath for item in cache.list_resources()}
        ids_path = [rpath for rname, rpath in cache_items.items() if 'ID.gz' in rname][0]
        return pl.read_csv(ids_path)['rail_id'].cast(str).to_list()

    def _id_matrix_market(self, mm_dataframe: pl.DataFrame, cache: BiocFileCacheType) -> pl.DataFrame:
        ''' set column names to sample identifiers '''
        ids: List[str] = self._get_jxn_ids(cache)
        return mm_dataframe.rename(dict(zip(mm_dataframe.columns, ids)))

    def _metadata(self, cache_resources: List[models.Resource]):
        join_cols = ['rail_id', 'external_id', 'study']

        cache_dataframe: pl.DataFrame = self.metadata

        for resource in cache_resources:
            if self.dbase not in resource.rname:
                continue

            current_dataframe = pl.read_csv(resource.rpath, separator='\t')

            try:
                cache_dataframe = cache_dataframe.join(current_dataframe, on=join_cols, suffix=resource.rid)
            except Exception as e:
                logging.error(f'Error joining resource: {e}')

        cache_dataframe = cache_dataframe.drop([col for col in cache_dataframe.columns if 'BFC' in col])

        return replace_organism(cache_dataframe).unique()

    def _jxn(self, cache: BiocFileCacheType, cache_resources: List[models.Resource]) -> Tuple[pl.DataFrame, pl.DataFrame]:
        join_cols = ['rail_id']

        # XXX: decide mapping from sample id to jxn db entry
        # only set here to avoid schema conflict in downstream join
        #cache_dataframe: pl.DataFrame = None

        mm_dataframe = None
        for resource in cache_resources:
            if self.dbase not in resource.rname:
                continue
            if 'ID' in resource.rname:
                continue

            if 'MM' in resource.rname:
                mm_array = mmread(resource.rpath).toarray()
                mm_dataframe = self._id_matrix_market(pl.from_numpy(mm_array), cache)
            else:
                # The samples in the MM jxn table are not in the same order as the metadata
                # XXX: decide mapping from sample id to jxn db entry
                #ids = self._get_jxn_ids(cache)
                # skip copy of metadata table
                current_dataframe = pl.read_csv(resource.rpath, separator='\t')
                if all(col in current_dataframe.columns for col in join_cols):
                    continue
                else:
                    # XXX: decide mapping from sample id to jxn db entry
                    #cache_dataframe = pl.concat([pl.DataFrame(ids), current_dataframe], how='horizontal')
                    cache_dataframe = current_dataframe

        return mm_dataframe, cache_dataframe

class Metadata():
    def __init__(
            self,
            organism: str,
            cache_location: Optional[str],
            root_url: str = 'http://duffel.rail.bio/recount3/'
    ):
        self.organism: str = organism
        self.cache_location: str = cache_location
        self.root_url: str = root_url

    def cache(self) -> None:
        try:
            endpoints = EndpointConnector(root_url=self.root_url, organism=self.organism)
            metadata = MetadataLocator(root_organism_url=endpoints.root_organism_url, data_sources=endpoints.data_sources)

            QCache(
                fpaths=metadata.fpaths,
                cache_location=self.cache_location
            # not spawning threads for so few resources
            ).biocache_serial()

        except Exception as e:
            log.error(e)

    def load(self) -> pl.DataFrame:
        # TODO: extract extension + pattern from .sqlite db
        ext = 'MD.gz'
        pattern = 'recount_project'
        sep = '\t'

        fpaths = glob(path.join(self.cache_location, f'*{pattern}*{ext}'))

        meta_dataframe = pl.read_csv(fpaths[0], separator=sep)

        for table_fpath in fpaths[1:]:
            current_table = pl.read_csv(table_fpath, separator='\t')
            meta_dataframe = pl.concat([meta_dataframe, current_table])

        return replace_organism(meta_dataframe).unique()
