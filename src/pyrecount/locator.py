#! /usr/bin/env python3
from os import path
from itertools import product
from .models import Dtype, Annotation, Tags
from typing import Union, Optional, List, Dict

class ProjectLocator:
    def __init__(
            self,
            root_organism_url: str,
            data_sources: Dict[str, str],
            dbase: str,
            dtype: Dtype,
            annotation: Annotation,
            project: Union[str, List[str]],
            sample: Optional[Union[str, List[str]]] = None,
            jxn_format: Optional[str] = None
    ):

        self.root_organism_url: str = root_organism_url
        self.data_sources: Dict[str, str] = data_sources
        self.dbase: str = dbase
        self.dtype: Dtype = dtype
        self.annotation: Annotation = annotation
        self.project: List[str] = [project] if isinstance(project, str) else project 
        self.sample: List[str] = [sample] if isinstance(sample, str) else sample
        self.jxn_format: Optional[str] = jxn_format


    @property
    def _extensions(self) -> List[str]:
        match self.dtype:
            case Dtype.METADATA:
                return ['MD.gz']
            case Dtype.GENE | Dtype.EXON:
                return ['.gz']
            case Dtype.JXN:
                return ['MM.gz', 'RR.gz', 'ID.gz']
            case Dtype.BW:
                return ['ALL.bw']
            case _:
                raise ValueError(f'Invalid dtype: {self.dtype}')

    @property
    def _tags(self) -> List[str]:
        match self.dtype:
            case Dtype.METADATA:
                return [self.dbase] + Tags.METADATA.value
            case Dtype.GENE | Dtype.EXON | Dtype.JXN | Dtype.BW:
                return [self.dtype.value]
            case _:
                raise ValueError(f'Invalid dtype: {self.dtype}')

    def _get_indices(self, attribute_name):
        return {value: value[-2:] for value in getattr(self, attribute_name)}

    @property
    def _project_indices(self) -> str:
        return self._get_indices(attribute_name='project')

    @property
    def _sample_indices(self) -> str:
        return self._get_indices(attribute_name='sample')

    @property
    def fpaths(self) -> List[str]:
        paths = list()
        base = path.join(self.root_organism_url, self.data_sources[self.dbase], self.dtype.value)

        for project_id, project_index in self._project_indices.items():

            project_base = path.join(base, project_index, project_id)
            tag_extension_prod = product(self._tags, self._extensions)

            match self.dtype:
                case Dtype.METADATA:
                    file_names = [f'{self.dbase}.{tag}.{project_id}.{ext}' for tag, ext in tag_extension_prod]
                case Dtype.JXN:
                    file_names = [f'{self.dbase}.{tag}.{project_id}.{self.jxn_format.upper()}.{ext}' for tag, ext in tag_extension_prod]
                case Dtype.BW:
                    file_names = list()
                    for sample_id, sample_index in self._sample_indices.items():
                        file_names.extend(
                            [f'{sample_index}/{self.dbase}.{tag}.{project_id}_{sample_id}.{ext}' for tag, ext in tag_extension_prod]
                        )
                case _:
                    raise ValueError(f'Invalid dtype: {self.dtype}')

            paths.extend([path.join(project_base, fn) for fn in file_names])

        return paths

class MetadataLocator():
    def __init__(
            self,
            root_organism_url: str,
            data_sources: Dict[str, str],
    ):

        self.root_organism_url: str = root_organism_url
        self.data_sources: Dict[str, str] = data_sources

    @property
    def fpaths(self) -> List[str]:
        dtype, ext = 'metadata', '.recount_project.MD.gz'
        return [path.join(self.root_organism_url, dsource, dtype, path.basename(dsource) + ext) for dsource in self.data_sources.values()]
