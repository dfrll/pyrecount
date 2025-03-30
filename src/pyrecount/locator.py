#! /usr/bin/env python3
from os import path
from itertools import product
from .models import Dtype, Annotation, Tags, Extensions
from typing import Union, Optional, List, Dict

from pprint import pprint


class ProjectLocator:
    def __init__(
        self,
        root_organism_url: str,
        data_sources: Dict[str, str],
        dbase: str,
        dtype: List[Dtype],
        annotation: Annotation,
        project_ids: List[str],
        sample: Optional[Union[str, List[str]]] = None,
        jxn_format: Optional[str] = None,
    ):
        if dtype == Dtype.BW and sample is None:
            raise ValueError(
                f"Parameter `sample` is required when `dtype` is {Dtype.BW}"
            )
        # if not isinstance(dtype, List):
        # raise TypeError('Parameter value for `dtype` must be a list.')
        if not isinstance(project_ids, List):
            raise TypeError("Parameter value for `project_ids` must be a list.")

        self.root_organism_url: str = root_organism_url
        self.data_sources: Dict[str, str] = data_sources
        self.dbase: str = dbase
        self.dtype: List[Dtype] = [dtype] if isinstance(dtype, str) else dtype
        self.annotation: Annotation = annotation
        self.project_ids: List[str] = project_ids
        self.sample: List[str] = [sample] if isinstance(sample, str) else sample
        self.jxn_format: Optional[str] = jxn_format

    # TODO: replace with .models.extensions
    @property
    def _extensions(self) -> List[str]:
        match self.dtype:
            case Dtype.METADATA:
                return ["MD.gz"]
            case Dtype.GENE | Dtype.EXON:
                return ["gtf.gz"]
            case Dtype.JXN:
                return ["MM.gz", "RR.gz", "ID.gz"]
            case Dtype.BW:
                return ["ALL.bw"]
            case _:
                raise ValueError(f"Invalid dtype: {self.dtype}")

    @property
    def _tags(self) -> List[str]:
        match self.dtype:
            case Dtype.METADATA:
                return [self.dbase] + Tags.METADATA.value
            case Dtype.GENE | Dtype.EXON | Dtype.JXN | Dtype.BW:
                return [self.dtype.value]
            case _:
                raise ValueError(f"Invalid dtype: {self.dtype}")

    def _get_indices(self, attribute_name) -> Dict[str, str]:
        return {value: value[-2:] for value in getattr(self, attribute_name)}

    @property
    def _project_indices(self) -> Dict[str, str]:
        return self._get_indices(attribute_name="project_ids")

    @property
    def _sample_indices(self) -> Dict[str, str]:
        return self._get_indices(attribute_name="sample")

    @property
    def urls(self) -> List[str]:
        paths = list()
        base = path.join(self.root_organism_url)

        tag_extension_prod = list(product(self._tags, self._extensions))

        match self.dtype:
            case Dtype.METADATA:
                project_base = path.join(
                    base, self.data_sources[self.dbase], self.dtype.value
                )
                for project_id, project_index in self._project_indices.items():
                    file_names = [
                        f"{self.dbase}.{tag}.{project_id}.{ext}"
                        for tag, ext in tag_extension_prod
                    ]
                    paths.extend(
                        path.join(project_base, project_index, project_id, fn)
                        for fn in file_names
                    )

            case Dtype.JXN:
                project_base = path.join(
                    base, self.data_sources[self.dbase], self.dtype.value
                )
                for project_id, project_index in self._project_indices.items():
                    file_names = [
                        f"{self.dbase}.{tag}.{project_id}.{self.jxn_format.upper()}.{ext}"
                        for tag, ext in tag_extension_prod
                    ]
                    paths.extend(
                        path.join(project_base, project_index, project_id, fn)
                        for fn in file_names
                    )

            case Dtype.GENE | Dtype.EXON:
                # annotations (.gtf)
                organism = path.basename(base)
                annotation_base = path.join(base, "annotations", self.dtype.value)
                annotation_files = [
                    f"{organism}.{self.dtype.value}.{self.annotation.value}.{ext}"
                    for ext in self._extensions
                ]
                paths.extend(path.join(annotation_base, fn) for fn in annotation_files)

                # raw coverage counts
                project_base = path.join(
                    base, self.data_sources[self.dbase], self.dtype.value
                )
                for project_id, project_index in self._project_indices.items():
                    current_base = path.join(project_base, project_index, project_id)
                    project_files = [
                        f"{self.dbase}.{tag}.{project_id}.{self.annotation.value}.gz"
                        for tag in self._tags
                    ]
                    paths.extend(path.join(current_base, fn) for fn in project_files)

            # TODO:
            # case Dtype.BW:
            # project_base = path.join(base, self.data_sources[self.dbase], self.dtype.value)

            # for project_id, project_index in self._project_indices.items():
            # for sample_id, sample_index in self._sample_indices.items():
            # file_names = [
            # path.join(sample_index, f'{self.dbase}.{tag}.{project_id}_{sample_id}.{ext}')
            # for tag, ext in tag_extension_prod
            # ]
            # paths.extend(path.join(project_base, project_index, project_id, fn) for fn in file_names)

            # for project_id, project_index in self._project_indices.items():
            # file_names = list()
            # project_base = path.join(base, self.data_sources[self.dbase], self.dtype.value, project_index, project_id)
            # for sample_id, sample_index in self._sample_indices.items():
            # file_names.extend(
            # [f'{sample_index}/{self.dbase}.{tag}.{project_id}_{sample_id}.{ext}' for tag, ext in tag_extension_prod]
            # )
            # paths.extend([path.join(project_base, fn) for fn in file_names])
            case _:
                raise ValueError(f"Invalid dtype: {self.dtype}")

        return paths


class MetadataLocator:
    def __init__(
        self,
        root_organism_url: str,
        data_sources: Dict[str, str],
    ):
        self.root_organism_url: str = root_organism_url
        self.data_sources: Dict[str, str] = data_sources

    @property
    def urls(self) -> List[str]:
        dtype, ext = Dtype.METADATA, ".recount_project.MD.gz"
        return [
            path.join(
                self.root_organism_url,
                dsource,
                dtype.value,
                path.basename(dsource) + ext,
            )
            for dsource in self.data_sources.values()
        ]
