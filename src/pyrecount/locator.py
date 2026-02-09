#! /usr/bin/env python3
from os import path
from itertools import product
from .models import Dtype, Annotation, Tags


class ProjectLocator:
    def __init__(
        self,
        root_organism_url: str,
        data_sources: dict[str, str],
        dbase: str,
        dtype: Dtype,
        annotation: Annotation | None,
        project_ids: list[str],
        sample: str | list[str] | None = None,
        jxn_format: str | None = None,
    ):
        if dtype in (Dtype.GENE, Dtype.EXON) and annotation is None:
            raise ValueError(f"`annotation` is required when `dtype` is {dtype}")
        if dtype == Dtype.BW and sample is None:
            raise ValueError(f"`sample` is required when `dtype` is {Dtype.BW}")
        if not isinstance(project_ids, list):
            raise TypeError("`project_ids` must be a list.")

        self.root_organism_url: str = root_organism_url
        self.data_sources: dict[str, str] = data_sources
        self.dbase: str = dbase
        self.dtype: Dtype = dtype
        self.annotation: Annotation | None = annotation
        self.project_ids: list[str] = project_ids
        self.sample: list[str] = (
            [sample]
            if isinstance(sample, str)
            else sample
            if sample is not None
            else []
        )
        self.jxn_format: str | None = jxn_format

    # TODO: replace with .models.extensions
    @property
    def _extensions(self) -> list[str]:
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
    def _tags(self) -> list[str]:
        match self.dtype:
            case Dtype.METADATA:
                return [self.dbase] + Tags.METADATA.value
            case Dtype.GENE | Dtype.EXON | Dtype.JXN | Dtype.BW:
                return [self.dtype.value]
            case _:
                raise ValueError(f"Invalid dtype: {self.dtype}")

    def _get_indices(self, attribute_name) -> dict[str, str]:
        values = getattr(self, attribute_name)
        if not values:
            return {}
        return {v: v[-2:] for v in values}

    @property
    def _project_indices(self) -> dict[str, str]:
        # return last last two characters of project id
        return self._get_indices(attribute_name="project_ids")

    @property
    def _sample_indices(self) -> dict[str, str]:
        # return last last two characters of sample id
        return self._get_indices(attribute_name="sample")

    @property
    def urls(self) -> list[str]:
        paths = list()
        base = path.join(self.root_organism_url)
        tag_extension_prod = list(product(self._tags, self._extensions))

        # TODO: condense to case Dtype.METADATA | Dtype.JXN
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

            case Dtype.BW:
                project_base = path.join(
                    base, self.data_sources[self.dbase], self.dtype.value
                )
                for project_id, project_index in self._project_indices.items():
                    for sample_id, sample_index in self._sample_indices.items():
                        for tag, ext in tag_extension_prod:
                            filename = (
                                f"{self.dbase}.{tag}.{project_id}_{sample_id}.{ext}"
                            )

                            paths.append(
                                path.join(
                                    project_base,
                                    project_index,
                                    project_id,
                                    sample_index,
                                    filename,
                                )
                            )

            case _:
                raise ValueError(f"Invalid dtype: {self.dtype}")

        return paths


class MetadataLocator:
    def __init__(
        self,
        root_organism_url: str,
        data_sources: dict[str, str],
    ):
        self.root_organism_url: str = root_organism_url
        self.data_sources: dict[str, str] = data_sources

    @property
    def urls(self) -> list[str]:
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
