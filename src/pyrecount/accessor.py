#! /usr/bin/env python3
import re
import asyncio
import logging
import polars as pl
from os import path, makedirs
from dataclasses import dataclass, field

from scipy.io import mmread
from urllib.parse import urlparse
from urllib.request import urlretrieve
from typing import Optional, Union, List, Tuple

from .api import EndpointConnector
from .models import Dtype, Annotation, Extensions
from .locator import ProjectLocator, MetadataLocator

from .utils import (
    replace_organism,
    download_url_to_path,
)

from pprint import pprint

log = logging.getLogger()


@dataclass
class Project:
    metadata: pl.DataFrame
    dbase: str
    organism: str
    dtype: List[Dtype]
    annotation: Optional[Annotation] = None
    jxn_format: str = None
    root_url: str = "http://duffel.rail.bio/recount3/"

    project_ids: List[str] = field(init=False)
    sample: List[str] = field(init=False)
    endpoints: "EndpointConnector" = field(init=False)

    def __post_init__(self):
        if not isinstance(self.metadata, pl.DataFrame):
            raise TypeError(
                "Parameter value for `metadata` must be a Polars DataFrame."
            )
        if not isinstance(self.dtype, list):
            raise TypeError("Parameter `dtype` must be a list.")

        if Dtype.JXN in self.dtype and self.jxn_format is None:
            raise ValueError(
                f"Parameter `jxn_format` is required when `dtype` is {Dtype.JXN}."
            )
        if (
            Dtype.GENE in self.dtype or Dtype.EXON in self.dtype
        ) and self.annotation is None:
            raise ValueError(
                f"Parameter `annotation` is required when `dtype` is {self.dtype}."
            )

        self.project_ids: List[str] = self.metadata["project"].unique().to_list()
        self.sample: List[str] = self.metadata["external_id"].unique().to_list()
        self.endpoints = EndpointConnector(
            root_url=self.root_url, organism=self.organism
        )

    def _get_project_urls(self, dtype) -> List[str]:
        project = ProjectLocator(
            root_organism_url=self.endpoints.root_organism_url,
            data_sources=self.endpoints.data_sources,
            dbase=self.dbase,
            dtype=dtype,
            annotation=self.annotation,
            project_ids=self.project_ids,
            sample=self.sample,
            jxn_format=self.jxn_format,
        )

        return project.urls

    async def cache(self) -> None:
        tasks = list()
        url_list = list()
        for dtype in self.dtype:
            url_list.extend(self._get_project_urls(dtype))
            for url in url_list:
                fpath = urlparse(url).path.lstrip("/")
                if path.exists(fpath):
                    continue
                makedirs(path.dirname(fpath), exist_ok=True)
                tasks.append(download_url_to_path(url=url, fpath=fpath))

        await asyncio.gather(*tasks)

    def load(self, dtype) -> Union[pl.DataFrame, Tuple[pl.DataFrame, pl.DataFrame]]:
        match dtype:
            case Dtype.METADATA:
                return self._metadata_load()
            case Dtype.JXN:
                return self._jxn_load()
            case Dtype.GENE:
                return self._gene_load()
            case Dtype.EXON:
                return self._exon_load()
            case Dtype.BW:
                return self._bw_load()
            case _:
                raise ValueError(f"Invalid dtype: {self.dtype}")

    def _metadata_load(self) -> pl.DataFrame:
        cache_dataframe = None
        join_cols = ["rail_id", "external_id", "study"]

        for project_id in self.project_ids:
            project_dataframe = None
            for url in self._get_project_urls(Dtype.METADATA):
                if project_id not in url:
                    continue

                if self.dbase not in url:
                    continue
                if Dtype.METADATA.value not in url:
                    continue

                # skip metadata pred if loading gtex or tcga
                if self.dbase in ["gtex", "tcga"] and "pred" in url:
                    continue

                fpath = urlparse(url).path.lstrip("/")

                current_dataframe = pl.read_csv(
                    fpath, separator="\t", infer_schema=False
                )

                try:
                    if project_dataframe is None:
                        project_dataframe = current_dataframe
                    else:
                        project_dataframe = project_dataframe.join(
                            current_dataframe, on=join_cols
                        )
                except Exception as e:
                    logging.error(f"Error joining resource {fpath} to metadata: {e}")
                    return

            if cache_dataframe is None:
                cache_dataframe = project_dataframe
            else:
                cache_dataframe, project_dataframe = self._add_missing_columns(
                    cache_dataframe, project_dataframe
                )
                try:
                    cache_dataframe = pl.concat(
                        [cache_dataframe, project_dataframe], how="vertical"
                    )
                except Exception as e:
                    logging.error(f"Error concatenating {fpath} to metadata: {e}")
                    return

        return replace_organism(cache_dataframe).unique()

    def _add_missing_columns(
        self, df1: pl.DataFrame, df2: pl.DataFrame
    ) -> Tuple[pl.DataFrame, pl.DataFrame]:
        missing_in_df1 = [col for col in df2.columns if col not in df1.columns]
        missing_in_df2 = [col for col in df1.columns if col not in df2.columns]

        df1_types = {col: df1.schema[col] for col in missing_in_df2}
        df2_types = {col: df2.schema[col] for col in missing_in_df1}

        for col in missing_in_df1:
            df1 = df1.with_columns(pl.lit(None).cast(df2_types[col]).alias(col))
        for col in missing_in_df2:
            df2 = df2.with_columns(pl.lit(None).cast(df1_types[col]).alias(col))
        return df1, df2

    def _jxn_load(self) -> Tuple[pl.DataFrame, pl.DataFrame]:
        cache_mm_dataframe: Optional[pl.DataFrame] = None
        cache_dataframe: Optional[pl.DataFrame] = None

        for project_id in self.project_ids:
            url_list = self._get_project_urls(Dtype.JXN)

            # rename junction dataframe columns using project ID + rail IDs
            for url in url_list:
                if project_id not in url:
                    continue
                if self.dbase not in url:
                    continue
                if Dtype.JXN.value not in url:
                    continue
                fpath = urlparse(url).path.lstrip("/")
                if "ID" in fpath:
                    ids = [
                        f"{project_id}_{rail_id}"
                        for rail_id in pl.read_csv(fpath)["rail_id"]
                        .cast(pl.String)
                        .to_list()
                    ]

            for url in url_list:
                if project_id not in url:
                    continue
                if self.dbase not in url:
                    continue
                if Dtype.JXN.value not in url:
                    continue
                if "ID" in url:
                    continue

                fpath = urlparse(url).path.lstrip("/")

                if "MM" in url:
                    # the samples in the MM jxn table are not in the same order as the metadata.
                    # this is the reason for renaming using rail IDs.
                    project_mm_dataframe = pl.from_numpy(mmread(fpath).toarray())
                    project_mm_dataframe = project_mm_dataframe.rename(
                        dict(zip(project_mm_dataframe.columns, ids))
                    )

                    if cache_mm_dataframe is None:
                        cache_mm_dataframe = project_mm_dataframe
                    else:
                        cache_mm_dataframe = pl.concat(
                            [cache_mm_dataframe, project_mm_dataframe], how="horizontal"
                        )

                else:
                    current_dataframe = pl.read_csv(
                        fpath, separator="\t", infer_schema=False
                    )
                    current_dataframe = current_dataframe.rename(
                        dict(
                            zip(
                                current_dataframe.columns,
                                [
                                    f"{project_id}_{colname}"
                                    for colname in current_dataframe.columns
                                ],
                            )
                        )
                    )
                    if cache_dataframe is None:
                        cache_dataframe = current_dataframe

                    else:
                        cache_dataframe = pl.concat(
                            [cache_dataframe, current_dataframe], how="horizontal"
                        )

        return cache_mm_dataframe, cache_dataframe

    def _bw_load(self) -> pl.DataFrame:
        for resource in self.cache_resources:
            if self.dbase not in resource.rname:
                continue
        return pl.DataFrame()

    def _read_gtf(self, rpath: str) -> pl.DataFrame:
        annotation_dataframe = pl.read_csv(
            rpath,
            comment_prefix="#",
            separator="\t",
            new_columns=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
        )

        fields = [
            "gene_id",
            "transcript_id",
            "exon_number",
            "gene_name",
            "gene_source",
            "gene_biotype",
            "transcript_name",
            "transcript_source",
            "transcript_biotype",
            "protein_id",
            "exon_id",
            "tag",
        ]

        return annotation_dataframe.with_columns(
            [
                annotation_dataframe["attribute"]
                .map_elements(
                    lambda x: re.findall(rf'{field} "([^"]*)"', x)[0]
                    if rf'{field} "' in x
                    else "",
                    return_dtype=pl.Utf8,
                )
                .alias(field)
                for field in fields
            ]
        )

    # TODO: extract first column (chromosome|start_1base|end_1ba…)
    def _read_counts(self, rname: str):
        counts_dataframe = pl.read_csv(
            rname,
            comment_prefix="#",
            separator="\t",
        )
        return counts_dataframe

    # def _gene_exon_load(self, dtype: Dtype):
    # for url in self._get_project_urls(dtype):
    # fpath = urlparse(url).path.lstrip("/")
    # if self.annotation.value in fpath:
    # if any(fpath.endswith(ext) for ext in Extensions.GENE.value):
    # annotation = self._read_gtf(fpath)
    # if fpath.endswith(f"{self.annotation.value}.gz"):
    # counts = self._read_counts(fpath)
    # return annotation, counts
    # return

    def _gene_load(self) -> pl.DataFrame:
        for url in self._get_project_urls(Dtype.GENE):
            fpath = urlparse(url).path.lstrip("/")
            if self.annotation.value in fpath:
                if any(fpath.endswith(ext) for ext in Extensions.GENE.value):
                    annotation = self._read_gtf(fpath)
                if fpath.endswith(f"{self.annotation.value}.gz"):
                    counts = self._read_counts(fpath)
        return annotation, counts

    def _exon_load(self) -> pl.DataFrame:
        for url in self._get_project_urls(Dtype.EXON):
            fpath = urlparse(url).path.lstrip("/")
            if self.annotation.value in fpath:
                if any(url.endswith(ext) for ext in Extensions.EXON.value):
                    annotation = self._read_gtf(fpath)
                if url.endswith(f"{self.annotation.value}.gz"):
                    counts = self._read_counts(fpath)
        return annotation, counts


class Metadata:
    def __init__(
        self,
        organism: str,
        root_url: str = "http://duffel.rail.bio/recount3/",
    ):
        self.organism = organism
        self.eps = EndpointConnector(organism, root_url)

    def cache(self) -> None:
        self.urls = MetadataLocator(
            self.eps.root_organism_url, self.eps.data_sources
        ).urls

        for url in self.urls:
            fpath = urlparse(url).path.lstrip("/")
            if path.exists(fpath):
                continue
            makedirs(path.dirname(fpath), exist_ok=True)
            urlretrieve(url, fpath)

    def load(self) -> pl.DataFrame:
        sep = "\t"
        pattern = ".recount_project."

        dataframes = list()
        for url in self.urls:
            if self.organism not in url:
                continue
            if pattern not in url:
                continue
            try:
                fpath = urlparse(url).path.lstrip("/")
                df = pl.read_csv(path.join(fpath), separator=sep, infer_schema=None)
            except Exception as e:
                logging.error(f"Error reading file {fpath}: {e}")
                return

            dataframes.append(df)

        if not dataframes:
            logging.error("All file reads failed.")
            return pl.DataFrame()

        meta_dataframe = pl.concat(dataframes, how="vertical")

        return replace_organism(meta_dataframe).unique()
