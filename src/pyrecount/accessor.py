#! /usr/bin/env python3
import re
import asyncio
import logging
import polars as pl
from functools import reduce
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

# caching policy
CACHEABLE_DTYPES = {
    Dtype.METADATA,
    Dtype.JXN,
    Dtype.GENE,
    Dtype.EXON,
    # Dtype.BW
}


@dataclass
class Project:
    metadata: pl.DataFrame
    dbase: str
    organism: str
    dtype: List[Dtype]
    annotation: Optional[Annotation] = None
    jxn_format: Optional[str] = None
    root_url: str = "http://duffel.rail.bio/recount3/"

    project_ids: List[str] = field(init=False)
    sample: List[str] = field(init=False)
    endpoints: "EndpointConnector" = field(init=False)

    def __post_init__(self):
        if not isinstance(self.metadata, pl.DataFrame):
            raise TypeError(
                "Parameter value for `metadata` must be a Polars DataFrame."
            )

        if not isinstance(self.dtype, list) or not all(
            isinstance(d, Dtype) for d in self.dtype
        ):
            raise TypeError("Parameter value for `dtype` must be a list of Dtype.")

        if Dtype.JXN in self.dtype and self.jxn_format is None:
            raise ValueError(
                "Parameter value for `jxn_format` is required when `Dtype.JXN` is included in `dtype`."
            )

        if (
            Dtype.GENE in self.dtype or Dtype.EXON in self.dtype
        ) and self.annotation is None:
            raise ValueError(
                f"Parameter value for `annotation` is required when `dtype` is {self.dtype}."
            )

        self.project_ids = self.metadata["project"].unique().to_list()
        self.sample = self.metadata["external_id"].unique().to_list()
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
        tasks: List[asyncio.Task] = list()

        for dtype in self.dtype:
            if dtype not in CACHEABLE_DTYPES:
                continue

            for url in self._get_project_urls(dtype):
                fpath = urlparse(url).path.lstrip("/")
                if path.exists(fpath):
                    continue

                makedirs(path.dirname(fpath), exist_ok=True)
                tasks.append(download_url_to_path(url=url, fpath=fpath))

        if tasks:
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

    def _valid_metadata_url(self, url: str, project_id: str) -> bool:
        if project_id not in url:
            return False
        if self.dbase not in url:
            return False
        if Dtype.METADATA.value not in url:
            return False
        if self.dbase in {"gtex", "tcga"} and "pred" in url:
            return False
        return True

    def _metadata_load(self) -> pl.DataFrame:
        cache_meta: list[pl.DataFrame] = []
        join_cols = ["rail_id", "external_id", "study"]
        urls = self._get_project_urls(Dtype.METADATA)

        for project_id in self.project_ids:
            dfs_for_project: list[pl.DataFrame] = []

            for url in urls:
                if not self._valid_metadata_url(url, project_id):
                    continue

                fpath = urlparse(url).path.lstrip("/")
                df = pl.read_csv(fpath, separator="\t", infer_schema=False)
                dfs_for_project.append(df)

            if not dfs_for_project:
                continue

            # join all metadata files for the project
            project_dataframe = (
                dfs_for_project[0]
                if len(dfs_for_project) == 1
                else reduce(
                    lambda left, right: left.join(right, on=join_cols, how="inner"),
                    dfs_for_project,
                )
            )

            cache_meta.append(project_dataframe)

        if not cache_meta:
            raise RuntimeError("No metadata loaded.")

        # concatenate all project dataframes once
        cache_dataframe = cache_meta[0]
        for df in cache_meta[1:]:
            cache_dataframe, df = self._add_missing_columns(cache_dataframe, df)
            cache_dataframe = pl.concat([cache_dataframe, df], how="vertical")

        return replace_organism(cache_dataframe).unique()

    def _add_missing_columns(
        self, df1: pl.DataFrame, df2: pl.DataFrame
    ) -> Tuple[pl.DataFrame, pl.DataFrame]:
        missing_in_df1 = [col for col in df2.columns if col not in df1.columns]
        missing_in_df2 = [col for col in df1.columns if col not in df2.columns]

        if missing_in_df1:
            df1 = df1.with_columns(
                [
                    pl.lit(None, dtype=df2.schema[col]).alias(col)
                    for col in missing_in_df1
                ]
            )

        if missing_in_df2:
            df2 = df2.with_columns(
                [
                    pl.lit(None, dtype=df1.schema[col]).alias(col)
                    for col in missing_in_df2
                ]
            )

        common_order = sorted(set(df1.columns) | set(df2.columns))
        df1 = df1.select(common_order)
        df2 = df2.select(common_order)

        return df1, df2

    def _jxn_load(self) -> Tuple[pl.DataFrame, pl.DataFrame]:
        cache_mm: List[pl.DataFrame] = []
        cache_meta: List[pl.DataFrame] = []

        urls = self._get_project_urls(Dtype.JXN)

        for project_id in self.project_ids:
            project_urls = [
                u
                for u in urls
                if project_id in u and self.dbase in u and Dtype.JXN.value in u
            ]

            ids: list[str] | None = None

            for url in project_urls:
                if "ID" in url:
                    fpath = urlparse(url).path.lstrip("/")
                    ids = pl.read_csv(fpath)["rail_id"].cast(pl.String).to_list()

            for url in project_urls:
                if "ID" in url:
                    continue

                fpath = urlparse(url).path.lstrip("/")

                if "MM" in url:
                    if ids is None:
                        raise RuntimeError(f"No ID file found for {project_id}")

                    mm = mmread(fpath).toarray()
                    mm_df = pl.from_numpy(mm)

                    if len(ids) != mm_df.width:
                        raise ValueError("Mismatch between MM columns and IDs")

                    mm_df = mm_df.rename(dict(zip(mm_df.columns, ids)))
                    cache_mm.append(mm_df)

                else:
                    df = pl.read_csv(
                        fpath, separator="\t", infer_schema=False
                    ).with_columns(pl.lit(project_id).alias("project_id"))

                    cache_meta.append(df)

        if not cache_mm or not cache_meta:
            raise RuntimeError("No junction data loaded.")

        return (
            pl.concat(cache_mm, how="horizontal"),
            pl.concat(cache_meta, how="vertical"),
        )

    def _bw_load(self) -> pl.DataFrame:
        urls = self._get_project_urls(Dtype.BW)

        for project_id in self.project_ids:
            project_urls = [
                u
                for u in urls
                if project_id in u and self.dbase in u and Dtype.BW.value in u
            ]

        return pl.DataFrame(project_urls, schema=["url"])

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

    def _read_counts(self, rname: str):
        counts_dataframe = pl.read_csv(
            rname,
            comment_prefix="#",
            separator="\t",
        )
        return counts_dataframe

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
                    # TODO: extract first column (chromosome|start_1base|end_1ba…)
                    exon_colname = counts.columns[0]
                    exon_fields = ["chrom", "start", "end", "strand"]
                    counts = (
                        counts.with_columns(
                            pl.col(exon_colname)
                            .str.split_exact("|", 4)
                            .struct.rename_fields(exon_fields)
                            .alias("exon_parts")
                        )
                        .unnest("exon_parts")
                        .drop(exon_colname)
                    )
                    reorder_cols = exon_fields + [
                        col for col in counts.columns if col not in set(exon_fields)
                    ]
                    counts = counts.select(reorder_cols)

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
