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
from typing import Protocol, Any
from collections.abc import Iterable


from .api import EndpointConnector
from .models import Dtype, Annotation, Extensions
from .locator import ProjectLocator, MetadataLocator

from .utils import (
    replace_organism,
    download_url_to_path,
)


log = logging.getLogger()


class Loader(Protocol):
    dtype: Dtype

    async def cache(self) -> None: ...
    def load(self) -> Any: ...


@dataclass
class Project:
    metadata: pl.DataFrame
    dbase: str
    organism: str
    annotation: Annotation | None = None
    jxn_format: str | None = None
    root_url: str = "http://duffel.rail.bio/recount3/"

    project_ids: list[str] = field(init=False)
    sample: list[str] = field(init=False)
    endpoints: EndpointConnector = field(init=False)

    _cached_metadata: pl.DataFrame | None = field(init=False, default=None)

    def __post_init__(self):
        if not isinstance(self.metadata, pl.DataFrame):
            raise TypeError("`metadata` must be a Polars DataFrame.")

        self.project_ids = self.metadata["project"].unique().to_list()
        self.sample = self.metadata["external_id"].unique().to_list()
        self.endpoints = EndpointConnector(
            root_url=self.root_url,
            organism=self.organism,
        )

    def _get_loader(self, dtype: Dtype) -> Loader:
        loaders: dict[Dtype, type[Loader]] = {
            Dtype.METADATA: MetadataLoader,
            Dtype.GENE: GeneLoader,
            Dtype.EXON: ExonLoader,
            Dtype.JXN: JunctionLoader,
            Dtype.BW: BigWigLoader,
        }
        try:
            return loaders[dtype](self)
        except KeyError:
            raise ValueError(f"No loader registered for {dtype}")

    async def cache(self, dtypes: Dtype | Iterable[Dtype]) -> None:
        if isinstance(dtypes, Dtype):
            dtypes = (dtypes,)

        tasks: list[asyncio.Task] = []

        for dtype in dtypes:
            loader = self._get_loader(dtype)
            tasks.append(asyncio.create_task(loader.cache()))

        if tasks:
            await asyncio.gather(*tasks)

    def load(self, dtype):
        loader = self._get_loader(dtype)
        return loader.load()

    def get_project_urls(self, dtype: Dtype) -> list[str]:
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

    async def _cache_urls(self, urls: list[str]) -> None:
        tasks: list[asyncio.Task] = []

        for url in urls:
            fpath = urlparse(url).path.lstrip("/")
            if path.exists(fpath):
                continue

            makedirs(path.dirname(fpath), exist_ok=True)
            tasks.append(download_url_to_path(url=url, fpath=fpath))

        # launches all tasks, without rate limiting
        if tasks:
            await asyncio.gather(*tasks)

    def scale_mapped_reads(
        self, counts: pl.DataFrame, target_size: float, L: int
    ) -> pl.DataFrame:
        md = self.load_metadata()

        mapped_reads = pl.col("star.all_mapped_reads").cast(pl.Float64)
        avg_mapped_len = pl.col("star.average_mapped_length").cast(pl.Float64)
        avg_read_len = pl.col("avg_len").cast(pl.Float64)

        # paired-end detection
        ratio = (avg_mapped_len / avg_read_len).round(0)
        paired_end = ratio == 2
        paired_factor = pl.when(paired_end).then(2).otherwise(1)

        # scale factors
        sf = md.select(
            [
                pl.col("external_id"),
                (
                    target_size * L * paired_factor / (mapped_reads * avg_mapped_len**2)
                ).alias("sf"),
            ]
        )

        sf_map = dict(zip(sf["external_id"], sf["sf"]))

        return counts.with_columns(
            [
                (pl.col(c) * sf_map[c])
                for c in counts.select(pl.selectors.numeric()).columns
            ]
        )

    def scale_auc(self, counts: pl.DataFrame, target_size: float) -> pl.DataFrame:
        md = self.load_metadata()

        auc = pl.col("bc_auc.all_reads_all_bases").cast(pl.Float64)

        # scale factors
        sf = md.select(
            "external_id",
            (target_size / auc).alias("sf"),
        )

        sf_map = dict(zip(sf["external_id"], sf["sf"]))

        return counts.with_columns(
            [
                pl.col(c).mul(sf_map[c]).round(0).cast(pl.Int64).alias(c)
                for c in counts.columns
                if c != "gene_id"
            ]
        )

    def load_metadata(self) -> pl.DataFrame:
        if self._cached_metadata is None:
            self._cached_metadata = self._get_loader(Dtype.METADATA).load()
        return self._cached_metadata

    def _add_missing_columns(
        self, df1: pl.DataFrame, df2: pl.DataFrame
    ) -> tuple[pl.DataFrame, pl.DataFrame]:
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

    def _gtf_read(self, rpath: str) -> pl.DataFrame:
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
                    lambda x: (
                        re.findall(rf'{field} "([^"]*)"', x)[0]
                        if rf'{field} "' in x
                        else ""
                    ),
                    return_dtype=pl.Utf8,
                )
                .alias(field)
                for field in fields
            ]
        )

    def _counts_read(self, rname: str):
        df = pl.read_csv(
            rname,
            comment_prefix="#",
            separator="\t",
        )

        first_col = df.columns[0]

        if not self.sample:
            return df

        keep = [first_col] + self.sample
        missing = set(keep) - set(df.columns)
        if missing:
            raise KeyError(f"Missing columns in counts file: {missing}")
        return df.select(keep)


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


class GeneLoader:
    dtype = Dtype.GENE

    def __init__(self, project: Project):
        self.project = project

        if project.annotation is None:
            raise ValueError("GeneLoader requires project.annotation to be set")

    def _urls(self) -> list[str]:
        return self.project.get_project_urls(self.dtype)

    async def cache(self) -> None:
        await self.project._cache_urls(self._urls())

    def load(self) -> tuple[pl.DataFrame, pl.DataFrame]:
        annotation: pl.DataFrame | None = None
        counts: pl.DataFrame | None = None

        for url in self._urls():
            fpath = urlparse(url).path.lstrip("/")

            if self.project.annotation.value in fpath:
                if any(fpath.endswith(ext) for ext in Extensions.GENE.value):
                    annotation = self.project._gtf_read(fpath)
                elif fpath.endswith(f"{self.project.annotation.value}.gz"):
                    counts = self.project._counts_read(fpath)

        if annotation is None or counts is None:
            raise RuntimeError("Missing gene annotation or counts file")

        return annotation, counts


class JunctionLoader:
    dtype = Dtype.JXN

    def __init__(self, project: Project):
        self.project = project

        if project.jxn_format is None:
            raise ValueError("JunctionLoader requires project.jxn_format")

    def _urls(self) -> list[str]:
        return self.project.get_project_urls(self.dtype)

    async def cache(self) -> None:
        await self.project._cache_urls(self._urls())

    def load(self) -> tuple[pl.DataFrame, pl.DataFrame]:
        cache_mm: list[pl.DataFrame] = []
        cache_meta: list[pl.DataFrame] = []

        for project_id in self.project.project_ids:
            project_urls = [u for u in self._urls() if project_id in u]

            ids: list[str] | None = None

            for url in project_urls:
                fpath = urlparse(url).path.lstrip("/")

                if "ID" in url:
                    ids = pl.read_csv(fpath)["rail_id"].cast(pl.String).to_list()

            for url in project_urls:
                fpath = urlparse(url).path.lstrip("/")

                if "ID" in url:
                    continue

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


class MetadataLoader:
    dtype = Dtype.METADATA

    def __init__(self, project: Project):
        self.project = project

    def _urls(self) -> list[str]:
        return self.project.get_project_urls(self.dtype)

    async def cache(self) -> None:
        await self.project._cache_urls(self._urls())

    def load(self) -> pl.DataFrame:
        cache_meta: list[pl.DataFrame] = []
        join_cols = ["rail_id", "external_id", "study"]

        for project_id in self.project.project_ids:
            dfs_for_project: list[pl.DataFrame] = []

            for url in self._urls():
                if project_id not in url:
                    continue

                fpath = urlparse(url).path.lstrip("/")
                df = pl.read_csv(fpath, separator="\t", infer_schema=False)

                if self.project.sample:
                    df = df.filter(
                        pl.col("external_id").is_in(self.project.sample),
                    )
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
            cache_dataframe, df = self.project._add_missing_columns(cache_dataframe, df)
            cache_dataframe = pl.concat([cache_dataframe, df], how="vertical")

        return replace_organism(cache_dataframe).unique()


class ExonLoader:
    dtype = Dtype.EXON

    def __init__(self, project: Project):
        self.project = project

        if project.annotation is None:
            raise ValueError("ExonLoader requires project.annotation to be set")

    def _urls(self) -> list[str]:
        return self.project.get_project_urls(self.dtype)

    async def cache(self) -> None:
        await self.project._cache_urls(self._urls())

    def load(self) -> pl.DataFrame:
        annotation: pl.DataFrame | None = None
        counts: pl.DataFrame | None = None

        for url in self._urls():
            fpath = urlparse(url).path.lstrip("/")

            if self.project.annotation.value not in fpath:
                continue

            # annotation GTF
            if any(url.endswith(ext) for ext in Extensions.EXON.value):
                annotation = self.project._gtf_read(fpath)

            # exon counts
            elif fpath.endswith(f"{self.project.annotation.value}.gz"):
                counts = self.project._counts_read(fpath)

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

        if annotation is None or counts is None:
            raise RuntimeError("Missing exon annotation or counts file")

        return annotation, counts


class BigWigLoader:
    dtype = Dtype.BW

    def __init__(self, project: Project):
        self.project = project

    def _urls(self) -> list[str]:
        return self.project.get_project_urls(self.dtype)

    async def cache(self) -> None:
        await self.project._cache_urls(self._urls())

    def load(self) -> pl.DataFrame:
        rows: list[dict[str, str]] = []

        for url in self._urls():
            fpath = urlparse(url).path.lstrip("/")

            project_id = next(
                (pid for pid in self.project.project_ids if pid in url),
                None,
            )

            if project_id is None:
                continue

            rows.append(
                {
                    "project_id": project_id,
                    "url": url,
                    "path": fpath,
                }
            )

        if not rows:
            raise RuntimeError("No BigWig files found")

        return pl.DataFrame(rows)
