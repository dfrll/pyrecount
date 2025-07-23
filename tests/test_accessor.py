#! /usr/bin/env python3
import pytest
import polars as pl

from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation

# TODO: test dbs other than sra for jxn, exon, gene dtypes
# TODO: handle BigWig
# TODO: transform raw counts
# TODO: multi-project support for exon, gene dtypes
# TODO: define dataframe schemas
# TODO: expose Lazyframes


@pytest.mark.parametrize(
    "organism, expected_shape", [("human", (347005, 8)), ("mouse", (416859, 8))]
)
def test_recount_metadata_accessor(organism, expected_shape):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()
    recount_meta_dataframe = recount_metadata.load()

    assert recount_meta_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, dbase, project, expected_shape",
    [
        ("human", "sra", ["SRP009615", "SRP075759"], (43, 174)),
        ["human", "gtex", ["FALLOPIAN_TUBE", "CERVIX_UTERI"], (28, 197)],
        ("human", "tcga", ["CHOL", "DLBC"], (93, 936)),
    ],
)
async def test_multi_project_metadata_accessor(
    organism, dbase, project, expected_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project))
    )

    dtype = Dtype.METADATA

    project = Project(
        metadata=project_meta_dataframe,
        dbase=dbase,
        organism=organism,
        dtype=[dtype],
        annotation=Annotation,
        jxn_format=None,
        root_url=root_url,
    )

    await project.cache()
    project_dataframe = project.load(dtype)

    assert project_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, project, expected_mm_shape, expected_shape",
    [
        ("human", ["SRP009615", "SRP075759"], (436480, 43), (717928, 11)),
        ("mouse", ["SRP111354", "SRP200978"], (325976, 27), (634751, 11)),
    ],
)
async def test_multi_project_jxn_accessor(
    organism, project, expected_mm_shape, expected_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project))
    )

    dtype = Dtype.JXN

    jxn = Project(
        metadata=project_meta_dataframe,
        dbase="sra",
        organism=organism,
        dtype=[dtype],
        annotation=Annotation,
        jxn_format="all",
        root_url=root_url,
    )

    await jxn.cache()
    jxn_mm_dataframe, jxn_dataframe = jxn.load(dtype)

    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, project, expected_shape",
    [("human", ["SRP009615"], (12, 174)), ("mouse", ["SRP111354"], (15, 176))],
)
async def test_project_metadata_accessor(organism, project, expected_shape):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project))
    )

    dtype = Dtype.METADATA

    project = Project(
        metadata=project_meta_dataframe,
        dbase="sra",
        organism=organism,
        dtype=[dtype],
        annotation=Annotation,
        jxn_format=None,
        root_url=root_url,
    )

    await project.cache()
    meta_dataframe = project.load(dtype)

    assert meta_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, project, annotation, expected_annotation_shape, expected_counts_shape",
    [
        ("human", ["SRP009615"], Annotation.GENCODE_V29, (1377600, 21), (1377601, 16)),
        ("mouse", ["SRP111354"], Annotation.GENCODE_V23, (841915, 21), (841916, 19)),
    ],
)
async def test_project_exon_accessor(
    organism, annotation, project, expected_annotation_shape, expected_counts_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project))
    )

    dtype = Dtype.EXON

    exon = Project(
        metadata=project_meta_dataframe,
        dbase="sra",
        organism=organism,
        dtype=[dtype],
        annotation=annotation,
        jxn_format=None,
        root_url=root_url,
    )

    await exon.cache()
    exon_annotation, exon_counts = exon.load(dtype)

    assert exon_annotation.shape == expected_annotation_shape
    assert exon_counts.shape == expected_counts_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, project, annotation, expected_annotation_shape, expected_counts_shape",
    [
        ("human", ["SRP009615"], Annotation.GENCODE_V29, (64836, 21), (64837, 13)),
        ("mouse", ["SRP111354"], Annotation.GENCODE_V23, (55420, 21), (55421, 16)),
    ],
)
async def test_project_gene_accessor(
    organism, project, annotation, expected_annotation_shape, expected_counts_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project))
    )

    dtype = Dtype.GENE

    gene = Project(
        metadata=project_meta_dataframe,
        dbase="sra",
        organism=organism,
        dtype=[dtype],
        annotation=annotation,
        jxn_format=None,
        root_url=root_url,
    )

    await gene.cache()
    gene_annotation, gene_counts = gene.load(dtype)

    assert gene_annotation.shape == expected_annotation_shape
    assert gene_counts.shape == expected_counts_shape
