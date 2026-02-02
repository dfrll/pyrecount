#! /usr/bin/env python3
import pytest
import polars as pl

from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation

# TODO: transform raw counts
# TODO: materialize requested columns: e.g., async polars for study id column
# TODO: multi-project support for exon, gene dtypes
# TODO: expand sra attributes
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
    "organism, dbase, project_ids, expected_shape",
    [
        ("human", "sra", ["SRP009615", "SRP075759"], (43, 174)),
        ("human", "gtex", ["FALLOPIAN_TUBE", "CERVIX_UTERI"], (28, 197)),
        ("human", "tcga", ["CHOL", "DLBC"], (93, 936)),
    ],
)
async def test_multi_project_metadata_accessor(
    organism, dbase, project_ids, expected_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project_ids))
    )

    dtype = Dtype.METADATA

    project = Project(
        metadata=project_meta_dataframe,
        dbase=dbase,
        organism=organism,
        dtype=[dtype],
        annotation=None,
        jxn_format=None,
        root_url=root_url,
    )

    await project.cache()
    project_dataframe = project.load(dtype)

    assert project_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, dbase, project_ids, expected_mm_shape, expected_shape",
    [
        ("human", "sra", ["SRP009615", "SRP075759"], (436480, 43), (717928, 11)),
        (
            "human",
            "gtex",
            ["FALLOPIAN_TUBE", "CERVIX_UTERI"],
            (807104, 28),
            (1424738, 11),
        ),
        ("human", "tcga", ["CHOL", "DLBC"], (824438, 93), (1645794, 11)),
        ("mouse", "sra", ["SRP111354", "SRP200978"], (325976, 27), (634751, 11)),
    ],
)
async def test_multi_project_jxn_accessor(
    organism, dbase, project_ids, expected_mm_shape, expected_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project_ids))
    )

    dtype = Dtype.JXN

    jxn = Project(
        metadata=project_meta_dataframe,
        dbase=dbase,
        organism=organism,
        dtype=[dtype],
        annotation=None,
        jxn_format="all",
        root_url=root_url,
    )

    await jxn.cache()
    jxn_mm_dataframe, jxn_dataframe = jxn.load(dtype)

    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, dbase, project_ids, annotation, expected_annotation_shape, expected_counts_shape",
    [
        (
            "human",
            "sra",
            ["SRP009615"],
            Annotation.GENCODE_V29,
            (1377600, 21),
            (1377601, 16),
        ),
        (
            "human",
            "gtex",
            ["FALLOPIAN_TUBE"],
            Annotation.GENCODE_V29,
            (1377600, 21),
            (1377601, 13),
        ),
        (
            "human",
            "tcga",
            ["CHOL"],
            Annotation.GENCODE_V29,
            (1377600, 21),
            (1377601, 49),
        ),
        (
            "mouse",
            "sra",
            ["SRP111354"],
            Annotation.GENCODE_V23,
            (841915, 21),
            (841916, 19),
        ),
    ],
)
async def test_project_exon_accessor(
    organism,
    dbase,
    annotation,
    project_ids,
    expected_annotation_shape,
    expected_counts_shape,
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project_ids))
    )

    dtype = Dtype.EXON

    exon = Project(
        metadata=project_meta_dataframe,
        dbase=dbase,
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
    "organism, dbase, project_ids, annotation, expected_annotation_shape, expected_counts_shape",
    [
        (
            "human",
            "sra",
            ["SRP009615"],
            Annotation.GENCODE_V29,
            (64836, 21),
            (64837, 13),
        ),
        (
            "human",
            "gtex",
            ["FALLOPIAN_TUBE"],
            Annotation.GENCODE_V29,
            (64836, 21),
            (64837, 10),
        ),
        ("human", "tcga", ["CHOL"], Annotation.GENCODE_V29, (64836, 21), (64837, 46)),
        (
            "mouse",
            "sra",
            ["SRP111354"],
            Annotation.GENCODE_V23,
            (55420, 21),
            (55421, 16),
        ),
    ],
)
async def test_project_gene_accessor(
    organism,
    dbase,
    project_ids,
    annotation,
    expected_annotation_shape,
    expected_counts_shape,
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project_ids))
    )

    dtype = Dtype.GENE

    gene = Project(
        metadata=project_meta_dataframe,
        dbase=dbase,
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


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "organism, dbase, project_ids, annotation, expected_shape",
    [
        ("human", "sra", ["SRP009615"], Annotation.GENCODE_V29, (12, 1)),
        (
            "human",
            "gtex",
            ["FALLOPIAN_TUBE"],
            Annotation.GENCODE_V29,
            (9, 1),
        ),
        ("human", "tcga", ["CHOL"], Annotation.GENCODE_V29, (45, 1)),
        ("mouse", "sra", ["SRP111354"], Annotation.GENCODE_V23, (15, 1)),
    ],
)
def test_project_bigwig_accessor(
    organism, dbase, project_ids, annotation, expected_shape
):
    root_url = "http://duffel.rail.bio/recount3"
    recount_metadata = Metadata(organism=organism, root_url=root_url)

    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col("project").is_in(project_ids))
    )

    dtype = Dtype.BW
    bw = Project(
        metadata=project_meta_dataframe,
        dbase=dbase,
        organism=organism,
        dtype=[dtype],
        annotation=None,
        jxn_format=None,
        root_url=root_url,
    )

    bw = bw.load(dtype)

    assert bw.shape == expected_shape
