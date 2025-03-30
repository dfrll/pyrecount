#! /usr/bin/env python3
import pytest
import polars as pl

from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation


@pytest.mark.asyncio
@pytest.mark.xfail(reason="feature not implemented: multi-project jxn load support.")
@pytest.mark.parametrize(
    "organism, project, expected_mm_shape, expected_shape",
    [
        ("human", ["SRP009615", "SRP075759"], (562896, 43), (562896, 10)),
    ],
)
async def test_multi_project_jxn_accessor(
    organism, project, expected_shape, expected_mm_shape
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
    )

    await jxn.cache()

    jxn_mm_dataframe, jxn_dataframe = jxn.load(dtype)

    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape


@pytest.mark.asyncio
@pytest.mark.xfail(
    reason="Project method _read_counts() needs to split first column to multiple columns."
    "This method passes after .load() method because polars ignores the column type in reading."
    "polars.exceptions.InvalidOperationError: cannot cast List type (inner: 'Int64', to: 'String')"
)
@pytest.mark.parametrize(
    "organism, project, annotation, expected_counts_shape",
    [
        ("human", ["SRP009615"], Annotation.GENCODE_V29, (1377601, 13)),
        ("mouse", ["SRP111354"], Annotation.GENCODE_V23, (841916, 16)),
    ],
)
async def test_project_exon_accessor(
    organism, annotation, project, expected_counts_shape
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
    _, exon_counts = exon.load(dtype)

    assert exon_counts == expected_counts_shape
