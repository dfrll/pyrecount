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
