#! /usr/bin/env python3
import pytest
import polars as pl
from os import path

from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation

outpath = path.dirname(__file__)


@pytest.mark.xfail(reason='bug in multi-project support.')
@pytest.mark.parametrize('organism, project, expected_shape', [
    ('human', ['SRP009615', 'SRP075759'], (43, 179)),
])
def test_project_metadata_accessor(organism, project, expected_shape):

    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_metadata_accessor')

    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location)
    metadata.cache()

    meta_dataframe = metadata.load().filter(
        (pl.col('project').is_in(project))
    )

    project = Project(
        metadata = meta_dataframe,
        dbase = 'sra',
        dtype = Dtype.METADATA,
        annotation = Annotation,
        cache_location = cache_location,
        jxn_format = None,
        root_url = root_url
    )

    project.cache()
    meta_dataframe = project.load()

    assert meta_dataframe.shape == expected_shape
