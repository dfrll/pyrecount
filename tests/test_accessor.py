#! /usr/bin/env python3
import pytest
from os import path
from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation
import polars as pl

outpath = path.join('output', path.dirname(__file__))

@pytest.mark.parametrize('organism, project, expected_mm_shape, expected_shape', [
    ('human', 'SRP009615', (281448, 12), (281448, 10))
])
def test_project_jxn_accessor(organism, project, expected_shape, expected_mm_shape):
    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_jxn_accessor')
    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    metadata.cache()

    meta_dataframe = metadata.load().filter(pl.col('project') == project)

    jxn = Project(
        metadata = meta_dataframe,
        dbase = 'sra',
        dtype = Dtype.JXN,
        annotation = Annotation,
        cache_location = cache_location,
        jxn_format = 'ALL',
        root_url = root_url
    )

    jxn.cache()
    jxn_mm_dataframe, jxn_dataframe = jxn.load()

    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape

@pytest.mark.parametrize('organism, expected_shape', [
    ('human', (347005, 8))
])
def test_metadata_accessor(organism, expected_shape):
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_metadata_accessor')

    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 

    metadata.cache()
    meta_dataframe = metadata.load()

    assert meta_dataframe.shape == expected_shape

@pytest.mark.parametrize('organism, project, expected_shape', [
    ('human', 'SRP009615', (12, 173)),
])
def test_project_metadata_accessor(organism, project, expected_shape):
    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_metadata_accessor')
    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    metadata.cache()

    meta_dataframe = metadata.load().filter(pl.col('project') == project)

    project_metadata = Project(
        metadata = meta_dataframe,
        dbase = 'sra',
        dtype = Dtype.METADATA,
        annotation = Annotation,
        cache_location = cache_location,
        jxn_format = None,
        root_url = root_url
    )

    project_metadata.cache()
    project_meta_dataframe = project_metadata.load()

    assert project_meta_dataframe.shape == expected_shape
