#! /usr/bin/env python3
from os import path
from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation
import polars as pl

outpath = path.join('output', path.dirname(__file__))

def test_project_jxn_accessor():
    expected_shape = (281448, 10)
    expected_mm_shape = (281448, 12)

    # XXX: mock dataframe instead
    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'
    cache_location = path.join(outpath, 'test_project_jxn_accessor')
    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    metadata.cache()

    meta_dataframe = metadata.load().filter(pl.col('project') == 'SRP009615')

    jxn = Project(
        metadata = meta_dataframe,
        dbase = 'sra',
        dtype = Dtype.JXN,
        annotation = Annotation,
        cache_location = cache_location,
        jxn_format = 'ALL',
        root_url = 'http://duffel.rail.bio/recount3/'
    )

    jxn.cache()
    jxn_mm_dataframe, jxn_dataframe = jxn.load()

    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape

def test_metadata_accessor():

    expected_shape = (347005, 8)

    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'
    cache_location = path.join(outpath, 'test_metadata_accessor')

    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 

    metadata.cache()
    meta_dataframe = metadata.load()

    assert meta_dataframe.shape == expected_shape

def test_project_metadata_accessor():
    expected_shape = (12, 173)

    # XXX: mock dataframe instead
    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'
    cache_location = path.join(outpath, 'test_project_metadata_accessor')
    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    metadata.cache()

    meta_dataframe = metadata.load().filter(pl.col('project') == 'SRP009615')

    project_metadata = Project(
        metadata = meta_dataframe,
        dbase = 'sra',
        dtype = Dtype.METADATA,
        annotation = Annotation,
        cache_location = cache_location,
        jxn_format = None,
        root_url = 'http://duffel.rail.bio/recount3/'
    )

    project_metadata.cache()
    project_meta_dataframe = project_metadata.load()

    assert project_meta_dataframe.shape == expected_shape
