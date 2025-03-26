#! /usr/bin/env python3
import pytest
from os import path
from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation
import polars as pl

outpath = path.join('output', path.dirname(__file__))

@pytest.mark.parametrize('organism, project, expected_shape', [
    ('human', ['SRP009615', 'ERP110066'], (1024, 173))
])
def test_project_metadata_accessor(organism, project, expected_shape):
    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_metadata_accessor')
    metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    metadata.cache()

    data_of_interest = metadata.load().filter(pl.col('project').is_in(project))

    project = Project(
        metadata = data_of_interest,
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

#@pytest.mark.parametrize('organism, project, expected_shape', [
    #('human', 'SRP009615', 0)
#])
#def test_project_jxn_accessor(organism, project, expected_shape):
    ## XXX: mock dataframe instead
    #root_url = 'http://duffel.rail.bio/recount3'
    #cache_location = path.join(outpath, 'test_project_jxn_accessor')
    #metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    #metadata.cache()

    #data_of_interest = metadata.load().filter(pl.col('project') == project)

    #gene = Project(
        #metadata = data_of_interest,
        #dbase = 'sra',
        #dtype = Dtype.GENE,
        #cache_location = cache_location,
        #annotation = Annotation.GENCODE_V29
    #)

    #gene.cache()
    #df = gene.load()

    ##assert

#@pytest.mark.parametrize('organism, project, expected_mm_shape, expected_shape', [
    #('human', ['SRP009615', 'ERP110066'], (281448, 12), (281448, 10))
#])
#def test_project_jxn_accessor(organism, project, expected_shape, expected_mm_shape):
    ## XXX: mock dataframe instead
    #root_url = 'http://duffel.rail.bio/recount3'
    #cache_location = path.join(outpath, 'test_project_jxn_accessor')
    #metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 
    #metadata.cache()

    #data_of_interest = metadata.load().filter(pl.col('project') == project)

    #jxn = Project(
        #metadata = data_of_interest,
        #dbase = 'sra',
        #dtype = Dtype.JXN,
        #annotation = Annotation,
        #cache_location = cache_location,
        #jxn_format = 'ALL',
        #root_url = root_url
    #)

    #jxn.cache()
    #jxn_mm_dataframe, jxn_dataframe = jxn.load()

    #assert jxn_mm_dataframe.shape == expected_mm_shape
    #assert jxn_dataframe.shape == expected_shape

#@pytest.mark.parametrize('organism, expected_shape', [
    #('human', (347005, 8))
#])
#def test_metadata_accessor(organism, expected_shape):
    #root_url = 'http://duffel.rail.bio/recount3'
    #cache_location = path.join(outpath, 'test_metadata_accessor')

    #metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 

    #metadata.cache()
    #meta_dataframe = metadata.load()

    #assert meta_dataframe.shape == expected_shape
