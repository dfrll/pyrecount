#! /usr/bin/env python3
import pytest
import polars as pl
from os import path

from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation

outpath = path.dirname(__file__)


@pytest.mark.xfail(reason='feature not implemented: multi-project jxn load support.')
@pytest.mark.parametrize('organism, project, expected_mm_shape, expected_shape', [
    ('human', ['SRP009615', 'SRP075759'], (562896, 43), (562896, 10)),
])
def test_multi_project_jxn_accessor(organism, project, expected_shape, expected_mm_shape):

    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_multi_project_jxn_accessor')

    recount_metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location)
    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col('project').is_in(project))
    )

    jxn = Project(
        metadata = project_meta_dataframe,
        dbase = 'sra',
        organism = organism,
        dtype = [Dtype.JXN],
        annotation = Annotation,
        cache_location = cache_location,
        jxn_format = 'all',
    )

    jxn.cache()
    jxn_mm_dataframe, jxn_dataframe = jxn.load()
    
    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape


@pytest.mark.xfail(reason='Project method _read_counts() needs to split first column to multiple columns.'
                   'This method passes after .load() method because polars ignores the column type in reading.')
@pytest.mark.parametrize('organism, project, annotation, expected_counts_shape', [
    ('human', ['SRP009615'], Annotation.GENCODE_V29, (64837, 13)),
    ('mouse', ['SRP111354'], Annotation.GENCODE_V23, (55421, 16))
])
def test_project_gene_accessor(organism, project, annotation, expected_counts_shape):

    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_gene_accessor')

    recount_metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location)
    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col('project').is_in(project))
    )

    gene = Project(
        metadata = project_meta_dataframe,
        dbase = 'sra',
        organism = organism,
        dtype = [Dtype.GENE],
        annotation = annotation,
        cache_location = cache_location,
        jxn_format = None,
        root_url = root_url
    )

    gene.cache()
    _, gene_counts = gene.load()

    # TODO: extract first column, assert doesn't work on list types:
    # polars.exceptions.InvalidOperationError: cannot cast List type (inner: 'Int64', to: 'String')
    assert gene_counts.shape == expected_counts_shape


@pytest.mark.xfail(reason='Project method _read_counts() needs to split first column to multiple columns.'
                   'This method passes after .load() method because polars ignores the column type in reading.')
@pytest.mark.parametrize('organism, project, annotation, expected_counts_shape', [
    ('human', ['SRP009615'], Annotation.GENCODE_V29, (1377601, 13)),
    ('mouse', ['SRP111354'], Annotation.GENCODE_V23, (841916, 16))
])
def test_project_exon_accessor(organism, annotation, project, expected_counts_shape):

    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_exon_accessor')

    recount_metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location)
    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col('project').is_in(project))
    )

    exon = Project(
        metadata = project_meta_dataframe,
        dbase = 'sra',
        organism = organism,
        dtype = [Dtype.EXON],
        annotation = annotation,
        cache_location = cache_location,
        jxn_format = None,
        root_url = root_url
    )

    exon.cache()
    _, exon_counts = exon.load()
    
    # TODO: extract first column, assert doesn't work on list types:
    # polars.exceptions.InvalidOperationError: cannot cast List type (inner: 'Int64', to: 'String')
    assert exon_counts == expected_counts_shape
