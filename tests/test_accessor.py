#! /usr/bin/env python3
import pytest
import polars as pl
from os import path

from pyrecount.accessor import Metadata, Project
from pyrecount.models import Dtype, Annotation

from pprint import pprint

outpath = path.dirname(__file__)

# TODO: add gene.get_read_counts to testing
# TODO: multi-project support for jxn, exon, gene
# TODO: transform raw counts
# TODO: profile transpose vs pivot in counts transform methods
# TODO: test dbs other than sra
# TODO: handle BigWig


@pytest.mark.parametrize('organism, expected_shape', [
    ('human', (347005, 8)),
    ('mouse', (416859, 8))
])
def test_recount_metadata_accessor(organism, expected_shape):

    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_recount_metadata_accessor')

    recount_metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location) 

    recount_metadata.cache()
    recount_meta_dataframe = recount_metadata.load()

    assert recount_meta_dataframe.shape == expected_shape


@pytest.mark.parametrize('organism, project, expected_mm_shape, expected_shape', [
    ('human', ['SRP009615'], (281448, 12), (281448, 10)),
    ('human', ['SRP075759'], (436480, 31), (436480, 10)),
    ('mouse', ['SRP111354'], (308775, 15), (308775, 10))
])
def test_project_jxn_accessor(organism, project, expected_shape, expected_mm_shape):

    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_jxn_accessor')

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
        cache_location = cache_location,
        annotation = Annotation,
        jxn_format = 'all',
        root_url = root_url
    )

    jxn.cache()
    jxn_mm_dataframe, jxn_dataframe = jxn.load()

    assert jxn_mm_dataframe.shape == expected_mm_shape
    assert jxn_dataframe.shape == expected_shape


@pytest.mark.parametrize('organism, project, expected_shape', [
    ('human', ['SRP009615'], (12, 173)),
    ('mouse', ['SRP111354'], (15, 175))
])
def test_project_metadata_accessor(organism, project, expected_shape):

    # XXX: mock dataframe instead
    root_url = 'http://duffel.rail.bio/recount3'
    cache_location = path.join(outpath, 'test_project_metadata_accessor')

    recount_metadata = Metadata(organism=organism, root_url=root_url, cache_location=cache_location)
    recount_metadata.cache()

    recount_meta_dataframe = recount_metadata.load()

    project_meta_dataframe = recount_meta_dataframe.filter(
        (pl.col('project').is_in(project))
    )

    project = Project(
        metadata = project_meta_dataframe,
        dbase = 'sra',
        organism = organism,
        dtype = [Dtype.METADATA],
        cache_location = cache_location,
        annotation = Annotation,
        jxn_format = None,
        root_url = root_url
    )

    project.cache()
    meta_dataframe = project.load()

    assert meta_dataframe.shape == expected_shape


@pytest.mark.parametrize('organism, project, annotation, expected_annotation_shape, expected_counts_shape', [
    ('human', ['SRP009615'], Annotation.GENCODE_V29, (1377600, 21), (1377601, 13)),
    ('mouse', ['SRP111354'], Annotation.GENCODE_V23, (841915, 21), (841916, 16))
])
def test_project_exon_accessor(organism, annotation, project, expected_annotation_shape, expected_counts_shape):

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
        cache_location = cache_location,
        annotation = annotation,
        jxn_format = None,
        root_url = root_url
    )

    exon.cache()
    exon_annotation, exon_counts = exon.load()
    
    assert exon_annotation.shape == expected_annotation_shape
    # see test_known_failure.py
    #assert exon_counts == expected_counts_shape


@pytest.mark.parametrize('organism, project, annotation, expected_annotation_shape, expected_counts_shape', [
    ('human', ['SRP009615'], Annotation.GENCODE_V29, (64836, 21), (64837, 13)),
    ('mouse', ['SRP111354'], Annotation.GENCODE_V23, (55420, 21), (55421, 16))
])
def test_project_gene_accessor(organism, project, annotation, expected_annotation_shape, expected_counts_shape):

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
        cache_location = cache_location,
        annotation = annotation,
        jxn_format = None,
        root_url = root_url
    )

    gene.cache()
    gene_annotation, gene_counts = gene.load()

    assert gene_annotation.shape == expected_annotation_shape
    # see test_known_failure.py
    #assert gene_counts.shape == expected_counts_shape
