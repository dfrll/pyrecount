#! /usr/bin/env python3
import pytest
import logging
from os import path, getcwd
from pyrecount.locator import ProjectLocator, MetadataLocator
from pyrecount.models import Dtype, Annotation
from pyrecount.cache import QCache

outpath = path.join('output', path.dirname(__file__))

def test_cache_metadata():

    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'

    cache_location = path.join(outpath, 'test_cache_metadata')

    root_organism_url = path.join(root_url, organism)
    data_sources = {'sra': 'data_sources/sra', 'gtex': 'data_sources/gtex', 'tcga': 'data_sources/tcga'}

    m = MetadataLocator(root_organism_url=root_organism_url, data_sources=data_sources)

    cache = QCache(fpaths=m.fpaths, cache_location=cache_location)
    cache.biocache()

def test_cache_project_metadata():

    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'

    cache_location = path.join(outpath, 'test_cache_metadata')

    root_organism_url = path.join(root_url, organism)
    data_sources = {'sra': 'data_sources/sra', 'gtex': 'data_sources/gtex', 'tcga': 'data_sources/tcga'}

    project_metadata=ProjectLocator(
        root_organism_url = root_organism_url,
        data_sources = data_sources,
        dbase='sra',
        project=['SRP009615', 'ERP110066'],
        dtype=Dtype.METADATA,
        annotation=Annotation.GENCODE_V29
    )

    cache = QCache(fpaths = project_metadata.fpaths, cache_location=cache_location)
    cache.biocache()
 