#! /usr/bin/env python3
import pytest
from os import path
from pyrecount.locator import ProjectLocator, MetadataLocator
from pyrecount.models import Dtype, Annotation, Tags

def test_locator_metadata():

    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'

    root_organism_url = path.join(root_url, organism)
    data_sources = {'sra': 'data_sources/sra', 'gtex': 'data_sources/gtex', 'tcga': 'data_sources/tcga'}

    metadata = MetadataLocator(root_organism_url=root_organism_url, data_sources=data_sources)

    assert len(metadata.fpaths) == 3

def test_project_locator_metadata():

    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'

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

    assert len(project_metadata.fpaths) == 10

def test_project_locator_jxn():

    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'

    root_organism_url = path.join(root_url, organism)
    data_sources = {'sra': 'data_sources/sra', 'gtex': 'data_sources/gtex', 'tcga': 'data_sources/tcga'}

    jxn=ProjectLocator(
        root_organism_url = root_organism_url,
        data_sources = data_sources,
        dbase='sra',
        project=['SRP009615', 'ERP110066'],
        dtype=Dtype.JXN,
        annotation=Annotation.GENCODE_V29,
        jxn_format='all'
    )

    assert len(jxn.fpaths) == 6
