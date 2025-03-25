##! /usr/bin/env python3
#import pytest
#from os import path, getcwd
#from pyrecount.locator import ProjectLocator, MetadataLocator
#from pyrecount.models import Dtype, Annotation
#from pyrecount.cache import QCache

#outpath = path.join('output', path.dirname(__file__))


#@pytest.mark.parametrize('organism, expected_shape', [
    #('human', 3),
    ##('mouse', None)
#])
#def test_cache_metadata(organism, expected_shape):
    #root_url = 'http://duffel.rail.bio/recount3'
    #root_organism_url = path.join(root_url, organism)
    #cache_location = path.join(outpath, 'test_cache_metadata')

    #data_sources = {
        #'sra': 'data_sources/sra',
        #'gtex': 'data_sources/gtex',
        #'tcga': 'data_sources/tcga'
    #}

    #metadata = MetadataLocator(root_organism_url=root_organism_url, data_sources=data_sources)

    #assert len(metadata.fpaths) == expected_shape

    #cache = QCache(fpaths=metadata.fpaths, cache_location=cache_location)
    #cache.biocache()

#@pytest.mark.parametrize('organism, project, expected_shape', [
    #('human', ['SRP009615', 'ERP110066'], 10)
#])
#def test_cache_project_metadata(organism, project, expected_shape):
    #root_url = 'http://duffel.rail.bio/recount3'
    #root_organism_url = path.join(root_url, organism)
    #cache_location = path.join(outpath, 'test_cache_metadata')

    #data_sources = {
        #'sra': 'data_sources/sra',
        #'gtex': 'data_sources/gtex',
        #'tcga': 'data_sources/tcga'
    #}

    #project_metadata=ProjectLocator(
        #root_organism_url = root_organism_url,
        #data_sources = data_sources,
        #dbase='sra',
        #project=project,
        #dtype=Dtype.METADATA,
        #annotation=Annotation.GENCODE_V29
    #)

    #print(project_metadata.fpaths)

    #assert len(project_metadata.fpaths) == expected_shape

    #cache = QCache(fpaths = project_metadata.fpaths, cache_location=cache_location)
    #cache.biocache()
 