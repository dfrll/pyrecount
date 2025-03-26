#! /usr/bin/env python3
import polars as pl
from pyrecount.models import Dtype, Annotation, Extensions
from pyrecount.accessor import Metadata, Project

# TODO: extract counts data
# TODO: testing dbs other than sra
# TODO: handle BigWig

cache_location = 'myproject'
organism = 'human'
dbase = 'sra'

# recount3 metadata
recount_metadata = Metadata(organism=organism, cache_location=cache_location)
recount_metadata.cache()

# load dataframe
recount_meta_dataframe = recount_metadata.load()

print(recount_meta_dataframe)

# number of samples per project
n_sample_project = recount_meta_dataframe.group_by('project').len()

print(n_sample_project)

# subset
project_dataframe = recount_meta_dataframe.filter(
    (pl.col('project').is_in(['SRP009615']))
)

print(project_dataframe)


# project metadata
dtype = Dtype.METADATA

project_metadata = Project(
    metadata = project_dataframe,
    dbase = dbase,
    dtype = dtype,
    cache_location = cache_location
)

project_metadata.cache()
project_meta_dataframe = project_metadata.load()

print(project_meta_dataframe)

# project junctions
dtype = Dtype.JXN

project_jxn = Project(
    metadata = project_dataframe,
    dbase = dbase,
    dtype = dtype,
    cache_location = cache_location,
    jxn_format = 'ALL'
)

project_jxn.cache()
jxn_mm_dataframe, jxn_dataframe = project_jxn.load()

print(jxn_mm_dataframe)
print(jxn_dataframe)


# project gene annotation
dtype = Dtype.GENE
annotation = Annotation.GENCODE_V29

gene = Project(
    metadata = project_dataframe,
    dbase = dbase,
    dtype = dtype,
    cache_location = cache_location,
    annotation = annotation
)

gene.cache()
gene_dataframe = gene.load()

print(gene_dataframe)
