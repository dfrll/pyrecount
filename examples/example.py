#! /usr/bin/env python3
import polars as pl
from pyrecount.models import Dtype, Annotation
from pyrecount.accessor import Metadata, Project

cache_location = 'myproject'
organism = 'human'
dbase = 'sra'

# recount3 metadata
metadata = Metadata(organism=organism, cache_location=cache_location)
metadata.cache()

# load dataframe
meta_dataframe = metadata.load()

print(meta_dataframe)

# number of samples per project
n_sample_project = meta_dataframe.group_by('project').len()

print(n_sample_project)

# subset to avoid caching all recount3 data
data_of_interest = meta_dataframe.filter(
    # TODO: fix multi-project case
    #(pl.col('project') == 'SRP009615') | (pl.col('project') == 'SRP069357')
    (pl.col('project') == 'SRP009615')
)

print(data_of_interest)

# project metadata
####################################
dtype = Dtype.METADATA

project_metadata = Project(
    metadata = data_of_interest,
    dbase = dbase,
    dtype = dtype,
    cache_location = cache_location
)

project_metadata.cache()
project_meta_dataframe = project_metadata.load()

print(project_meta_dataframe)


## project junctions
#######################################
dtype = Dtype.JXN

jxn = Project(
    metadata = data_of_interest,
    dbase = dbase,
    dtype = dtype,
    cache_location = cache_location,
    jxn_format = 'ALL'
)

jxn.cache()
jxn_mm_dataframe, jxn_dataframe = jxn.load()

print(jxn_mm_dataframe)
print(jxn_dataframe)

# TODO: extract gene and exon counts
