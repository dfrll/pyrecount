#! /usr/bin/env python3
import polars as pl
from pyrecount.models import Dtype, Annotation
from pyrecount.accessor import Metadata, Project

# TODO: extract gene and exon counts

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

#print(n_sample_project)

# subset to avoid caching all recount3 data
data_of_interest = meta_dataframe.filter(
    (pl.col('project').is_in(['SRP009615', 'ERP110066']))
)

print(data_of_interest)

# project metadata
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

### project junctions
#dtype = Dtype.JXN

#jxn = Project(
    #metadata = data_of_interest,
    #dbase = dbase,
    #dtype = dtype,
    #cache_location = cache_location,
    #jxn_format = 'ALL'
#)

#jxn.cache()
#jxn_mm_dataframe, jxn_dataframe = jxn.load()

#print(jxn_mm_dataframe)
#print(jxn_dataframe)

#dtype = Dtype.BW

#bw = Project(
    #metadata = data_of_interest,
    #dbase = dbase,
    #dtype = dtype,
    #cache_location = cache_location,
#)

#bw.cache()
#df = bw.load()

#print(df)

#dtype = Dtype.GENE
#annotation = Annotation.GENCODE_V29

#print(dtype)

#gene = Project(
    #metadata = data_of_interest,
    #dbase = dbase,
    #dtype = dtype,
    #cache_location = cache_location,
    #annotation = annotation
#)

#gene.cache()
#df = gene.load()

#print(df)
