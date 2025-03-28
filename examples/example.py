#! /usr/bin/env python3
import polars as pl
from pyrecount.models import Dtype, Annotation
from pyrecount.accessor import Metadata, Project


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
project_meta_dataframe = recount_meta_dataframe.filter(
    (pl.col('project').is_in(['SRP009615']))
)

print(project_meta_dataframe)


# project metadata
dtype = [Dtype.METADATA]

project = Project(
    metadata = project_meta_dataframe,
    dbase = dbase,
    organism = organism,
    dtype = dtype,
    cache_location = cache_location
)

project.cache()
project_dataframe = project.load()

print(project_dataframe)


# project junctions
dtype = [Dtype.JXN]

project_jxn = Project(
    metadata = project_dataframe,
    dbase = dbase,
    organism = organism,
    dtype = dtype,
    cache_location = cache_location,
    jxn_format = 'all'
)

project_jxn.cache()
jxn_mm_dataframe, jxn_dataframe = project_jxn.load()

print(jxn_mm_dataframe)
print(jxn_dataframe)


# project exon annotation, raw coverage counts
dtype = [Dtype.EXON]
annotation = Annotation.GENCODE_V29

exon = Project(
    metadata = project_dataframe,
    dbase = dbase,
    organism = organism,
    dtype = dtype,
    cache_location = cache_location,
    annotation = annotation
)

exon.cache()
exon_annotation, exon_raw_counts = exon.load()

print(exon_annotation)
print(exon_raw_counts)


# project gene annotation, raw coverage counts
dtype = [Dtype.GENE]
annotation = Annotation.GENCODE_V29

gene = Project(
    metadata = project_dataframe,
    dbase = dbase,
    organism = organism,
    dtype = dtype,
    cache_location = cache_location,
    annotation = annotation
)

gene.cache()
gene_annotation, gene_raw_counts = gene.load()

print(gene_annotation)
print(gene_raw_counts)
