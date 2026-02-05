#! /usr/bin/env python3
import polars as pl
import asyncio
from pyrecount.models import Dtype, Annotation
from pyrecount.accessor import Metadata, Project

organism = "human"
dbase = "sra"

# recount3 metadata
recount_metadata = Metadata(organism=organism)

recount_metadata.cache()

# load dataframe
recount_meta_dataframe = recount_metadata.load()

print(recount_meta_dataframe)

# number of samples per project
n_sample_project = (
    recount_meta_dataframe.group_by("project").len().sort(by="len", descending=True)
)

print(n_sample_project)

# subset
project_meta_dataframe = recount_meta_dataframe.filter(
    pl.col("project").is_in(["SRP009615"]) & pl.col("external_id").is_in(["SRR389077"])
)

print(project_meta_dataframe)


# project metadata
dtype = [Dtype.METADATA, Dtype.JXN, Dtype.GENE, Dtype.EXON, Dtype.BW]
annotation = Annotation.GENCODE_V29

project = Project(
    metadata=project_meta_dataframe,
    dbase=dbase,
    organism=organism,
    dtype=dtype,
    jxn_format="all",
    annotation=annotation,
)

asyncio.run(project.cache())

project_metadata = project.load(Dtype.METADATA)

print(project_metadata)

jxn_mm_dataframe, jxn_dataframe = project.load(Dtype.JXN)

print(jxn_dataframe)
print(jxn_mm_dataframe)

gene_annotation, gene_counts = project.load(Dtype.GENE)

print(gene_annotation)
print(gene_counts)

exon_annotation, exon_counts = project.load(Dtype.EXON)

print(exon_annotation)
print(exon_counts)

bw = project.load(Dtype.BW)
print(bw)
