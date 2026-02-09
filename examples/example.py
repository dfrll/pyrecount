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
    pl.col("project").is_in(["SRP009615"])
)

print(project_meta_dataframe)


# project metadata
annotation = Annotation.GENCODE_V29

project = Project(
    metadata=project_meta_dataframe,
    dbase=dbase,
    organism=organism,
    jxn_format="all",
    annotation=annotation,
)

dtypes = [Dtype.METADATA, Dtype.JXN, Dtype.GENE, Dtype.EXON, Dtype.BW]

asyncio.run(project.cache(dtypes))

project_metadata = project.load(Dtype.METADATA)

print(project_metadata)

gene_annotation, gene_counts = project.load(Dtype.GENE)

print(gene_annotation)
print(gene_counts)

scaled_counts = project.scale_mapped_reads(
    gene_counts,
    target_size=4e7,
    L=100,
)

print(scaled_counts)

scaled_counts = project.scale_auc(
    gene_counts,
    target_size=4e7,
)

print(scaled_counts)

jxn_mm_dataframe, jxn_dataframe = project.load(Dtype.JXN)

print(jxn_dataframe)
print(jxn_mm_dataframe)

exon_annotation, exon_counts = project.load(Dtype.EXON)

print(exon_annotation)
print(exon_counts)

bw = project.load(Dtype.BW)
print(bw)
