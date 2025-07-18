# pyrecount

## Installation

You can install `pyrecount` either with with **pip** (recommended for general usage) or **poetry** (recommended for development).

### Option 1: Install with pip

```
pip install git+https://github.com/dfrll/pyrecount.git
```

### Option 2: Install with poetry

```
pip clone https://github.com/dfrll/pyrecount.git
cd pyrecount
poetry install
```

## Examples

``` python
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
    (pl.col("project").is_in(["SRP009615"]))
)

print(project_meta_dataframe)


# project metadata
dtype = [Dtype.METADATA, Dtype.JXN, Dtype.GENE, Dtype.EXON]
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
```

## Citations

This project uses data from the following publication:

Collado-Torres L (2025). _Explore and download data from the recount3
project_. doi:10.18129/B9.bioc.recount3
<https://doi.org/10.18129/B9.bioc.recount3>,
https://github.com/LieberInstitute/recount3 - R package version
1.16.0, <http://www.bioconductor.org/packages/recount3>.
