#! /usr/bin/env python3
from enum import Enum

HOMES_INDEX = "homes_index"


class Annotation(Enum):
    GENCODE_V29 = "G029"
    GENCODE_V26 = "G026"
    FANTOM6_CAT = "F006"
    REFSEQ = "R109"
    ERCC = "ERCC"
    SIRV = "SIRV"
    GENCODE_V23 = "M023"


class Dtype(Enum):
    METADATA = "metadata"
    GENE = "gene_sums"
    EXON = "exon_sums"
    JXN = "junctions"
    BW = "base_sums"


class Tags(Enum):
    RECOUNT_METADATA = "recount_project"
    METADATA = [
        "recount_project",
        "recount_qc",
        "recount_seq_qc",
        "recount_pred",
    ]


class Extensions(Enum):
    METADATA = ["MD.gz"]
    GENE = ["gtf.gz"]
    EXON = ["gtf.gz"]
    JXN = ["MM.gz", "RR.gz", "ID.gz"]
    BW = ["ALL.bw"]
