#! /usr/bin/env python3
from enum import Enum

HOMES_INDEX = 'homes_index'

class Annotation(Enum):
    GENCODE_V26 = 'G026'
    GENCODE_V29 = 'G029'
    ERCC = 'ERCC'
    FANTOM6_CAT = 'F006'
    REFSEQ = 'R109'
    SIRV = 'SIRV'

class Dtype(Enum):
    METADATA = 'metadata'
    GENE = 'gene_sums'
    EXON = 'exon_sums'
    JXN = 'junctions'
    BW = 'base_sums'

class Tags(Enum):
    METADATA = [
        'recount_project',
        'recount_qc',
        'recount_seq_qc',
        'recount_pred',
    ]
