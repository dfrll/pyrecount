#! /usr/bin/env python3
import polars as pl

def replace_organism(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col('organism').replace(['Homo sapiens', 'Mus musculus'], ['human', 'mouse'])
    )
