#! /usr/bin/env python3
import aiohttp
import polars as pl


def replace_organism(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col("organism").replace(["Homo sapiens", "Mus musculus"], ["human", "mouse"])
    )


async def download_url_to_path(url: str, fpath: str):
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            with open(fpath, "wb") as fh:
                while True:
                    chunk = await response.content.read()
                    if not chunk:
                        break
                    fh.write(chunk)
