#! /usr/bin/env python3
import logging
from typing import List
from .models import Dtype
from pybiocfilecache.cache import BiocFileCache
from concurrent.futures import ThreadPoolExecutor

log = logging.getLogger()

class QCache():
    def __init__(self, organism_cache_location: str, fpaths: List[str]):
        self.fpaths: List[str] = fpaths
        self.max_workers: int = len(self.fpaths)
        self.cache: BiocFileCache = BiocFileCache(organism_cache_location)

    def _add_to_bfc(self, fpath: str) -> None:
        # fpath: path to the source file.
        try:
            self.cache.add(rname=fpath, fpath=fpath, rtype='web', download=True)
        except Exception as e:
            log.error(f'Error caching {fpath}: {e}')

    def biocache(self) -> None:
        log.info(f'Caching {Dtype.METADATA.value} with ThreadPoolExecutor.')
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for fpath in self.fpaths:
                if self.cache.get(fpath) is not None:
                    continue
                log.info(f'Caching {fpath}')
                executor.submit(self._add_to_bfc, fpath)

    def biocache_serial(self) -> None:
        log.info(f'Caching {Dtype.METADATA.value} serially.')
        for fpath in self.fpaths:
            self._add_to_bfc(fpath)
