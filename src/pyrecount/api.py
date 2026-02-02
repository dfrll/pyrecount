#! /usr/bin/env python3
import time
import logging
from os import path
from typing import Optional, Dict
from .models import HOMES_INDEX
from requests import Response, get
from dataclasses import dataclass, field
from requests.exceptions import RequestException

log = logging.getLogger()


@dataclass
class EndpointConnector:
    organism: str
    root_url: str = "http://duffel.rail.bio/recount3"

    root_organism_url: str = field(init=False)
    data_sources: Dict[str, str] = field(init=False)

    def __post_init__(self):
        self.root_organism_url = path.join(self.root_url, self.organism)

        index = path.join(self.root_organism_url, HOMES_INDEX)
        resp = self._validate_endpoint(endpoint=index)
        if resp:
            self._set_data_sources(resp)

    def _set_data_sources(self, resp: Response):
        self.data_sources = {
            path.basename(dsource): dsource
            for dsource in resp.text.strip().splitlines()
        }

    def _validate_endpoint(self, endpoint: str) -> Optional[Response]:
        attempts = 3
        backoff = 2

        for attempt in range(1, attempts + 1):
            try:
                resp = get(endpoint, timeout=10)
                resp.raise_for_status()
                return resp

            except RequestException as e:
                log.warning(f"Attempt {attempt}/{attempts} failed for {endpoint}: {e}")

                if attempt < attempts:
                    time.sleep(backoff**attempt)
                else:
                    log.error(f"All retries failed for {endpoint}")
                    return None
