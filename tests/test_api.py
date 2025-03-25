#! /usr/bin/env python3
import pytest
from os import path
from pyrecount.api import EndpointConnector

def test_endpoint_connector():
    root_url, organism = 'http://duffel.rail.bio/recount3', 'human'

    expected = path.join(root_url, organism)

    endpoints = EndpointConnector(root_url=root_url, organism=organism)
    assert endpoints.root_organism_url == expected
