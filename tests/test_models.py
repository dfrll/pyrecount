#! /usr/bin/env python3
import pytest
from pyrecount.models import Dtype, Annotation, Tags

@pytest.mark.parametrize('test_dtype, expected', [
    (Dtype.METADATA, 'metadata'),
    (Dtype.JXN, 'junctions')
])
def test_dtypes(test_dtype, expected):
    assert test_dtype.value == expected
