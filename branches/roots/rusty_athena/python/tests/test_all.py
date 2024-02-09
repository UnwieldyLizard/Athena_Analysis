import pytest
import rusty_athena


def test_sum_as_string():
    assert rusty_athena.sum_as_string(1, 1) == "2"
