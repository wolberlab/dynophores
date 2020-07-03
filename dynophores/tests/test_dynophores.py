"""
Unit and regression test for the dynophores package.
"""

# Import package, test suite, and other packages as needed
import dynophores
import pytest
import sys

def test_dynophores_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "dynophores" in sys.modules
