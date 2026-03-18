"""Shared test fixtures for mofforge test suite."""

from pathlib import Path

import pytest

# Paths to test data
TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"
CRYSTAL_DIR = DATA_DIR / "crystals"
MOIETY_DIR = DATA_DIR / "moieties"


@pytest.fixture
def crystal_dir():
    """Path to test crystal CIF files."""
    return CRYSTAL_DIR


@pytest.fixture
def moiety_dir():
    """Path to test moiety XYZ files."""
    return MOIETY_DIR
