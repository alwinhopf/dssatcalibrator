"""Shared test fixtures / path constants.

Tests that need the local DSSAT install or the smoke-run outputs are skipped
automatically when those paths are absent, so the suite still runs on CI.
"""
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SMOKE = REPO / "_smoke"
HEMP = Path("C:/Users/alwin/Documents/GitHub/DSSAT48Hemp/Hemp")

TARGET_EXPERIMENTS = [
    "CNKU2101", "YUFE2101", "YUFE2201", "YUBA2201",
    "YUBA2101", "YUKU2101", "YUKU2201",
]

HEMP_REFERENCE_FILES = [
    "CNKU2101.HMT",
    "YUFE2101.HMA",
    "YUFE2201.HMA",
    "YUBA2101.HMA",
    "YUKU2101.HMA",
    "YUKU2101.HMT",
    "YUKU2101.HMX",
    "YUKU2201.HMA",
]


@pytest.fixture
def smoke_dir():
    smoke = SMOKE if (SMOKE / "PlantGro.OUT").exists() else REPO / "tests" / "fixtures"
    if not (smoke / "PlantGro.OUT").exists():
        pytest.skip("smoke outputs not present (run the hemp smoke test first)")
    return smoke


@pytest.fixture
def hemp_dir():
    if not HEMP.exists():
        pytest.skip("local DSSAT hemp install not present")
    missing = [name for name in HEMP_REFERENCE_FILES if not (HEMP / name).exists()]
    if missing:
        pytest.skip("local DSSAT hemp reference files not present: " + ", ".join(missing))
    return HEMP


@pytest.fixture
def target_experiments():
    return list(TARGET_EXPERIMENTS)
