"""Shared test fixtures / path constants.

Tests that need the local DSSAT install or the smoke-run outputs are skipped
automatically when those paths are absent, so the suite still runs on CI.
"""
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SMOKE = REPO / "_smoke"
HEMP = Path("C:/DSSAT48/Hemp")

TARGET_EXPERIMENTS = [
    "CNKU2101", "YUFE2101", "YUFE2201", "YUBA2201",
    "YUBA2101", "YUKU2101", "YUKU2201",
]


@pytest.fixture
def smoke_dir():
    if not (SMOKE / "PlantGro.OUT").exists():
        pytest.skip("smoke outputs not present (run the hemp smoke test first)")
    return SMOKE


@pytest.fixture
def hemp_dir():
    if not HEMP.exists():
        pytest.skip("local DSSAT hemp install not present")
    return HEMP


@pytest.fixture
def target_experiments():
    return list(TARGET_EXPERIMENTS)
