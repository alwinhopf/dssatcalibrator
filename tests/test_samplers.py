"""Tests for the parameter space and samplers."""
from pathlib import Path

import numpy as np
import pytest

from dssatcalibrator.config import load_config
from dssatcalibrator.spaces import ParameterSpace
from dssatcalibrator.samplers import sample

REPO = Path(__file__).resolve().parents[1]
HEMP_CFG = REPO / "config_hemp.yaml"


@pytest.fixture
def space():
    return ParameterSpace.from_config(load_config(HEMP_CFG))


def test_space_from_config(space):
    assert space.ndim > 0
    assert space.ndim == len(space.low) == len(space.high)
    assert space.ndim == len(space.start) == len(space.specs)
    assert space.names == [spec["name"] for spec in space.specs]
    assert (space.low < space.high).all()
    assert (space.low <= space.start).all() and (space.start <= space.high).all()


@pytest.mark.parametrize("engine", ["lhs", "sobol", "montecarlo", "grid"])
def test_sampler_bounds_and_start(space, engine):
    df = sample(space, n=32, engine=engine, seed=1, include_start=True)
    assert list(df.columns) == space.names
    # all draws within bounds
    assert (df.to_numpy() >= space.low - 1e-9).all()
    assert (df.to_numpy() <= space.high + 1e-9).all()
    # first row is the configured start point
    np.testing.assert_allclose(df.iloc[0].to_numpy(), space.start, rtol=1e-6)


def test_to_theta(space):
    df = sample(space, n=4, engine="lhs")
    theta = space.to_theta(df.iloc[1].to_numpy())
    assert set(theta) == set(space.names)
