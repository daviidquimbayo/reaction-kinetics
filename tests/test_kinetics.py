import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from kinetics import (activation_energy, arrhenius_k, concentration, consecutive, half_life,
                      intermediate_analytical, intermediate_max_time)


@pytest.mark.parametrize("order", [0, 1, 2])
def test_half_life_halves_concentration(order):
    a0, k = 2.0, 0.3
    assert math.isclose(concentration(order, a0, k, half_life(order, a0, k)), a0 / 2)


def test_first_order_half_life_independent_of_a0():
    assert half_life(1, 1.0, 0.2) == half_life(1, 5.0, 0.2)


def test_activation_energy_roundtrip():
    a, ea = 1e13, 75_000
    k1, k2 = arrhenius_k(a, ea, 300), arrhenius_k(a, ea, 350)
    assert math.isclose(activation_energy(k1, 300, k2, 350), ea)


def test_rate_increases_with_temperature():
    assert arrhenius_k(1e10, 50_000, 320) > arrhenius_k(1e10, 50_000, 300)


def test_consecutive_matches_analytical():
    ts, _, B, _ = consecutive(1.0, 0.5, 1.0, t_final=5)
    for t, b in zip(ts[::500], B[::500]):
        assert math.isclose(b, intermediate_analytical(1.0, 0.5, 1.0, t), abs_tol=1e-9)


def test_mass_conservation():
    _, A, B, C = consecutive(2.0, 0.3, 1.5, t_final=8)
    assert all(math.isclose(a + b + c, 1.5) for a, b, c in zip(A, B, C))


def test_intermediate_maximum():
    ts, _, B, _ = consecutive(1.0, 0.5, 1.0, t_final=5)
    t_max = ts[B.index(max(B))]
    assert math.isclose(t_max, intermediate_max_time(1.0, 0.5), abs_tol=1e-3)
