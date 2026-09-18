"""The protonation assignment is easy to get subtly wrong, so it is tested
from several directions: the two examples in the spec, totals, determinism,
and the site order inside a molecule."""
import pytest

from lipidmd.config import Config
from lipidmd.protonation import assign_protonation, n_protons, summarize

ECO_ORDER = ["amine_primary", "amine_tertiary"]
MC3_ORDER = ["amine_tertiary"]
N = 64


# ---- the two cases the spec calls out explicitly ---------------------------
def test_eco_50_percent_every_molecule_plus_one():
    states = assign_protonation(N, ECO_ORDER, 0.5, seed=1)
    assert len(states) == N
    assert all(s.charge == 1 for s in states)
    assert summarize(states) == {1: 64}


def test_mc3_50_percent_half_the_molecules_plus_one():
    states = assign_protonation(N, MC3_ORDER, 0.5, seed=1)
    assert summarize(states) == {0: 32, 1: 32}


# ---- the full matrix --------------------------------------------------------
@pytest.mark.parametrize("order,fraction,expected", [
    (ECO_ORDER, 0.0, {0: 64}),
    (ECO_ORDER, 0.5, {1: 64}),
    (ECO_ORDER, 1.0, {2: 64}),
    (MC3_ORDER, 0.0, {0: 64}),
    (MC3_ORDER, 0.5, {0: 32, 1: 32}),
    (MC3_ORDER, 1.0, {1: 64}),
])
def test_matrix(order, fraction, expected):
    assert summarize(assign_protonation(N, order, fraction, seed=7)) == expected


@pytest.mark.parametrize("order", [ECO_ORDER, MC3_ORDER])
@pytest.mark.parametrize("fraction", [0.0, 0.25, 0.5, 0.75, 1.0])
def test_fraction_of_amines_is_exact(order, fraction):
    """The definition: protonated amines / all amines == fraction."""
    states = assign_protonation(N, order, fraction, seed=3)
    protonated_amines = sum(len(s.protonated_sites) for s in states)
    assert protonated_amines == fraction * N * len(order)
    assert sum(s.charge for s in states) == protonated_amines


def test_even_spread_never_differs_by_more_than_one():
    states = assign_protonation(N, ECO_ORDER, 0.75, seed=3)  # 96 protons over 64 molecules
    assert summarize(states) == {1: 32, 2: 32}


# ---- which site / which molecule -------------------------------------------
def test_eco_plus_one_uses_first_site_in_order():
    states = assign_protonation(N, ECO_ORDER, 0.5, seed=1)
    assert all(s.protonated_sites == ("amine_primary",) for s in states)
    states = assign_protonation(N, ECO_ORDER[::-1], 0.5, seed=1)
    assert all(s.protonated_sites == ("amine_tertiary",) for s in states)


def test_eco_plus_two_has_both_sites():
    states = assign_protonation(N, ECO_ORDER, 1.0, seed=1)
    assert all(set(s.protonated_sites) == set(ECO_ORDER) for s in states)


def test_same_seed_same_molecules():
    a = assign_protonation(N, MC3_ORDER, 0.5, seed=42)
    b = assign_protonation(N, MC3_ORDER, 0.5, seed=42)
    assert a == b


def test_different_seed_different_molecules():
    a = [s.charge for s in assign_protonation(N, MC3_ORDER, 0.5, seed=1)]
    b = [s.charge for s in assign_protonation(N, MC3_ORDER, 0.5, seed=2)]
    assert a != b


def test_mc3_choice_is_not_trivially_first_half():
    charges = [s.charge for s in assign_protonation(N, MC3_ORDER, 0.5, seed=1)]
    assert charges != [1] * 32 + [0] * 32 and charges != [0] * 32 + [1] * 32


def test_indices_are_positions():
    states = assign_protonation(N, MC3_ORDER, 0.5, seed=1)
    assert [s.index for s in states] == list(range(N))


# ---- refusing bad input -----------------------------------------------------
def test_non_integer_proton_count_rejected():
    with pytest.raises(ValueError, match="whole number"):
        n_protons(63, 1, 0.5)          # 31.5 protons


def test_fraction_out_of_range():
    with pytest.raises(ValueError):
        n_protons(64, 1, 1.2)


def test_partial_protonation_needs_order():
    with pytest.raises(ValueError, match="protonation_order"):
        assign_protonation(N, None, 0.5, seed=1, n_sites=2)       # ECO at +1: order matters
    # ...but not when every molecule is fully (de)protonated
    assert summarize(assign_protonation(N, None, 1.0, seed=1, n_sites=2)) == {2: 64}


def test_matches_real_config_site_counts():
    cfg = Config()
    for lip in cfg.systems["lipids"]:
        n_sites = cfg.n_sites(lip)
        states = assign_protonation(cfg.systems["n_molecules"], [f"s{i}" for i in range(n_sites)],
                                    0.5, cfg.systems["seeds"]["protonation"])
        expected = {1: 64} if n_sites == 2 else {0: 32, 1: 32}
        assert summarize(states) == expected, lip
