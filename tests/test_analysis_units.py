"""Pure-numpy pieces of the analysis that don't need a trajectory."""
import numpy as np
import pytest

pytest.importorskip("scipy")

from lipidmd.analysis.contacts import EventTracker
from lipidmd.analysis.convergence import block_standard_errors, running_average
from lipidmd.analysis.sasa import shrake_rupley
from lipidmd.analysis.water import PermeationCounter


def test_event_tracker_gap_tolerance():
    frames = [{"a"}, {"a"}, set(), {"a"}, set(), set(), set(), {"a"}]
    strict, tolerant = EventTracker(0), EventTracker(1)
    for k, f in enumerate(frames):
        strict.update(k, f)
        tolerant.update(k, f)
    strict.close()
    tolerant.close()
    assert [(s, e, c) for _, s, e, c in strict.events] == [(0, 1, False), (3, 3, False), (7, 7, True)]
    # a one-frame break is bridged; the three-frame break is not
    assert [(s, e, c) for _, s, e, c in tolerant.events] == [(0, 3, False), (7, 7, True)]


def test_permeation_counts_only_full_crossings():
    c = PermeationCounter(4, core_half_width=0.5, bulk_boundary=2.0)
    traj = np.array([
        # w0 crosses down, w1 enters core and returns, w2 wraps via PBC, w3 crosses up
        [3.0, 3.0, 3.0, -3.0],
        [1.0, 1.0, 3.5, -1.0],
        [0.0, 0.0, -3.5, 0.0],
        [-1.0, 1.0, -3.0, 1.0],
        [-3.0, 3.0, -3.0, 3.0],
    ])
    for t, dz in enumerate(traj):
        c.update(float(t), dz)
    assert sorted((i, d) for _, i, d in c.events) == [(0, -1), (3, 1)]


def test_sasa_isolated_sphere():
    r, p = 0.17, 0.14
    sasa = shrake_rupley(np.zeros((1, 3)), np.array([r]), p, n_points=500)
    assert sasa[0] == pytest.approx(4 * np.pi * (r + p) ** 2, rel=1e-6)


def test_sasa_two_overlapping_spheres():
    R = 0.3  # r + probe
    d = 0.3
    xyz = np.array([[0, 0, 0], [d, 0, 0]], float)
    sasa = shrake_rupley(xyz, np.array([0.16, 0.16]), 0.14, n_points=2000)
    cap = 2 * np.pi * R * (R - d / 2)            # area of each buried cap
    exact = 4 * np.pi * R**2 - cap
    assert sasa == pytest.approx([exact, exact], rel=0.02)


def test_sasa_occluder_and_periodic():
    xyz = np.array([[0.05, 5, 5]])
    occ = np.array([[9.95, 5, 5]])               # 0.1 nm away through the x boundary
    r = np.array([0.17])
    alone = shrake_rupley(xyz, r, 0.14, box=np.array([10.0, 10, 10]))
    buried = shrake_rupley(xyz, r, 0.14, box=np.array([10.0, 10, 10]), occluder_xyz=occ, occluder_radii=r)
    assert buried[0] < 0.7 * alone[0]


def test_running_average_and_blocking():
    assert running_average([1, 2, 3]).tolist() == [1, 1.5, 2]
    rng = np.random.default_rng(0)
    x = rng.normal(size=1024)
    b = block_standard_errors(x)
    assert b.block_frames.tolist()[:3] == [1, 2, 4]
    assert b["sem"].iloc[0] == pytest.approx(1 / 32, rel=0.2)   # white noise: sem = 1/sqrt(1024)
