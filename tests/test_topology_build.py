import numpy as np
import pytest

from lipidmd import topology as top
from lipidmd.build import BuildError, add_proton, headgroup_planes, place_aggregate, reorder_to


def test_gro_roundtrip(tmp_path):
    res = [top.Residue("POPC", ["P", "C1"], np.array([[1.0, 2.0, 3.0], [1.1, 2.1, 3.1]])),
           top.Residue("TIP3", ["OH2", "H1", "H2"], np.zeros((3, 3)))]
    p = tmp_path / "x.gro"
    top.write_gro(p, res, np.array([5.0, 5.0, 8.0]))
    back, box = top.read_gro(p)
    assert [r.resname for r in back] == ["POPC", "TIP3"]
    assert back[0].names == ["P", "C1"]
    np.testing.assert_allclose(back[0].xyz, res[0].xyz, atol=1e-3)
    np.testing.assert_allclose(box, [5, 5, 8])


def test_run_length_keeps_order():
    res = [top.Residue(n, ["X"], np.zeros((1, 3))) for n in ["MC3", "MC3", "MC3H", "MC3", "TIP3", "TIP3"]]
    rle = top.run_length(res, {"MC3": "MC3", "MC3H": "MC3H", "TIP3": "TIP3"})
    assert rle == [("MC3", 2), ("MC3H", 1), ("MC3", 1), ("TIP3", 2)]


def test_parse_itp(tmp_path):
    p = tmp_path / "m.itp"
    p.write_text("""[ moleculetype ]
MC3H 3
[ atoms ]
1 NG3P1 1 MC3H N1 1 -0.3 14.007
2 HGP2  1 MC3H HN1 2 0.3 1.008 ; comment
3 CG331 1 MC3H C1 3 1.0 12.011
[ bonds ]
1 2 1
1 3 1
""")
    mt = top.parse_itp(p)["MC3H"]
    assert mt.atom_names == ["N1", "HN1", "C1"]
    assert mt.charge == pytest.approx(1.0)
    assert sorted(mt.neighbors("N1")) == ["C1", "HN1"]


def _prm(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_merge_prm_dedups_including_reversed(tmp_path):
    a = _prm(tmp_path, "a.prm", "[ bondtypes ]\nCG321 NG301 1 0.147 200000\n")
    b = _prm(tmp_path, "b.prm", "[ bondtypes ]\nNG301 CG321 1 0.147 200000\nCG321 NG3P1 1 0.148 190000\n")
    merged = top.merge_prm([a, b])
    assert merged.count("NG301") == 1 and "NG3P1" in merged


def test_merge_prm_conflict_raises(tmp_path):
    a = _prm(tmp_path, "a.prm", "[ angletypes ]\nA B C 5 110.0 300\n")
    b = _prm(tmp_path, "b.prm", "[ angletypes ]\nA B C 5 111.0 300\n")
    with pytest.raises(top.ParameterConflict):
        top.merge_prm([a, b])


def test_merge_prm_multiterm_dihedral(tmp_path):
    text = "[ dihedraltypes ]\nA B C D 9 0 1.0 1\nA B C D 9 180 0.5 2\n"
    merged = top.merge_prm([_prm(tmp_path, "a.prm", text), _prm(tmp_path, "b.prm", text)])
    assert merged.count("A  B  C  D") == 2  # both terms kept once


def test_add_proton_geometry():
    # tertiary amine: N at origin, three C in a pyramid below it (lone pair points +z)
    xyz = np.array([[0, 0, 0], [0.14, 0, -0.05], [-0.07, 0.12, -0.05], [-0.07, -0.12, -0.05]], float)
    res = top.Residue("MC3", ["N1", "C1", "C2", "C3"], xyz)
    out = add_proton(res, "N1", ["C1", "C2", "C3"], "HN1")
    h = out.xyz[out.index("HN1")]
    assert np.linalg.norm(h) == pytest.approx(0.101, abs=1e-6)
    assert h[2] > 0.09                     # points away from the carbons
    assert len(res.names) == 4             # original untouched


def test_add_proton_needs_three_neighbors():
    res = top.Residue("X", ["N1", "C1"], np.array([[0, 0, 0], [0.1, 0, 0]], float))
    with pytest.raises(BuildError):
        add_proton(res, "N1", ["C1"], "H")


def test_reorder_to_matches_topology_order():
    res = top.Residue("MC3", ["C1", "HN1", "N1"], np.arange(9, dtype=float).reshape(3, 3))
    out = reorder_to(res, ["N1", "HN1", "C1"], "MC3H")
    assert out.names == ["N1", "HN1", "C1"] and out.resname == "MC3H"
    np.testing.assert_allclose(out.xyz[0], [6, 7, 8])
    with pytest.raises(BuildError, match="don't match"):
        reorder_to(res, ["N1", "C1", "H9"], "MC3H")


def test_planes_and_placement():
    memb = [top.Residue("POPC", ["P"], np.array([[0, 0, z]], float)) for z in (2.0, 2.0, 6.0, 6.0)]
    planes = headgroup_planes(memb, {"POPC": "P"})
    assert planes["lower"] == 2.0 and planes["upper"] == 6.0
    agg = [top.Residue("MC3", ["A", "B"], np.array([[0, 0, 0], [1, 1, 2]], float))]
    placed = place_aggregate(agg, np.array([10.0, 10.0]), planes["upper"], gap=1.0)
    xyz = placed[0].xyz
    assert xyz[:, 2].min() == pytest.approx(7.0)
    np.testing.assert_allclose(xyz[:, :2].mean(0), [5, 5])
