from pathlib import Path

import pytest

from lipidmd.params import (atoms_within, check_penalties, expected_formula, format_penalty_report,
                            hill_formula, parse_formula, parse_str)

STR = Path(__file__).parent / "fixtures" / "toy_amine.str"


@pytest.fixture
def mol():
    return parse_str(STR)


def test_parse_header_atoms_bonds(mol):
    assert mol.resname == "TOYH"
    assert mol.net_charge == 1.0
    assert mol.param_penalty == 70.0 and mol.charge_penalty == 25.0
    assert len(mol.atom_names) == 12
    assert mol.atom_types["N1"] == "NG3P1"
    assert mol.atom_charge_penalties["C7"] == 15.0
    assert ("N1", "C10") in mol.bonds and ("C7", "C8") in mol.bonds
    assert len(mol.bonds) == 11


def test_parse_parameters(mol):
    kinds = [(p.kind, p.penalty) for p in mol.parameters]
    assert kinds == [("bond", 12.0), ("bond", 3.0), ("angle", 55.0), ("dihedral", 61.0), ("dihedral", 70.0)]
    assert mol.parameters[3].types == ("CG321", "CG321", "CG321", "CG331")


def test_atoms_within(mol):
    assert atoms_within(mol, ["N1"], 0) == {"N1"}
    assert atoms_within(mol, ["N1"], 1) == {"N1", "HN1", "C1", "C9", "C10"}
    assert "C3" in atoms_within(mol, ["N1"], 3) and "C4" not in atoms_within(mol, ["N1"], 3)


def test_threshold_and_ordering(mol):
    items = check_penalties(mol, ["N1"], threshold=10)
    assert [it.penalty for it in items] == [70.0, 61.0, 55.0, 25.0, 15.0, 12.0]
    assert all(it.penalty > 10 for it in items)          # the 3.0 bond and 2.0 charge are dropped
    assert len(check_penalties(mol, ["N1"], threshold=50)) == 3


def test_severity(mol):
    items = {it.penalty: it for it in check_penalties(mol, ["N1"])}
    assert items[70.0].severity == "ATTENTION"
    assert items[55.0].severity == "ATTENTION"
    assert items[25.0].severity == "review"
    assert items[12.0].severity == "review"


def test_near_amine_flags(mol):
    items = {it.penalty: it for it in check_penalties(mol, ["N1"], near_bonds=3)}
    assert items[70.0].near_amine          # C9-N1-C1-C2 (matched in reverse order)
    assert not items[61.0].near_amine      # C5-C6-C7-C8, far down the tail
    assert items[55.0].near_amine          # includes C1-C2-C3
    assert items[25.0].near_amine          # the N1 charge itself
    assert not items[15.0].near_amine      # C7 charge
    assert items[12.0].near_amine          # C1-N1 bond


def test_instances_resolved_to_atoms(mol):
    items = {it.penalty: it for it in check_penalties(mol, ["N1"])}
    assert items[61.0].instances == [("C5", "C6", "C7", "C8")]
    assert {frozenset(i) for i in items[70.0].instances} == {frozenset({"C9", "N1", "C1", "C2"}),
                                                             frozenset({"C10", "N1", "C1", "C2"})}
    assert len(items[55.0].instances) == 5  # C1-C2-C3 ... C5-C6-C7


def test_no_amine_names_means_no_near_flags(mol):
    items = check_penalties(mol, [])
    assert not any(it.near_amine for it in items)
    text = format_penalty_report(mol, items, [])
    assert "NOT checked" in text


def test_report_text(mol):
    items = check_penalties(mol, ["N1"])
    text = format_penalty_report(mol, items, ["N1"])
    assert "ATTENTION" in text and ">>>" in text
    assert "6 above threshold, 3 need manual attention, 4 near an ionizable N" in text


def test_missing_amine_warned(mol):
    text = format_penalty_report(mol, check_penalties(mol, ["NX"]), ["NX"])
    assert "not found in this .str" in text


def test_not_a_str_file(tmp_path):
    p = tmp_path / "x.str"
    p.write_text("* nothing\nEND\n")
    with pytest.raises(ValueError, match="RESI"):
        parse_str(p)


def test_formula_helpers():
    assert hill_formula(parse_formula("C43H79NO2")) == "C43H79NO2"
    assert expected_formula("C43H79NO2", 1) == "C43H80NO2"
    assert hill_formula({"O": 1, "H": 2}) == "H2O"
