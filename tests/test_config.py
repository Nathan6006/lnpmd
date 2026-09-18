import pytest
import yaml

from lipidmd.config import Config, ConfigError


def edit(cfg_dir, name, fn):
    p = cfg_dir / f"{name}.yaml"
    data = yaml.safe_load(p.read_text())
    fn(data)
    p.write_text(yaml.safe_dump(data))


def test_real_config_loads():
    cfg = Config()
    names = [s.name for s in cfg.systems_list()]
    assert names == ["ECO_p000", "ECO_p050", "ECO_p100", "MC3_p000", "MC3_p050", "MC3_p100",
                     "SM102_p000", "SM102_p050", "SM102_p100"]
    assert len(cfg.runs()) == 27


def test_lipid_site_counts():
    cfg = Config()
    assert cfg.n_sites("ECO") == 2
    assert cfg.n_sites("MC3") == 1
    assert cfg.n_sites("SM102") == 1


def test_runs_seeds_and_filtering():
    cfg = Config()
    runs = cfg.runs(systems=["MC3_p050"], replicas=[2])
    assert len(runs) == 1 and runs[0].name == "MC3_p050/rep2"
    assert runs[0].velocity_seed == cfg.systems["seeds"]["velocities"][1]


def test_unknown_system_rejected():
    with pytest.raises(ConfigError):
        Config().runs(systems=["MC3_p075"])


def test_membrane_counts_refuse_to_round():
    cfg = Config()
    assert cfg.membrane_counts("popc") == {"POPC": 256}
    if cfg.membrane("popc_pops_chol")["counts_per_leaflet"] is None:
        with pytest.raises(ConfigError, match="non-integer"):
            cfg.membrane_counts("popc_pops_chol")


def test_counts_must_sum(config_copy):
    edit(config_copy, "membranes", lambda d: d["membranes"]["popc_pops_chol"].update(
        counts_per_leaflet={"POPC": 128, "POPS": 77, "CHL1": 50}))
    with pytest.raises(ConfigError, match="sum"):
        Config(config_copy).membrane_counts("popc_pops_chol")


def test_target_membrane_is_one_line(config_copy):
    edit(config_copy, "membranes", lambda d: d.update(target="popc"))
    assert Config(config_copy).target_membrane() == "popc"


def test_bad_fraction_rejected(config_copy):
    edit(config_copy, "systems", lambda d: d.update(protonation_levels=[0.0, 1.5]))
    with pytest.raises(ConfigError, match="fraction"):
        Config(config_copy)


def test_seed_count_must_match_replicas(config_copy):
    edit(config_copy, "systems", lambda d: d.update(replicas=4))
    with pytest.raises(ConfigError, match="velocities"):
        Config(config_copy)


def test_states_must_cover_every_charge(config_copy):
    edit(config_copy, "lipids", lambda d: d["ECO"]["states"].pop(2))
    with pytest.raises(ConfigError, match="states"):
        Config(config_copy)


def test_bad_protonation_order(config_copy):
    edit(config_copy, "lipids", lambda d: d["ECO"].update(protonation_order=["amine_primary"]))
    with pytest.raises(ConfigError, match="protonation_order"):
        Config(config_copy)


def test_replica_variation_changes_protonation_seed(config_copy):
    edit(config_copy, "systems", lambda d: d.update(replica_variation="velocities"))
    seeds = {r.protonation_seed for r in Config(config_copy).runs(systems=["MC3_p050"])}
    assert len(seeds) == 1
    edit(config_copy, "systems", lambda d: d.update(replica_variation="velocities+protonation"))
    seeds = {r.protonation_seed for r in Config(config_copy).runs(systems=["MC3_p050"])}
    assert len(seeds) == 3


def test_unresolved_lists_nulls_and_skips_optional():
    open_items = Config().unresolved()
    assert all(not p.startswith("hpc.account") for p in open_items)
    assert all(".pka_apparent." not in p for p in open_items)
    if Config().md["thermostat"] is None:
        assert "md.thermostat" in open_items


def test_require_names_the_file():
    cfg = Config()
    if cfg.md["thermostat"] is None:
        with pytest.raises(ConfigError, match="config/md.yaml"):
            cfg.require("md.thermostat")


def test_smoke_fallbacks_fill_only_nulls():
    cfg = Config()
    smoke = cfg.with_smoke_fallbacks()
    assert smoke.md["thermostat"] is not None
    assert smoke.md["temperature_K"] == cfg.md["temperature_K"]      # set values untouched
    if cfg.md["thermostat"] is None:
        assert "md.thermostat" in smoke.smoke_fallbacks_used
        assert cfg.md["thermostat"] is None                          # original not mutated
