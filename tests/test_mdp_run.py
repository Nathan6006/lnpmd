import pytest

from lipidmd.config import Config, ConfigError
from lipidmd.mdp import SMOKE_STEPS, render_mdp, stages
from lipidmd.run import last_log_time_ps, render_job_script, walltime_hours


@pytest.fixture
def smoke_cfg():
    return Config().with_smoke_fallbacks()


def test_charmm36_nonbonded_always_present(smoke_cfg):
    for st in stages(smoke_cfg, smoke=True):
        text = render_mdp(smoke_cfg, st, gen_seed=1, smoke=True)
        for line in ("coulombtype             = PME", "rcoulomb                = 1.2",
                     "vdw-modifier            = Force-switch", "rvdw-switch             = 1.0",
                     "rvdw                    = 1.2", "DispCorr                = no"):
            assert line in text, (st.name, line)


def test_stage_sequence_and_velocities(smoke_cfg):
    sts = stages(smoke_cfg, smoke=True)
    assert sts[0].kind == "em" and sts[-1].kind == "prod"
    assert all(s.nsteps == SMOKE_STEPS for s in sts)
    texts = {s.name: render_mdp(smoke_cfg, s, gen_seed=1234, smoke=True) for s in sts}
    assert "gen-seed                = 1234" in texts[sts[1].name]
    assert all("gen-vel                 = no" in texts[s.name] for s in sts[2:])


def test_restraint_defines(smoke_cfg):
    sts = stages(smoke_cfg, smoke=True)
    eq1 = render_mdp(smoke_cfg, sts[1], 1, smoke=True)
    assert "-DPOSRES -DPOSRES_FC_LIPID=1000" in eq1 and "-DPOSRES_AGG" in eq1
    prod = render_mdp(smoke_cfg, sts[-1], 1, smoke=True)
    assert "define                  =    ;" in prod  # no restraints in production


def test_semiisotropic_npt(smoke_cfg):
    prod = render_mdp(smoke_cfg, stages(smoke_cfg, smoke=True)[-1], 1)
    assert "pcoupltype              = semiisotropic" in prod
    assert "ref-t                   = 310 310" in prod


def test_real_config_refuses_open_decisions():
    cfg = Config()
    if cfg.md["equilibration"] is None:
        with pytest.raises(ConfigError):
            stages(cfg)


def test_walltime():
    assert walltime_hours("24:00:00") == 24
    assert walltime_hours("1-12:30:00") == pytest.approx(36.5)


def test_job_script_placeholders_and_chain(smoke_cfg, tmp_path):
    text = render_job_script(smoke_cfg, tmp_path, "job", stages(smoke_cfg, smoke=True), placeholders=True)
    if smoke_cfg.hpc["partition"] is None:
        assert "TODO_PARTITION" in text
    assert "-cpi $name.cpt -maxh $MAXH" in text
    assert "run_stage em em.mdp start.gro em" in text
    assert "run_stage prod prod.mdp eq2.gro md" in text
    if smoke_cfg.hpc["partition"] is None:
        with pytest.raises(ConfigError):
            render_job_script(smoke_cfg, tmp_path, "job", stages(smoke_cfg, smoke=True), placeholders=False)


def test_log_time_parser(tmp_path):
    log = tmp_path / "prod.log"
    log.write_text("""
           Step           Time
              0        0.00000

           Step           Time
         500000     1000.00000

   Energies (kJ/mol)
""")
    assert last_log_time_ps(log) == 1000.0
    assert last_log_time_ps(tmp_path / "missing.log") is None
