import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def config_copy(tmp_path):
    """A writable copy of config/ so tests can break it on purpose."""
    dst = tmp_path / "config"
    shutil.copytree(ROOT / "config", dst)
    return dst
