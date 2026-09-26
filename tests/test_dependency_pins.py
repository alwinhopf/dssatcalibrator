import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "dssatutils": "a4202fbc6377a62b391340c317c92c96d157e031",
    "dssatengine": "2f20fd8c4afcf5f116fa5c6ed8aea086069bfa57",
}


def test_python_and_r_shared_dependency_pins_match():
    python_manifest = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    r_manifest = (ROOT / "DESCRIPTION").read_text(encoding="utf-8")
    for package, revision in EXPECTED.items():
        pattern = rf"alwinhopf/{package}(?:[.]git)?@([0-9a-f]{{40}})"
        python_match = re.search(pattern, python_manifest)
        r_match = re.search(pattern, r_manifest)
        assert python_match and r_match, f"{package} must use immutable pins"
        assert python_match.group(1) == revision
        assert r_match.group(1) == revision

