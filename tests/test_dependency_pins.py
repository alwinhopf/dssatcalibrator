import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "dssatutils": "f728cd810923465360e4730d8893bb566dc0773c",
    "dssatengine": "14871db233af97a378fd42fb0c9749b5ff9b9453",
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

