import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "dssatutils": "e9c859fa1d915623df23e2eb13084cb085dbfe3e",
    "dssatengine": "2280b11977ad373b9ae19d2d4497e8f276f7b133",
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

