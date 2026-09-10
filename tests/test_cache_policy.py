from types import SimpleNamespace

import pytest

from openbio_singlecell.cache_policy import check_file_artifact_cache


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (SimpleNamespace(cache_classic=True, cache_none=False, cache_lru=0), "classic"),
        (SimpleNamespace(cache_classic=False, cache_none=True, cache_lru=0), "none"),
        (SimpleNamespace(cache_classic=False, cache_none=False, cache_lru=12), "lru"),
    ],
)
def test_file_artifact_cache_accepts_supported_comfy_modes(args, expected):
    assert check_file_artifact_cache(args) == expected


def test_file_artifact_cache_allows_ram_pressure_default_with_launch_hint(caplog):
    args = SimpleNamespace(cache_classic=False, cache_none=False, cache_lru=0, cache_ram=[])

    assert check_file_artifact_cache(args) == "ram_pressure"
    assert "RAM-pressure" in caplog.text
    assert "disk space" in caplog.text
    assert "Consider --cache-classic, --cache-none, or --cache-lru N" in caplog.text


def test_file_artifact_cache_allows_explicit_ram_pressure_with_launch_hint(caplog):
    args = SimpleNamespace(cache_classic=False, cache_none=False, cache_lru=0, cache_ram=[2.0])

    assert check_file_artifact_cache(args) == "ram_pressure"
    assert "RAM-pressure" in caplog.text
    assert "Consider --cache-classic, --cache-none, or --cache-lru N" in caplog.text
