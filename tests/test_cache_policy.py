from types import SimpleNamespace

import pytest

from openbio_singlecell.cache_policy import require_file_artifact_cache


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (SimpleNamespace(cache_classic=True, cache_none=False, cache_lru=0), "classic"),
        (SimpleNamespace(cache_classic=False, cache_none=True, cache_lru=0), "none"),
        (SimpleNamespace(cache_classic=False, cache_none=False, cache_lru=12), "lru"),
    ],
)
def test_file_artifact_cache_accepts_supported_comfy_modes(args, expected):
    assert require_file_artifact_cache(args) == expected


def test_file_artifact_cache_rejects_ram_pressure_default_with_launch_hint():
    args = SimpleNamespace(cache_classic=False, cache_none=False, cache_lru=0, cache_ram=[])

    with pytest.raises(RuntimeError, match=r"RAM-pressure.*--cache-classic"):
        require_file_artifact_cache(args)


def test_file_artifact_cache_rejects_explicit_ram_pressure_too():
    args = SimpleNamespace(cache_classic=False, cache_none=False, cache_lru=0, cache_ram=[2.0])

    with pytest.raises(RuntimeError, match=r"RAM-pressure.*--cache-classic"):
        require_file_artifact_cache(args)
