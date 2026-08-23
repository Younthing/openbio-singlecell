from __future__ import annotations

import pytest

from install import main as manager_install
from openbio_singlecell import dependencies
from scripts.generate_demo import KNOWN_MARKERS, build_demo, validate_existing

pytestmark = pytest.mark.skipif(
    not dependencies.AVAILABLE,
    reason="Scientific dependencies are unavailable.",
)


def test_manager_install_generates_demo(tmp_path, monkeypatch):
    comfy_root = tmp_path / "ComfyUI"
    comfy_root.mkdir()
    (comfy_root / "main.py").touch()
    monkeypatch.delenv("OPENBIO_COMFYUI_ROOT", raising=False)
    monkeypatch.setenv("COMFYUI_FOLDERS_BASE_PATH", str(comfy_root))

    assert manager_install() == 0

    output = comfy_root / "input" / "openbio-singlecell" / "openbio_singlecell_demo.h5ad"
    valid, reason = validate_existing(output, dependencies.ad, dependencies.sparse)
    assert valid, reason


def test_demo_validation_rejects_forged_same_shape_anndata(tmp_path):
    demo = build_demo(dependencies.np, dependencies.pd, dependencies.sparse, dependencies.ad)
    valid_path = tmp_path / "valid.h5ad"
    demo.write_h5ad(valid_path)
    assert demo.uns["openbio_singlecell"]["display_name"] == "OpenBio single-cell demo"
    assert demo.uns["openbio_singlecell"]["analysis_history"] == {}
    assert validate_existing(valid_path, dependencies.ad, dependencies.sparse) == (True, "valid")

    cases = [
        ("markers", "known marker genes"),
        ("groups", "exactly 200"),
        ("samples", "samples do not match"),
        ("batches", "batches do not match"),
        ("fractional", "integer counts"),
        ("negative", "negative counts"),
    ]
    for issue, expected_reason in cases:
        forged = demo.copy()
        if issue == "markers":
            names = list(forged.var_names)
            names[names.index(next(iter(KNOWN_MARKERS["T cell"])))] = "FORGED_MARKER"
            forged.var_names = names
        elif issue == "groups":
            groups = forged.obs["cell_type"].astype(str).to_numpy().copy()
            groups[0] = next(group for group in KNOWN_MARKERS if group != groups[0])
            forged.obs["cell_type"] = dependencies.pd.Categorical(groups, categories=list(KNOWN_MARKERS))
        elif issue == "samples":
            forged.obs["sample"] = dependencies.pd.Categorical(["sample_1"] * forged.n_obs)
        elif issue == "batches":
            forged.obs["batch"] = dependencies.pd.Categorical(["batch_1"] * forged.n_obs)
        elif issue == "fractional":
            forged.X = forged.X.astype(float)
            forged.X.data[0] = 0.5
        elif issue == "negative":
            forged.X = forged.X.copy()
            forged.X.data[0] = -1

        path = tmp_path / f"forged_{issue}.h5ad"
        forged.write_h5ad(path)
        valid, reason = validate_existing(path, dependencies.ad, dependencies.sparse)
        assert not valid
        assert expected_reason in reason
