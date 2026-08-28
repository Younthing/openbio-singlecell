# Filter Predicted Doublets: official usage research

Researched: 2026-08-28

## Upstream official contract

Scanpy documents `obs['predicted_doublet']` as the boolean output of `scanpy.pp.scrublet`; `True` indicates a predicted doublet and `False` a retained singlet candidate.

- Scanpy Scrublet API: https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pp.scrublet.html
- Scrublet project best practices: https://github.com/swolock/scrublet
- Wolock SL, Lopez R, Klein AM. *Cell Systems*. 2019;8(4):281-291.e9. https://doi.org/10.1016/j.cels.2018.11.005

The original project explicitly recommends inspecting threshold separation and embedding localization. Therefore filtering should not be silently fused into detection.

## Scientific and data-contract implications

- Subset with the inverse of an actual boolean prediction column.
- Never coerce strings or general objects with `astype(bool)`: values such as `"False"` are non-empty strings and become `True`, removing the wrong cells.
- Reject missing predictions. Nullable booleans with missing values require an upstream decision rather than implicit retention/removal.
- Report input, retained and removed cells, predicted rate, prediction column, unchanged gene count, warnings, and upstream-method limitations.
- Allow an all-empty result but warn prominently; also warn when no cells are predicted, so either extreme is visible rather than silently presented as routine filtering.

AnnData supports a zero-observation slice. Keeping at least one cell is a downstream recommendation, not a requirement of this atomic boolean-mask operation.
