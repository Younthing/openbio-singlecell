# Third-party notices

`openbio-singlecell` is based on and designed to run with the following upstream projects. Their copyright and license terms remain in effect. The complete GNU GPL version 3 text is included in `LICENSE`.

## ComfyUI

- Project: https://github.com/comfyanonymous/ComfyUI
- License: GNU General Public License version 3
- The paired source baseline is recorded in `release_manifest.json`.

## ComfyUI_frontend

- Project: https://github.com/Comfy-Org/ComfyUI_frontend
- License: GNU General Public License version 3 only
- The paired source baseline is recorded in `release_manifest.json`.
- A prebuilt paired frontend must include `dist/LICENSE` and the generated `dist/THIRD_PARTY_NOTICES.md`. Generate and verify them with `scripts/generate_frontend_notices.mjs` after every build. The generated notice inventories the frozen production dependency graph and deduplicates identical license texts.

## Scanpy

- Project: https://github.com/scverse/scanpy
- License: BSD 3-Clause License
- Required version: `>=1.12.3,<1.13`, including its `leiden` extra

## AnnData

- Project: https://github.com/scverse/anndata
- License: BSD 3-Clause License
- Required version: `>=0.13.2,<0.14`

Python packages installed as transitive dependencies are not vendored by this repository and retain their own notices and license terms. Redistributors must include the applicable upstream license texts and notices with the packaged artifacts. This source-package notice is not a substitute for the generated notice shipped with a prebuilt frontend.
