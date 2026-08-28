# Preview Result: official usage research

Researched: 2026-08-28

## Official framework usage

ComfyUI V3 output nodes declare `is_output_node=True`. They are terminal side-effect/UI
adapters and may legally declare no graph outputs. Output nodes are always considered for
execution, while preview data is returned through the node's UI payload rather than through
a scientific data socket.

- ComfyUI custom-node server overview and output-node contract:
  https://docs.comfy.org/custom-nodes/backend/server_overview
- ComfyUI V3 migration/schema guide:
  https://docs.comfy.org/custom-nodes/v3_migration

For plot previews, Pillow's official `Image.verify()` API checks whether an image file is
broken without decoding pixels; the documentation requires reopening the image if it must
then be loaded. A robust PNG check therefore uses `verify()`, reopens the bytes, calls
`load()`, and confirms `format == "PNG"`.

- Pillow image verification:
  https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.verify

Python documents `os.replace(src, dst)` as a same-filesystem replacement operation that is
atomic when successful. A temporary file used for that operation must be created in the
destination directory, flushed, and closed before replacement, which is especially
important on Windows.

- `os.replace`: https://docs.python.org/3/library/os.html#os.replace
- `tempfile.NamedTemporaryFile`:
  https://docs.python.org/3/library/tempfile.html#tempfile.NamedTemporaryFile

## OpenBio preview contract

- A summary or table preview is serialized into the bounded OpenBio UI payload and creates
  no persistent output file.
- A plot preview uses one deterministic temporary filename per workflow/node identity, so a
  rerun replaces that node's previous preview without colliding with another workflow.
- Plot bytes must be a genuinely decodable PNG, not merely bytes with a PNG signature.
- Replacement must be atomic. A validation or write failure must leave the preceding preview
  intact and must remove the staging file.
- The UI path remains under ComfyUI's temporary output root.

This node is a framework/UI adapter. It does not perform a scientific method and therefore
does not require a method citation or analysis-package citation.
