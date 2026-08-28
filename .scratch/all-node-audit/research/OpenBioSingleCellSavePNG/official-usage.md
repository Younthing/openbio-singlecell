# Save PNG: official usage research

Researched: 2026-08-28

## Official image-validation usage

Pillow documents `Image.verify()` as a structural check for a broken image. It does not
decode pixels, and the image must be reopened before subsequent use. A durable PNG adapter
should therefore verify, reopen, fully `load()`, and confirm the detected format.

```python
from io import BytesIO
from PIL import Image

with Image.open(BytesIO(png_bytes)) as image:
    image.verify()
with Image.open(BytesIO(png_bytes)) as image:
    image.load()
    assert image.format == "PNG"
```

- Pillow `Image.verify`:
  https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.verify
- Pillow `Image.open`:
  https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.open

Python documents `os.replace` as atomic when successful on the same filesystem. Secure
staging files should be created in the destination directory, closed before replacement on
Windows, and removed after failure.

- `os.replace`: https://docs.python.org/3/library/os.html#os.replace
- `NamedTemporaryFile`:
  https://docs.python.org/3/library/tempfile.html#tempfile.NamedTemporaryFile

## OpenBio durable-image contract

- Accept only a typed `PlotResult` and verify that its bytes are a complete, decodable PNG.
- Preserve the original bytes exactly; validation must not re-encode or alter metadata.
- Write and flush a same-directory staging file, then commit atomically.
- A failed validation/write/commit must retain any previous destination byte-for-byte and
  remove staging material.
- `overwrite=False` must not replace a concurrently created file.

This adapter performs no statistical method and needs no scientific method citation.
