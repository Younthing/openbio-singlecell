from __future__ import annotations

import logging


def check_file_artifact_cache(args: object) -> str:
    """Report the cache mode and advise on RAM-pressure cache disk retention."""
    if bool(getattr(args, "cache_classic", False)):
        return "classic"
    if int(getattr(args, "cache_lru", 0)) > 0:
        return "lru"
    if bool(getattr(args, "cache_none", False)):
        return "none"
    logging.getLogger(__name__).warning(
        "OpenBio is using ComfyUI's RAM-pressure cache. "
        "RAM pressure does not account for disk space retained by file artifacts. "
        "Consider --cache-classic, --cache-none, or --cache-lru N to limit artifact retention."
    )
    return "ram_pressure"


__all__ = ["check_file_artifact_cache"]
