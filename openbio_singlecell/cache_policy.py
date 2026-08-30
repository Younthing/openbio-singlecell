from __future__ import annotations


def require_file_artifact_cache(args: object) -> str:
    """Reject ComfyUI's RAM-pressure cache, which cannot size file tickets."""
    if bool(getattr(args, "cache_classic", False)):
        return "classic"
    if int(getattr(args, "cache_lru", 0)) > 0:
        return "lru"
    if bool(getattr(args, "cache_none", False)):
        return "none"
    raise RuntimeError(
        "OpenBio file artifacts do not support ComfyUI's RAM-pressure cache. "
        "Restart with --cache-classic, --cache-none, or --cache-lru N."
    )


__all__ = ["require_file_artifact_cache"]
