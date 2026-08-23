from __future__ import annotations

from comfy_api.latest import ComfyExtension, io

NODE_CLASSES: list[type[io.ComfyNode]] = []


class OpenBioSingleCellExtension(ComfyExtension):
    async def get_node_list(self) -> list[type[io.ComfyNode]]:
        return NODE_CLASSES


async def comfy_entrypoint() -> OpenBioSingleCellExtension:
    return OpenBioSingleCellExtension()


__all__ = ["NODE_CLASSES", "OpenBioSingleCellExtension", "comfy_entrypoint"]
