if __package__:
    from .openbio_singlecell.extension import comfy_entrypoint
else:
    from openbio_singlecell.extension import comfy_entrypoint

WEB_DIRECTORY = "./web"

__all__ = ["WEB_DIRECTORY", "comfy_entrypoint"]
