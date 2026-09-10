from __future__ import annotations

from . import operations_abundance as _operations_abundance  # noqa: F401
from . import operations_annotation as _operations_annotation  # noqa: F401
from . import operations_cnv as _operations_cnv  # noqa: F401
from . import operations_communication as _operations_communication  # noqa: F401
from . import operations_correction as _operations_correction  # noqa: F401
from . import operations_data as _operations_data  # noqa: F401
from . import operations_differential as _operations_differential  # noqa: F401
from . import operations_embedding as _operations_embedding  # noqa: F401
from . import operations_enrichment as _operations_enrichment  # noqa: F401
from . import operations_factorization as _operations_factorization  # noqa: F401
from . import operations_integration as _operations_integration  # noqa: F401
from . import operations_lineage as _operations_lineage  # noqa: F401
from . import operations_monocle2 as _operations_monocle2  # noqa: F401
from . import operations_population as _operations_population  # noqa: F401
from . import operations_preprocess as _operations_preprocess  # noqa: F401
from . import operations_qc as _operations_qc  # noqa: F401
from . import operations_regulatory as _operations_regulatory  # noqa: F401
from . import operations_results as _operations_results  # noqa: F401
from . import operations_schist as _operations_schist  # noqa: F401
from . import operations_trajectory as _operations_trajectory  # noqa: F401
from . import operations_velocity as _operations_velocity  # noqa: F401
from .operations_input import load_h5ad
from .worker_protocol import register_operation

register_operation("openbio.node.loadh5ad", load_h5ad)

__all__: list[str] = []
