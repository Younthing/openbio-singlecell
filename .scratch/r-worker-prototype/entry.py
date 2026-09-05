"""THROWAWAY: register one R-backed operation before the ordinary Python entry."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1]), str(HERE)]

import operations_r  # noqa: E402,F401
from openbio_singlecell.worker_entry import main  # noqa: E402

raise SystemExit(main())
