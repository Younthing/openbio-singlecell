from __future__ import annotations

from openbio_singlecell.analysis_utils import matrix_totals_and_nonzero


def test_sparse_matrix_totals_and_nonzero_operates_on_worker_owned_input(science):
    class WorkerOwnedCSR(science.sparse.csr_matrix):
        def copy(self, order="C"):
            raise AssertionError("worker-owned sparse input must not be copied")

    matrix = WorkerOwnedCSR(
        (
            science.np.array([1.0, 2.0, 0.0, 4.0]),
            science.np.array([0, 0, 1, 2]),
            science.np.array([0, 3, 4]),
        ),
        shape=(2, 3),
    )

    totals, nonzero = matrix_totals_and_nonzero(matrix, axis=1)

    assert totals.tolist() == [3.0, 4.0]
    assert nonzero.tolist() == [1, 1]
