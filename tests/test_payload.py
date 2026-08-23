from __future__ import annotations

from openbio_singlecell.contracts import SingleCellResult
from openbio_singlecell.payload import MAX_COLUMNS, MAX_ROWS, result_to_payload


class FakeFrame:
    def __init__(self, rows, columns):
        self._rows = rows
        self.columns = columns
        self.iloc = FakeIndexer(self)

    def __len__(self):
        return len(self._rows)

    def itertuples(self, index=False, name=None):
        return iter(tuple(row) for row in self._rows)


class FakeIndexer:
    def __init__(self, frame):
        self.frame = frame

    def __getitem__(self, key):
        row_slice, column_slice = key
        return FakeFrame(
            [row[column_slice] for row in self.frame._rows[row_slice]],
            self.frame.columns[column_slice],
        )


def make_result(**kwargs):
    values = dict(
        kind="summary",
        title="result",
        parameters={},
        description="description",
        warnings=[],
        input_cells=1,
        input_genes=1,
        random_seed=0,
        elapsed_seconds=0.25,
        source={},
    )
    values.update(kwargs)
    return SingleCellResult(**values)


def test_table_payload_is_bounded_and_normalizes_nonfinite_values():
    rows = [[row * 1000 + column for column in range(70)] for row in range(105)]
    rows[0][0] = float("nan")
    rows[0][1] = float("inf")
    table = FakeFrame(rows, [f"column_{index}" for index in range(70)])

    payload = result_to_payload(make_result(kind="table", table=table))

    assert len(payload["columns"]) == MAX_COLUMNS
    assert len(payload["rows"]) == MAX_ROWS
    assert payload["total_rows"] == 105
    assert payload["rows"][0][:2] == [None, None]
    assert any("Export CSV" in warning for warning in payload["warnings"])


def test_summary_collections_and_strings_are_bounded():
    summary = {
        "values": list(range(100)),
        "long": "x" * 5000,
        **{f"key_{index}": index for index in range(100)},
    }
    payload = result_to_payload(make_result(summary=summary))

    assert len(payload["summary"]["values"]) == 65
    assert len(payload["summary"]["long"]) <= 2000
    assert payload["summary"]["…"] == "truncated"
