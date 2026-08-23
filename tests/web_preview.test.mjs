import assert from "node:assert/strict";
import test from "node:test";

import { normalizeSingleCellPayload } from "../web/single_cell_result_renderer.mjs";

test("accepts each discriminated result kind", () => {
    for (const kind of ["summary", "table", "plot"]) {
        const payload = normalizeSingleCellPayload({
            schema_version: 1,
            kind,
            title: "result",
            columns: ["gene"],
            rows: [["CD3D"]],
            total_rows: 1,
            summary: { cells: 10 },
            warnings: [],
            elapsed_seconds: 0.25,
        });

        assert.equal(payload.kind, kind);
        assert.equal(payload.schema_version, 1);
    }
});

test("rejects unknown schemas and kinds", () => {
    assert.equal(normalizeSingleCellPayload(null), null);
    assert.equal(normalizeSingleCellPayload({ schema_version: 2, kind: "summary" }), null);
    assert.equal(normalizeSingleCellPayload({ schema_version: 1, kind: "unknown" }), null);
});

test("normalizes untrusted collection fields", () => {
    const payload = normalizeSingleCellPayload({
        schema_version: 1,
        kind: "table",
        title: null,
        columns: "gene",
        rows: ["not-a-row", ["CD3D", null]],
        total_rows: -1,
        warnings: "warning",
        elapsed_seconds: Number.NaN,
    });

    assert.deepEqual(payload.columns, []);
    assert.deepEqual(payload.rows, [["CD3D", ""]]);
    assert.equal(payload.total_rows, 0);
    assert.deepEqual(payload.warnings, []);
    assert.equal(payload.elapsed_seconds, null);
});
