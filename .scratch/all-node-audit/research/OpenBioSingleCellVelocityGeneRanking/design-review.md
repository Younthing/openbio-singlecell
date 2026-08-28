# RNA Velocity Gene Ranking: module design review

## Decision

**Retain as a thin, honest view over a typed dynamics artifact and rename the display to “Recovered Dynamics Fit Ranking.”** It remains separate from recovery because it produces a bounded tabular reporting view; do not expand it into group-specific `rank_velocity_genes` behind a mode switch.

Inputs: `velocity_state` at `dynamics_recovered` or later stage with the same fit fingerprint; visible `top_n=50`; advanced `include_failed=False`; advanced `max_output_rows=100000`. Remove free-form `likelihood_column`; the artifact defines canonical `fit_likelihood`.

Outputs: canonical `table`, strict JSON `summary`, and `code`. Table columns are `rank`, `gene`, `fit_status`, `fit_likelihood`, `fit_r2`, `velocity_gene`, followed only by a fixed audited subset of available kinetic diagnostics. Sort finite successful fits by descending likelihood and original feature position; failed genes, if requested, follow deterministically and never receive a numeric rank.

`summary` reports eligible/fitted/failed/returned counts, likelihood and R2 quantiles, truncation, upstream model/fit fingerprint, references, versions and the non-inferential limitation. `code` accepts portable fitted AnnData and returns `(dataframe, summary_dict)` with identical schema/order.

Migration ignores no old custom column silently: only the canonical old default maps; a nondefault column blocks with guidance to use a generic table selector. Tests cover artifact stage/fingerprint, ties, NaN/Inf/failures, truncation/guards, exact schema/types, strict JSON, compilation and equivalence.
