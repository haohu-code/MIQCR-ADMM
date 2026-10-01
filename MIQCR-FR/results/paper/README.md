# Archived paper measurements

Source study: direct_miqp_20260925_r3, 431 instances and 1,293 fresh method runs.
`records.csv` preserves the collected records, including raw bounds, native
statuses and original validation fields. It is not a new optimization run.
`reference_incumbents.json` contains historical verified objectives used only
for validation; its provenance is recorded separately. `checksums.json`
protects archive inputs; expected_tables_sha256.json identifies the eight
historical LaTeX tables. Their `ADMM` method label is now displayed as
`MIQCR–FR` in the paper. The records and historical checksums remain unchanged.

From the MIQCR-FR directory run:

    python3 scripts/reproduce_paper.py outputs/paper_tables

The renderer verifies the historical table bytes first, then changes only the
displayed method label to `MIQCR–FR`. CSV and validation outputs retain `ADMM`.

The script recomputes validation and requires precisely the recorded task 131
(nug25, CB) bound contradiction. It then excludes that bound, gap and
certification from publication summaries, retaining measured costs and the
verified incumbent. This expected rejection is not a successful solver
certificate. Eleven CB runs have no Phase-2 model; finite-gap and P2 averages
use the smaller denominators. Certification counts are 296/156/187 for
MIQCR–FR/MIQCR-CB/direct Gurobi.

This compact archive supports table regeneration and record-level checks.
Native MAT workspaces and solver logs are not included; consequently it does
not independently revalidate solution vectors or rerun solvers. For that use
the full benchmark instructions. Timing differences across machines are expected.
