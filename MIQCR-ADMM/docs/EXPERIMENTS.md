# Unified six-family comparison: direct MIQP protocol

## Inputs and methods

Use all 431 instances in data/manifest.csv, in its hashed order: 61 QAP,
102 cycle-cover, 30 semi-assignment, 90 k-cluster, 72 partitioning and76 QSPP.
Task3*(k-1)+[1,2,3] selects MIQCR–FR, MIQCR-CB, and direct Gurobi for
manifest rowk. The saved method identifiers remain `ADMM`, `CB`, `GUROBI`.
MIQCR–FR is the paper name for the ADMM implementation in this protocol.
All methods
start from the same structural feasible point. No previous incumbents are
used as starts, cutoffs, objective bounds or stopping criteria.

MIQCR–FR uses ADMM to solve the equality-face DNN relaxation, then constructs
its convex MIQP by continuous-bound optimization. QAP retains the box/gangster treatment;
other families use the arrow/nonnegative block. No product pruning or extra
relaxation strengthening is introduced. The binary-exact Hessian guard is
recorded. CB uses SMIQP's existing reformulation with coefficient-omission
accounting. Direct Gurobi receives the original binary quadratic model.

## One Phase-2 solve

Each obtained model is submitted directly to Gurobi once. There is no separate
continuous-relaxation solve, root-bound diagnostic, preliminary optimization,
or root-QP time allowance. Gurobi manages its own root node and search.

The published comparison used Gurobi13.0.2 for all three methods, MATLAB
R2025b, one AMD EPYC9654 core and64GiB per task. MATLAB, BLAS and Gurobi use
one thread. Seed=0; MIPGap=MIPGapAbs=1e-6. Presolve, cuts, heuristics and the root algorithm use automatic defaults.
For both convex reformulation methods, PreQLinearize=0 preserves their
quadratic relaxation; direct Gurobi retains automatic prelinearization. NonConvex=0 checks the convex reformulations; NonConvex=2
allows the original direct formulation. There is no case-specific tuning.
Local users must record their
actual solver version; changing it is not an exact replication of the paper.

## Budget and accounting

Total algorithmic budget:3600seconds. Phase-1 solver cap:900seconds, clipped
to the total budget. ADMM uses tolerance1e-3, at most10000 iterations and
eight penalty changes. Its face construction is included in its solver time.

CB input writing, the external solve and export validation are included in
Phase1. After completed finite oracle evaluations, a child process constructs
an export in isolation; the caller waits for it, so checkpoint work is charged
to the same single-core P1 budget. Atomic rename protects complete checkpoints.
The external process group is killed at the remaining P1 deadline with no
shutdown grace. Only the latest fully completed export/checkpoint is read;
no SDP is solved after the deadline. No completed checkpoint means PHASE1_LIMIT.
The ConicBundle and CSDP libraries are unchanged; the calling/export interface
contains the previously introduced checkpoint and oracle-fallback code.
An accepted CSDP status is not an exact-arithmetic SDP dual certificate.

P1=SDPSeconds+Post. Post includes reformulation and model-specific start checks.
P2 includes the single integer solve and returned-solution validation. Its
solver time allowance is3600-P1, less measured P2 setup overhead. Total=P1+P2
when P2 exists. Shared input loading and common original-model/start preparation,
MATLAB startup, and output serialization are excluded for all methods.
Noninterruptible numerical work/termination may slightly overrun; retain actual
measured times. Cluster watchdogs are safety limits, not extra algorithm time.

## Numerical validation

For returned points, check integrality and original feasibility; QSPP also
checks the original network and a single simple source-to-sink path.
Reconstruct every retained product as x_i*x_j after rounding the binary
variables, check all model constraints and bounds, and compare the original
objective with the reconstructed model objective. Allow only the recorded CB
coefficient-omission bound plus1e-5+1e-9*abs(original objective).
Independently evaluate the model at the raw solver vector and compare with
its reported objval. Record the raw-to-reconstructed objective change;
small variable deviations are not themselves a failed binary identity.
No tolerance is enlarged to hide the earlier discrepancies.

Final L is the integer solver's objbound minus the CB omission allowance.
It must not exceed the run's verified feasible U. The collector also checks
against other methods' verified incumbents and, for this local study, a frozen
list of historical verified feasible objectives used ONLY for validation.
This list is never passed to optimization. A contradicted bound remains an
explicit validation failure even if the solver reports OPTIMAL.

Numerical certification requires U-L<1-1e-4 for integer-valued objectives,
or U-L<max(1e-6,1e-6*abs(U)) otherwise. These are floating-point tests, not
interval proofs. Preserve pipeline Status and native SolverStatus separately.
Save model, solver outputs, effective parameters, validation measurements,
input hashes and version metadata. Missing/error records cannot count as solved.

## Reproduction and history

run_experiment and run_batch execute the protocol. collect_results.py validates
and summarizes records; its output directory must be new. A failed collector
still writes its audit and returns nonzero. Local cluster launchers and additional
regression tests are outside the release. The original cluster concurrency cap was 800; it is a scheduling
limit, not an algorithm parameter or a workstation default.

The completed paper experiment is direct_miqp_20260925_r3. Its compact
records and validation references are included in results/paper/. Run
`python3 scripts/reproduce_paper.py outputs/paper_tables` to reproduce its tables.
New optimization runs must use fresh directories and are not substituted into
that archive. The public collector averages finite accepted gaps and records
the denominators; measured times retain flagged runs.

The cluster build uses ConicBundle's standard OPTI mode (`-O3 -DNDEBUG`)
and builds CSDP with `-O3`, without fast-math assumptions. These are build
options; third-party optimization sources are unchanged.
