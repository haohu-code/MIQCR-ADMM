# MIQCR–FR: ADMM implementation

MATLAB implementation accompanying **Exact Convex Reformulations of Binary
Quadratic Programs from Facially Reduced Semidefinite Relaxations** by Hao Hu
and Mingming Xu. MIQCR–FR combines facial reduction with exact convex
reformulation for equality-constrained binary quadratic problems. This
implementation uses ADMM to solve the reduced relaxation and continuous-bound
optimization to construct the convex MIQP.

The paper and generated LaTeX tables use the name **MIQCR–FR**. Runner options
and CSV/MAT records retain the method identifier `ADMM` to match the archived
measurements. The repository and implementation directory are named `MIQCR-FR`.
The saved identifier and paper name refer to the same tested implementation.

## Reproduce the paper tables (no solver required)

From this directory, with Python 3.9 or newer:

```sh
python3 scripts/reproduce_paper.py outputs/paper_tables
```

Choose a new output directory. This checks the archived records, regenerates
all eight LaTeX tables, and checks them byte for byte against the historical
tables before replacing the displayed `ADMM` label with `MIQCR–FR`.
It also writes per-instance CSV results and validation/summary.csv,
including P1/P2 times and finite-gap denominators. The expected counts are
296 MIQCR–FR, 156 CB and 187 direct-Gurobi certifications. The known nug25 CB
bound failure remains explicitly flagged; its time and verified incumbent
are retained. See [the archive description](results/paper/README.md).
This reproduces reporting from saved measurements, not new solver timings.

## First solver run

Install MATLAB and a licensed Gurobi with its MATLAB interface. Set
`GUROBI_HOME` to its installation directory, or configure the MATLAB path
using Gurobi's setup instructions. Paper environment: MATLAB R2025b and
Gurobi 13.0.2 on Linux; a short macOS smoke check also passed with MATLAB
R2026a and Gurobi 11.0.2. Other versions may produce different results.

Start MATLAB from this directory with one computational thread:

```sh
matlab -singleCompThread
```

Then run:

```matlab
setup;
run_batch([1,3], 'outputs/quickstart', 30, 10);
```

These tasks run MIQCR–FR and direct Gurobi on chr12a (144 binary variables), with
30 seconds total and at most 10 seconds for the Phase-1 solver. Each creates
CSV, MAT and log files. The known optimum is 9552; a short run may stop before
proving it. Collect the results with:

```sh
python3 scripts/collect_results.py outputs/quickstart outputs/quickstart_report --tasks 1,3
```

Existing result files are never overwritten. Missing or invalid records make
the collector exit nonzero and are listed in audit.json. Gap means use finite,
validated bounds; all available measured times are retained.

## Optional MIQCR-CB comparison

MIQCR–FR and direct Gurobi do not require CB. To run CB, first follow
[the Linux build instructions](docs/BUILD_CB.md); Python 3 is also required
for its timeout wrapper. Set `MIQCR_EXTERNAL` before starting MATLAB to the
absolute directory containing the built `MIQCR-CB/` folder. Then:

```matlab
run_experiment(2,30,10,'outputs/cb_quickstart');
```

A short CB run may finish without a usable Phase-2 model (PHASE1_LIMIT).
The source and notices are in [external/MIQCR-CB](../external/MIQCR-CB/README.md).
The CB build is tested on x86-64 Linux; other platforms are not validated.

## Rerun the full benchmark

The [hashed manifest](data/manifest.csv) fixes 431 instances from six families.
For manifest row k, tasks `3*(k-1)+[1,2,3]` select MIQCR–FR, MIQCR-CB and
direct Gurobi (saved identifiers `ADMM`, `CB`, `GUROBI`).
All three use the same saved initial feasible solution.

```matlab
run_batch([], 'outputs/full_comparison', 3600, 900);
```

This is sequential and allows up to 1,293 solver-hours plus overhead.
For independent tasks on a workstation or allocated compute node, the Python
launcher starts separate single-thread MATLAB processes:

```sh
python3 scripts/run_parallel.py --matlab matlab --jobs 2 --tasks 1,3 --output outputs/parallel_example --total 30 --p1 10
# Full run, after building CB; choose jobs to fit memory and license capacity:
python3 scripts/run_parallel.py --matlab matlab --jobs 2 --output outputs/full_parallel
```

Each paper task had 64 GiB memory and one AMD EPYC 9654 core. Parallel jobs
must fit the available memory; run computation on compute nodes when using
a cluster. The launcher is for Linux/macOS; it is not a scheduler submission
script. For Slurm, set an array of 1–1293 with one CPU and 64 GiB per task,
export the one-thread settings below, and run from this directory:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
matlab -singleCompThread -batch "setup; run_experiment(${SLURM_ARRAY_TASK_ID},3600,900,'outputs/full_array');"
```

See [the protocol](docs/EXPERIMENTS.md) for budgets, solver settings and
validation. Changing versions or hardware is not a timing replication.
Collect a fresh run using:

```sh
python3 scripts/collect_results.py outputs/full_comparison outputs/full_report --reference-incumbents results/paper/reference_incumbents.json
```

Reference incumbents are used only to detect invalid bounds, never by the
optimizer. The fixed paper-table renderer applies only to the archived study;
new runs are summarized independently by collect_results.py.

## Reading the code

To follow the mathematics, read these MATLAB routines in order:

1. [equality_admm.m](src/equality_admm.m): form the reduced relaxation and
   perform the ADMM projections.
2. [continuous_bound_model.m](src/continuous_bound_model.m): extract the
   multipliers, choose the eigenvalue shift, and complete the quadratic matrix.
3. [phase2_model.m](src/phase2_model.m): add product variables and their
   McCormick constraints.
4. [run_experiment.m](scripts/run_experiment.m): put the two phases together,
   apply the time budget, and save the measurements.

The problem structure uses `p.Q`, `p.c`, and `p.constant` for the objective
`x'*p.Q*x + p.c'*x + p.constant`. Its equalities `p.G*x = p.h` correspond
to `Ax = b` in the paper; `p.A*x <= p.b` holds additional inequalities.
The mathematical routines retain symbols such as `N`, `W`, `mu`, and `Gamma`
so the calculations can be read alongside the paper.

For reporting, [collect_results.py](scripts/collect_results.py) checks new
runs and computes family averages.
[reproduce_paper.py](scripts/reproduce_paper.py) uses the saved paper records:
its `main()` lists the steps, followed through functions for validation,
summary tables, and detailed instance tables. All Python scripts use only
the standard library.

## Layout and scope

- `src/`: ADMM, reformulation, adapters and incumbent validation.
- `scripts/`: solver runners, collection and paper-table reproduction.
- `data/`: selected inputs, manifest and provenance.
- `results/paper/`: compact archived measurements and reporting checksums.
- `docs/`: numerical protocol and optional CB build instructions.

The MATLAB algorithms and settings are those of the final direct-MIQP
experiment `direct_miqp_20260925_r3`. Subsequent release changes improve
readability, reporting and packaging; they do not change the optimization methods. Phase 2 submits each
model directly to Gurobi once, without a separate continuous-relaxation solve.

Original code and documentation are MIT licensed; see [LICENSE](LICENSE) and [license scope](LICENSING.md).
Benchmark inputs and third-party code are excluded from that grant. Read
[data/PROVENANCE.md](data/PROVENANCE.md) for upstream sources and outstanding
data-distribution checks before publishing the complete data bundle.
