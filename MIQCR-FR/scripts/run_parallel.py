#!/usr/bin/env python3
"""Run independent benchmark tasks on an already allocated Linux/macOS machine."""

import argparse
import concurrent.futures
import os
from pathlib import Path
import shutil
import subprocess
import sys
import csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matlab", default="matlab")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--tasks", help="comma-separated task IDs; default all 1293")
    parser.add_argument("--output", type=Path, required=True, help="new directory")
    parser.add_argument("--total", type=float, default=3600)
    parser.add_argument("--p1", type=float, default=900)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    count = 3 * len(list(csv.DictReader((root / "data/manifest.csv").open())))
    try:
        tasks = (
            [int(t) for t in args.tasks.split(",")]
            if args.tasks
            else list(range(1, count + 1))
        )
    except ValueError:
        parser.error("--tasks must contain integers")
    if (
        args.jobs < 1
        or not 0 < args.total < float("inf")
        or not 0 < args.p1 < float("inf")
        or len(set(tasks)) != len(tasks)
        or not tasks
        or any(t < 1 or t > count for t in tasks)
    ):
        parser.error("invalid jobs, budgets, or task IDs")
    matlab = shutil.which(args.matlab)
    if matlab is None:
        parser.error("MATLAB executable not found: " + args.matlab)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(
        os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    folder = str(output).replace("'", "''")

    def run(task):
        command = f"setup; run_experiment({task},{args.total},{args.p1},'{folder}');"
        with (output / f"task_{task}.stdout.log").open("w") as log:
            result = subprocess.run(
                [matlab, "-singleCompThread", "-batch", command],
                cwd=root,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        return task, result.returncode

    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            task, code = future.result()
            print(f"Task {task}: exit {code}", flush=True)
            if code:
                failed.append(task)
    print("Failed tasks:", sorted(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
