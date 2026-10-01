#!/usr/bin/env python3
"""Collect a run without changing raw records. Uses only the Python standard library."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

METHODS = ("ADMM", "CB", "GUROBI")


def number(row, key):
    try:
        return float(row[key])
    except (ValueError, KeyError, TypeError):
        return float("nan")


def collect(
    manifest, runs, output, tasks=None, require_native=True, reference_incumbents=None
):
    """Check saved run records and report family means without altering inputs.

    A task is one instance-method pair. Reference incumbents, when supplied,
    check lower bounds only; they do not replace the objectives being reported.
    """
    # 1. Read the expected tasks and check each saved record.
    cases = list(csv.DictReader(manifest.open()))
    expected = {
        3 * k + j + 1: (case, method)
        for k, case in enumerate(cases)
        for j, method in enumerate(METHODS)
    }
    wanted = set(expected if tasks is None else tasks)
    if not wanted <= set(expected):
        raise ValueError("Task outside manifest")
    output.mkdir(parents=True, exist_ok=False)
    rows, flags, issues = {}, defaultdict(list), []
    for path in sorted(runs.glob("*.csv")):
        try:
            data = list(csv.DictReader(path.open()))
            if len(data) != 1:
                raise ValueError("expected one record")
            row = data[0]
            task = int(row["Task"])
            if task not in wanted or task in rows:
                raise ValueError("unexpected or duplicate task {}".format(task))
            rows[task] = row
            case, method = expected[task]
            if (row["Instance"], row["Method"], row["Family"]) != (
                case["Instance"],
                method,
                case["Family"],
            ):
                flags[task].append("manifest mismatch")
            if path.stem != case["Instance"] + "_" + method:
                flags[task].append("filename mismatch")
            if require_native and not path.with_suffix(".mat").is_file():
                flags[task].append("missing native MAT record")
            if row["Status"] == "ERROR":
                flags[task].append("runner ERROR")
            if row["Status"] not in (
                "OPTIMAL",
                "TIME_LIMIT",
                "SUBOPTIMAL",
                "INTERRUPTED",
                "NODE_LIMIT",
                "ITERATION_LIMIT",
                "SOLUTION_LIMIT",
                "WORK_LIMIT",
                "MEM_LIMIT",
                "TOTAL_LIMIT",
                "PHASE1_LIMIT",
                "ERROR",
            ):
                flags[task].append("unexpected solver status")
            if not math.isfinite(number(row, "Objective")):
                flags[task].append("missing verified incumbent")
            lower = number(row, "Bound")
            if math.isnan(lower) or lower == float("inf"):
                flags[task].append("invalid bound")
            for key in ("P1", "Total"):
                if not math.isfinite(number(row, key)) or number(row, key) < 0:
                    flags[task].append("invalid " + key)
            p2 = number(row, "P2")
            if not math.isfinite(p2) and row["Status"] not in ("PHASE1_LIMIT", "ERROR"):
                flags[task].append("missing P2")
            total = number(row, "P1") + (p2 if math.isfinite(p2) else 0)
            if abs(total - number(row, "Total")) > 1e-5 * max(1, total):
                flags[task].append("timing inconsistency")
        except (ValueError, KeyError) as err:
            issues.append("{}: {}".format(path.name, err))

    # 2. Check lower bounds against verified objectives, then apply the
    #    paper's optimality test to each run's own U and L.
    best = dict(reference_incumbents or {})
    if any(not math.isfinite(float(v)) for v in best.values()):
        raise ValueError("Nonfinite validation-only reference objective")
    for task, row in rows.items():
        if not flags[task]:
            key = row["Instance"]
            best[key] = min(best.get(key, float("inf")), number(row, "Objective"))
    for task, row in rows.items():
        upper = best.get(row["Instance"], number(row, "Objective"))
        if number(row, "Bound") > upper + 1e-5 + 1e-9 * abs(upper):
            flags[task].append("bound exceeds verified incumbent")
        upper, lower = number(row, "Objective"), number(row, "Bound")
        integer = str(row.get("IntegerObjective", "")).lower() in ("1", "true")
        threshold = 1 - 1e-4 if integer else max(1e-6, 1e-6 * abs(upper))
        row["ValidatedProved"] = int(
            not flags[task] and math.isfinite(lower) and upper - lower < threshold
        )
        row["ValidatedGap"] = (
            100 * max(0, upper - lower) / max(1, abs(upper))
            if not flags[task]
            else float("nan")
        )
        row["Validation"] = "; ".join(flags[task]) or "OK"
    if rows:
        with (output / "results.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(next(iter(rows.values()))))
            writer.writeheader()
            writer.writerows(rows[k] for k in sorted(rows))

    # 3. Summarize each family and method. Keep measured times even when a
    #    bound fails validation. Gap means use only finite, accepted gaps;
    #    record each denominator because some runs never reach Phase 2.
    summaries = []
    for family in sorted({expected[t][0]["Family"] for t in wanted}):
        for method in METHODS:
            ids = [
                t
                for t in wanted
                if expected[t][0]["Family"] == family and expected[t][1] == method
            ]
            if not ids:
                continue
            observed = [rows[t] for t in ids if t in rows]
            good = [rows[t] for t in ids if t in rows and not flags[t]]
            times = [
                number(r, "Total")
                for r in observed
                if math.isfinite(number(r, "Total"))
            ]
            gaps = [r["ValidatedGap"] for r in good if math.isfinite(r["ValidatedGap"])]
            p1 = [number(r, "P1") for r in observed if math.isfinite(number(r, "P1"))]
            p2 = [number(r, "P2") for r in observed if math.isfinite(number(r, "P2"))]
            summaries.append(
                dict(
                    Family=family,
                    Method=method,
                    Expected=len(ids),
                    Present=len(observed),
                    Flagged=len(observed) - len(good),
                    Certified=sum(r["ValidatedProved"] for r in good),
                    MeanTotal=sum(times) / len(times) if times else float("nan"),
                    TimeN=len(times),
                    MeanGap=sum(gaps) / len(gaps) if gaps else float("nan"),
                    QualityN=len(good),
                    GapN=len(gaps),
                    MeanP1=sum(p1) / len(p1) if p1 else float("nan"),
                    P1N=len(p1),
                    MeanP2=sum(p2) / len(p2) if p2 else float("nan"),
                    P2N=len(p2),
                )
            )
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    audit = dict(
        expected=len(wanted),
        present=len(rows),
        missing=sorted(wanted - set(rows)),
        validation_only_reference_count=len(reference_incumbents or {}),
        flags={str(k): v for k, v in flags.items() if v},
        issues=issues,
    )
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    lines = [
        "# Run summary",
        "",
        "Times retain available flagged records; quality excludes flagged records.",
        "Gap means include finite accepted gaps only; GapN and P2N are in summary.csv.",
        "Missing records are counted explicitly. See results.csv for U, L and native statuses.",
        "",
        "| Family | Method | Present/expected | Flagged | Certified | Mean total (s) | Mean gap (%) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in summaries:
        lines.append(
            "| {Family} | {Method} | {Present}/{Expected} | {Flagged} | "
            "{Certified} | {MeanTotal:.2f} | {MeanGap:.2f} |".format(**r)
        )
    (output / "report.md").write_text("\n".join(lines) + "\n")
    return audit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path)
    parser.add_argument("output", type=Path, help="new directory for generated tables")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data/manifest.csv",
    )
    parser.add_argument("--tasks", help="comma-separated task numbers; default: all")
    parser.add_argument(
        "--published-records",
        action="store_true",
        help="inspect curated CSV/log records without native MAT workspaces",
    )
    parser.add_argument(
        "--reference-incumbents",
        type=Path,
        help="JSON of verified objectives used only for bound validation",
    )
    args = parser.parse_args()
    audit = collect(
        args.manifest,
        args.runs,
        args.output,
        [int(t) for t in args.tasks.split(",")] if args.tasks else None,
        require_native=not args.published_records,
        reference_incumbents=(
            json.loads(args.reference_incumbents.read_text())
            if args.reference_incumbents
            else None
        ),
    )
    print(json.dumps(audit, indent=2))
    raise SystemExit(2 if audit["missing"] or audit["flags"] or audit["issues"] else 0)
