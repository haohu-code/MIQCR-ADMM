#!/usr/bin/env python3
"""Regenerate the paper tables from saved measurements, without running solvers.

Usage (from MIQCR-FR):
    python3 scripts/reproduce_paper.py outputs/paper_tables

The script validates the saved records, excludes the known invalid CB bound,
and writes the summary and per-instance tables. It checks those tables against
the archived paper versions before updating the display name from ADMM to
MIQCR--FR. Saved method identifiers are unchanged. This script describes one
fixed experiment; use collect_results.py to summarize new solver runs.
"""

import argparse
import csv
import hashlib
import json
import math
import statistics
import tempfile
from pathlib import Path

from collect_results import collect

# Keep the order used in the manuscript.
FAMILIES = {
    "qap": "QAP",
    "cycle_cover": "Cycle cover",
    "semi_assignment": "Semi-assignment",
    "kcluster": "$k$-cluster",
    "partition": "Partitioning",
    "qspp": "QSPP",
}
METHODS = ["ADMM", "CB", "GUROBI"]
STATUS_LABELS = {
    "OPTIMAL": "O",
    "TIME_LIMIT": "T",
    "PHASE1_LIMIT": "P1",
    "VALIDATION_FAILURE": "V",
}


def format_number(value, decimals=2):
    """Format finite values and missing/infinite bounds for LaTeX."""
    value = float(value)
    if math.isnan(value):
        return "--"
    if math.isinf(value):
        if value < 0:
            return r"$-\infty$"
        return r"$\infty$"
    return f"{value:,.{decimals}f}"


def finite_mean(rows, column):
    """Average the available finite measurements, as in the paper."""
    values = [float(row[column]) for row in rows]
    values = [value for value in values if math.isfinite(value)]
    if not values:
        return math.nan
    return statistics.mean(values)


def latex_row(entries):
    return " & ".join(str(entry) for entry in entries) + r" \\" + "\n"


def check_file_hashes(folder, expected_hashes):
    """Check that the supplied files agree with the saved reference copies."""
    for name, expected in expected_hashes.items():
        actual = hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit(f"Checksum mismatch: {folder / name}")


def validate_archive(root, archive, output):
    """Apply the same record-level bound checks used when collecting new runs."""
    expected_hashes = json.loads((archive / "checksums.json").read_text())
    check_file_hashes(archive, expected_hashes)
    references = json.loads((archive / "reference_incumbents.json").read_text())

    # The archive stores one CSV; the collector reads one CSV per task.
    # Split it temporarily so both paths use exactly the same validation rules.
    with tempfile.TemporaryDirectory() as temporary_folder:
        runs = Path(temporary_folder) / "runs"
        runs.mkdir()
        with (archive / "records.csv").open() as stream:
            for row in csv.DictReader(stream):
                filename = row["Instance"] + "_" + row["Method"] + ".csv"
                with (runs / filename).open("w", newline="") as task_file:
                    writer = csv.DictWriter(task_file, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)

        validation = collect(
            root / "data/manifest.csv",
            runs,
            output / "validation",
            require_native=False,
            reference_incumbents=references,
        )

    # This failure is retained, not reclassified as a valid bound.
    expected_flags = {"131": ["bound exceeds verified incumbent"]}
    if (
        validation["missing"]
        or validation["issues"]
        or validation["flags"] != expected_flags
    ):
        raise SystemExit("Unexpected validation outcome; inspect validation/audit.json")

    source = output / "validation/results.csv"
    with source.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1293
    assert len({row["Task"] for row in rows}) == 1293
    return rows, references, source


def prepare_publication_rows(rows, references):
    """Keep original measurements, but remove rejected bounds from reported results."""
    for row in rows:
        row["SourceRun"] = "direct_miqp_20260925_r3"
        for column in ("Status", "Bound", "Gap", "Proved"):
            row["Source" + column] = row[column]

        # These runs retain the common initial feasible solution.
        initial_only = row["Status"] == "PHASE1_LIMIT"
        rejected_cb_case = row["Instance"] == "nug25" and row["Method"] == "CB"
        row["InitialOnly"] = str(int(initial_only or rejected_cb_case))
        row["PublicationValidation"] = row["Validation"]
        row["Gap"] = row["ValidatedGap"]
        row["Proved"] = row["ValidatedProved"]
        if row["Validation"] != "OK":
            row["Bound"] = "nan"
            row["Gap"] = "nan"
            row["Proved"] = "0"
            row["Status"] = "VALIDATION_FAILURE"

    excluded = [
        (row["Instance"], row["Method"]) for row in rows if row["Validation"] != "OK"
    ]
    assert excluded == [("nug25", "CB")]

    # U is the verified feasible objective and L is the adjusted lower bound.
    # Historical reference objectives are used only to check L, never as starts.
    for row in rows:
        if row["Validation"] != "OK":
            continue
        U = float(row["Objective"])
        L = float(row["Bound"])
        observed_best = min(
            float(other["Objective"])
            for other in rows
            if other["Instance"] == row["Instance"]
        )
        reference_best = min(references[row["Instance"]], observed_best)
        assert L <= reference_best + 1e-5 + 1e-9 * abs(reference_best)
        if row["IntegerObjective"] == "1":
            tolerance = 1 - 1e-4
        else:
            tolerance = max(1e-6, 1e-6 * abs(U))
        certified = math.isfinite(L) and U - L < tolerance
        assert int(row["Proved"]) == int(certified)


def write_summary(rows, output):
    """Write family means and certification counts for the main paper."""
    text = (
        r"\begin{table}[ht]\centering\small\setlength{\tabcolsep}{3pt}" + "\n"
        r"\caption{Results for all 431 instances. Times are in seconds; reported "
        r"times and gaps are arithmetic means.}\label{tab:final-summary}" + "\n"
        r"\begin{tabular}{llrrrrr}\toprule" + "\n"
    )
    text += latex_row(
        ["Family", "Method", "Certified", "P1", "P2", "Total", r"Gap (\%)"]
    )
    text += r"\midrule" + "\n"
    for family, family_name in FAMILIES.items():
        for method in METHODS:
            group = [
                row
                for row in rows
                if row["Family"] == family and row["Method"] == method
            ]
            certified = sum(int(row["Proved"]) for row in group)
            text += latex_row(
                [
                    family_name if method == "ADMM" else "",
                    method.replace("GUROBI", "Gurobi"),
                    f"{certified}/{len(group)}",
                    format_number(finite_mean(group, "P1")),
                    format_number(finite_mean(group, "P2")),
                    format_number(finite_mean(group, "Total")),
                    format_number(finite_mean(group, "Gap")),
                ]
            )
        text += r"\addlinespace" + "\n"
    text += r"""\bottomrule\end{tabular}
\par\smallskip
\begin{minipage}{\linewidth}\footnotesize\raggedright
P2 and gap means include all instances except for CB: four cycle-cover and
seven QSPP instances have no Phase~2 model or finite gap; one QAP gap is
excluded by validation. Total-time means include all runs.
\end{minipage}
\end{table}
"""
    (output / "summary.tex").write_text(text)


def write_instance_tables(rows, output):
    """Write the six detailed tables, retaining native outcomes and star/dagger marks."""
    for family, family_name in FAMILIES.items():
        text = (
            r"\begingroup\scriptsize\setlength{\tabcolsep}{2pt}" + "\n"
            r"\begin{longtable}{ll l rrrrrr l}"
            + "\n"
            + rf"\caption{{{family_name}: all instances; times in seconds.}}"
            + rf"\label{{tab:instances-{family}}}"
            + r"\\"
            + "\n"
        )
        header = r"\toprule" + "\n"
        header += latex_row(
            [
                "Instance",
                "$n$",
                "Method",
                "P1",
                "P2",
                "Total",
                "$U$",
                "$L$",
                r"Gap (\%)",
                "Status",
            ]
        )
        header += r"\midrule" + "\n"
        text += header + r"\endfirsthead" + "\n"
        text += header + r"\endhead" + "\n"
        text += r"\midrule\multicolumn{10}{r}{Continued on next page}\\\endfoot" + "\n"
        text += r"\bottomrule\endlastfoot" + "\n"

        # ADMM rows select each instance once, in the original manifest order.
        instances = [
            row["Instance"]
            for row in rows
            if row["Family"] == family and row["Method"] == "ADMM"
        ]
        for instance in instances:
            group = [row for row in rows if row["Instance"] == instance]
            best = min(float(row["Objective"]) for row in group)
            for method in METHODS:
                row = next(row for row in group if row["Method"] == method)
                U = float(row["Objective"])
                decimals = 4 if family == "semi_assignment" else 0
                objective_text = format_number(U, decimals)
                if abs(U - best) <= 1e-6 + 1e-9 * abs(best):
                    objective_text = r"\textbf{" + objective_text + "}"
                if row["InitialOnly"] == "1":
                    objective_text += r"$^{\dagger}$"
                status = STATUS_LABELS[row["Status"]]
                if row["Proved"] == "1":
                    status += r"$^{*}$"
                line = latex_row(
                    [
                        instance.replace("_", r"\_") if method == "ADMM" else "",
                        row["N"] if method == "ADMM" else "",
                        method.replace("GUROBI", "Gurobi"),
                        format_number(row["P1"]),
                        format_number(row["P2"]),
                        format_number(row["Total"]),
                        objective_text,
                        format_number(row["Bound"]),
                        format_number(row["Gap"]),
                        status,
                    ]
                )
                # Keep the three methods for an instance on the same PDF page.
                if method != "GUROBI":
                    line = line.rstrip() + "*\n"
                text += line
            text += r"\addlinespace[2pt]" + "\n"
        text += r"\end{longtable}\endgroup" + "\n"
        (output / f"instances_{family}.tex").write_text(text)


def write_incumbent_table(rows, output):
    """Count matches to the best observed objective, allowing numerical ties."""
    text = (
        r"\begin{table}[ht]\centering\small\caption{Matches to the best observed "
        r"feasible objective, with ties allowed.}\label{tab:final-incumbents}"
        r"\begin{tabular}{lrrr}\toprule" + "\n"
    )
    text += latex_row(["Family", "ADMM", "CB", "Gurobi"]) + r"\midrule" + "\n"
    for family, family_name in FAMILIES.items():
        counts = []
        for method in METHODS:
            group = [
                row
                for row in rows
                if row["Family"] == family and row["Method"] == method
            ]
            matches = 0
            for row in group:
                best = min(
                    float(other["Objective"])
                    for other in rows
                    if other["Instance"] == row["Instance"]
                )
                if abs(float(row["Objective"]) - best) <= 1e-6 + 1e-9 * abs(best):
                    matches += 1
            counts.append(f"{matches}/{len(group)}")
        text += latex_row([family_name] + counts)
    text += r"\bottomrule\end{tabular}\end{table}" + "\n"
    (output / "incumbents.tex").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", type=Path, help="new directory for regenerated tables"
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    archive = root / "results/paper"
    output = args.output
    output.mkdir(parents=True, exist_ok=False)

    # 1. Validate saved bounds and prepare the rows used in the paper.
    rows, references, source = validate_archive(root, archive, output)
    prepare_publication_rows(rows, references)
    with (output / "results.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    certified = {
        method: sum(int(row["Proved"]) for row in rows if row["Method"] == method)
        for method in METHODS
    }
    assert certified == {"ADMM": 296, "CB": 156, "GUROBI": 187}
    audit = {
        "records": len(rows),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source": "validation/results.csv",
        "certified": certified,
        "excluded_bounds": [
            {
                "task": row["Task"],
                "instance": row["Instance"],
                "method": row["Method"],
                "reason": row["Validation"],
            }
            for row in rows
            if row["Validation"] != "OK"
        ],
        "rule": (
            "Use all fresh production runs only. Retain all measured costs and "
            "verified feasible objectives; exclude nug25 CB bound, gap and "
            "certification. Do not replace production with regression or diagnostic "
            "results."
        ),
        "quality_denominators": (
            "Finite accepted gaps only. All verified incumbents are eligible for "
            "best-observed matches. PHASE1_LIMIT and nug25 CB retain only their "
            "structural starts."
        ),
        "qualification": (
            "Floating-point validation, not interval certificates. The original "
            "collector correctly failed on task 131; this reporting exclusion does "
            "not clear that audit."
        ),
    }
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")

    # 2. Write the two summary tables and the six per-family tables.
    write_summary(rows, output)
    write_instance_tables(rows, output)
    write_incumbent_table(rows, output)

    # 3. Check that the generated LaTeX agrees with the archived paper tables.
    expected = json.loads((archive / "expected_tables_sha256.json").read_text())
    check_file_hashes(output, expected)

    # 4. Use the current paper name after checking every historical table byte.
    # ADMM remains the method identifier in the archived CSV and validation files.
    # Only the displayed method cell changes; all measurements stay unchanged.
    for filename in expected:
        table = output / filename
        original_text = table.read_text()
        table.write_text(original_text.replace(" & ADMM & ", " & MIQCR--FR & "))

    print(json.dumps(audit, indent=2))
    print("All 8 historical tables verified byte for byte before relabeling.")
    print("Output tables use MIQCR--FR; archived ADMM identifiers are unchanged.")
    print("Known nug25 CB bound rejection retained; see validation/audit.json.")


if __name__ == "__main__":
    main()
