# Selected comparison inputs

Inputs are organized by problem class:

| Folder | Instances |
|---|---:|
| `qap/` | 61 |
| `cycle_cover/` | 102 |
| `semi_assignment/` | 30 |
| `kcluster/` | 90 |
| `partition/` | 72 |
| `qspp/` | 76 |
| **Total** | **431** |

`manifest.csv` records each instance's study, family, name, current path
(relative to the release root), and SHA-256 hash. QSPP network groups are retained in the Subfamily column. Input contents and shared
starting solutions are unchanged.

`metadata/paper_equality_cases.csv` preserves the original 355-case ordering
and source provenance; `metadata/paper_equality_validation.csv` retains its
validation records. `qspp/cases.csv` preserves the 76-case QSPP ordering.
These tables retain historical ordering and provenance. The unified runner uses
the row order of manifest.csv for all six classes.

These are the selected inputs for the completed 431-instance paper experiment,
not the full upstream archives. Historical study names identify provenance
only. See [PROVENANCE.md](PROVENANCE.md) for sources and distribution status.
The manifest paths are relative to MIQCR-ADMM; no private checkout is needed
to load the supplied MAT files.
