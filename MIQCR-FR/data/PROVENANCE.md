# Benchmark provenance and distribution status

The MIT code license does **not** license these benchmark inputs. Availability
for download is not treated as proof of permission to redistribute. Preserve
upstream notices. The converted files and shared feasible starts are hashed
in manifest.csv; metadata/paper_equality_cases.csv and qspp/cases.csv record
selection and source names. Historical paths in metadata are provenance, not
runtime dependencies.

| Family | Upstream source / retrieval route | Distribution evidence |
|---|---|---|
| QAP (61) | [QAPLIB](https://coral.ise.lehigh.edu/data-sets/qaplib/), obtained through [QAP R package mirror](https://github.com/mhahsler/qap/tree/a0e8cbacf270e7f8927afe79e65bf8fe89eb4412/inst/qaplib) | Mirror GPL-3.0 notice retained in licenses/; confirm its coverage of upstream benchmark files. |
| Cycle cover (102) | [de Meijer–Sotirov repository](https://github.com/frankdemeijer/SDPforQCCP), snapshot c90f32b25c934ba9725ca7543083e7b134d66ff2 | No explicit data redistribution grant located in the retained source. |
| Semi-assignment (30) | [The Optimization Firm, QSAP1](https://minlp.com/optimization-test-problems) | Distribution permission remains to be confirmed. |
| k-cluster (90) | [Lambert repository](https://github.com/amelie-lambert/k-cluster), snapshot 0c42bc7a0d20f3bad6eb838d86706fd1b500cbf2 | Repository EPL-2.0 notice retained in licenses/. |
| Partitioning (72) | [Trick's DIMACS/Stanford GraphBase collection](https://mat.tepper.cmu.edu/COLOR/instances.html); derived three/four-part cases | Confirm upstream graph redistribution terms; derived models do not remove this question. |
| QSPP (76) | Author-supplied Hu–Sotirov archive accompanying [On solving the quadratic shortest path problem](https://doi.org/10.1287/ijoc.2018.0861) | Explicit release authorization for the selected archive remains to be recorded. These are not the six unrelated QPLIB cases. |

Before public distribution of the full MAT bundle, resolve the items above
with the data owners or replace affected files with an upstream retrieval and
conversion procedure. The URLs provide upstream retrieval routes; they are
not yet an automated reconstruction of the hashed MAT inputs. No permission
has been inferred or requested on the authors' behalf by this packaging update.
