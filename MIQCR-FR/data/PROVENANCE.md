# Benchmark sources and licenses

This directory contains the 431 selected inputs used in the paper. The MIT
license for our implementation does **not** apply to third-party benchmark
data. The source references and available license information are listed below;
the absence of a dataset-specific license is not a grant of additional rights
from this repository.

## Sources

| Family | Instances | Source |
|---|---:|---|
| QAP | 61 | [QAPLIB](https://coral.ise.lehigh.edu/data-sets/qaplib/), originally retrieved from the [QAP R package mirror](https://github.com/mhahsler/qap/tree/a0e8cbacf270e7f8927afe79e65bf8fe89eb4412/inst/qaplib). |
| Cycle cover | 102 | [F. de Meijer and R. Sotirov's collection](https://github.com/frankdemeijer/SDPforQCCP/tree/c90f32b25c934ba9725ca7543083e7b134d66ff2). |
| Semi-assignment | 30 | QSAP1 from [The Optimization Firm's benchmark collection](https://minlp.com/optimization-test-problems), accompanying Nohra, Raghunathan, and Sahinidis, *Spectral relaxations and branching strategies for global optimization of mixed-integer quadratic programs* (2021). |
| k-cluster | 90 | [Amélie Lambert's collection](https://github.com/amelie-lambert/k-cluster/tree/0c42bc7a0d20f3bad6eb838d86706fd1b500cbf2). |
| Partitioning | 72 | Three- and four-part instances constructed from 36 graphs in [Michael Trick's DIMACS/Stanford GraphBase collection](https://mat.tepper.cmu.edu/COLOR/instances.html). |
| QSPP | 76 | Archive supplied by Hao Hu, accompanying Hu and Sotirov, [*On solving the quadratic shortest path problem*](https://doi.org/10.1287/ijoc.2018.0861). |

## License information

**QAP.** The [University of Edinburgh QAPLIB deposit](https://doi.org/10.7488/ds/3428)
is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/).
Its authors are R. E. Burkard, E. Çela, S. E. Karisch, F. Rendl, M. Anjos,
and P. Hahn (2022). All 61 selected raw instance files match this deposit's
data after whitespace normalization. We converted the selected instances to
the MATLAB formulation used by this solver. The R mirror's separate
[GPL-3.0 notice](licenses/qap-mirror-GPL-3.0.txt) is also retained.

**k-cluster.** The upstream repository is distributed under
[EPL-2.0](licenses/kcluster-EPL-2.0.txt). The original data are available at
the linked source snapshot; this release converts the selected instances to
the solver's MATLAB format. The upstream license and notices continue to
apply to this material.

**Partitioning.** Knuth identifies the programs and data of
[Stanford GraphBase as public domain](https://www-cs-faculty.stanford.edu/~knuth/sgb.html).
The selected GraphBase-derived graphs are `anna`, `david`, `huck`, `jean`,
`games120`, the five `miles` graphs, and the thirteen `queen` graphs. These
23 graphs yield 46 partitioning instances. The remaining graphs are Gary
Lewandowski's eight `mulsol`/`zeroin` register-allocation graphs and Michael
Trick's five `myciel` graphs. Their source headers give attribution but do
not state a dataset-specific license.

**Cycle cover.** The linked repository provides the benchmark archives and
format documentation. No dataset-specific license is stated in the checked
repository or retained archives.

**Semi-assignment.** The QSAP1 archive does not state a dataset-specific
license. The provider's [Terms of Service](https://minlp.com/terms-service)
contain restrictions on redistribution of its materials without authorization;
this repository does not grant a separate license for these inputs.

**QSPP.** The supplied Hu–Sotirov archive identifies the instance structure
and source study but contains no separate data-license statement. These
inputs are not licensed under this repository's MIT code license.

Source and license information was checked on September 30–October 1, 2026.

## Conversion and instance identification

The release stores the selected problems in a common MATLAB representation
and includes shared feasible starting solutions. Partitioning models are
derived from the source graphs; the other families are converted from the
linked benchmark collections. These conversions do not imply endorsement
by the original data authors.

[manifest.csv](manifest.csv) records the input paths and SHA-256 hashes.
[metadata/paper_equality_cases.csv](metadata/paper_equality_cases.csv) and
[qspp/cases.csv](qspp/cases.csv) record instance selection and source names.
Historical paths in the metadata document provenance and are not runtime
dependencies. The upstream links identify the original collections; they
are not an automated reconstruction procedure for the supplied MAT files.
