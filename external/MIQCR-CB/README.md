# MIQCR-CB dependency

`source/Smiqp-1.0/` contains the SMIQP 1.0 source, including its license and
bundled solver dependencies. The MIQCR–FR experiment runner invokes
`source/Smiqp-1.0/src/Alg/smiqp` for MIQCR-CB Phase 1 and reads the opt-in
`SMIQP_EXPORT_PHASE2` output. Build the executable from this source before
running MIQCR-CB tasks; compiled objects and executables are not distributed
in this release.

The source was downloaded from the authors' SMIQP distribution at
<https://cedric.cnam.fr/~lamberta/smiqp/Smiqp-1.0.zip> on 2026-07-20.
The original archive has SHA-256
`0994a8a50a9885bb40afc4cf7d4ccf0475e7d4418c1fd189e14a9e4293168262`.
The included source contains local portability and Phase-2 export changes,
followed by the stopping/checkpoint changes below. All bundled third-party
license notices are retained. The current Linux build is documented in
[BUILD_CB.md](../../MIQCR-FR/docs/BUILD_CB.md). Historical Python experiment
settings are not used by this release.

## Graceful binary Phase-1 stopping (September 24, 2026)

`src/SdpSolver/solver_sdp.c` has a local calling-routine modification:
`SMIQP_CB_WALL_LIMIT` sets a positive soft wall deadline checked after each
bundle descent step. On expiration the routine leaves the working set intact,
retrieves the current center, and performs the existing final SDP evaluation.
A CSDP return code other than 0 (success) or 3 (reduced accuracy) now aborts the executable rather than exporting
unsuccessful/interrupted oracle data. Rebuild after updating this source.
ConicBundle and CSDP library sources are unchanged. This is a soft limit;
a single step and final evaluation may overrun it. It does not by itself
certify dual feasibility in exact arithmetic.

Reduced-accuracy oracle returns are recorded in the log and saved export metadata; they are not exact dual-feasibility certificates.

### Oracle-failure fallback
The graceful-stopping revision now saves the last completed finite oracle
evaluation (status0 or3), including expanded multipliers and its matching SDP
solution. On a subsequent oracle failure it restores that snapshot and marks
ORACLE_FAILURE_WITH_SNAPSHOT; without a snapshot it fails explicitly. Failed
oracle data are never passed as successful subgradients. Normal model checks
and coefficient-omission corrections still apply; this is not an exact SDP
dual-feasibility certificate. This supersedes the abort-on-any-oracle-failure
behavior above and requires a fresh frozen comparison.

## Hard Phase-1 deadline with disk checkpoints
This revision supersedes the soft-time-limit protocol. After completed finite
oracle evaluations, a child process constructs the existing CB export and
atomically renames it into a separate checkpoint file. Checkpoints are attempted
at the first successful evaluation and at most once every30 seconds thereafter.
Child isolation leaves the ongoing solver state unchanged; all checkpoint work
counts toward P1. The P1 watchdog kills the process group at900 seconds without
a shutdown grace. The adapter uses the latest complete checkpoint, validates
it, and constructs P2 without another SDP solve. It records
TIME_LIMIT_WITH_CHECKPOINT. If no evaluation/export completed before the
deadline, PHASE1_LIMIT remains explicit. Model construction/validation still
consume the3600-second total budget. Never use incomplete .tmp files.
