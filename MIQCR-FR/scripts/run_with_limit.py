#!/usr/bin/env python3
"""Run an external solver with a wall limit, terminating its process group."""

import os
import signal
import subprocess
import sys


def main():
    hard = sys.argv[1] == "--hard"
    if hard:
        del sys.argv[1]
    seconds = float(sys.argv[1])
    if seconds <= 0:
        raise ValueError("Time limit must be positive")
    # A separate process group lets the limit stop solver subprocesses as well.
    process = subprocess.Popen(sys.argv[2:], start_new_session=True)
    try:
        return process.wait(timeout=seconds)
    except subprocess.TimeoutExpired:
        if hard:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            return 124
        # Allow a short shutdown period before forcing termination.
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass  # The solver can finish between the wait and the signal.
        try:
            process.wait(timeout=30)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            return 137
        return 124


if __name__ == "__main__":
    sys.exit(main())
