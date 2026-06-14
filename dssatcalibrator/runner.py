"""Parallel execution of spawns.

The unit of work is one ``spawn_and_run`` call (an isolated DSSAT subprocess).
Spawns are independent and I/O/CPU-bound on the external binary, so a thread
pool launching subprocesses is sufficient and avoids Windows pickling issues.

Mirrors the gridded tutorial's warm-up-then-parallel idea: callers should run a
single spawn per experiment first (to let DSSAT resolve/cache weather+soil)
before fanning out — exposed here via ``warmup`` experiments run serially.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

from .spawn import SpawnResult, spawn_and_run


def resolve_cores(num_cores: int) -> int:
    if num_cores and num_cores > 0:
        return num_cores
    return max(1, (os.cpu_count() or 2) - 2)


def run_many(
    jobs: list[dict],
    *,
    n_workers: int,
    on_done: Callable[[SpawnResult], None] | None = None,
) -> list[SpawnResult]:
    """Run a list of spawn jobs in parallel.

    Each job dict is forwarded as kwargs to :func:`spawn_and_run`. Results are
    returned in job order. ``on_done`` (if given) is called as each completes,
    for progress reporting.
    """
    results: list[SpawnResult | None] = [None] * len(jobs)
    n_workers = max(1, n_workers)
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(spawn_and_run, job.pop("theta"), **job): i
                for i, job in enumerate(jobs)}
        for fut in as_completed(futs):
            i = futs[fut]
            try:
                res = fut.result()
            except Exception as e:  # noqa: BLE001 — a failed spawn must not kill the batch
                res = SpawnResult(status="error", run_dir=Path("."), theta={}, message=repr(e))
            results[i] = res
            if on_done:
                on_done(res)
    return results  # type: ignore[return-value]
