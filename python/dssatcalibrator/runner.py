"""Parallel execution of spawns.

The unit of work is one :func:`spawn.spawn_and_run` call — an isolated DSSAT
subprocess. Because the heavy lifting happens in that external binary, a *thread*
pool that launches subprocesses gives true multi-core parallelism (Python releases
the GIL while waiting on the process) while sidestepping the pickling/spawn cost
that a process pool incurs on Windows. So: one OS process per DSSAT run, scheduled
across all cores; the Python threads just wait.

Warm-up then parallel
---------------------
``run_many(..., warmup=k)`` runs the first ``k`` jobs **serially** before fanning
the rest across the pool. This mirrors the gridded tutorial's pattern: materialise
any shared, lazily-created inputs (e.g. a site's downloaded+cached weather/soil)
exactly once, so the parallel workers only ever read them and never race to create
them. In the current prototype DSSAT resolves weather/soil centrally from the
install, so there is nothing to warm up and ``warmup`` defaults to 0 — but the hook
is here for when on-demand ``dssatutils`` downloading is wired in.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable

from .spawn import SpawnResult, spawn_and_run


def resolve_cores(num_cores: int) -> int:
    """How many workers to use. ``num_cores > 0`` is taken as-is; ``0`` (or less)
    means "leave 2 logical cores free for the rest of the machine"."""
    if num_cores and num_cores > 0:
        return num_cores
    return max(1, (os.cpu_count() or 2) - 2)


def run_many(
    jobs: list[dict],
    *,
    n_workers: int,
    on_done: Callable[[SpawnResult], None] | None = None,
    warmup: int = 0,
) -> list[SpawnResult]:
    """Run a list of spawn jobs and return their results in job order.

    Each job dict is forwarded as kwargs to :func:`spawn.spawn_and_run` (it must
    contain ``theta`` plus the spawn keyword arguments). A spawn that raises is
    captured as an ``error`` result so one bad run never kills the batch.
    ``on_done`` (if given) is called as each job finishes, for progress reporting.
    ``warmup`` runs the first N jobs serially before parallelising the remainder.
    """
    results: list[SpawnResult | None] = [None] * len(jobs)
    n_workers = max(1, n_workers)

    def _run(i: int) -> SpawnResult:
        job = jobs[i]
        theta = job.get("theta")
        kwargs = {k: v for k, v in job.items() if k != "theta"}
        try:
            return spawn_and_run(theta, **kwargs)
        except Exception as e:  # noqa: BLE001 — a failed spawn must not kill the batch
            return SpawnResult(status="error", run_dir=Path("."), theta=theta or {}, message=repr(e))

    warm = min(max(warmup, 0), len(jobs))
    for i in range(warm):                       # serial warm-up pass
        results[i] = _run(i)
        if on_done:
            on_done(results[i])

    with ThreadPoolExecutor(max_workers=n_workers) as ex:   # parallel remainder
        futs = {ex.submit(_run, i): i for i in range(warm, len(jobs))}
        for fut in as_completed(futs):
            i = futs[fut]
            results[i] = fut.result()
            if on_done:
                on_done(results[i])
    return results  # type: ignore[return-value]
