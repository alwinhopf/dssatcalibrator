param(
    [string]$Config = "calibration_global_hemp\stage1_phenology_local_refine_around_n600_best.yaml",
    [string]$Parent = "results\global_hemp_calibration\global_anthesis_lhs5000_sharded_16targets_ukabaug9",
    [int]$FirstSeed = 45000,
    [int]$Shards = 10,
    [int]$ShardSize = 500
)

$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "python"

New-Item -ItemType Directory -Force -Path $Parent | Out-Null
$manifest = Join-Path $Parent "shard_manifest.csv"
if (-not (Test-Path $manifest)) {
    "shard,seed,n,include_start,outdir,status,started,ended" | Set-Content -Path $manifest
}

for ($i = 0; $i -lt $Shards; $i++) {
    $seed = $FirstSeed + $i
    $includeStart = ($i -eq 0)
    $n = if ($includeStart) { $ShardSize - 1 } else { $ShardSize }
    $outdir = Join-Path $Parent ("shard_{0:D3}" -f $i)
    $done = Join-Path $outdir "design.csv"
    $stdout = Join-Path $outdir "run.stdout.log"
    $stderr = Join-Path $outdir "run.stderr.log"

    New-Item -ItemType Directory -Force -Path $outdir | Out-Null

    if (Test-Path $done) {
        $now = Get-Date -Format s
        "$i,$seed,$n,$includeStart,$outdir,skipped,$now,$now" | Add-Content -Path $manifest
        continue
    }

    $started = Get-Date -Format s
    $args = @("run_calibration.py", $Config, "--n", "$n", "--seed", "$seed", "--outdir", $outdir)
    if (-not $includeStart) {
        $args += "--no-include-start"
    }

    $proc = Start-Process -FilePath "python" -ArgumentList $args -Wait -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr -WindowStyle Hidden
    $exitCode = $proc.ExitCode
    $ended = Get-Date -Format s
    $status = if ($exitCode -eq 0) { "complete" } else { "failed_$exitCode" }
    "$i,$seed,$n,$includeStart,$outdir,$status,$started,$ended" | Add-Content -Path $manifest
    if ($exitCode -ne 0) {
        throw "Shard $i failed with exit code $exitCode. See $stderr"
    }
}
