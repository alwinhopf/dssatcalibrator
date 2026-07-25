param(
    [Parameter(Mandatory = $true)][string]$Config,
    [Parameter(Mandatory = $true)][string]$Parent,
    [int]$ChunkSize = 500,
    [string]$Python = "python",
    [switch]$NoCombine
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

New-Item -ItemType Directory -Force -Path $Parent | Out-Null
$design = Join-Path $Parent "full_design.csv"
$manifest = Join-Path $Parent "chunk_manifest.csv"

if (-not (Test-Path $design)) {
    & $Python run_calibration.py $Config --generate-design $design
    if ($LASTEXITCODE -ne 0) { throw "Failed to generate fixed design." }
}

$rowCount = @(Import-Csv $design).Count
if ($rowCount -lt 1) { throw "Fixed design is empty: $design" }
if (-not (Test-Path $manifest)) {
    "chunk,start,stop,outdir,status,started,ended" | Set-Content $manifest
}

$chunk = 0
for ($start = 0; $start -lt $rowCount; $start += $ChunkSize) {
    $stop = [Math]::Min($start + $ChunkSize, $rowCount)
    $outdir = Join-Path $Parent ("chunk_{0:D3}" -f $chunk)
    $resultDesign = Join-Path $outdir "design.csv"
    if (Test-Path $resultDesign) {
        $chunk += 1
        continue
    }

    New-Item -ItemType Directory -Force -Path $outdir | Out-Null
    $stdout = Join-Path $outdir "run.stdout.log"
    $stderr = Join-Path $outdir "run.stderr.log"
    $started = (Get-Date).ToString("s")
    $args = @(
        "run_calibration.py", $Config,
        "--design-csv", $design,
        "--design-start", $start,
        "--design-stop", $stop,
        "--outdir", $outdir,
        "--design-only"
    )
    $proc = Start-Process -FilePath $Python -ArgumentList $args -Wait -PassThru `
        -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $ended = (Get-Date).ToString("s")
    $status = if ($proc.ExitCode -eq 0 -and (Test-Path $resultDesign)) { "complete" } else { "failed" }
    if ($status -eq "complete") {
        $evaluated = @(Import-Csv $resultDesign)
        $expectedRows = $stop - $start
        if ($evaluated.Count -ne $expectedRows) {
            $status = "failed_row_count"
        } else {
            $nObs = @($evaluated | ForEach-Object { [int]$_.n_obs } | Sort-Object -Unique)
            if ($nObs.Count -ne 1 -or $nObs[0] -lt 1) {
                $status = "failed_observation_count"
            }
        }
    }
    "$chunk,$start,$stop,$outdir,$status,$started,$ended" | Add-Content $manifest
    if ($status -ne "complete") {
        throw "Chunk $chunk failed validation ($status). See $stdout and $stderr."
    }
    $chunk += 1
}

if (-not $NoCombine) {
    $chunkDirs = Get-ChildItem $Parent -Directory -Filter "chunk_*" |
        Where-Object { Test-Path (Join-Path $_.FullName "design.csv") } |
        Sort-Object Name |
        ForEach-Object { $_.FullName }
    $combined = "${Parent}_combined"
    & $Python run_calibration.py $Config --combine $chunkDirs --outdir $combined
    if ($LASTEXITCODE -ne 0) { throw "Combining fixed-design chunks failed." }
}
