param(
    [Parameter(Mandatory = $true)][string]$Config,
    [Parameter(Mandatory = $true)][string]$Parent,
    [int]$ChunkSize = 500,
    [string]$Python = "python",
    [int]$ExpectedNObs = 0,
    [string]$DeadlineUtc = "",
    [int]$MinRemainingMinutes = 35,
    [switch]$NoCombine
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

New-Item -ItemType Directory -Force -Path $Parent | Out-Null
$design = Join-Path $Parent "full_design.csv"
$manifest = Join-Path $Parent "chunk_manifest.csv"
$contextPath = Join-Path $Parent "run_context.json"

if (-not (Test-Path $design)) {
    & $Python run_calibration.py $Config --generate-design $design
    if ($LASTEXITCODE -ne 0) { throw "Failed to generate fixed design." }
}

$rowCount = @(Import-Csv $design).Count
if ($rowCount -lt 1) { throw "Fixed design is empty: $design" }
$context = [ordered]@{
    config_path = (Resolve-Path $Config).Path
    config_sha256 = (Get-FileHash -Algorithm SHA256 $Config).Hash
    design_sha256 = (Get-FileHash -Algorithm SHA256 $design).Hash
    row_count = $rowCount
    chunk_size = $ChunkSize
    expected_n_obs = $ExpectedNObs
}
if (Test-Path $contextPath) {
    $saved = Get-Content $contextPath -Raw | ConvertFrom-Json
    foreach ($field in @("config_sha256", "design_sha256", "row_count", "chunk_size", "expected_n_obs")) {
        if ([string]$saved.$field -ne [string]$context[$field]) {
            throw "Fixed-design resume context changed at '$field'. Use a new Parent directory."
        }
    }
} else {
    $context | ConvertTo-Json | Set-Content $contextPath
}
if (-not (Test-Path $manifest)) {
    "chunk,start,stop,outdir,status,started,ended" | Set-Content $manifest
}

$chunk = 0
for ($start = 0; $start -lt $rowCount; $start += $ChunkSize) {
    if ($DeadlineUtc) {
        $deadline = [datetime]::Parse($DeadlineUtc).ToUniversalTime()
        if ([datetime]::UtcNow.AddMinutes($MinRemainingMinutes) -ge $deadline) {
            Write-Warning "Pausing before chunk $chunk to respect the wall-clock budget."
            exit 75
        }
    }
    $stop = [Math]::Min($start + $ChunkSize, $rowCount)
    $outdir = Join-Path $Parent ("chunk_{0:D3}" -f $chunk)
    $resultDesign = Join-Path $outdir "design.csv"
    if (Test-Path $resultDesign) {
        $evaluated = @(Import-Csv $resultDesign)
        $expectedRows = $stop - $start
        $nObs = @($evaluated | ForEach-Object { [int]$_.n_obs } | Sort-Object -Unique)
        $scoresFinite = @($evaluated | Where-Object {
            $value = 0.0
            -not [double]::TryParse([string]$_.score, [ref]$value) -or
            [double]::IsNaN($value) -or [double]::IsInfinity($value)
        }).Count -eq 0
        $obsValid = $nObs.Count -eq 1 -and $nObs[0] -gt 0 -and
            ($ExpectedNObs -le 0 -or $nObs[0] -eq $ExpectedNObs)
        if ($evaluated.Count -eq $expectedRows -and $obsValid -and $scoresFinite) {
            $chunk += 1
            continue
        }
        Write-Warning "Chunk $chunk checkpoint is incomplete or invalid; rerunning it."
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
            if ($nObs.Count -ne 1 -or $nObs[0] -lt 1 -or
                ($ExpectedNObs -gt 0 -and $nObs[0] -ne $ExpectedNObs)) {
                $status = "failed_observation_count"
            }
            $badScores = @($evaluated | Where-Object {
                $value = 0.0
                -not [double]::TryParse([string]$_.score, [ref]$value) -or
                [double]::IsNaN($value) -or [double]::IsInfinity($value)
            })
            if ($badScores.Count -gt 0) {
                $status = "failed_nonfinite_score"
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
    $combinedRows = @(Import-Csv (Join-Path $combined "design.csv")).Count
    if ($combinedRows -ne $rowCount) {
        throw "Combined design has $combinedRows rows; expected $rowCount."
    }
}
