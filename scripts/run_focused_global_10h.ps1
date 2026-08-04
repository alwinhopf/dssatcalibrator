param(
    [string]$RunTag = "20260726",
    [string]$Python = "python",
    [double]$MaxHours = 10.0,
    [int]$ChunkSize = 250,
    [int]$PhenologySamples = 1200,
    [int]$GrowthSamples = 1800
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

$root = "results\global_hemp_calibration\focused_global_${RunTag}"
$configDir = Join-Path $root "configs"
$statusPath = Join-Path $root "workflow_status.json"
$logPath = Join-Path $root "workflow.log"
$startedPath = Join-Path $root "workflow_started_utc.txt"
New-Item -ItemType Directory -Force -Path $configDir | Out-Null

if (Test-Path $startedPath) {
    $started = [datetime]::Parse((Get-Content $startedPath -Raw).Trim()).ToUniversalTime()
} else {
    $started = [datetime]::UtcNow
    $started.ToString("o") | Set-Content $startedPath
}
$deadline = $started.AddHours($MaxHours)

$base = "calibration_global_hemp\stage2_global_16_joint_growth_from_phen_de_lhs5000.yaml"
$reference = $base
$best = "results\global_hemp_calibration\global_16_best_local_sensitivity_20260726\best_theta_writer_effective.json"
$recommendations = "results\global_hemp_calibration\global_16_best_local_sensitivity_20260726\focused_parameter_recommendations.csv"

function Write-Status {
    param([string]$Phase, [string]$State, [string]$Message)
    $record = [ordered]@{
        run_tag = $RunTag
        phase = $Phase
        state = $State
        message = $Message
        started_utc = $started.ToString("o")
        deadline_utc = $deadline.ToString("o")
        updated_utc = [datetime]::UtcNow.ToString("o")
    }
    $record | ConvertTo-Json | Set-Content $statusPath
    "$(Get-Date -Format s) [$Phase] $State - $Message" | Tee-Object -FilePath $logPath -Append
}

function Assert-Time {
    param([string]$Phase, [double]$ReserveHours = 0.0)
    if ([datetime]::UtcNow.AddHours($ReserveHours) -ge $deadline) {
        Write-Status $Phase "paused_budget" "Ten-hour budget reached; checkpoints are resumable."
        throw "Wall-clock budget reached before $Phase."
    }
}

function Invoke-Python {
    param([string[]]$Arguments)
    & $Python @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Python failed ($LASTEXITCODE): $($Arguments -join ' ')"
    }
}

function Invoke-Chunks {
    param([string]$Config, [string]$Parent, [string]$Phase, [int]$ExpectedNObs)
    $combined = "${Parent}_combined"
    if (Test-Path (Join-Path $combined "best_theta.json")) {
        Write-Status $Phase "skipped" "Validated combined result already exists."
        return
    }
    Assert-Time $Phase 0.6
    Write-Status $Phase "running" "Starting or resuming fixed-design chunks."
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File `
        "scripts\run_fixed_design_chunks.ps1" -Config $Config `
        -Parent $Parent -ChunkSize $ChunkSize -Python $Python `
        -ExpectedNObs $ExpectedNObs -DeadlineUtc $deadline.ToString("o")
    if ($LASTEXITCODE -eq 75) {
        Write-Status $Phase "paused_budget" "Paused cleanly between validated chunks."
        throw "$Phase paused at the wall-clock budget."
    }
    if ($LASTEXITCODE -ne 0) { throw "$Phase fixed-design run failed." }
    if (-not (Test-Path (Join-Path $combined "best_theta.json"))) {
        throw "$Phase did not produce a combined best theta."
    }
    Write-Status $Phase "complete" "All chunks validated and combined."
}

function Invoke-Calibration {
    param([string]$Config, [string]$Outdir, [string]$Phase)
    if (Test-Path (Join-Path $Outdir "best_theta.json")) {
        Write-Status $Phase "skipped" "Completed result already exists."
        return
    }
    Assert-Time $Phase 0.4
    Write-Status $Phase "running" "Starting calibration."
    Invoke-Python @("run_calibration.py", $Config, "--outdir", $Outdir)
    if (-not (Test-Path (Join-Path $Outdir "best_theta.json"))) {
        throw "$Phase did not produce best_theta.json."
    }
    Write-Status $Phase "complete" "Calibration completed."
}

try {
    $phenConfig = Join-Path $configDir "phenology_lhs${PhenologySamples}.yaml"
    if (-not (Test-Path $phenConfig)) {
        Invoke-Python @(
            "scripts\build_focused_global_config.py", $base, $reference, $best,
            $recommendations, $phenConfig, "--stage", "phenology",
            "--name", "focused_global_${RunTag}_phen_lhs${PhenologySamples}",
            "--samples", "$PhenologySamples"
        )
    }
    $phenParent = Join-Path $root "phenology_lhs${PhenologySamples}"
    Invoke-Chunks $phenConfig $phenParent "phenology_lhs" 16
    $phenCombined = "${phenParent}_combined"

    $phenDeConfig = Join-Path $configDir "phenology_de.yaml"
    if (-not (Test-Path $phenDeConfig)) {
        Invoke-Python @(
            "scripts\build_local_calibration_config.py", $phenConfig,
            (Join-Path $phenCombined "best_theta.json"), $phenDeConfig,
            "--name", "focused_global_${RunTag}_phen_de",
            "--fraction", "0.18", "--engine", "de",
            "--de-popsize", "1", "--de-maxiter", "5"
        )
    }
    $phenDeOut = Join-Path $root "phenology_de"
    Invoke-Calibration $phenDeConfig $phenDeOut "phenology_de"

    $growthConfig = Join-Path $configDir "growth_lhs${GrowthSamples}.yaml"
    if (-not (Test-Path $growthConfig)) {
        Invoke-Python @(
            "scripts\build_focused_global_config.py", $base, $reference, $best,
            $recommendations, $growthConfig, "--stage", "growth",
            "--name", "focused_global_${RunTag}_growth_lhs${GrowthSamples}",
            "--samples", "$GrowthSamples",
            "--override", (Join-Path $phenDeOut "best_theta.json")
        )
    }
    $growthParent = Join-Path $root "growth_lhs${GrowthSamples}"
    Invoke-Chunks $growthConfig $growthParent "growth_lhs" 456
    $growthCombined = "${growthParent}_combined"

    $growthDeConfig = Join-Path $configDir "growth_de.yaml"
    if (-not (Test-Path $growthDeConfig)) {
        Invoke-Python @(
            "scripts\build_local_calibration_config.py", $growthConfig,
            (Join-Path $growthCombined "best_theta.json"), $growthDeConfig,
            "--name", "focused_global_${RunTag}_growth_de",
            "--fraction", "0.15", "--engine", "de",
            "--de-popsize", "1", "--de-maxiter", "5"
        )
    }
    $growthDeOut = Join-Path $root "growth_de"
    Invoke-Calibration $growthDeConfig $growthDeOut "growth_de"

    $finalConfig = Join-Path $configDir "final_verification.yaml"
    if (-not (Test-Path $finalConfig)) {
        Invoke-Python @(
            "scripts\build_focused_global_config.py", $base, $reference, $best,
            $recommendations, $finalConfig, "--stage", "final",
            "--name", "focused_global_${RunTag}_final", "--samples", "1",
            "--override", (Join-Path $phenDeOut "best_theta.json"),
            "--override", (Join-Path $growthDeOut "best_theta.json")
        )
    }
    Invoke-Calibration $finalConfig (Join-Path $root "final") "final_verification"
    Write-Status "workflow" "complete" "End-to-end focused calibration completed."
} catch {
    Write-Status "workflow" "failed_or_paused" $_.Exception.Message
    throw
}
