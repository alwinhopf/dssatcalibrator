param(
    [string]$Python = "python",
    [int]$PollSeconds = 60,
    [int]$ChunkSize = 500
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

$phenBase = "calibration_global_hemp\stage1_global_16_from_1348_maxbias_lhs5000.yaml"
$phenLhsParent = "results\global_hemp_calibration\global_16_from_1348_rmse_maxbias_lhs5000_fixed"
$phenLhsCombined = "${phenLhsParent}_combined"
$phenDeConfig = "calibration_global_hemp\stage1_global_16_from_1348_rmse_maxbias_de.yaml"
$phenDeOut = "results\global_hemp_calibration\global_16_from_1348_rmse_maxbias_de"
$jointLhsConfig = "calibration_global_hemp\stage2_global_16_joint_growth_from_phen_de_lhs5000.yaml"
$jointLhsParent = "results\global_hemp_calibration\global_16_joint_growth_from_phen_de_lhs5000_fixed"
$jointLhsCombined = "${jointLhsParent}_combined"
$jointDeConfig = "calibration_global_hemp\stage2_global_16_joint_growth_de.yaml"
$jointDeOut = "results\global_hemp_calibration\global_16_joint_growth_de"

function Invoke-CheckedPython {
    param([string[]]$Arguments)
    & $Python @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Python command failed with exit code $LASTEXITCODE`: $($Arguments -join ' ')"
    }
}

function Wait-ForResult {
    param([string]$Path, [string]$Label)
    while (-not (Test-Path $Path)) {
        Write-Output "$(Get-Date -Format s) waiting for $Label"
        Start-Sleep -Seconds $PollSeconds
    }
    Write-Output "$(Get-Date -Format s) found $Label`: $Path"
}

function Invoke-Calibration {
    param([string]$Config, [string]$Outdir, [string]$Label)
    $best = Join-Path $Outdir "best_theta.json"
    if (Test-Path $best) {
        Write-Output "$(Get-Date -Format s) skipping completed $Label"
        return
    }
    Write-Output "$(Get-Date -Format s) starting $Label"
    Invoke-CheckedPython @("run_calibration.py", $Config, "--outdir", $Outdir)
    if (-not (Test-Path $best)) {
        throw "$Label finished without producing $best"
    }
}

Wait-ForResult (Join-Path $phenLhsCombined "best_theta.json") "phenology LHS result"

if (-not (Test-Path $phenDeConfig)) {
    Invoke-CheckedPython @(
        "scripts\build_local_calibration_config.py", $phenBase,
        (Join-Path $phenLhsCombined "best_theta.json"), $phenDeConfig,
        "--name", "global_16_from_1348_rmse_maxbias_de",
        "--fraction", "0.12", "--engine", "de",
        "--de-popsize", "2", "--de-maxiter", "5"
    )
}
Invoke-Calibration $phenDeConfig $phenDeOut "phenology DE refinement"

if (-not (Test-Path $jointLhsConfig)) {
    Invoke-CheckedPython @(
        "scripts\build_joint_growth_config.py", $phenBase,
        (Join-Path $phenDeOut "best_theta.json"), $jointLhsConfig,
        "--name", "global_16_joint_growth_from_phen_de_lhs5000",
        "--samples", "5000", "--phenology-fraction", "0.10",
        "--growth-fraction", "0.25"
    )
}

if (-not (Test-Path (Join-Path $jointLhsCombined "best_theta.json"))) {
    Write-Output "$(Get-Date -Format s) starting joint phenology/growth LHS"
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File `
        "scripts\run_fixed_design_chunks.ps1" -Config $jointLhsConfig `
        -Parent $jointLhsParent -ChunkSize $ChunkSize -Python $Python
    if ($LASTEXITCODE -ne 0) {
        throw "Joint fixed-design LHS failed with exit code $LASTEXITCODE"
    }
}
Wait-ForResult (Join-Path $jointLhsCombined "best_theta.json") "joint LHS result"

if (-not (Test-Path $jointDeConfig)) {
    Invoke-CheckedPython @(
        "scripts\build_local_calibration_config.py", $jointLhsConfig,
        (Join-Path $jointLhsCombined "best_theta.json"), $jointDeConfig,
        "--name", "global_16_joint_growth_de",
        "--fraction", "0.10", "--engine", "de",
        "--de-popsize", "1", "--de-maxiter", "5"
    )
}
Invoke-Calibration $jointDeConfig $jointDeOut "joint phenology/growth DE refinement"

Write-Output "$(Get-Date -Format s) global hemp calibration sequence complete"
