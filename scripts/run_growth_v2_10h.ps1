param(
    [string]$RunTag = "20260727",
    [string]$Python = "python",
    [double]$MaxHours = 10.0,
    [int]$LhsSamples = 800,
    [int]$ChunkSize = 200
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

$root = "results\global_hemp_calibration\growth_v2_${RunTag}"
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

$base = "results\global_hemp_calibration\focused_global_20260726_v2\configs\final_verification.yaml"
$reference = "calibration_global_hemp\stage2_global_16_joint_growth_from_phen_de_lhs5000.yaml"
$baseline = "results\global_hemp_calibration\focused_global_20260726_v2\final\best_theta.json"

function Write-Status {
    param([string]$Phase, [string]$State, [string]$Message)
    [ordered]@{
        run_tag = $RunTag; phase = $Phase; state = $State; message = $Message
        started_utc = $started.ToString("o")
        deadline_utc = $deadline.ToString("o")
        updated_utc = [datetime]::UtcNow.ToString("o")
    } | ConvertTo-Json | Set-Content $statusPath
    "$(Get-Date -Format s) [$Phase] $State - $Message" |
        Tee-Object -FilePath $logPath -Append
}

function Assert-Time {
    param([string]$Phase, [double]$ReserveHours = 0.25)
    if ([datetime]::UtcNow.AddHours($ReserveHours) -ge $deadline) {
        Write-Status $Phase "paused_budget" "Paused between resumable phases."
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

function Invoke-Calibration {
    param([string]$Config, [string]$Outdir, [string]$Phase)
    if (Test-Path (Join-Path $Outdir "best_theta.json")) {
        Write-Status $Phase "skipped" "Completed result exists."
        return
    }
    Assert-Time $Phase
    Write-Status $Phase "running" "Starting calibration."
    Invoke-Python @("run_calibration.py", $Config, "--outdir", $Outdir)
    if (-not (Test-Path (Join-Path $Outdir "best_theta.json"))) {
        throw "$Phase did not produce best_theta.json."
    }
    Write-Status $Phase "complete" "Calibration completed."
}

try {
    foreach ($stress in @("source", "potential")) {
        $cfg = Join-Path $configDir "diagnostic_${stress}.yaml"
        if (-not (Test-Path $cfg)) {
            Invoke-Python @(
                "scripts\build_growth_v2_config.py", $base, $reference, $baseline, $cfg,
                "--mode", "diagnostic", "--stress", $stress,
                "--name", "growth_v2_${RunTag}_diagnostic_${stress}", "--samples", "1"
            )
        }
        Invoke-Calibration $cfg (Join-Path $root "diagnostic_${stress}") "diagnostic_${stress}"
    }

    $sourceScore = [double](Import-Csv (Join-Path $root "diagnostic_source\design.csv"))[0].score
    $potentialScore = [double](Import-Csv (Join-Path $root "diagnostic_potential\design.csv"))[0].score
    $stress = if ($sourceScore -le $potentialScore) { "source" } else { "potential" }
    [ordered]@{
        selected = $stress; source_score = $sourceScore; potential_score = $potentialScore
    } | ConvertTo-Json | Set-Content (Join-Path $root "stress_mode_selection.json")
    Write-Status "stress_selection" "complete" "Selected $stress controls."

    $searchConfig = Join-Path $configDir "targeted_lhs${LhsSamples}_${stress}.yaml"
    if (-not (Test-Path $searchConfig)) {
        Invoke-Python @(
            "scripts\build_growth_v2_config.py", $base, $reference, $baseline,
            $searchConfig, "--mode", "search", "--stress", $stress,
            "--name", "growth_v2_${RunTag}_lhs${LhsSamples}_${stress}",
            "--samples", "$LhsSamples"
        )
    }

    $sensitivityOut = Join-Path $root "targeted_local_sensitivity"
    if (-not (Test-Path (Join-Path $sensitivityOut "local_sensitivity.csv"))) {
        Assert-Time "targeted_sensitivity"
        Write-Status "targeted_sensitivity" "running" "Running real-DSSAT OAT screen."
        Invoke-Python @(
            "scripts\run_best_local_sensitivity.py", $searchConfig, $baseline,
            $searchConfig, $sensitivityOut, "--fraction", "0.08", "--cores", "0"
        )
        Write-Status "targeted_sensitivity" "complete" "OAT screen completed."
    }

    $lhsParent = Join-Path $root "targeted_lhs${LhsSamples}"
    $lhsCombined = "${lhsParent}_combined"
    if (-not (Test-Path (Join-Path $lhsCombined "best_theta.json"))) {
        Assert-Time "targeted_lhs" 0.6
        Write-Status "targeted_lhs" "running" "Starting or resuming LHS chunks."
        & powershell.exe -NoProfile -ExecutionPolicy Bypass -File `
            "scripts\run_fixed_design_chunks.ps1" -Config $searchConfig `
            -Parent $lhsParent -ChunkSize $ChunkSize -Python $Python `
            -ExpectedNObs 456 -DeadlineUtc $deadline.ToString("o")
        if ($LASTEXITCODE -eq 75) { throw "Targeted LHS paused at budget." }
        if ($LASTEXITCODE -ne 0) { throw "Targeted LHS failed." }
        Write-Status "targeted_lhs" "complete" "All LHS chunks validated."
    }

    $starts = Join-Path $root "de_starts"
    if (-not (Test-Path (Join-Path $starts "selected_starts.csv"))) {
        Invoke-Python @(
            "scripts\select_diverse_lhs_starts.py", $searchConfig,
            (Join-Path $lhsCombined "design.csv"), $starts, "--count", "3", "--pool", "50"
        )
    }

    $deRows = @()
    foreach ($index in 1..3) {
        $seed = 2101 + 1000 * ($index - 1)
        $deConfig = Join-Path $configDir "de_${index}.yaml"
        $deOut = Join-Path $root "de_${index}"
        if (-not (Test-Path $deConfig)) {
            Invoke-Python @(
                "scripts\build_local_calibration_config.py", $searchConfig,
                (Join-Path $starts "start_${index}.json"), $deConfig,
                "--name", "growth_v2_${RunTag}_de_${index}", "--fraction", "0.18",
                "--engine", "de", "--de-popsize", "1", "--de-maxiter", "3",
                "--seed", "$seed"
            )
        }
        Invoke-Calibration $deConfig $deOut "de_${index}"
        $score = [double](Import-Csv (Join-Path $deOut "design.csv"))[0].score
        $deRows += [pscustomobject]@{
            run = $index; seed = $seed; score = $score
            theta = (Join-Path $deOut "best_theta.json")
        }
    }
    $deRows | Export-Csv (Join-Path $root "de_comparison.csv") -NoTypeInformation
    $bestDe = $deRows | Sort-Object score | Select-Object -First 1

    $finalConfig = Join-Path $configDir "final.yaml"
    if (-not (Test-Path $finalConfig)) {
        Invoke-Python @(
            "scripts\build_growth_v2_config.py", $base, $reference, $baseline,
            $finalConfig, "--mode", "diagnostic", "--stress", $stress,
            "--name", "growth_v2_${RunTag}_final", "--samples", "1",
            "--override", $bestDe.theta
        )
    }
    $finalOut = Join-Path $root "final"
    Invoke-Calibration $finalConfig $finalOut "final_verification"

    $cvConfig = Join-Path $configDir "loeo.yaml"
    if (-not (Test-Path $cvConfig)) {
        Invoke-Python @(
            "scripts\build_growth_v2_config.py", $base, $reference, $baseline,
            $cvConfig, "--mode", "search", "--stress", $stress,
            "--name", "growth_v2_${RunTag}_loeo", "--samples", "25",
            "--override", (Join-Path $finalOut "best_theta.json")
        )
    }
    $cvOut = Join-Path $root "loeo"
    if (-not (Test-Path (Join-Path $cvOut "cv_summary.json"))) {
        Assert-Time "loeo" 1.0
        Write-Status "loeo" "running" "Running 12 leave-one-environment-out refits."
        Invoke-Python @("scripts\run_targeted_loeo.py", $cvConfig, $cvOut)
        Write-Status "loeo" "complete" "Cross-validation completed."
    }

    Write-Status "workflow" "complete" "Growth v2 workflow completed end to end."
} catch {
    Write-Status "workflow" "failed_or_paused" $_.Exception.Message
    throw
}
