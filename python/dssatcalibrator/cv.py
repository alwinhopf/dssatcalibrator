"""Cross-validation framework for dssatcalibrator.

Optional layer that runs calibrate() on train folds and evaluate_thetas() on
held-out folds, producing overfit diagnostics and a cv_report.

Activated by calling :func:`run_cross_validation` with an enabled
``cross_validation`` config block.
"""

from dataclasses import dataclass
from copy import deepcopy
from pathlib import Path
import json

import pandas as pd
import numpy as np

from . import orchestrator
from .orchestrator import _make_folds, _setup, _year_key
from .objective import ObjectiveResult

@dataclass
class FoldResult:
    fold: int
    label: str
    train_exps: list[str]
    test_exps: list[str]
    best_theta: dict
    cal_obj: ObjectiveResult
    val_obj: ObjectiveResult

@dataclass
class CVResult:
    strategy: str
    folds: list[FoldResult]
    summary: dict
    report_df: pd.DataFrame
    final_theta: dict
    final_theta_method: str

def parse_cv_config(cfg: dict) -> dict:
    """Extract and validate cross_validation block from config."""
    cv_cfg = cfg.get('cross_validation', {}) or {}
    strategy = str(cv_cfg.get('strategy', 'leave_one_out'))
    valid_strategies = {'leave_one_out', 'k_fold', 'leave_site_out', 'temporal_forward'}
    if strategy not in valid_strategies:
        raise ValueError(
            f"Unknown cross-validation strategy '{strategy}'. Expected one of: "
            f"{', '.join(sorted(valid_strategies))}."
        )
    k = int(cv_cfg.get('k', 5))
    if strategy == 'k_fold' and k < 2:
        raise ValueError("cross_validation.k must be at least 2 for k_fold.")
    final_theta = str(cv_cfg.get('final_theta', cv_cfg.get('final_model', 'full_refit')))
    if final_theta not in {'full_refit', 'best_fold', 'ensemble_mean'}:
        raise ValueError(
            "cross_validation.final_theta must be full_refit, best_fold, or ensemble_mean."
        )

    overfit_threshold = float(cv_cfg.get('overfit_threshold', 1.5))
    if not np.isfinite(overfit_threshold) or overfit_threshold < 1.2:
        raise ValueError("cross_validation.overfit_threshold must be finite and at least 1.2.")

    parsed = {
        'enabled': bool(cv_cfg.get('enabled', False)),
        'strategy': strategy,
        'k': k,
        'seed': int(cv_cfg.get('seed', cfg.get('calibrator', {}).get('seed', 42))),
        'report': bool(cv_cfg.get('report', True)),
        'overfit_threshold': overfit_threshold,
        'final_theta': final_theta,
    }
    return parsed

def _overfit_ratio(cal_rmse: float, val_rmse: float) -> float:
    """Calculate ratio of validation to calibration RMSE safely."""
    if pd.isna(cal_rmse) or cal_rmse == 0:
        return float('nan')
    return float(val_rmse) / float(cal_rmse)

def _aggregate_results(folds: list[FoldResult], cv_cfg: dict) -> tuple[dict, pd.DataFrame]:
    """Aggregate per-fold results into a report dataframe and summary dict."""
    rows = []

    all_vars = set()
    for f in folds:
        if f.cal_obj and getattr(f.cal_obj, 'per_var', None):
            all_vars.update(f.cal_obj.per_var.keys())
        if f.val_obj and getattr(f.val_obj, 'per_var', None):
            all_vars.update(f.val_obj.per_var.keys())
    all_vars = sorted(all_vars)

    for f in folds:
        row = {
            'fold': f.fold,
            'strategy': cv_cfg['strategy'],
            'train_experiments': ','.join(f.train_exps),
            'test_experiments': ','.join(f.test_exps),
            'cal_score': f.cal_obj.score if f.cal_obj else float('nan'),
            'val_score': f.val_obj.score if f.val_obj else float('nan')
        }

        fold_overfit_ratios = []
        for v in all_vars:
            cal_rmse = f.cal_obj.per_var.get(v, {}).get('RMSE', float('nan')) if f.cal_obj and f.cal_obj.per_var else float('nan')
            val_rmse = f.val_obj.per_var.get(v, {}).get('RMSE', float('nan')) if f.val_obj and f.val_obj.per_var else float('nan')
            cal_d = f.cal_obj.per_var.get(v, {}).get('d', float('nan')) if f.cal_obj and f.cal_obj.per_var else float('nan')
            val_d = f.val_obj.per_var.get(v, {}).get('d', float('nan')) if f.val_obj and f.val_obj.per_var else float('nan')

            row[f'cal_RMSE_{v}'] = cal_rmse
            row[f'val_RMSE_{v}'] = val_rmse
            row[f'cal_d_{v}'] = cal_d
            row[f'val_d_{v}'] = val_d

            ratio = _overfit_ratio(cal_rmse, val_rmse)
            if not pd.isna(ratio):
                fold_overfit_ratios.append(ratio)

        if fold_overfit_ratios:
            row['overfit_ratio'] = max(fold_overfit_ratios)
        else:
            row['overfit_ratio'] = float('nan')

        rows.append(row)

    report_df = pd.DataFrame(rows)

    cal_score_mean = report_df['cal_score'].mean() if not report_df.empty else float('nan')
    val_score_mean = report_df['val_score'].mean() if not report_df.empty else float('nan')
    overfit_ratio_mean = report_df['overfit_ratio'].mean() if not report_df.empty else float('nan')
    overfit_ratio_max = report_df['overfit_ratio'].max() if not report_df.empty else float('nan')

    per_var = {}
    for v in all_vars:
        if f'cal_RMSE_{v}' in report_df:
            v_cal = report_df[f'cal_RMSE_{v}'].mean()
            v_val = report_df[f'val_RMSE_{v}'].mean()
            ratio = _overfit_ratio(v_cal, v_val)
            per_var[v] = {
                'cal_RMSE_mean': float(v_cal) if not pd.isna(v_cal) else None,
                'val_RMSE_mean': float(v_val) if not pd.isna(v_val) else None,
                'overfit_ratio': float(ratio) if not pd.isna(ratio) else None,
            }

    if pd.isna(overfit_ratio_max):
        rec = 'good'
    elif overfit_ratio_max <= 1.2:
        rec = 'good'
    elif overfit_ratio_max <= float(cv_cfg['overfit_threshold']):
        rec = 'mild_overfit'
    else:
        rec = 'overfit'

    if not report_df.empty and not pd.isna(report_df['val_score'].max()):
        worst_idx = report_df['val_score'].idxmax()
        worst_fold = int(report_df.loc[worst_idx, 'fold'])
    else:
        worst_fold = None

    summary = {
        'cal_score_mean': float(cal_score_mean) if not pd.isna(cal_score_mean) else None,
        'val_score_mean': float(val_score_mean) if not pd.isna(val_score_mean) else None,
        'overfit_ratio_mean': float(overfit_ratio_mean) if not pd.isna(overfit_ratio_mean) else None,
        'overfit_ratio_max': float(overfit_ratio_max) if not pd.isna(overfit_ratio_max) else None,
        'per_variable': per_var,
        'recommendation': rec,
        'worst_fold': worst_fold
    }

    return summary, report_df

def _select_final_theta(cfg: dict, folds: list[FoldResult], progress: bool) -> tuple[dict, str]:
    """Select the final recommended parameter set."""
    cv_cfg = parse_cv_config(cfg)
    method = cv_cfg['final_theta']

    if method == 'full_refit':
        if progress:
            print("Running full refit for final theta...", flush=True)
        res = orchestrator.calibrate(cfg, progress=progress)
        return res.best_theta, method
    elif method == 'best_fold':
        if not folds:
            return {}, method
        valid = [f for f in folds if f.val_obj and np.isfinite(f.val_obj.score)]
        if not valid:
            raise ValueError("Cross-validation produced no finite validation scores.")
        best_fold = min(valid, key=lambda f: f.val_obj.score)
        return best_fold.best_theta, method
    elif method == 'ensemble_mean':
        if not folds:
            return {}, method
        first_theta = folds[0].best_theta
        keys = list(first_theta)
        expected = set(keys)
        if any(set(f.best_theta) != expected for f in folds[1:]):
            raise ValueError("Cross-validation folds have inconsistent parameter sets.")
        final_theta = {}
        for k in keys:
            final_theta[k] = sum(f.best_theta[k] for f in folds) / len(folds)
        return final_theta, method
    else:
        # Fallback to full refit
        res = orchestrator.calibrate(cfg, progress=progress)
        return res.best_theta, 'full_refit'

def run_cross_validation(cfg: dict, *, progress: bool = True) -> CVResult:
    """Main entry point for running cross validation."""
    cv_cfg = parse_cv_config(cfg)
    if not cv_cfg['enabled']:
        raise ValueError("Cross-validation is not enabled in config.")

    _, _, _, _, _, _, experiments, _ = _setup(cfg)

    strategy_map = {
        'leave_one_out': 'loeo',
        'k_fold': 'random',
        'leave_site_out': 'site',
        'temporal_forward': 'year'
    }
    scheme = strategy_map.get(cv_cfg['strategy'], 'loeo')
    seed = cv_cfg['seed']

    if len(experiments) < 2:
        raise ValueError("Cross-validation requires at least two active experiments.")
    fold_splits = _make_folds(experiments, scheme, seed, k=cv_cfg['k'])
    if cv_cfg['strategy'] == 'temporal_forward':
        # Expanding-window validation: a held-out year may only use experiments
        # from earlier years for training (never later years).
        fold_splits = sorted(
            fold_splits,
            key=lambda split: _year_key(split[1][0]) if split[1] else split[0],
        )
        accumulated = []
        fold_jobs = []
        for label, test_exps in fold_splits:
            if accumulated:
                fold_jobs.append((label, test_exps, list(accumulated)))
            accumulated.extend(test_exps)
    else:
        fold_jobs = [
            (label, test_exps, [e for e in experiments if e not in set(test_exps)])
            for label, test_exps in fold_splits
        ]

    folds = []
    for i, (label, test_exps, train_exps) in enumerate(fold_jobs):
        if progress:
            print(f"Running Fold {i+1}/{len(fold_jobs)}: {label}", flush=True)

        cfg_train = deepcopy(cfg)
        if not train_exps or not test_exps:
            continue
        cfg_train['experiments'] = train_exps

        cal_result = orchestrator.calibrate(cfg_train, progress=progress)

        cfg_test = deepcopy(cfg)
        cfg_test['experiments'] = test_exps

        val_objs, _ = orchestrator.evaluate_thetas(cfg_test, [cal_result.best_theta])
        if not val_objs:
            raise RuntimeError(f"Cross-validation fold '{label}' returned no validation result.")
        val_obj = val_objs[0]

        folds.append(FoldResult(
            fold=i+1,
            label=label,
            train_exps=train_exps,
            test_exps=test_exps,
            best_theta=cal_result.best_theta,
            cal_obj=cal_result.best,
            val_obj=val_obj
        ))

    if not folds:
        raise ValueError("Cross-validation produced no non-empty train/test folds.")
    final_theta, final_theta_method = _select_final_theta(cfg, folds, progress)
    summary, report_df = _aggregate_results(folds, cv_cfg)

    return CVResult(
        strategy=cv_cfg['strategy'],
        folds=folds,
        summary=summary,
        report_df=report_df,
        final_theta=final_theta,
        final_theta_method=final_theta_method
    )

def write_cv_report(result: CVResult, outdir: Path) -> list[Path]:
    """Write cross-validation summary and table to outdir."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths = []

    csv_path = outdir / 'cv_report.csv'
    result.report_df.to_csv(csv_path, index=False)
    paths.append(csv_path)

    json_path = outdir / 'cv_summary.json'
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(result.summary, f, indent=2, allow_nan=False)
    paths.append(json_path)

    print("\n--- Cross-Validation Summary ---")
    print(f"Strategy: {result.strategy}")
    print(f"Mean Cal Score: {result.summary.get('cal_score_mean')}")
    print(f"Mean Val Score: {result.summary.get('val_score_mean')}")
    print(f"Max Overfit Ratio: {result.summary.get('overfit_ratio_max')}")
    print(f"Recommendation: {result.summary.get('recommendation')}")
    print(f"Final Theta Method: {result.final_theta_method}")
    print("--------------------------------\n", flush=True)

    return paths
