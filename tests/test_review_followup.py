from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import time
import numpy as np
import pandas as pd
from dssatcalibrator.engines.mcmc import chain_diagnostics
from dssatcalibrator.eval_cache import EvaluationCache
import dssatcalibrator.spawn as sp

def test_diagnostics_detect_nonmixing_and_do_not_hide_undefined_parameters():
    rng=np.random.default_rng(8)
    d=pd.DataFrame([{'step':s,'walker':w,'moving':rng.normal()+w*10,'constant':1.} for w in range(4) for s in range(500)])
    ess,rhat=chain_diagnostics(d,['moving'])
    assert rhat > 2 and ess < 20
    assert np.isnan(chain_diagnostics(d,['constant','moving'])[1])
    assert np.isnan(chain_diagnostics(d,['moving','constant'])[1])

def test_identical_spawns_are_serialized(tmp_path,monkeypatch):
    paths={k:tmp_path/k for k in ('soil','weather','genotype')};paths['root']=tmp_path
    for p in paths.values():p.mkdir(exist_ok=True)
    (tmp_path/'E1.MZX').write_text('*TREATMENTS\n@N R O C TNAME\n 1 1 1 0 first\n*CULTIVARS\n')
    exe=tmp_path/'fake.exe';exe.write_text('fake')
    cfg={'source':{'hemp_dir':str(tmp_path)},'calibrator':{'cache_spawns':True}}
    crop={'code':'MZ','genotype_stem':'MZCER048','filex_ext':'MZX','model':'MZCER048'}
    monkeypatch.setattr(sp,'resolve_dssat_paths',lambda cfg:paths)
    calls=[]
    def run(run_dir,*args,**kwargs):
        calls.append(1);time.sleep(.05);(run_dir/'PlantGro.OUT').write_text('new');return ''
    monkeypatch.setattr(sp,'_run_native_dssat',run)
    monkeypatch.setattr(sp.dssat_io,'collect_run_outputs',lambda p:{'plantgro':pd.DataFrame({'treatment':[1],'LAID':[2.]})})
    def job():return sp.spawn_and_run({'x':1},exp_id='E1',cfg=cfg,crop=crop,param_specs=[{'name':'x','group':'unused'}],run_root=tmp_path/'runs',exe=exe)
    with ThreadPoolExecutor(2) as pool:results=list(pool.map(lambda _:job(),range(2)))
    assert len(calls)==1
    assert sorted(x.status for x in results)==['cached','success']

def test_objective_cache_tracks_nested_weather_and_implementation(tmp_path,monkeypatch):
    import dssatcalibrator.eval_cache as ec
    paths={k:tmp_path/k for k in ('soil','weather','genotype')};paths['root']=tmp_path
    for p in paths.values():p.mkdir(exist_ok=True)
    monkeypatch.setattr(ec,'resolve_dssat_paths',lambda cfg:paths)
    w=paths['weather']/'nested'/'secondary.wth';w.parent.mkdir();w.write_text('old')
    cfg={'source':{'hemp_dir':str(tmp_path)},'calibrator':{'workdir':str(tmp_path)}}
    def context():
        c=EvaluationCache.from_setup(cfg,crop={'genotype_stem':'MZCER048','filex_ext':'MZX','code':'MZ'},specs=[],experiments=['E1'],treatments={'E1':[1]},obs_table=pd.DataFrame(),exe=str(tmp_path/'fake'))
        return c.context
    before=context();w.write_text('changed');assert context()!=before
