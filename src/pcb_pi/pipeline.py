"""Stage caching, coupled solves and auditable source preservation."""
from pathlib import Path
import json
import csv
import os
import shutil
import numpy as np
from .common import ROOT,WorkflowError,digest,fingerprint,write_json
from .config import load
from . import extract,mesh as geometry,fem,circuit,report

def implementation_hashes():
    return {str(p.relative_to(ROOT)):digest(p) for p in (ROOT/'src').rglob('*') if p.suffix in ('.py','.F90')}

def identity(case,board):
    return fingerprint({'case':case,'board_sha256':digest(board),'image':os.environ.get('PCB_PI_IMAGE_ID'),
                        'implementation':implementation_hashes()})

def seal(folder,key):
    folder=Path(folder)
    # Mesh symlinks refer to directories inside this stage; hash their files via base/.
    files={str(p.relative_to(folder)):digest(p) for p in folder.rglob('*') if p.is_file() and p.name!='stage.json'}
    write_json(folder/'stage.json',{'key':key,'artifacts':files})

def cached(folder,key):
    p=Path(folder)/'stage.json'
    if not p.exists():return False
    try:
        record=json.loads(p.read_text())
        return record['key']==key and bool(record['artifacts']) and all((p.parent/k).is_file() and digest(p.parent/k)==v for k,v in record['artifacts'].items())
    except (OSError,ValueError,KeyError):return False

def reset(folder):
    folder=Path(folder)
    if folder.exists():shutil.rmtree(folder)
    folder.mkdir(parents=True)

def conductivity(case,mesh,T):
    c=case['materials']['copper'];temperature=T[mesh['tets'][mesh['electric_full']]].mean(axis=1)
    if np.max(T)>c['maximum_temperature_c'] or np.min(T)<=-273.15:
        raise WorkflowError('Temperature exceeds the declared copper model range')
    rho=c['resistivity_ohm_m']*(1+c['temperature_coefficient_per_k']*(temperature-c['reference_temperature_c']))
    if np.any(rho<=0):raise WorkflowError('Nonpositive temperature-dependent resistivity')
    return 1/rho

def dc_solve(case,mesh,T,out,module):
    sigma=conductivity(case,mesh,T)
    network=fem.electrical(mesh,sigma,out/'extraction',module,case['numerics']['solver_timeout_s'])
    terminals,loads=circuit.dc(case,network,out/'circuit')
    fields=fem.recover(network,terminals)
    source=case['source'];drops={}
    for n,l in case['loads'].items():
        drops[n]={'supply_drop_v':terminals[source['positive']]-terminals[l['positive']],
                  'return_rise_v':terminals[l['negative']]-terminals[source['negative']]}
    result={k:v for k,v in fields.items() if k not in ('potential','current_density','q')}
    total_current=sum(l['current_a'] for l in case['loads'].values())
    result['operating_point_loss_equivalent_resistance_ohm']=fields['copper_loss_w']/total_current**2 if total_current>0 else None
    result.update(load_voltages_v=loads,terminal_voltages_v=terminals,path_drops=drops,temperature_min_c=float(T.min()),temperature_max_c=float(T.max()),
                  design_limits_passed=all(case['limits']['minimum_load_v']<=v<=case['limits']['maximum_load_v'] for v in loads.values()))
    write_json(out/'summary.json',result);np.savez_compressed(out/'state.npz',temperature=T,potential=fields['potential'],q=fields['q'],current_density=fields['current_density'])
    with (out/'terminals.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['terminal','voltage_v','current_into_copper_a'])
        w.writerows((name,terminals[name],fields['terminal_currents_a'][name]) for name in network['names'])
    with (out/'loads.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['load','voltage_v','current_a','supply_drop_v','return_rise_v'])
        w.writerows((name,v,case['loads'][name]['current_a'],drops[name]['supply_drop_v'],drops[name]['return_rise_v']) for name,v in loads.items())
    return network,fields,result

def restore(mesh,out):
    a=dict(np.load(out/'extraction/basis.npz'));s=dict(np.load(out/'state.npz'))
    names=json.loads((out/'extraction/conductance.json').read_text())['terminals']
    K,vol,grad=fem.matrix(mesh['nodes'][mesh['electric_nodes']],mesh['electric_tets'],a['sigma'])
    n={'names':names,'G':a['conductance'],'U':a['potentials'],'sigma':a['sigma'],'K':K,'vol':vol,'grad':grad,'mesh':mesh}
    r=json.loads((out/'summary.json').read_text())
    fields=fem.recover(n,r['terminal_voltages_v'])
    return n,fields,r,s['temperature']

def run(path,out,analysis='all',force=False):
    case,project,board,_=load(path);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    key=identity(case,board);before=digest(board);case_before=digest(path)
    write_json(out/'case.resolved.json',case)
    provenance={'status':'running','key':key,'board_sha256':before,'image_id':os.environ.get('PCB_PI_IMAGE_ID'),
                'case':str(path.relative_to(project)),'case_sha256':case_before,'implementation_sha256':implementation_hashes(),
                'analysis':analysis,'temperature_policy':'steady electrothermal bias; fixed during fast electrical transients'}
    write_json(out/'provenance.json',provenance)
    try:
        module=fem.native_module(out/'native')
        seal(out/'native',key)
        for stage in ('extract','mesh'):
            folder=out/stage
            if force or not cached(folder,key):
                reset(folder)
                if stage=='extract':extract.extract(board,folder/'board.json',case['mesh']['chord_mm'])
                else:geometry.build(json.loads((out/'extract/board.json').read_text()),case,folder)
                seal(folder,key)
        mesh=geometry.read(out/'mesh');T=np.full(len(mesh['nodes']),case['cooling']['ambient_c'])
        if force or not cached(out/'dc',key):
            reset(out/'dc');network,fields,dc=dc_solve(case,mesh,T,out/'dc',module);seal(out/'dc',key)
        else:network,fields,dc,T=restore(mesh,out/'dc')
        summaries={'dc':dc};checks={'dc':{'solver_completed':True,'numerically_valid':True,'design_limits_passed':dc['design_limits_passed']}}
        if analysis!='dc':
            folder=out/'electrothermal'
            if force or not cached(folder,key):
                folder.mkdir(exist_ok=True);history=[];consecutive=0;previous=None;n=case['numerics']
                for iteration in range(int(n['max_iterations'])):
                    it=folder/f'iteration-{iteration:03d}';itkey=fingerprint([key,T.tolist()])
                    print(f'Electrothermal iteration {iteration+1}',flush=True)
                    if not force and cached(it,itkey):
                        network,fields,result,usedT=restore(mesh,it)
                        newT=np.load(it/'thermal/temperature.npy');balance=json.loads((it/'thermal/balance.json').read_text())
                    else:
                        reset(it);network,fields,result=dc_solve(case,mesh,T,it,module)
                        q=np.zeros(len(mesh['tets']));q[mesh['electric_full']]=fields['q']
                        newT,balance=fem.thermal(mesh,case,q,it/'thermal',n['solver_timeout_s'],module=module);np.save(it/'thermal/temperature.npy',newT);seal(it,itkey)
                    delta=float(np.max(abs(newT-T)))
                    metrics=np.array([*result['load_voltages_v'].values(),result['copper_loss_w']])
                    relative=None if previous is None else float(np.max(abs(metrics-previous)/np.maximum(abs(previous),1e-12)))
                    passed=delta<n['temperature_tolerance_c'] and relative is not None and relative<n['relative_tolerance']
                    consecutive=consecutive+1 if passed else 0
                    row={'iteration':iteration+1,'temperature_residual_c':delta,'electrical_relative_change':relative,'copper_loss_w':result['copper_loss_w'],**balance}
                    history.append(row);write_json(folder/'history.json',history)
                    with (folder/'history.csv').open('w',newline='') as f:
                        w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerows(history)
                    write_json(out/'checkpoint.json',{'stage':'electrothermal','iteration':iteration+1,'key':key,'temperature_residual_c':delta})
                    if np.max(newT)>case['materials']['copper']['maximum_temperature_c']:raise WorkflowError('Thermal result outside material validity')
                    if consecutive>=2:break
                    T=T+n['relaxation']*(newT-T);previous=metrics
                else:raise WorkflowError('Electrothermal iteration did not converge')
                # T is the last temperature actually used for the electrical solve.
                write_json(folder/'final.json',{'iteration_directory':it.name,'electrical':result,'thermal':balance,'temperature_mismatch_c':delta,'iterations':len(history)})
                seal(folder,key)
            final=json.loads((folder/'final.json').read_text())
            network,fields,result,T=restore(mesh,folder/final['iteration_directory'])
            summaries['electrothermal']=final;checks['electrothermal']={'solver_completed':True,'numerically_valid':True,'design_limits_passed':result['design_limits_passed']}
        for stage,func in [('ac',circuit.ac),('transient',circuit.transient)]:
            if analysis not in (stage,'all'):continue
            folder=out/stage
            if force or not cached(folder,key):reset(folder);func(case,network,folder);seal(folder,key)
            result=json.loads((folder/'summary.json').read_text());summaries[stage]=result
            if stage=='ac':valid=True;passed=all(l['within_target'] for p in result.values() for l in p.values())
            else:
                valid=all(v['passed'] for p in result.values() for v in p['timestep_verification'].values())
                passed=all(v['within_limits'] for p in result.values() for v in p['loads'].values())
            checks[stage]={'solver_completed':True,'numerically_valid':valid,'design_limits_passed':passed}
        viewT=T if analysis=='dc' else np.load(out/'electrothermal'/final['iteration_directory']/'thermal/temperature.npy')
        report.fields(mesh,network,fields,viewT,out/'fields');seal(out/'fields',key)
        report.publish(out,case,summaries,checks)
        if digest(board)!=before or digest(path)!=case_before or identity(case,board)!=key:
            raise WorkflowError('Inputs or implementation changed during execution')
        provenance['status']='completed' if all(c['numerically_valid'] for c in checks.values()) else 'verification_failed'
        write_json(out/'provenance.json',provenance);write_json(out/'checkpoint.json',{'stage':provenance['status'],'key':key,'analyses':list(checks)})
        if provenance['status']!='completed':raise WorkflowError('Numerical verification failed; inspect checks.json')
        publications=['case.resolved.json','checks.json','report.md','provenance.json','checkpoint.json']
        if 'ac' in summaries:publications.append('impedance.png')
        if 'transient' in summaries:publications.append('transient.png')
        write_json(out/'publication.json',{'key':key,'artifacts':{name:digest(out/name) for name in publications}})
        print(json.dumps(checks,indent=2),flush=True)
        return checks
    except Exception as error:
        provenance.update(status='failed',error=str(error));write_json(out/'provenance.json',provenance)
        raise

def verify(path,out):
    case,project,board,_=load(path);out=Path(out);key=identity(case,board)
    if not (out/'checks.json').exists():raise WorkflowError('No completed run to verify')
    provenance=json.loads((out/'provenance.json').read_text())
    if provenance.get('status')!='completed' or provenance.get('key')!=key:
        raise WorkflowError('Run is incomplete, failed or stale; inspect provenance.json')
    checks=json.loads((out/'checks.json').read_text())
    publication=json.loads((out/'publication.json').read_text())
    if publication['key']!=key or not all((out/name).is_file() and digest(out/name)==value for name,value in publication['artifacts'].items()):
        raise WorkflowError('Altered reports or run metadata')
    stages=['native','extract','mesh','fields']+list(checks)
    invalid=[x for x in stages if not cached(out/x,key)]
    if invalid:raise WorkflowError('Stale inputs or altered artifacts: '+', '.join(invalid))
    if not all(x['numerically_valid'] for x in checks.values()):raise WorkflowError('Numerical verification failed')
    write_json(out/'verification.json',{'artifact_hashes_passed':True,'source_hash_passed':True,'stages':stages,'mesh_convergence':'Run scripts/verify-convergence.py separately; not inferred from solver success'})
    print('Source and artifact verification passed; design limits remain separate.')
