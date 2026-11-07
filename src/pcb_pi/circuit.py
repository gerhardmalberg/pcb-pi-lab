"""Explicit passive PDN circuits; no extracted or implied AC parasitics."""
from pathlib import Path
import copy
import subprocess
import numpy as np
from .common import WorkflowError, write_json

def base(case,network,population):
    names=network['names'];G=network['G']
    lines=['* PCB PI Lab: extracted DC conductance plus declared circuit parasitics',
           '.options reltol=1e-8 abstol=1e-12 vntol=1e-10 gmin=1e-15']
    # VCCS matrix retains the complete reciprocal conductance network, including
    # floating references. Do not collapse a net into an ideal circuit node.
    for i,a in enumerate(names):
        for j,b in enumerate(names):
            if G[i,j]!=0:lines.append(f'Gpcb_{i}_{j} {a} 0 {b} 0 {G[i,j]:.17g}')
    source=case['source'];p=source['positive'];n=source['negative']
    lines += [f'Vreference {n} 0 0',f'Vsupply supply_ideal {n} DC {source["voltage_v"]:.17g}',
              f'Rsource supply_ideal supply_series {source["resistance_ohm"]:.17g}']
    if source['inductance_h']>0:lines.append(f'Lsource supply_series {p} {source["inductance_h"]:.17g}')
    else:lines.append(f'Vsource_link supply_series {p} 0')
    for name,load in case['loads'].items():
        if load['inductance_h']>0:lines.append(f'Lload_{name} {load["positive"]} load_{name} {load["inductance_h"]:.17g}')
        else:lines.append(f'Vload_{name} {load["positive"]} load_{name} 0')
    for i,c in enumerate(case['populations'][population]):
        p=c['positive'];n=c['negative'];L=c['esl_h']+c['mounting_inductance_h'];r=c['esr_ohm']
        if L>0:lines.append(f'Lcap{i} {p} cap_l{i} {L:.17g}')
        else:lines.append(f'Vcap_l{i} {p} cap_l{i} 0')
        if r>0:lines.append(f'Rcap{i} cap_l{i} cap_c{i} {r:.17g}')
        else:lines.append(f'Vcap_r{i} cap_l{i} cap_c{i} 0')
        lines.append(f'Ccap{i} cap_c{i} {n} {c["capacitance_f"]:.17g}')
    return lines

def solve(lines,commands,expressions,out,timeout=120):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    deck=lines+['.control','set noaskquit','set wr_singlescale','set wr_vecnames','set numdgt=17']+commands+[
        'wrdata values.dat '+' '.join(expressions),'quit','.endc','.end']
    (out/'circuit.cir').write_text('\n'.join(deck)+'\n')
    with (out/'ngspice.log').open('w') as log:
        try:r=subprocess.run(['ngspice','-b','circuit.cir'],cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        except subprocess.TimeoutExpired as e:raise WorkflowError(f'ngspice timed out: {out}') from e
    log=(out/'ngspice.log').read_text(errors='replace').lower()
    if r.returncode or any(x in log for x in ('error:', 'singular matrix','timestep too small','no convergence','fatal')) or not (out/'values.dat').exists():
        raise WorkflowError(f'ngspice failed; inspect {out}/ngspice.log')
    a=np.loadtxt(out/'values.dat',skiprows=1,ndmin=2)
    if not np.isfinite(a).all():raise WorkflowError('Nonfinite circuit results')
    return a

def loads(case,transient=False):
    lines=[]
    for name,l in case['loads'].items():
        text=f'Iload_{name} load_{name} {l["negative"]} DC {l["current_a"]:.17g}'
        if transient and name==case['transient']['load']:
            t=case['transient'];a=l['current_a'];b=a+t['delta_a']
            if t.get('kind','pulse')=='step':
                text+=f' PWL(0 {a:.17g} {t["delay_s"]} {a:.17g} {t["delay_s"]+t["rise_s"]} {b:.17g} {t["stop_s"]} {b:.17g})'
            else:text+=f' PULSE({a:.17g} {b:.17g} {t["delay_s"]} {t["rise_s"]} {t["fall_s"]} {t["width_s"]} {2*t["stop_s"]})'
        lines.append(text)
    return lines

def dc(case,network,out):
    names=network['names'];loadnames=list(case['loads'])
    expressions=[f'v({n})' for n in names]+[f'v(load_{n},{case["loads"][n]["negative"]})' for n in loadnames]
    a=solve(base(case,network,case['baseline_population'])+loads(case),['op'],expressions,out,case['numerics']['solver_timeout_s'])
    values=dict(zip(names,map(float,a[0,1:1+len(names)])))
    voltages=dict(zip(loadnames,map(float,a[0,1+len(names):])))
    return values,voltages

def target(case,f):
    t=case['target_impedance']
    if 'envelope' in t:
        a=np.array(t['envelope']);return np.interp(np.log10(f),np.log10(a[:,0]),a[:,1])
    return np.full_like(f,t['allowed_dynamic_v']/t['delta_current_a'])

def ac(case,network,out):
    out=Path(out);results={}
    for population in case['populations']:
        rows={}
        for excited,l in case['loads'].items():
            lines=base(case,network,population)+loads(case)+[f'Iprobe {l["negative"]} load_{excited} AC 1']
            names=list(case['loads']);expressions=[f'v(load_{n},{case["loads"][n]["negative"]})' for n in names]
            ac=case['ac'];folder=out/population/excited
            data=solve(lines,[f'ac dec {ac["points_per_decade"]} {ac["start_hz"]} {ac["stop_hz"]}'],expressions,folder,case['numerics']['solver_timeout_s'])
            f=data[:,0];Z=data[:,1::2]+1j*data[:,2::2]
            limit=target(case,f)
            for j,name in enumerate(names):
                np.savetxt(folder/f'{name}.csv',np.column_stack((f,Z[:,j].real,Z[:,j].imag,abs(Z[:,j]),limit)),delimiter=',',header='frequency_hz,real_ohm,imag_ohm,magnitude_ohm,target_ohm',comments='')
            driving=abs(Z[:,names.index(excited)]);worst=int(np.argmax(driving/limit))
            rows[excited]={'maximum_impedance_ohm':float(driving.max()),'worst_ratio':float((driving/limit)[worst]),'worst_frequency_hz':float(f[worst]),'within_target':bool(np.all(driving<=limit))}
        results[population]=rows
    write_json(out/'summary.json',results)
    return results

def transient_metrics(time,voltage,t,limits,steady_v=None):
    baseline=float(np.mean(voltage[time<t['delay_s']*.8]))
    is_step=t.get('kind','pulse')=='step'
    if is_step and steady_v is None:raise WorkflowError('Step settling requires a solved final DC operating point')
    target=steady_v if is_step else baseline
    pulse_end=sum(t[k] for k in (('delay_s','rise_s') if is_step else ('delay_s','rise_s','width_s','fall_s')))
    after=time>=pulse_end;tailtime=time[after];tail=voltage[after]
    if len(tail)<3:raise WorkflowError('Insufficient transient recovery samples')
    bad=np.flatnonzero(abs(tail-target)>limits['settling_band_v'])
    if not len(bad):settling=0.
    elif bad[-1]==len(tail)-1:settling=None
    else:
        i=bad[-1];a=abs(tail[i]-target);b=abs(tail[i+1]-target)
        crossing=tailtime[i]+(a-limits['settling_band_v'])/(a-b)*(tailtime[i+1]-tailtime[i])
        settling=float(crossing-pulse_end)
    sufficient=t['stop_s']-pulse_end>=limits['settling_s']
    result={'baseline_v':baseline,'minimum_v':float(voltage.min()),'maximum_v':float(voltage.max()),
            'droop_v':max(0.,float(baseline-voltage.min())),'overshoot_v':max(0.,float(voltage.max()-baseline)),
            'settling_s':settling,'settling_target_v':target,'recovery_window_sufficient':sufficient}
    result['within_limits']=bool(result['droop_v']<=limits['maximum_droop_v'] and result['overshoot_v']<=limits['maximum_overshoot_v'] and result['minimum_v']>=limits['minimum_load_v'] and result['maximum_v']<=limits['maximum_load_v'] and sufficient and settling is not None and settling<=limits['settling_s'])
    return result

def transient(case,network,out):
    out=Path(out);results={};t=case['transient'];names=list(case['loads'])
    final={}
    if t.get('kind','pulse')=='step':
        final_case=copy.deepcopy(case);final_case['loads'][t['load']]['current_a']+=t['delta_a']
        _,final=dc(final_case,network,out/'final-operating-point')
    for population in case['populations']:
        expressions=[f'v(load_{n},{case["loads"][n]["negative"]})' for n in names]
        runs=[]
        for label,factor in [('nominal',1),('half-step',.5)]:
            dt=t['step_s']*factor
            data=solve(base(case,network,population)+loads(case,True),[f'tran {dt:.17g} {t["stop_s"]:.17g} 0 {dt:.17g}'],expressions,out/population/label,case['numerics']['solver_timeout_s'])
            np.savetxt(out/population/label/'voltages.csv',data,delimiter=',',header='time_s,'+','.join(n+'_v' for n in names),comments='')
            runs.append({n:transient_metrics(data[:,0],data[:,i+1],t,case['limits'],final.get(n)) for i,n in enumerate(names)})
        comparisons={}
        for n in names:
            differences={}
            for metric in ('droop_v','overshoot_v','settling_s'):
                x,y=runs[0][n][metric],runs[1][n][metric]
                atol=1e-12 if metric=='settling_s' else 1e-9
                differences[metric]=None if x is None or y is None else abs(x-y)/max(abs(y),atol/.01)
            comparisons[n]={'relative_changes':differences,'passed':all(v is not None and v<=.01 for v in differences.values())}
        results[population]={'loads':runs[0],'timestep_verification':comparisons}
    write_json(out/'summary.json',results)
    return results
