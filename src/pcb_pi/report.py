"""Compact numerical evidence and portable visualization artifacts."""
from pathlib import Path
import html
import json
import subprocess
import numpy as np
from .common import write_json

def vtu(path,nodes,tets,point_data,cell_data):
    def array(name,a,typ='Float64'):
        a=np.asarray(a);components=a.shape[1] if a.ndim>1 else 1
        return f'<DataArray type="{typ}" Name="{html.escape(name)}" NumberOfComponents="{components}" format="ascii">'+' '.join(map(str,a.ravel()))+'</DataArray>'
    text=f'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="{len(nodes)}" NumberOfCells="{len(tets)}">'
    text+='<Points>'+array('Points',nodes)+'</Points><Cells>'+array('connectivity',tets.ravel(),'Int32')+array('offsets',np.arange(1,len(tets)+1)*4,'Int32')+array('types',np.full(len(tets),10),'UInt8')+'</Cells>'
    text+='<PointData>'+''.join(array(k,v) for k,v in point_data.items())+'</PointData><CellData>'+''.join(array(k,v) for k,v in cell_data.items())+'</CellData></Piece></UnstructuredGrid></VTKFile>'
    Path(path).write_text(text)

def fields(mesh,network,recovery,T,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    vtu(out/'electrical.vtu',mesh['nodes'][mesh['electric_nodes']],mesh['electric_tets'],{'Potential_V':recovery['potential']},
        {'Current_density_A_m2':recovery['current_density'],'Joule_heating_W_m3':recovery['q']})
    vtu(out/'thermal.vtu',mesh['nodes'],mesh['tets'],{'Temperature_C':T},{})
    script='''from paraview.simple import *
from pathlib import Path
root=Path(__file__).resolve().parent
e=XMLUnstructuredGridReader(FileName=[str(root/'electrical.vtu')])
t=XMLUnstructuredGridReader(FileName=[str(root/'thermal.vtu')])
Show(t)
ColorBy(GetDisplayProperties(t), ('POINTS', 'Temperature_C'))
ResetCamera()
SaveState(str(root/'view.pvsm'))
'''
    (out/'open-view.py').write_text(script)
    with (out/'paraview.log').open('w') as log:
        subprocess.run(['xvfb-run','-a','pvpython','--force-offscreen-rendering',str(out/'open-view.py')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)

def publish(out,case,summaries,checks):
    out=Path(out)
    lines=['# '+case['name'],'','CPU power-integrity model. Circuit inductances and capacitor parasitics are declared inputs; no broadband layout extraction.','',
           '| Analysis | Solver completed | Numerical checks passed | Design limits passed |','|---|---|---|---|']
    for name,s in checks.items():lines.append(f'| {name} | {s["solver_completed"]} | {s["numerically_valid"]} | {s.get("design_limits_passed","not applicable")} |')
    lines+=['','## Electrical operating points','','| State | Copper loss [W] | Load voltages [V] |','|---|---:|---|']
    for label,result in [('Cold',summaries['dc'])]+([('Heated',summaries['electrothermal']['electrical'])] if 'electrothermal' in summaries else []):
        loads=', '.join(f'{name}: {v:.6g}' for name,v in result['load_voltages_v'].items())
        lines.append(f'| {label} | {result["copper_loss_w"]:.6g} | {loads} |')
    if 'electrothermal' in summaries:
        final=summaries['electrothermal']
        lines+=['',f'Maximum board temperature: {final["thermal"]["maximum_temperature_c"]:.6g} C after {final["iterations"]} iterations. Final temperature residual: {final["temperature_mismatch_c"]:.6g} C.']
    lines+=['',f'Declared circuit-model frequency range: {case["model_frequency_range_hz"][0]:g}–{case["model_frequency_range_hz"][1]:g} Hz.']
    if 'ac' in summaries:lines+=['','![Driving-point impedance](impedance.png)']
    if 'transient' in summaries:lines+=['','![Load voltages](transient.png)']
    lines+=['','## Assumptions','']+['- '+x for x in case.get('assumptions',[])]+['','Electrical transients hold the recorded temperature field fixed. Copper current-density peaks at sharp corners are mesh-dependent. Component heat is included only through explicit extra_heat inputs.','',
      'Numerical tolerances do not establish hardware correlation. Mesh refinement is a separate verification; see verification.json when available.','',
      '## Results','', 'Machine-readable summaries are in each analysis directory. Native inputs, logs, meshes and circuit decks are retained. `fields/` contains voltage, current, loss and temperature data for ParaView.']
    (out/'report.md').write_text('\n'.join(lines)+'\n');write_json(out/'checks.json',checks)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if 'ac' in summaries:
        fig,ax=plt.subplots(figsize=(7,4))
        for pop in case['populations']:
            for load in case['loads']:
                p=out/'ac'/pop/load/(load+'.csv')
                a=np.loadtxt(p,delimiter=',',skiprows=1);ax.loglog(a[:,0],a[:,3],label=f'{pop}: {load}')
        ax.loglog(a[:,0],a[:,4],'k--',label='Target');ax.set(xlabel='Frequency [Hz]',ylabel='Driving-point impedance [ohm]');ax.legend();ax.grid(True,which='both',alpha=.2);fig.tight_layout();fig.savefig(out/'impedance.png',dpi=140);plt.close(fig)
    if 'transient' in summaries:
        fig,ax=plt.subplots(figsize=(7,4))
        for pop in case['populations']:
            a=np.loadtxt(out/'transient'/pop/'nominal/voltages.csv',delimiter=',',skiprows=1)
            for i,load in enumerate(case['loads']):ax.plot(a[:,0]*1e6,a[:,i+1],label=f'{pop}: {load}')
        ax.set(xlabel='Time [us]',ylabel='Load voltage [V]');ax.legend();ax.grid(True,alpha=.2);fig.tight_layout();fig.savefig(out/'transient.png',dpi=140);plt.close(fig)
