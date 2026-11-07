"""Native Elmer solves and independent tetrahedral energy/reaction audits."""
import json
import os
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from .common import ROOT, WorkflowError, digest, write_json
from .mesh import metrics

def native_module(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    source=Path('/opt/elmer-source/fem/src/modules/StatCurrentSolve.F90')
    target=out/'StatCurrentSolve.so';record=out/'module.json'
    helper=ROOT/'src/pcb_pi/PiFields.F90'
    key={'source_sha256':digest(source),'coefficient_source_sha256':digest(helper),'image_id':os.environ.get('PCB_PI_IMAGE_ID','unrecorded')}
    if record.exists() and target.exists():
        old=json.loads(record.read_text())
        if all(old.get(k)==v for k,v in key.items()) and old.get('binary_sha256')==digest(target) and (out/'PiFields.so').exists() and old.get('coefficient_binary_sha256')==digest(out/'PiFields.so'):return target
    subprocess.run(['elmerf90','-O2','-shared','-fPIC',str(source),'-o',str(target)],check=True,capture_output=True)
    subprocess.run(['elmerf90','-O2','-shared','-fPIC',str(helper),'-o',str(out/'PiFields.so')],check=True,capture_output=True)
    write_json(record,{**key,'binary_sha256':digest(target),'coefficient_binary_sha256':digest(out/'PiFields.so'),'upstream':'Elmer 26.1','license':'Static-current module LGPL-2.1-or-later; PiFields Apache-2.0'})
    return target

def matrix(nodes,tets,coefficient):
    vol,grad=metrics(nodes,tets)
    c=np.asarray(coefficient)
    if c.ndim==1:local=np.einsum('eik,ejk,e,e->eij',grad,grad,vol,c)
    else:local=np.einsum('eik,ek,ejk,e->eij',grad,c,grad,vol)
    i=np.repeat(tets,4,axis=1).ravel();j=np.tile(tets,(1,4)).ravel()
    return coo_matrix((local.ravel(),(i,j)),shape=(len(nodes),len(nodes))).tocsr(),vol,grad

def write_native(folder,nodes,tets,bodies,boundaries):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    face_map={tuple(sorted(face)):i+1 for i,t in enumerate(tets) for face in t[[[0,1,2],[0,1,3],[0,2,3],[1,2,3]]]}
    with (folder/'mesh.header').open('w') as f:f.write(f'{len(nodes)} {len(tets)} {sum(len(x) for x in boundaries.values())}\n2\n504 {len(tets)}\n303 {sum(len(x) for x in boundaries.values())}\n')
    with (folder/'mesh.nodes').open('w') as f:
        for i,p in enumerate(nodes,1):f.write(f'{i} -1 '+' '.join(f'{x:.17g}' for x in p)+'\n')
    with (folder/'mesh.elements').open('w') as f:
        for i,(t,b) in enumerate(zip(tets,bodies),1):f.write(f'{i} {b} 504 '+' '.join(str(int(x)+1) for x in t)+'\n')
    with (folder/'mesh.boundary').open('w') as f:
        index=0
        for tag,faces in boundaries.items():
            for face in faces:
                index+=1
                f.write(f'{index} {tag} {face_map[tuple(sorted(face))]} 0 303 '+' '.join(str(int(x)+1) for x in face)+'\n')

def header(variable,procedure):
    return ['Check Keywords Warn','Header','  Mesh DB "." "mesh"','End','Simulation',
      '  Coordinate System = Cartesian 3D','  Simulation Type = Steady State','  Steady State Max Iterations = 1','End',
      'Solver 1',f'  Equation = "{variable}"',f'  Procedure = {procedure}',f'  Variable = {variable}',
      '  Variable DOFs = 1','  Linear System Solver = Direct','  Linear System Direct Method = Cholmod',
      '  Linear System Symmetric = Logical True','  Linear System Abort Not Converged = Logical True',
      '  Linear System Convergence Tolerance = 1.0e-12','  Nonlinear System Max Iterations = 1',
      '  Steady State Convergence Tolerance = 1.0e-12','  Stabilize = Logical False','End',
      'Solver 2','  Exec Solver = After Simulation','  Equation = OutputFields',
      '  Procedure = "ResultOutputSolve" "ResultOutputSolver"','  Output File Name = "fields"',
      '  Vtu Format = Logical True','  Binary Output = Logical True','  Single Precision = Logical False',
      '  Save Geometry Ids = Logical True','End','Equation 1','  Active Solvers(2) = 1 2','End']

def execute(folder,lines,nodes,variable,timeout):
    folder=Path(folder)
    prefix=os.path.relpath(folder.resolve()/'fields',(folder/'mesh').resolve())
    lines=[f'  Output File Name = "{prefix}"' if 'Output File Name =' in line else line for line in lines]
    (folder/'case.sif').write_text('\n'.join(lines)+'\n')
    with (folder/'solver.log').open('w') as log:
        try:r=subprocess.run(['ElmerSolver','case.sif'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        except subprocess.TimeoutExpired as e:raise WorkflowError(f'Elmer timed out: {folder}') from e
    text=(folder/'solver.log').read_text(errors='replace')
    bad=('fatal:','did not converge','failed to converge','solution not converged','singular matrix')
    if r.returncode or 'ALL DONE' not in text.upper() or any(x in text.lower() for x in bad) or re.search(r'\bnan\b',text,re.I):
        raise WorkflowError(f'Elmer failed; inspect {folder}/solver.log')
    paths=sorted(folder.rglob('*.vtu'))
    if len(paths)!=1:raise WorkflowError(f'Expected one Elmer VTU result: {folder}')
    raw=paths[0].read_bytes();start=raw.find(b'<AppendedData')
    if start>=0:
        marker=raw.index(b'_',raw.index(b'>',start))+1
        root=ET.fromstring(raw[:start]+b'</VTKFile>')
        block=raw[marker:]
    else:root=ET.fromstring(raw);block=b''
    def read_array(a):
        if a.attrib.get('format')=='ascii':return np.fromstring(a.text,sep=' ')
        import struct
        offset=int(a.attrib['offset']);fmt='<Q' if root.attrib.get('header_type')=='UInt64' else '<I'
        size=struct.calcsize(fmt);length=struct.unpack_from(fmt,block,offset)[0]
        return np.frombuffer(block[offset+size:offset+size+length],dtype={'Float64':'<f8','Float32':'<f4'}[a.attrib['type']]).copy()
    piece=root.find('.//Piece')
    coords=read_array(piece.find('Points/DataArray')).reshape(-1,3)
    array=next((x for x in piece.findall('PointData/DataArray') if x.attrib.get('Name','').lower()==variable.lower()),None)
    if array is None:raise WorkflowError(f'Missing Elmer {variable} output')
    values=read_array(array)
    from scipy.spatial import cKDTree
    distances,indices=cKDTree(coords).query(nodes)
    if distances.max()>1e-10 or len(set(indices))!=len(nodes) or not np.isfinite(values).all():
        raise WorkflowError('Elmer node mapping or output values invalid')
    return values[indices]

def electrical(mesh,sigma,out,module,timeout=120):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    used=mesh['electric_nodes'];nodes=mesh['nodes'][used];tets=mesh['electric_tets']
    names=list(mesh['terminals']);K,vol,grad=matrix(nodes,tets,sigma)
    boundaries={i+1:np.array(mesh['terminals'][name]['faces']) for i,name in enumerate(names)}
    base=out/'base';write_native(base,nodes,tets,np.ones(len(tets),int),boundaries)
    potentials=[]
    for col,name in enumerate(names):
        folder=out/name;folder.mkdir(exist_ok=True)
        link=folder/'mesh'
        if not link.exists():link.symlink_to('../base',target_is_directory=True)
        np.savetxt(folder/'sigma.dat',sigma,header=str(len(sigma)),comments='',fmt='%.17g')
        lines=header('Potential',f'"{module.with_suffix("")}" "StatCurrentSolver"')
        lines += ['Material 1',f'  Electric Conductivity = Variable Coordinate 1',f'    Real Procedure "{module.parent}/PiFields" "PiSigma"','End',
                  'Body 1','  Target Bodies(1) = 1','  Equation = 1','  Material = 1','End']
        for i in range(len(names)):
            lines += [f'Boundary Condition {i+1}',f'  Target Boundaries(1) = {i+1}',f'  Potential = Real {int(col==i)}','End']
        potentials.append(execute(folder,lines,nodes,'Potential',timeout))
    U=np.column_stack(potentials);reaction=K@U
    boundary=np.unique(np.concatenate([mesh['terminals'][k]['nodes'] for k in names]))
    interior=np.setdiff1d(np.arange(len(nodes)),boundary)
    residual=float(np.max(np.abs(reaction[interior]))/max(np.max(np.abs(reaction[boundary])),1e-30)) if len(interior) else 0.
    G=U.T@reaction
    G=(G+G.T)/2
    raw=G.copy()
    # Remove roundoff in the floating-net gauge; never bridge separate nets.
    for i in range(len(names)):G[i,i]-=G[i].sum()
    scale=max(float(np.abs(G).max()),1e-30)
    if residual>1e-7 or np.linalg.eigvalsh(G).min() < -1e-8*scale:
        raise WorkflowError('Conductance extraction failed residual/passivity check')
    reactions=np.array([reaction[mesh['terminals'][k]['nodes']].sum(axis=0) for k in names])
    disagreement=float(np.max(np.abs(reactions-G))/scale)
    if disagreement>.005:raise WorkflowError('Terminal reaction and energy conductance disagree')
    np.savez_compressed(out/'basis.npz',potentials=U,conductance=G,sigma=sigma)
    write_json(out/'conductance.json',{'terminals':names,'siemens':G.tolist(),'interior_residual_relative':residual,
                                     'reaction_energy_disagreement':disagreement,'gauge_correction_relative':float(abs(raw-G).max()/scale)})
    return {'names':names,'G':G,'U':U,'sigma':np.asarray(sigma),'K':K,'vol':vol,'grad':grad,'mesh':mesh}

def recover(network,terminal_voltages):
    v=np.array([terminal_voltages[k] for k in network['names']]);u=network['U']@v
    t=network['mesh']['electric_tets'];gradient=np.einsum('ei,eik->ek',u[t],network['grad'])
    current=-network['sigma'][:,None]*gradient
    q=network['sigma']*np.sum(gradient**2,axis=1)
    powers=q*network['vol'];currents=network['G']@v
    terminal_power=float(v@currents);loss=float(powers.sum())
    error=abs(terminal_power-loss)/max(loss,1e-8)
    balance=abs(currents.sum())/max(np.abs(currents).sum()/2,1e-8)
    if error>.005 or balance>.005:raise WorkflowError('Electrical conservation failed')
    return {'potential':u,'current_density':current,'q':q,'copper_loss_w':loss,
            'terminal_currents_a':dict(zip(network['names'],map(float,currents))),
            'power_imbalance_relative':error,'current_imbalance_relative':balance,
            'maximum_current_density_a_m2':float(np.linalg.norm(current,axis=1).max())}

def thermal(mesh,case,q,out,timeout=120,module=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    nodes=mesh['nodes'];tets=mesh['tets'];vol=mesh['volumes']
    coeff=np.array([case['materials'][m]['k_w_mk'] for m in mesh['materials']])
    rho=np.array([case['materials'][m]['density_kg_m3'] for m in mesh['materials']])
    q=np.array(q,dtype=float)
    for heat in case.get('extra_heat',[]):
        ids=np.array(mesh['terminals'][heat['terminal']]['nodes'])
        full_nodes=mesh['electric_nodes'][ids]
        selected=np.sum(np.isin(tets,full_nodes),axis=1)>=3
        if not selected.any():raise WorkflowError('Extra heat has no adjacent volume')
        q[selected]+=heat['watts']/vol[selected].sum()
    properties=np.column_stack((coeff,rho))
    unique,groups=np.unique(properties,axis=0,return_inverse=True)
    write_native(out/'mesh',nodes,tets,groups+1,{1:mesh['exterior']})
    lines=header('Temperature','"HeatSolve" "HeatSolver"')
    module=module or native_module('/tmp/pcb-pi-native')
    np.savetxt(out/'heat.dat',q/rho,header=str(len(q)),comments='',fmt='%.17g')
    for i,(kx,ky,kz,density) in enumerate(unique,1):
        lines += [f'Material {i}',f'  Heat Conductivity(3) = {kx:.17g} {ky:.17g} {kz:.17g}',f'  Density = {density:.17g}','End',
          f'Body Force {i}','  Heat Source = Variable Coordinate 1',f'    Real Procedure "{module.parent}/PiFields" "PiHeat"','End',f'Body {i}',f'  Target Bodies(1) = {i}',
          '  Equation = 1',f'  Material = {i}',f'  Body Force = {i}','End']
    cooling=case['cooling'];lines += ['Boundary Condition 1','  Target Boundaries(1) = 1']
    if cooling['kind']=='convection':lines += [f'  Heat Transfer Coefficient = Real {cooling["h_w_m2k"]}',f'  External Temperature = Real {cooling["ambient_c"]+273.15}']
    else:lines += [f'  Temperature = Real {cooling["temperature_c"]+273.15}']
    lines += ['End']
    T=execute(out,lines,nodes,'Temperature',timeout)-273.15
    total=float(q@vol);faces=mesh['exterior'];xyz=nodes[faces]
    area=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
    if cooling['kind']=='convection':removed=float(np.sum(area*cooling['h_w_m2k']*(T[faces].mean(axis=1)-cooling['ambient_c'])))
    else:
        K,_,_=matrix(nodes,tets,coeff)
        source=np.bincount(tets.ravel(),weights=np.repeat(q*vol/4,4),minlength=len(nodes))
        removed=float(-np.sum((K@T-source)[np.unique(faces)]))
    imbalance=abs(removed-total)/max(abs(total),1e-8)
    if imbalance>.01:raise WorkflowError(f'Thermal energy imbalance {imbalance:.3%}')
    write_json(out/'balance.json',{'input_w':total,'removed_w':removed,'imbalance_relative':imbalance,'maximum_temperature_c':float(T.max())})
    return T,{'input_w':total,'removed_w':removed,'imbalance_relative':imbalance,'maximum_temperature_c':float(T.max())}
