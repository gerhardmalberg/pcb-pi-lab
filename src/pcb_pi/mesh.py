"""Net-preserving planar overlay, Gmsh triangulation and conformal tetrahedra.

The sweep follows the Apache-2.0 PCB Thermal Lab construction, but never
unions different electrical nets or replaces them with thermal mixtures.
"""
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union, polygonize
from shapely.affinity import rotate, translate
from .common import WorkflowError, write_json

FACES = np.array([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])

def polygons(items):
    return unary_union([Polygon(p['outer'],p['holes']) for p in items])

def surfaces(tets):
    faces=tets[:,FACES].reshape(-1,3)
    keys=np.sort(faces,axis=1)
    _, first, counts=np.unique(keys,axis=0,return_index=True,return_counts=True)
    if counts.max()>2: raise WorkflowError('Nonmanifold volume mesh')
    idx=first[counts==1]
    return faces[idx],idx//4

def metrics(nodes,tets):
    x=nodes[tets]
    b=x[:,1:]-x[:,:1]
    det=np.linalg.det(b)
    if np.any(np.abs(det)<1e-24): raise WorkflowError('Degenerate tetrahedron')
    grad=np.zeros((len(tets),4,3)); grad[:,1:]=np.linalg.inv(b).transpose(0,2,1)
    grad[:,0]=-grad[:,1:].sum(axis=1)
    return np.abs(det)/6,grad

def triangulate(contours,region,size,out):
    import gmsh
    arrangement=shapely.unary_union(contours,grid_size=1e-6)
    cells=[p for p in polygonize(arrangement) if p.area>1e-12 and region.covers(p.representative_point())]
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        gmsh.model.add('pi-planar')
        vertices={};lines={}
        def point(p):
            key=tuple(round(float(x),6) for x in p[:2])
            if key not in vertices:vertices[key]=gmsh.model.geo.addPoint(*key,0,size)
            return vertices[key]
        def loop(ring):
            ids=[point(p) for p in list(ring.coords)[:-1]]
            ids=[v for i,v in enumerate(ids) if v!=ids[i-1]]
            if len(ids)<3:raise WorkflowError('Overlay sliver; reduce chord error')
            edges=[]
            for a,b in zip(ids,ids[1:]+ids[:1]):
                key=tuple(sorted((a,b)))
                if key not in lines:lines[key]=gmsh.model.geo.addLine(*key)
                edges.append(lines[key] if a<b else -lines[key])
            return gmsh.model.geo.addCurveLoop(edges)
        for p in cells:gmsh.model.geo.addPlaneSurface([loop(p.exterior)]+[loop(h) for h in p.interiors])
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber('Mesh.MeshSizeMax',size)
        gmsh.option.setNumber('Mesh.MeshSizeMin',size/10)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(out/'planar.msh'))
        tags,xyz,_=gmsh.model.mesh.getNodes();order=np.argsort(tags)
        tags=tags[order];xy=np.array(xyz).reshape(-1,3)[order,:2]
        _,ids=gmsh.model.mesh.getElementsByType(2)
        tri=np.searchsorted(tags,np.array(ids).reshape(-1,3))
        return xy,tri
    finally:gmsh.finalize()

def build(board,case,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    clean=lambda x: shapely.set_precision(shapely.make_valid(x),1e-6)
    outline=clean(polygons(board['outline']))
    if outline.is_empty:raise WorkflowError('Empty board outline')
    copper=defaultdict(list);pads=defaultdict(list)
    for item in board['copper']:
        g=clean(polygons(item['polygons']).intersection(outline))
        copper[item['layer'],item['net']].append(g)
        if item['kind']=='pad':pads[item['ref'],item['number'],item['layer']].append(g)
    copper={k:clean(unary_union(v)) for k,v in copper.items()}
    pads={k:clean(unary_union(v)) for k,v in pads.items()}
    if any('No saved fill' in x for x in board['warnings']):
        raise WorkflowError('Missing saved zone fill; save a filled board in KiCad first')
    cu_layers=[l['name'] for l in board['stackup'] if l['kind']=='copper']
    holes=[];barrels=defaultdict(list)
    for h in board['holes']:
        a,b=h['drill']; radius=min(a,b)/2
        if radius<=0:continue
        half=abs(a-b)/2
        shape=Point(0,0).buffer(radius,quad_segs=24) if not half else LineString([(-half,0),(half,0)] if a>b else [(0,-half),(0,half)]).buffer(radius,quad_segs=24)
        hole=clean(translate(rotate(shape,-h['angle_deg'],origin=(0,0)),*h['center']))
        if not outline.covers(hole):raise WorkflowError('Board edge intersects a drilled hole')
        holes.append(hole)
        if h['plated']:
            if set(h['layers'])!=set(cu_layers):raise WorkflowError('Only through-hole plating is supported; blind/buried via found')
            barrels[h['net']].append(clean(hole.buffer(case['geometry']['plating_mm'],join_style=2).difference(hole)))
    void=clean(unary_union(holes));region=clean(outline.difference(void))
    barrels={k:clean(unary_union(v)) for k,v in barrels.items()}
    copper={k:clean(v.difference(void)) for k,v in copper.items()}
    counts=Counter(x['ref'] for x in board['components'])
    duplicate=sorted({t['ref'] for t in case['terminals'].values() if counts[t['ref']]>1})
    if duplicate:raise WorkflowError('Duplicate component reference makes terminal selection ambiguous: '+', '.join(duplicate))
    pad_shapes={};components={x['ref']:x for x in board['components']}
    for name,t in case['terminals'].items():
        comp=components.get(t['ref'])
        if comp is None or comp['dnp']:raise WorkflowError(f'Missing or DNP terminal component: {t["ref"]}')
        matches=[p for p in comp['pads'] if p['number']==t['pad']]
        if len(matches)!=1 or matches[0]['net']!=t['net'] or t['layer'] not in matches[0]['layers']:
            raise WorkflowError(f'Terminal pad connectivity does not match declaration: {name}')
        key=(t['ref'],t['pad'],t['layer'])
        if key not in pads:raise WorkflowError(f'Terminal has no copper: {name}')
        pad_shapes[name]=pads[key].difference(void)
    for layer in cu_layers:
        items={net:clean(unary_union([g,barrels.get(net,Polygon())])) for (l,net),g in copper.items() if l==layer}
        for net,g in barrels.items():items.setdefault(net,g)
        for n,g in items.items():
            for m,h in items.items():
                if n<m and g.distance(h)<1e-7:raise WorkflowError(f'Electrical short between {n} and {m} on {layer}')
    contours=[region.boundary]+[g.boundary for g in copper.values()]+[g.boundary for g in barrels.values()]+[g.boundary for g in pad_shapes.values()]
    xy,tri=triangulate(contours,region,case['mesh']['size_mm']*case['mesh']['scale'],out)
    points=shapely.points(xy[tri].mean(axis=1));n=len(xy)
    stack=board['stackup']
    # Refine through thickness as well as in plane. Physical interfaces stay exact.
    divisions=max(1,int(round(2/case['mesh']['scale'])))
    levels=sorted({float(z) for l in stack for z in np.linspace(l['z_bottom_mm'],l['z_top_mm'],divisions+1)})
    nodes=np.vstack([np.column_stack((xy,np.full(n,z))) for z in levels])*1e-3
    elements=[];materials=[];nets=[]
    barrel_masks={net:shapely.covers(g,points) for net,g in barrels.items()}
    for j,(lo,hi) in enumerate(zip(levels,levels[1:])):
        layer=next((l for l in stack if l['z_bottom_mm']<=(lo+hi)/2<l['z_top_mm']),None)
        if layer is None:raise WorkflowError('Stackup has a gap')
        mat=np.full(len(tri),'soldermask' if layer['name'].endswith('.Mask') else 'dielectric',dtype=object)
        net=np.full(len(tri),'',dtype=object)
        if layer['kind']=='copper':
            for (l,k),g in copper.items():
                if l==layer['name']:
                    mask=shapely.covers(g,points);mat[mask]='copper';net[mask]=k
        if min(l['z_bottom_mm'] for l in stack if l['kind']=='copper')<(lo+hi)/2<max(l['z_top_mm'] for l in stack if l['kind']=='copper'):
            for k,mask in barrel_masks.items():mat[mask]='copper';net[mask]=k
        # Sorting makes shared prism-side diagonals agree globally.
        t=np.sort(tri,axis=1)+j*n;a,b,c=t.T;A,B,C=(t+n).T
        pieces=[np.column_stack(q) for q in ((a,b,c,C),(a,b,B,C),(a,A,B,C))]
        elements.extend(pieces);materials.extend([mat]*3);nets.extend([net]*3)
    tets=np.vstack(elements);mats=np.concatenate(materials);enets=np.concatenate(nets)
    vol,grad=metrics(nodes,tets)
    selected=np.isin(enets,list(case['nets'].values())) & (mats=='copper')
    if not selected.any():raise WorkflowError('No selected conductive volume')
    electric_full=np.flatnonzero(selected);used=np.unique(tets[selected]);lookup=np.full(len(nodes),-1);lookup[used]=np.arange(len(used))
    et=lookup[tets[selected]];efaces,eowners=surfaces(et)
    terminal_nodes={};terminal_faces={}
    for name,t in case['terminals'].items():
        layer=next(l for l in stack if l['name']==t['layer'])
        z=layer['z_top_mm'] if t['layer']=='F.Cu' else layer['z_bottom_mm']
        face_xyz=nodes[used][efaces]
        mask=(np.abs(face_xyz[:,:,2]-z*1e-3).max(axis=1)<1e-10) & shapely.covers(pad_shapes[name],shapely.points(face_xyz[:,:,:2].mean(axis=1)*1e3))
        faces=efaces[mask]
        if not len(faces):raise WorkflowError(f'No resolved pad contact faces for {name}')
        ids=np.unique(faces)
        if any(np.intersect1d(ids,x).size for x in terminal_nodes.values()):raise WorkflowError('Overlapping terminal contacts')
        terminal_nodes[name]=ids;terminal_faces[name]=faces
    # Every electrically selected connected component must meet >=2 terminals.
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    edges=np.concatenate([et[:,[0,i]] for i in (1,2,3)])
    graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(used),len(used)))
    count,labels=connected_components(graph,directed=False)
    for label in range(count):
        contacts=[name for name,ids in terminal_nodes.items() if np.any(labels[ids]==label)]
        if len(contacts)<2:raise WorkflowError('Disconnected selected copper or terminal; each island needs at least two terminals')
        if len({case['terminals'][k]['net'] for k in contacts})!=1:raise WorkflowError('Mesh shorts distinct nets')
    for net in case['nets'].values():
        groups={int(labels[terminal_nodes[k][0]]) for k,t in case['terminals'].items() if t['net']==net}
        if len(groups)!=1:raise WorkflowError(f'Disconnected terminals on net {net}')
    ext,owners=surfaces(tets)
    np.savez_compressed(out/'mesh.npz',nodes=nodes,tets=tets,materials=mats.astype(str),nets=enets.astype(str),electric_full=electric_full,electric_nodes=used,electric_tets=et,exterior=ext,exterior_owners=owners)
    terminal_info={k:{'nodes':v.tolist(),'faces':terminal_faces[k].tolist()} for k,v in terminal_nodes.items()}
    write_json(out/'terminals.json',terminal_info)
    write_json(out/'geometry.json',{'nodes':len(nodes),'tetrahedra':len(tets),'electrical_tetrahedra':len(et),'volume_m3':float(vol.sum()),'minimum_tetrahedron_volume_m3':float(vol.min()),'units':'m','construction':'Conformal Gmsh planar mesh with three tetrahedra per swept prism','electrical_components':int(count)})
    write_msh(out/'board.msh',nodes,tets,np.array([{'copper':1,'dielectric':2,'soldermask':3}[x] for x in mats]))
    return read(out)

def read(out):
    import json
    data=dict(np.load(Path(out)/'mesh.npz',allow_pickle=False))
    data['terminals']=json.loads((Path(out)/'terminals.json').read_text())
    data['volumes'],data['gradients']=metrics(data['nodes'],data['tets'])
    return data

def write_msh(path,nodes,tets,bodies):
    with open(path,'w') as f:
        f.write(f'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n{len(nodes)}\n')
        for i,p in enumerate(nodes,1):f.write(f'{i} '+' '.join(f'{x:.15g}' for x in p)+'\n')
        f.write(f'$EndNodes\n$Elements\n{len(tets)}\n')
        for i,(t,b) in enumerate(zip(tets,bodies),1):f.write(f'{i} 4 2 {b} {b} '+' '.join(str(x+1) for x in t)+'\n')
        f.write('$EndElements\n')
