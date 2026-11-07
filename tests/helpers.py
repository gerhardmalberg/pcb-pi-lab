import numpy as np
from pcb_pi.mesh import metrics,surfaces

def swept(xy,tri,z):
    n=len(xy);nodes=np.vstack([np.column_stack((xy,np.full(n,h))) for h in z]);elements=[]
    for j in range(len(z)-1):
        t=np.sort(tri,axis=1)+j*n;a,b,c=t.T;A,B,C=(t+n).T
        elements.extend(np.column_stack(q) for q in ((a,b,c,C),(a,b,B,C),(a,A,B,C)))
    tets=np.vstack(elements);faces,owners=surfaces(tets);v,g=metrics(nodes,tets)
    terminals={}
    for name,h in [('a',z[0]),('b',z[-1])]:
        f=faces[np.max(abs(nodes[faces,2]-h),axis=1)<1e-12]
        terminals[name]={'nodes':np.unique(f).tolist(),'faces':f.tolist()}
    return {'nodes':nodes,'tets':tets,'materials':np.full(len(tets),'copper'),'nets':np.full(len(tets),'test'),
            'electric_full':np.arange(len(tets)),'electric_nodes':np.arange(len(nodes)),'electric_tets':tets,
            'terminals':terminals,'exterior':faces,'exterior_owners':owners,'volumes':v,'gradients':g}

def bar(length=.01,width=.001,thickness=.0001,n=10):
    xy=np.array([[0,0],[width,0],[width,thickness],[0,thickness]])
    return swept(xy,np.array([[0,1,2],[0,2,3]]),np.linspace(0,length,n+1))

def barrel(length=.001,inner=.0002,plating=.000025,n=96):
    angles=np.linspace(0,2*np.pi,n,endpoint=False)
    xy=np.vstack([np.column_stack((r*np.cos(angles),r*np.sin(angles))) for r in (inner,inner+plating)])
    tri=[]
    for i in range(n):
        j=(i+1)%n;tri.extend([[i,j,n+j],[i,n+j,n+i]])
    return swept(xy,np.array(tri),[0,length/2,length])
