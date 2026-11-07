#!/usr/bin/env python3
"""Sequential synthetic-case mesh refinement, retaining each native solution."""
import json
import argparse
import math
from pathlib import Path
import shutil
import subprocess
import sys
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from pcb_pi.config import load
from pcb_pi.common import digest,write_json

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('case',nargs='?',type=Path,default=root/'examples/synthetic-pdn/case.json')
p.add_argument('--scales',nargs='+',type=float,default=[1.,.75,.5,.375,.25,.1875,.125])
a=p.parse_args()
if len(a.scales)<3 or a.scales[0]!=1 or any(not math.isfinite(v) or v<=0 for v in a.scales) or any(x<=y for x,y in zip(a.scales,a.scales[1:])):
    p.error('scales must start at 1 and contain at least three positive, strictly decreasing values')
original=a.case.resolve()
c,project,board,base=load(original)
subprocess.run([str(root/'pcb-pi'),'verify',str(original)],check=True)
records=[];changes=[];accepted=False
for scale in a.scales:
    if scale==1.:
        out=base
    else:
        folder=root/'.cache/refinement'/f'{scale:g}';folder.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(board,folder/'board.kicad_pcb')
        case=json.loads(json.dumps(c));case['board']='board.kicad_pcb';case['project_root']='.';case['mesh']['scale']=c['mesh']['scale']*scale
        out=root/'results/refinement'/f'{scale:g}';case['output']=str(out)
        write_json(folder/'case.json',case)
        subprocess.run([str(root/'pcb-pi'),'run',str(folder/'case.json'),'--analysis','electrothermal'],check=True)
        subprocess.run([str(root/'pcb-pi'),'verify',str(folder/'case.json')],check=True)
    final=json.loads((out/'electrothermal/final.json').read_text());e=final['electrical'];t=final['thermal']
    metrics={'copper_loss_w':e['copper_loss_w'],'temperature_rise_c':t['maximum_temperature_c']-c['cooling']['ambient_c']}
    metrics.update({name+'_pcb_drop_v':drops['supply_drop_v']+drops['return_rise_v'] for name,drops in e['path_drops'].items()})
    records.append({'scale':scale,'output':str(out),'metrics':metrics,'final_sha256':digest(out/'electrothermal/final.json')})
    write_json(root/'results/refinement/checkpoint.json',records)
    if len(records)>1:
        old,new=records[-2:]
        relative={k:abs(new['metrics'][k]-v)/max(abs(new['metrics'][k]),1e-12) for k,v in old['metrics'].items()}
        changes.append({'from':old['scale'],'to':new['scale'],'relative_changes':relative,'passed':all(v<=.02 for v in relative.values())})
        print('Mesh comparison: '+json.dumps(changes[-1]),flush=True)
    if len(records)>=3 and changes[-1]['passed']:
        subprocess.run([str(root/'pcb-pi'),'run',str(folder/'case.json'),'--analysis','all'],check=True)
        subprocess.run([str(root/'pcb-pi'),'verify',str(folder/'case.json')],check=True)
        accepted=True
        break
result={'board_sha256':digest(board),'case_sha256':digest(original),'records':records,'comparisons':changes,'passed':accepted,'accepted_output':str(out) if accepted else None,'criterion':'.02 relative change in final refinement; at least three meshes'}
write_json(root/'results/refinement/summary.json',result)
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit('Mesh refinement did not meet 2% criterion')
