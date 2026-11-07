"""Container launcher: host Python standard library only."""
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path
from .common import ROOT, WorkflowError, write_json
from .config import load

def main():
    p = argparse.ArgumentParser(description='PCB Power Integrity Lab (CPU)')
    s = p.add_subparsers(dest='command', required=True)
    s.add_parser('doctor')
    b = s.add_parser('build'); b.add_argument('--from-scratch', action='store_true')
    for name in ('run','verify'):
        q = s.add_parser(name); q.add_argument('case', type=Path)
        if name == 'run':
            q.add_argument('--analysis', choices=('dc','electrothermal','ac','transient','all'), default='all')
            q.add_argument('--force', action='store_true')
    q = s.add_parser('test'); q.add_argument('--integration', action='store_true')
    a = p.parse_args()
    try:
        if a.command == 'build':
            subprocess.run(['python3',str(ROOT/'scripts/build.py')]+(['--from-scratch'] if a.from_scratch else []),check=True)
            return
        engine = os.environ.get('PCB_PI_ENGINE','podman')
        image = os.environ.get('PCB_PI_IMAGE','localhost/pcb-pi-cpu:0.1')
        ident = subprocess.check_output([engine,'image','inspect',image,'--format','{{.Id}}'], text=True).strip()
        cmd=[engine,'run','--rm','--init','--network=none','--read-only','--tmpfs','/tmp:rw,size=512m',
             '--env','HOME=/tmp/pcb-pi-home','--env','MPLCONFIGDIR=/tmp/mpl','--env','PYTHONPATH=/workflow/src',
             '--env','OPENBLAS_NUM_THREADS=1','--env','OMP_NUM_THREADS=1','--env',f'PCB_PI_IMAGE_ID={ident}',
             '--volume',f'{ROOT}:/workflow:ro','--workdir','/workflow','--entrypoint','python3']
        if engine == 'podman':cmd += ['--userns=keep-id','--user',f'{os.getuid()}:{os.getgid()}']
        if a.command in ('run','verify'):
            c, project, board, out = load(a.case)
            if shutil.disk_usage(out.parent if out.parent.exists() else project).free < 20*1024**3:
                raise WorkflowError('Less than 20 GiB free; stop disk-heavy work')
            out.mkdir(parents=True,exist_ok=True)
            cmd += ['--volume',f'{project}:/project:ro','--volume',f'{out}:/output:rw']
            args=[a.command,'/project/'+str(a.case.resolve().relative_to(project)),'--output','/output']
            if a.command=='run':args += ['--analysis',a.analysis]+(['--force'] if a.force else [])
        elif a.command=='test':
            out=ROOT/'results/validation';out.mkdir(parents=True,exist_ok=True)
            cmd += ['--volume',f'{out}:/output:rw']
            args=['test']+(['--integration'] if a.integration else [])
        else: args=['doctor']
        subprocess.run(cmd+[image,'-m','pcb_pi.worker']+args,check=True)
    except (WorkflowError,subprocess.CalledProcessError,FileNotFoundError,KeyError,ValueError) as e:
        p.exit(1, f'pcb-pi: {e}\n')
