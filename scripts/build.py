#!/usr/bin/env python3
"""Reuse installed layers or rebuild independently from verified upstream assets."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request

root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--from-scratch',action='store_true');a=p.parse_args()
engine=os.environ.get('PCB_PI_ENGINE','podman');pins=json.loads((root/'env/base-pins.json').read_text());runtime=json.loads((root/'env/runtime.json').read_text())
if shutil.disk_usage(root).free<20*1024**3:raise SystemExit('At least 20 GiB free is required')
def call(args):subprocess.run(args,cwd=root,check=True)
base=runtime['reuse_tools_image']
if a.from_scratch:
    context=root/'.cache/build';context.mkdir(parents=True,exist_ok=True)
    assets={'cacert.pem':'cacert.pem','elmer.tar.gz':'elmer-26.1.tar.gz','paraview.tar.gz':'ParaView-5.13.3.tar.gz'}
    for dest,source in assets.items():
        info=pins['assets'][source];path=context/dest
        def valid():
            if not path.exists():return False
            with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()==info['sha256']
        if not valid():urllib.request.urlretrieve(info['url'],path)
        if not valid():raise SystemExit('Asset hash mismatch: '+dest)
    call([engine,'build','--build-arg','BASE_IMAGE='+pins['base_image'],'--build-arg','APT_SNAPSHOT='+pins['ubuntu_snapshot'],'--build-arg','BUILD_JOBS=2','-t','localhost/pcb-pi-base:0.1','-f',str(root/'containers/Containerfile.base'),str(context)])
    call([engine,'build','--build-arg','CPU_IMAGE=localhost/pcb-pi-base:0.1','-t','localhost/pcb-pi-tools:0.1','-f',str(root/'containers/Containerfile.tools'),str(context)])
    base='localhost/pcb-pi-tools:0.1'
else:
    if subprocess.run([engine,'image','exists',base]).returncode:
        raise SystemExit('Reusable tools image missing. Run pcb-pi build --from-scratch.')
call([engine,'build','--pull=never','--build-arg','TOOLS_IMAGE='+base,'-t',runtime['image'],'-f',str(root/'containers/Containerfile'),str(root)])
ident=subprocess.check_output([engine,'image','inspect',runtime['image'],'--format','{{.Id}}'],text=True).strip()
(root/'.cache').mkdir(exist_ok=True)
(root/'.cache/image.json').write_text(json.dumps({'image':runtime['image'],'id':ident,'base':base,'runtime':runtime},indent=2)+'\n')
