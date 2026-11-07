"""Commands executed inside the pinned CPU image."""
import argparse
import json
import shutil
import sys
import unittest
from pathlib import Path
from .common import WorkflowError

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=('doctor','run','verify','test'));p.add_argument('case',nargs='?',type=Path)
    p.add_argument('--analysis',default='all');p.add_argument('--output',type=Path);p.add_argument('--force',action='store_true');p.add_argument('--integration',action='store_true');a=p.parse_args()
    try:
        if a.command=='doctor':
            import gmsh,pcbnew,numpy,scipy,shapely
            from .fem import native_module
            tools={n:shutil.which(n) for n in ('ElmerSolver','elmerf90','ngspice','pvpython')}
            if not all(tools.values()):raise WorkflowError('Missing tools: '+str(tools))
            import ctypes
            ctypes.CDLL(str(native_module('/tmp/pcb-pi-doctor')))
            print(json.dumps({'tools':tools,'kicad':pcbnew.GetBuildVersion(),'gmsh':gmsh.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'shapely':shapely.__version__,'cpu':True,'cuda':False},indent=2))
        elif a.command=='test':
            pattern='test_*.py' if a.integration else 'test_unit*.py'
            suite=unittest.defaultTestLoader.discover('/workflow/tests',pattern=pattern)
            if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():sys.exit(1)
        else:
            from .pipeline import run,verify
            if a.command=='run':run(a.case,a.output,a.analysis,a.force)
            else:verify(a.case,a.output)
    except (WorkflowError,KeyError,ValueError) as error:
        print(f'pcb-pi: {error}',file=sys.stderr);raise

if __name__=='__main__':main()
