"""Independent numerical benchmarks against analytic solutions, not snapshots."""
import copy
import json
import shutil
import unittest
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from pcb_pi.common import ROOT,WorkflowError,write_json
from pcb_pi import fem,circuit,mesh as geometry,extract,pipeline
from pcb_pi.config import validate
from helpers import bar,barrel

class NativeBenchmarks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=Path('/output/benchmarks');cls.out.mkdir(parents=True,exist_ok=True)
        cls.module=fem.native_module(cls.out/'native');cls.results={}
    def record(self,name,data):
        self.results[name]=data;write_json(self.out/'summary.json',self.results)
    def test_trace_resistance_and_joule_loss(self):
        mesh=bar();rho=1.724e-8;expected=rho*.01/(.001*.0001)
        net=fem.electrical(mesh,np.full(len(mesh['tets']),1/rho),self.out/'trace',self.module)
        actual=1/net['G'][0,0];err=abs(actual/expected-1);self.assertLess(err,.01)
        f=fem.recover(net,{'a':expected,'b':0.})
        self.assertLess(abs(f['copper_loss_w']/expected-1),.01)
        self.assertLess(f['power_imbalance_relative'],.005)
        self.record('trace',{'expected_ohm':expected,'computed_ohm':actual,'relative_error':err,'power_imbalance':f['power_imbalance_relative']})
    def test_plated_barrel(self):
        mesh=barrel();rho=1.724e-8;expected=rho*.001/(np.pi*((.000225)**2-(.0002)**2))
        net=fem.electrical(mesh,np.full(len(mesh['tets']),1/rho),self.out/'barrel',self.module)
        actual=1/net['G'][0,0];error=abs(actual/expected-1);self.assertLess(error,.01)
        self.record('barrel',{'expected_ohm':expected,'computed_ohm':actual,'relative_error':error})
    def test_uniform_temperature_and_zero_load(self):
        mesh=bar();rho=1.724e-8;factor=1+.00393*50
        net=fem.electrical(mesh,np.full(len(mesh['tets']),1/(rho*factor)),self.out/'uniform-hot',self.module)
        actual=1/net['G'][0,0];expected=rho*.01/(.001*.0001)*factor
        self.assertLess(abs(actual/expected-1),.01)
        fields=fem.recover(net,{'a':3.3,'b':3.3})
        self.assertLess(fields['copper_loss_w'],1e-12)
        c=validate(json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text()))
        T,balance=fem.thermal(mesh,c,np.zeros(len(mesh['tets'])),self.out/'zero-heat',module=self.module)
        np.testing.assert_allclose(T,c['cooling']['ambient_c'],atol=1e-6)
        self.record('uniform_temperature',{'expected_ohm':expected,'computed_ohm':actual,'zero_load_loss_w':fields['copper_loss_w'],'zero_heat_temperature_error_c':float(abs(T-c['cooling']['ambient_c']).max())})
    def test_parallel_current_sharing(self):
        a=bar();b=bar(width=.002);b['nodes'][:,0]+=.004;offset=len(a['nodes'])
        mesh={**a,'nodes':np.vstack((a['nodes'],b['nodes'])),'tets':np.vstack((a['tets'],b['tets']+offset))}
        mesh['electric_nodes']=np.arange(len(mesh['nodes']));mesh['electric_tets']=mesh['tets']
        mesh['terminals']={name:{'nodes':a['terminals'][name]['nodes']+(np.array(b['terminals'][name]['nodes'])+offset).tolist(),
                                'faces':a['terminals'][name]['faces']+(np.array(b['terminals'][name]['faces'])+offset).tolist()} for name in ('a','b')}
        net=fem.electrical(mesh,np.full(len(mesh['tets']),1/1.724e-8),self.out/'parallel',self.module)
        f=fem.recover(net,{'a':.001,'b':0.})
        powers=f['q']*net['vol'];ratio=float(powers[len(a['tets']):].sum()/powers[:len(a['tets'])].sum())
        self.assertAlmostEqual(ratio,2.,places=8)
        expected=1.724e-8*.01/(.003*.0001)
        self.assertLess(abs(1/net['G'][0,0]/expected-1),.01)
        self.record('parallel',{'current_ratio':ratio,'expected_ratio':2.,'equivalent_resistance_ohm':float(1/net['G'][0,0])})
    def test_heated_conductor_analytic(self):
        # Insulated lateral sides, fixed end temperatures. Uniform current gives
        # -k T'' = I^2 rho0/A^2 (1+alpha*(T-T0)); closed-form cosine solution.
        mesh=bar(length=.001,width=.0001,thickness=.0001,n=32)
        mesh['exterior']=np.vstack([mesh['terminals'][x]['faces'] for x in ('a','b')])
        c=validate(json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text()))
        c['materials']['copper']['k_w_mk']=[1,1,1];c['cooling']={'kind':'fixed','temperature_c':20.,'ambient_c':20.}
        T=np.full(len(mesh['nodes']),20.);rho0=1e-6;alpha=.01;current=.1;history=[]
        for iteration in range(30):
            rho=rho0*(1+alpha*(T[mesh['tets']].mean(axis=1)-20))
            folder=self.out/'heated'/str(iteration)
            net=fem.electrical(mesh,1/rho,folder/'dc',self.module)
            fields=fem.recover(net,{'a':current/net['G'][0,0],'b':0})
            solved,balance=fem.thermal(mesh,c,fields['q'],folder/'thermal')
            delta=float(abs(solved-T).max());history.append(delta)
            if delta<.0001:T=solved;break
            T=.5*(T+solved)
        else:self.fail('Analytic electrothermal benchmark did not converge')
        omega=current*np.sqrt(rho0*alpha)/1e-8;z=mesh['nodes'][:,2]
        exact=20+(np.cos(omega*(z-.0005))/np.cos(omega*.0005)-1)/alpha
        error=float(abs(T-exact).max()/max(exact-20))
        self.assertLess(error,.01);self.assertLess(balance['imbalance_relative'],.01)
        self.record('heated_conductor',{'relative_temperature_rise_error':error,'maximum_temperature_c':float(T.max()),'expected_maximum_c':float(exact.max()),'iterations':len(history),'heat_imbalance':balance['imbalance_relative']})
    def test_fixed_conductivity_and_extra_heat(self):
        mesh=bar(length=.001,width=.0001,thickness=.0001,n=24)
        mesh['exterior']=np.vstack([mesh['terminals'][x]['faces'] for x in ('a','b')])
        c=validate(json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text()));c['materials']['copper']['k_w_mk']=[1,1,1]
        c['cooling']={'kind':'fixed','temperature_c':20.,'ambient_c':20.}
        q=np.full(len(mesh['tets']),1e8);T,b=fem.thermal(mesh,c,q,self.out/'fixed-conductivity')
        expected=20+1e8*.001**2/8;self.assertLess(abs((T.max()-20)/(expected-20)-1),.01)
        c['extra_heat']=[{'terminal':'a','watts':.001}]
        _,extra=fem.thermal(mesh,c,q,self.out/'extra-heat')
        self.assertAlmostEqual(extra['input_w']-b['input_w'],.001,places=10)
        self.record('fixed_conductivity',{'temperature_c':float(T.max()),'expected_c':expected,'balance':b,'explicit_extra_heat_w':extra['input_w']-b['input_w']})
    def test_coupled_fixed_conductivity_limit(self):
        mesh=bar(length=.001,width=.0001,thickness=.0001,n=12)
        mesh['exterior']=np.vstack([mesh['terminals'][x]['faces'] for x in ('a','b')])
        c=validate(json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text()))
        c['materials']['copper'].update(k_w_mk=[1,1,1],temperature_coefficient_per_k=0.,resistivity_ohm_m=1e-6)
        c['cooling']={'kind':'fixed','ambient_c':20.,'temperature_c':20.}
        T=np.full(len(mesh['nodes']),20.);losses=[];temperatures=[]
        for i in range(2):
            folder=self.out/'fixed-coupling'/str(i)
            n=fem.electrical(mesh,pipeline.conductivity(c,mesh,T),folder/'dc',self.module)
            f=fem.recover(n,{'a':.1/n['G'][0,0],'b':0.})
            T,b=fem.thermal(mesh,c,f['q'],folder/'thermal',module=self.module)
            losses.append(f['copper_loss_w']);temperatures.append(float(T.max()))
        self.assertAlmostEqual(losses[0],losses[1],places=12)
        self.assertAlmostEqual(temperatures[0],temperatures[1],places=8)
        self.assertAlmostEqual(temperatures[1],32.5,places=6)
        self.record('fixed_coupling',{'copper_losses_w':losses,'maximum_temperatures_c':temperatures,'expected_maximum_c':32.5})

class CircuitBenchmarks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=Path('/output/circuit-benchmarks');cls.out.mkdir(parents=True,exist_ok=True);cls.results={}
    def setup_case(self):
        c=validate(json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text()))
        c['loads']={'load1':dict(c['loads']['load1'])};c['loads']['load1']['inductance_h']=0;c['loads']['load1']['current_a']=.2
        c['populations']={'baseline':[dict(c['populations']['baseline'][0])]}
        cap=c['populations']['baseline'][0];cap.update(positive='load1_p',negative='load1_n',capacitance_f=10e-6,esr_ohm=.03,esl_h=1e-9,mounting_inductance_h=1e-9)
        names=['src_p','src_n','load1_p','load1_n'];G=np.zeros((4,4))
        for i,j,r in [(0,2,.02),(1,3,.01)]:
            G[i,i]+=1/r;G[j,j]+=1/r;G[i,j]-=1/r;G[j,i]-=1/r
        return c,{'names':names,'G':G}
    def record(self,name,value):self.results[name]=value;write_json(self.out/'summary.json',self.results)
    def test_dc_return_path(self):
        c,n=self.setup_case();v,loads=circuit.dc(c,n,self.out/'dc')
        expected=3.3-.2*(.01+.02+.01)
        self.assertAlmostEqual(loads['load1'],expected,places=9);self.assertAlmostEqual(v['load1_n'],.2*.01,places=9)
        self.record('dc',{'load_v':loads['load1'],'expected_v':expected,'return_rise_v':v['load1_n']})
    def test_ac_series_parallel_rlc(self):
        c,n=self.setup_case();circuit.ac(c,n,self.out/'ac')
        a=np.loadtxt(self.out/'ac/baseline/load1/load1.csv',delimiter=',',skiprows=1);s=2j*np.pi*a[:,0]
        z1=.04+s*20e-9;z2=.03+s*2e-9+1/(s*10e-6);exact=1/(1/z1+1/z2);actual=a[:,1]+1j*a[:,2]
        error=float(np.max(abs(actual-exact)/abs(exact)));self.assertLess(error,.001)
        self.record('ac_rlc',{'maximum_relative_complex_error':error})
    def test_antiresonance_and_limit_failure(self):
        c,n=self.setup_case();cap=dict(c['populations']['baseline'][0]);cap.update(ref='C2',capacitance_f=.1e-6,esr_ohm=.001,esl_h=.3e-9,mounting_inductance_h=.2e-9)
        c['populations']['baseline'][0]['esr_ohm']=.001;c['populations']['baseline'].append(cap)
        c['target_impedance']={'allowed_dynamic_v':1e-6,'delta_current_a':1}
        r=circuit.ac(c,n,self.out/'antiresonance')
        a=np.loadtxt(self.out/'antiresonance/baseline/load1/load1.csv',delimiter=',',skiprows=1);s=2j*np.pi*a[:,0]
        exact=1/(1/(.04+s*20e-9)+1/(.001+s*2e-9+1/(s*10e-6))+1/(.001+s*.5e-9+1/(s*.1e-6)))
        err=float(np.max(abs((a[:,1]+1j*a[:,2])-exact)/abs(exact)))
        self.assertLess(err,.001);self.assertFalse(r['baseline']['load1']['within_target'])
        self.record('antiresonance',{'maximum_relative_error':err,'expected_limit_failure':True})
    def test_rc_transient_independent_ode(self):
        c,n=self.setup_case();c['source']['inductance_h']=0
        c['populations']['baseline'][0].update(esr_ohm=0,esl_h=0,mounting_inductance_h=0)
        t=c['transient'];t.update(step_s=1e-8,stop_s=200e-6)
        lines=circuit.base(c,n,'baseline')+circuit.loads(c,True)
        a=circuit.solve(lines,[f'tran {t["step_s"]} {t["stop_s"]} 0 {t["step_s"]}'],['v(load_load1,load1_n)'],self.out/'rc-transient')
        R=.04;C=10e-6;I=.2;delta=t['delta_a'];edges=[t['delay_s'],t['delay_s']+t['rise_s'],t['delay_s']+t['rise_s']+t['width_s'],t['delay_s']+t['rise_s']+t['width_s']+t['fall_s']]
        def load(time):return I+delta*np.interp(time,[0,*edges,t['stop_s']],[0,0,1,1,0,0])
        ref=solve_ivp(lambda time,y:[((3.3-y[0])/R-load(time))/C],(0,t['stop_s']),[3.3-I*R],rtol=1e-11,atol=1e-12,max_step=t['step_s'],dense_output=True)
        err=float(np.max(abs(a[:,1]-ref.sol(a[:,0])[0]))/ (delta*R))
        self.assertLess(err,.001);self.record('rc_transient',{'maximum_error_relative_to_load_step_droop':err})
    def test_parallel_current_sharing_and_transfer(self):
        # Independently assembled 3-terminal star, unequal return resistances.
        c,n=self.setup_case();c['loads']['load2']={'positive':'load2_p','negative':'load2_n','current_a':.3,'inductance_h':0}
        n['names']+=['load2_p','load2_n'];G=np.zeros((6,6));G[:4,:4]=n['G']
        for i,j,r in [(0,4,.03),(1,5,.02)]:G[i,i]+=1/r;G[j,j]+=1/r;G[i,j]-=1/r;G[j,i]-=1/r
        n['G']=G;c['populations']={'baseline':[]};c['source']['inductance_h']=0
        _,v=circuit.dc(c,n,self.out/'star-dc')
        self.assertAlmostEqual(v['load1'],3.3-.5*.01-.2*.03,places=9)
        self.assertAlmostEqual(v['load2'],3.3-.5*.01-.3*.05,places=9)
        circuit.ac(c,n,self.out/'star-ac');a=np.loadtxt(self.out/'star-ac/baseline/load1/load2.csv',delimiter=',',skiprows=1)
        np.testing.assert_allclose(a[:,1],.01,rtol=.001);np.testing.assert_allclose(a[:,2],0,atol=1e-12)
        self.record('multi_terminal',{'load_voltages_v':v,'transfer_impedance_ohm':float(a[0,1])})
    def test_rlc_transient_independent_ode(self):
        c,n=self.setup_case();c['source']['inductance_h']=20e-9
        c['populations']['baseline'][0].update(esr_ohm=0,esl_h=0,mounting_inductance_h=0)
        t=c['transient'];t.update(step_s=1e-8,stop_s=200e-6)
        a=circuit.solve(circuit.base(c,n,'baseline')+circuit.loads(c,True),
                        [f'tran {t["step_s"]} {t["stop_s"]} 0 {t["step_s"]}'],['v(load_load1,load1_n)'],self.out/'rlc-transient')
        R=.04;L=20e-9;C=10e-6;I=.2;delta=t['delta_a']
        edges=[t['delay_s'],t['delay_s']+t['rise_s'],t['delay_s']+t['rise_s']+t['width_s'],t['delay_s']+t['rise_s']+t['width_s']+t['fall_s']]
        def ode(time,y):
            current=I+delta*np.interp(time,[0,*edges,t['stop_s']],[0,0,1,1,0,0])
            return [(3.3-R*y[0]-y[1])/L,(y[0]-current)/C]
        ref=solve_ivp(ode,(0,t['stop_s']),[I,3.3-I*R],rtol=1e-11,atol=1e-12,max_step=t['step_s'],dense_output=True)
        exact=ref.sol(a[:,0])[1];err=float(np.max(abs(a[:,1]-exact))/np.ptp(exact))
        self.assertLess(err,.001);self.record('rlc_transient',{'maximum_error_relative_to_voltage_excursion':err})
    def test_finite_edge_step_and_final_settling(self):
        c,n=self.setup_case();c['source']['inductance_h']=0
        c['populations']['baseline'][0].update(ref='U1',esr_ohm=0,esl_h=0,mounting_inductance_h=0)
        t=c['transient'];t.update(kind='step',step_s=1e-8,stop_s=200e-6)
        del t['width_s'];del t['fall_s'];c=validate(c)
        result=circuit.transient(c,n,self.out/'step')['baseline']
        expected=3.3-(.2+t['delta_a'])*.04
        self.assertAlmostEqual(result['loads']['load1']['settling_target_v'],expected,places=10)
        self.assertTrue(result['loads']['load1']['within_limits']);self.assertTrue(result['timestep_verification']['load1']['passed'])
        a=np.loadtxt(self.out/'step/baseline/nominal/voltages.csv',delimiter=',',skiprows=1)
        tau=.04*10e-6;x=np.maximum(a[:,0]-t['delay_s'],0);rise=t['rise_s'];d=t['delta_a']
        response=lambda s:np.where(s>0,s-tau*(1-np.exp(-np.maximum(s,0)/tau)),0.)
        exact=3.3-.2*.04-.04*d/rise*(response(x)-response(x-rise))
        error=float(abs(a[:,1]-exact).max()/(d*.04));self.assertLess(error,.001)
        self.record('finite_edge_step',{'maximum_error_relative_to_voltage_excursion':error,'final_voltage_v':expected,'settling_s':result['loads']['load1']['settling_s']})
    def test_solver_failure_is_not_cached_success(self):
        with self.assertRaises(WorkflowError):circuit.solve(['* invalid circuit','NotADevice a b 1'],['op'],['v(a)'],self.out/'expected-failure')
        with self.assertRaises(WorkflowError):circuit.solve(['* inconsistent ideal sources','Vone a 0 1','Vtwo a 0 2'],['op'],['v(a)'],self.out/'singular')

class GeometryFailures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=Path('/output/geometry-checks');cls.out.mkdir(parents=True,exist_ok=True)
        cls.c=validate(json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text()))
        cls.board=extract.extract(ROOT/'examples/synthetic-pdn/board.kicad_pcb',cls.out/'extracted.json',.01)
    def test_missing_and_dnp_components(self):
        for kind in ('missing','dnp'):
            b=copy.deepcopy(self.board)
            if kind=='missing':b['components']=[x for x in b['components'] if x['ref']!='J1']
            else:next(x for x in b['components'] if x['ref']=='J1')['dnp']=True
            with self.subTest(kind=kind),self.assertRaisesRegex(WorkflowError,'Missing or DNP'):geometry.build(b,self.c,self.out/kind)
    def test_duplicate_component_reference(self):
        b=copy.deepcopy(self.board);b['components'].append(copy.deepcopy(b['components'][0]))
        with self.assertRaisesRegex(WorkflowError,'Duplicate component reference'):
            geometry.build(b,self.c,self.out/'duplicate-reference')
    def test_wrong_pad_connectivity(self):
        b=copy.deepcopy(self.board);b['components'][0]['pads'][0]['net']='WRONG'
        with self.assertRaisesRegex(WorkflowError,'connectivity'):geometry.build(b,self.c,self.out/'connectivity')
    def test_missing_fill(self):
        b=copy.deepcopy(self.board);b['warnings'].append('No saved fill for RETURN')
        with self.assertRaisesRegex(WorkflowError,'fill'):geometry.build(b,self.c,self.out/'fill')
    def test_blind_via(self):
        b=copy.deepcopy(self.board);b['holes'][0]['layers']=['F.Cu']
        with self.assertRaisesRegex(WorkflowError,'through-hole'):geometry.build(b,self.c,self.out/'blind')
    def test_short(self):
        b=copy.deepcopy(self.board);x=copy.deepcopy(next(x for x in b['copper'] if x['net']=='SUPPLY'));x['net']='RETURN';b['copper'].append(x)
        with self.assertRaisesRegex(WorkflowError,'short'):geometry.build(b,self.c,self.out/'short')
    def test_disconnected_load(self):
        b=copy.deepcopy(self.board)
        b['copper']=[x for x in b['copper'] if not(x['kind']=='track' and x['net']=='SUPPLY')]
        with self.assertRaisesRegex(WorkflowError,'Disconnected'):geometry.build(b,self.c,self.out/'disconnected')
    def test_electrothermal_nonconvergence_is_failed(self):
        folder=self.out/'nonconvergence';folder.mkdir(exist_ok=True)
        shutil.copyfile(ROOT/'examples/synthetic-pdn/board.kicad_pcb',folder/'board.kicad_pcb')
        c=copy.deepcopy(self.c);c['numerics']['max_iterations']=1;c['mesh']['scale']=2.;c['output']='results'
        write_json(folder/'case.json',c)
        with self.assertRaisesRegex(WorkflowError,'did not converge'):
            pipeline.run(folder/'case.json',folder/'results','electrothermal')
        self.assertEqual(json.loads((folder/'results/provenance.json').read_text())['status'],'failed')
        with self.assertRaises(WorkflowError):pipeline.verify(folder/'case.json',folder/'results')

if __name__=='__main__':unittest.main()
