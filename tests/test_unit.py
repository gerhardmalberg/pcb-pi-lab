import copy
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from pcb_pi.common import ROOT,WorkflowError
from pcb_pi.config import validate
from pcb_pi.mesh import metrics
from pcb_pi.fem import matrix
from pcb_pi.pipeline import seal,cached,conductivity
from pcb_pi.circuit import transient_metrics,target

def example():return json.loads((ROOT/'examples/synthetic-pdn/case.json').read_text())

class ConfigTests(unittest.TestCase):
    def test_explicit_case(self):self.assertEqual(validate(example())['numerics']['relaxation'],.5)
    def test_missing_parameters(self):
        for group,key in [('source','inductance_h'),('geometry','plating_mm'),('ac','stop_hz'),('materials','copper')]:
            c=example();del c[group][key]
            with self.subTest(group=group,key=key),self.assertRaises((KeyError,WorkflowError)):validate(c)
    def test_duplicate_terminal(self):
        c=example();c['terminals']['another']=dict(c['terminals']['src_p'])
        with self.assertRaisesRegex(WorkflowError,'Duplicate'):validate(c)
    def test_wrong_polarity(self):
        c=example();c['source']['positive']='src_n'
        with self.assertRaisesRegex(WorkflowError,'SUPPLY'):validate(c)
    def test_invalid_electrical_values(self):
        for v in (-1,0,float('nan'),float('inf'),True):
            c=example();c['source']['resistance_ohm']=v
            with self.subTest(v=v),self.assertRaises(WorkflowError):validate(c)
    def test_invalid_capacitors(self):
        for key in ('capacitance_f','esr_ohm','esl_h','mounting_inductance_h'):
            c=example();c['populations']['baseline'][0][key]=-1
            with self.subTest(key=key),self.assertRaises(WorkflowError):validate(c)
    def test_internal_node_collision(self):
        c=example();c['terminals']['supply_ideal']=c['terminals'].pop('src_p');c['source']['positive']='supply_ideal'
        with self.assertRaisesRegex(WorkflowError,'reserved'):validate(c)
    def test_insufficient_bandwidth(self):
        c=example();c['model_frequency_range_hz']=[1,1e5]
        with self.assertRaisesRegex(WorkflowError,'frequency range'):validate(c)
    def test_unresolved_edge(self):
        c=example();c['transient']['step_s']=1e-4
        with self.assertRaisesRegex(WorkflowError,'resolve'):validate(c)
    def test_invalid_envelope(self):
        c=example();c['target_impedance']={'envelope':[[10,.1],[1,.2]]}
        with self.assertRaisesRegex(WorkflowError,'ordered'):validate(c)
    def test_invalid_thermal_iteration(self):
        c=example();c['numerics']['relaxation']=2
        with self.assertRaisesRegex(WorkflowError,'iteration'):validate(c)
    def test_extra_heat_needs_location(self):
        c=example();c['extra_heat']=[{'terminal':'unknown','watts':1}]
        with self.assertRaisesRegex(WorkflowError,'named terminal'):validate(c)
    def test_temperature_model_rejects_invalid(self):
        mesh={'tets':np.array([[0,1,2,3]]),'electric_full':np.array([0])}
        with self.assertRaisesRegex(WorkflowError,'range'):conductivity(example(),mesh,np.full(4,200.))
    def test_uniform_temperature_resistance(self):
        mesh={'tets':np.array([[0,1,2,3]]),'electric_full':np.array([0])};c=example()
        cold=conductivity(c,mesh,np.full(4,20.));hot=conductivity(c,mesh,np.full(4,70.))
        self.assertAlmostEqual(float(cold[0]/hot[0]),1+.00393*50)

class NumericalTests(unittest.TestCase):
    def test_affine_tetrahedron_gradient(self):
        x=np.array([[2,3,4],[4,3,4],[3,6,4],[2,4,8]],float)*1e-3;t=np.array([[0,1,2,3]])
        v,g=metrics(x,t);u=x@np.array([2.,-3.,4.])+7
        np.testing.assert_allclose((u[t]@g[0])[0],[2,-3,4],atol=1e-10)
        self.assertAlmostEqual(v[0],4e-9)
    def test_stiffness_gauge_and_energy(self):
        x=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float);t=np.array([[0,1,2,3]])
        K,v,g=matrix(x,t,np.array([2.]))
        np.testing.assert_allclose(K@np.ones(4),0,atol=1e-14)
        self.assertAlmostEqual(float(x[:,0]@K@x[:,0]),1/3)
    def test_degenerate_mesh(self):
        with self.assertRaisesRegex(WorkflowError,'Degenerate'):metrics(np.zeros((4,3)),np.array([[0,1,2,3]]))
    def test_flat_and_piecewise_target(self):
        c=example();np.testing.assert_allclose(target(c,np.array([10.,100.])),.25)
        c['target_impedance']={'envelope':[[10,.1],[1000,.3]]}
        np.testing.assert_allclose(target(c,np.array([10.,100.,1000.])),[.1,.2,.3])
    def test_no_settling_is_not_passed(self):
        c=validate(example());t=np.linspace(0,c['transient']['stop_s'],1000);v=np.where(t<c['transient']['delay_s'],3.3,3.0)
        r=transient_metrics(t,v,c['transient'],c['limits']);self.assertIsNone(r['settling_s']);self.assertFalse(r['within_limits'])
    def test_recovery_window(self):
        c=validate(example());t=np.linspace(0,c['transient']['stop_s'],1000);v=np.full_like(t,3.3)
        c['limits']['settling_s']=1.
        r=transient_metrics(t,v,c['transient'],c['limits']);self.assertFalse(r['within_limits'])
    def test_settling_interpolates_band_crossing(self):
        c=validate(example());t=np.linspace(0,c['transient']['stop_s'],251)
        end=sum(c['transient'][k] for k in ('delay_s','rise_s','width_s','fall_s'));tau=20e-6
        v=3.3-np.where(t<c['transient']['delay_s'],0.,.01*np.exp(-np.maximum(t-end,0)/tau))
        r=transient_metrics(t,v,c['transient'],c['limits']);expected=tau*np.log(.01/c['limits']['settling_band_v'])
        self.assertLess(abs(r['settling_s']/expected-1),.001)
    def test_cache_detects_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'result').write_text('value');seal(p,'a');self.assertTrue(cached(p,'a'));self.assertFalse(cached(p,'b'))
            (p/'result').write_text('changed');self.assertFalse(cached(p,'a'))

if __name__=='__main__':unittest.main()
