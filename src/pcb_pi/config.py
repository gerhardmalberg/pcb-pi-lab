"""Validate explicit electrical interfaces before any solver is launched."""
import copy
import json
import re
from pathlib import Path
from .common import WorkflowError, number

def validate(case):
    c = copy.deepcopy(case)
    if c.get('schema') != 1:
        raise WorkflowError('case.schema must be 1')
    for key in ('board', 'project_root', 'output'):
        if not isinstance(c.get(key), str) or not c[key]:
            raise WorkflowError(f'{key} must be an explicit path')
    nets = c.get('nets', {})
    if set(nets) != {'supply', 'return'} or any(not isinstance(x, str) or not x for x in nets.values()) or nets['supply'] == nets['return']:
        raise WorkflowError('Declare distinct nets.supply and nets.return')
    terminals = c.get('terminals', {})
    if not isinstance(terminals, dict) or len(terminals) < 4:
        raise WorkflowError('At least four explicitly mapped terminals are required')
    pads = set()
    for name, t in terminals.items():
        if not re.fullmatch(r'[a-z][a-z0-9_]*', name):
            raise WorkflowError('Terminal names must use lowercase circuit identifiers')
        for k in ('ref', 'pad', 'net', 'layer'):
            if not isinstance(t.get(k), str) or not t[k]:
                raise WorkflowError(f'terminals.{name}.{k} is required')
        if t['net'] not in nets.values() or t['layer'] not in ('F.Cu', 'B.Cu'):
            raise WorkflowError(f'Invalid terminal net or exterior layer: {name}')
        pad = (t['ref'], t['pad'])
        if pad in pads:
            raise WorkflowError('Duplicate terminal pad; reuse its circuit node instead')
        pads.add(pad)
    def pair(item, name):
        for k, net in [('positive', nets['supply']), ('negative', nets['return'])]:
            if item.get(k) not in terminals or terminals[item[k]]['net'] != net:
                raise WorkflowError(f'{name}.{k} must name a terminal on {net}')
    s = c['source']; pair(s, 'source')
    for k in ('voltage_v', 'resistance_ohm'):
        number(s[k], 'source.' + k)
    number(s['inductance_h'], 'source.inductance_h', allow_zero=True)
    loads = c.get('loads', {})
    if not loads:
        raise WorkflowError('At least one load is required')
    for name, load in loads.items():
        if not re.fullmatch(r'[a-z][a-z0-9_]*', name): raise WorkflowError('Invalid load name')
        pair(load, name)
        for k in ('current_a', 'inductance_h'):
            number(load[k], name + '.' + k, allow_zero=True)
    populations = c.get('populations', {})
    if not populations or c.get('baseline_population') not in populations:
        raise WorkflowError('Declare capacitor populations and baseline_population')
    for name, population in populations.items():
        if not re.fullmatch(r'[a-z][a-z0-9_]*', name): raise WorkflowError('Invalid population name')
        seen = set()
        for cap in population:
            pair(cap, 'capacitor')
            if cap['ref'] in seen: raise WorkflowError('Duplicate capacitor reference')
            seen.add(cap['ref'])
            if any(terminals[cap[k]]['ref'] != cap['ref'] for k in ('positive','negative')):
                raise WorkflowError('Capacitor terminals must belong to its reference')
            for k in ('capacitance_f', 'esr_ohm', 'esl_h', 'mounting_inductance_h'):
                number(cap[k], 'capacitor.' + k, allow_zero=k != 'capacitance_f')
    internal={'supply_ideal','supply_series'} | {'load_'+name for name in loads}
    internal |= {prefix+str(i) for i in range(max(map(len,populations.values()))) for prefix in ('cap_l','cap_c')}
    if internal.intersection(terminals): raise WorkflowError('Terminal name collides with a reserved internal circuit node')
    number(c['geometry']['plating_mm'], 'plating_mm')
    for k in ('size_mm', 'chord_mm', 'scale'):
        number(c['mesh'][k], 'mesh.' + k)
    for name in ('copper', 'dielectric', 'soldermask'):
        m = c['materials'][name]
        if len(m['k_w_mk']) != 3: raise WorkflowError('Thermal conductivity requires three diagonal components')
        for k in m['k_w_mk']: number(k, name + '.k_w_mk')
        number(m['density_kg_m3'], name + '.density')
    copper = c['materials']['copper']
    number(copper['resistivity_ohm_m'], 'resistivity')
    number(copper['temperature_coefficient_per_k'], 'temperature coefficient', allow_zero=True)
    number(copper['reference_temperature_c'], 'reference temperature', -273.15)
    number(copper['maximum_temperature_c'], 'maximum temperature', copper['reference_temperature_c'])
    cooling = c['cooling']
    number(cooling['ambient_c'], 'ambient', -273.15)
    if cooling['kind'] == 'convection': number(cooling['h_w_m2k'], 'convection coefficient')
    elif cooling['kind'] == 'fixed': number(cooling['temperature_c'], 'boundary temperature', -273.15)
    else: raise WorkflowError('cooling.kind must be convection or fixed')
    for heat in c.get('extra_heat', []):
        if heat['terminal'] not in terminals: raise WorkflowError('Extra heat needs a named terminal')
        number(heat['watts'], 'extra heat', allow_zero=True)
    n = c.setdefault('numerics', {})
    for k, v in {'relaxation':.5, 'max_iterations':50, 'temperature_tolerance_c':.1,
                 'relative_tolerance':.001, 'solver_timeout_s':120}.items(): n.setdefault(k, v)
    for k, v in n.items(): number(v, 'numerics.' + k)
    if n['relaxation'] > 1 or int(n['max_iterations']) != n['max_iterations']:
        raise WorkflowError('Invalid electrothermal iteration settings')
    ac = c['ac']
    number(ac['start_hz'], 'ac.start_hz'); number(ac['stop_hz'], 'ac.stop_hz', ac['start_hz'])
    number(ac['points_per_decade'], 'ac.points_per_decade')
    if int(ac['points_per_decade']) != ac['points_per_decade']: raise WorkflowError('AC point count must be integer')
    valid = c['model_frequency_range_hz']
    if len(valid) == 2:
        number(valid[0], 'minimum model frequency'); number(valid[1], 'maximum model frequency', valid[0])
    if len(valid) != 2 or valid[0] > ac['start_hz'] or valid[1] < ac['stop_hz']:
        raise WorkflowError('AC sweep exceeds the declared circuit-model frequency range')
    t = c['transient']
    t.setdefault('kind','pulse')
    if t['kind'] not in ('pulse','step'):raise WorkflowError('transient.kind must be pulse or step')
    extra=('width_s','fall_s') if t['kind']=='pulse' else ()
    for k in ('step_s', 'stop_s', 'delay_s', 'rise_s')+extra:
        number(t[k], 'transient.' + k)
    if t['stop_s'] <= sum(t[k] for k in ('delay_s','rise_s')+extra):
        raise WorkflowError('Transient must include a settling interval after the last edge')
    edge=min(t['rise_s'],t['fall_s']) if t['kind']=='pulse' else t['rise_s']
    if t['step_s'] > edge / 5:
        raise WorkflowError('Transient step must resolve each edge with at least five intervals')
    if valid[1] < 1 / edge:
        raise WorkflowError('Circuit-model bandwidth is insufficient for the requested load edges')
    if t['load'] not in loads: raise WorkflowError('Unknown transient load')
    number(t['delta_a'], 'load change')
    limits = c['limits']
    for k in ('minimum_load_v','maximum_load_v','maximum_droop_v','maximum_overshoot_v','settling_band_v','settling_s'):
        number(limits[k], 'limits.' + k)
    if limits['maximum_load_v'] <= limits['minimum_load_v']: raise WorkflowError('Invalid load-voltage limits')
    target = c['target_impedance']
    if 'envelope' in target:
        points = target['envelope']
        if len(points) < 2: raise WorkflowError('Target envelope requires at least two [Hz, ohm] points')
        for f, z in points: number(f,'target frequency'); number(z,'target impedance')
        if any(a[0] >= b[0] for a,b in zip(points,points[1:])) or points[0][0] > ac['start_hz'] or points[-1][0] < ac['stop_hz']:
            raise WorkflowError('Target envelope must be ordered and cover the AC sweep')
    else:
        number(target['allowed_dynamic_v'],'allowed dynamic voltage')
        number(target['delta_current_a'],'target current change')
    return c

def load(path):
    path = Path(path).resolve()
    c = validate(json.loads(path.read_text()))
    project = (path.parent / c['project_root']).resolve()
    board = (path.parent / c['board']).resolve()
    out = (path.parent / c['output']).resolve()
    if not board.is_relative_to(project) or not board.is_file():
        raise WorkflowError('Board must be an existing file inside project_root')
    if out == project or project.is_relative_to(out) or out == board.parent:
        raise WorkflowError('Output must be a dedicated results directory')
    if not path.is_relative_to(project): raise WorkflowError('Case must be inside project_root')
    return c, project, board, out
