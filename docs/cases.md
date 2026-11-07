# Case configuration and physical assumptions

Use the supplied `examples/synthetic-pdn/case.json` as the complete schema-1
example. A board is read relative to its case file; `project_root` must contain
both files. `output` must point to a dedicated generated-results directory.

## Electrical interface

Declare distinct supply and return net names. Each named terminal identifies
`ref`, `pad`, `net` and an exterior `layer` (`F.Cu` or `B.Cu`). Names use lowercase
letters, digits and underscores. A pad can have only one terminal name; reuse
that name in the circuit. Missing, DNP, mismatched or overlapping contacts fail.
Each referenced footprint must have a unique component reference on the board;
duplicate references are rejected before meshing.
Names reserved by the generated source, load and capacitor circuit cannot be
reused as terminal names.

The electrical mesh preserves saved tracks, pads, filled polygons and through
barrels by net and layer. Conductivity applies to finite copper volumes, including
annular plating; bores are not solid copper. Blind/buried vias and holes crossing
the outline are rejected. Missing saved fills are rejected. Every selected net
must connect all its declared terminals. Unselected copper conducts heat only.
Drill dimensions are treated as finished bore dimensions; the declared plating
grows outward from that bore. Confirm this interpretation for an imported design.

A terminal is an equipotential surface contact. It does not short the whole
net. Separate terminal excitations generate a reciprocal conductance matrix.
The circuit contains that matrix as voltage-controlled current sources; source
and load nodes remain distinct on both supply and return. The independent
tetrahedral stiffness/reaction audit checks the Elmer solution and its energy.
This circuit representation is resistive and linear at the current temperature.

## Circuit inputs

The source declares differential `positive`/`negative` terminals, `voltage_v`,
positive `resistance_ohm` and nonnegative `inductance_h`. It represents an ideal
voltage source with passive series impedance, not a regulator control-loop model.

Each load declares its two terminals, `current_a`, and total series package/lead
`inductance_h`. Current flows from supply to return. This release uses constant
current loads plus the declared pulse; it does not silently infer constant-power
device behavior or heat from load power.

Each capacitor population is an explicit list of populated references, mapped
terminals, effective `capacitance_f`, `esr_ohm`, `esl_h`, and
`mounting_inductance_h`. The capacitor's two terminals must belong to its declared
reference. Supply and package inductances include only the portions the user
assigns to them. ESR excludes board resistance already extracted by FEM. No AC
parasitic is inferred from a footprint name, nominal capacitance or placement.
Declared pad ports remain equipotential contact regions in every population,
including when a capacitor is omitted from a population.

At DC capacitors are open and inductors short. AC uses a one-ampere small-signal
current injection at each load; the saved complex differential voltages are the
driving-point and transfer impedances in ohms. Source and other load perturbations
are zero. Results cover the declared circuit-model frequency range only.

`target_impedance` accepts either `allowed_dynamic_v` and `delta_current_a` for
a flat reference, or an ordered `envelope` of `[frequency_hz, impedance_ohm]`
pairs. Envelope interpolation is linear in impedance versus log frequency.
The reference is a design criterion, not a guaranteed time-domain voltage bound.
Envelope checks apply at the sampled frequencies. Increase sweep density around
sharp resonances; a finite sweep does not guarantee an unsampled peak stays below
the envelope.

Transients specify the excited load, positive current increment, delay, rise,
plateau width, fall, stop and maximum step, all in seconds/amperes. The DC
operating point initializes the run. The circuit-model bandwidth must cover at
least the inverse edge time, and the step must be at most one fifth of each edge.
Settling is measured after the falling edge relative to the initial voltage;
all remaining samples must remain within the declared band. An insufficient
recovery window or absent settling cannot pass the settling limit.

Set `transient.kind` to `pulse` (the default) for a finite pulse, or `step` for
a sustained load increase. A step omits `width_s` and `fall_s`. Its settling band
is centered on a separately solved DC operating point at the increased current,
and its settling time starts at the end of the rising edge. Pulse settling starts
at the end of the falling edge and targets the initial voltage. Overshoot is the
largest excursion above the initial voltage over the entire waveform.

## Thermal interface

Copper requires reference resistivity, reference temperature, temperature
coefficient and maximum model temperature. Conductivity is evaluated separately
for each tetrahedron using its mean nodal temperature. The full elementwise
coefficient array is passed unchanged through the compiled PiFields callback;
it is not quantized into a few temperature classes.

All materials require three diagonal thermal conductivities and density.
`cooling.kind` selects uniform exterior convection or a prescribed temperature
on the exterior. Convection requires `ambient_c` and `h_w_m2k`; prescribed
temperature requires `temperature_c`. Contact thermal resistance, radiation,
fluid flow, IC package bodies and thermal transients are outside this release.
All dielectric layers share the declared dielectric properties. Saved mask layers
are uniform sheets; mask apertures are not extracted. Non-copper portions of a
copper layer use the dielectric material. Declare these approximations among the
case assumptions when they matter to the intended thermal result.

Copper losses come from the reconstructed total potential field, so superposed
terminal fields retain their cross terms. The electrical/thermal grids share
the same tetrahedra. An optional `extra_heat` list maps explicit watts to named
terminal contact regions; heat is spread over adjacent tetrahedral volumes.
These inputs must describe actual deposited heat, not the electrical power of
an external load. Energy audits include all declared auxiliary heat.

Elmer's heat source is W/kg, so the implementation converts the W/m3 loss field
using each material's density. Reports and exported fields retain W/m3 and watts.

## Mesh and resistance metrics

The saved copper outlines and terminal boundaries are imprinted onto a common
Gmsh triangulation, then swept through the physical stack. Each prism is split
conformally into three tetrahedra. `mesh.scale` multiplies the planar size and
also sets the subdivisions per physical layer to `max(1, round(2 / scale))`.
The first three refinement levels therefore use two, three and four subdivisions
through each layer; additional levels use five, eight, eleven and sixteen. Polygon chord approximation is controlled separately by
`mesh.chord_mm`; changing a mesh size does not refine the original polygon arcs.

DC summaries include the operating-point loss-equivalent resistance
`copper_loss / total_load_current**2`. This includes supply and return copper and
depends on the load sharing. It is null for zero total current. The terminal
conductance matrix is the reusable multiport representation; do not interpret
the loss-equivalent scalar as a separate resistor to add to that matrix.

## Failure and validation behavior

Nonfinite or invalid parameters, shorts, disconnected terminals, missing source
data, solver errors, circuit singularities, material-range violations, energy
imbalance and electrothermal nonconvergence fail explicitly. A failing design
limit remains a valid numerical result and is reported separately.

Transient verification reruns at half the maximum timestep. Droop, overshoot
and settling use 1% relative comparison with absolute floors of 1 nV and 1 ps.
The final settling-band crossing is interpolated between samples rather than
rounded to the next timestep. Near-zero quantities must not
cause division by zero or disguise an unsettled response as a pass.

Electrical conservation uses 0.5% relative tolerance with 50 pW/50 pA absolute
floors. Thermal conservation uses 1% with a 100 pW absolute floor. Reported
normalized imbalances use those same floors at zero or very small loads.

Input and implementation changes invalidate cached stages. Artifact corruption
also invalidates them. Results are tied to the actual runtime image identity;
changing images requires revalidation. Provenance does not establish measured
component accuracy, fabrication tolerances or hardware correlation.
