# Numerical validation

Validation uses synthetic geometry and independently specified circuits only.
It establishes agreement with the stated models, not correlation with a physical
PCB, component datasheet, current rating or regulator control loop.

The verified tool identities are in
[validated-runtime.json](../env/locks/validated-runtime.json). The installed-image
build, native static-current module compilation/loading and CPU numerical runs
were exercised. The complete fresh-machine image build recipe was not executed.
KiCad's installed bindings emit nonfatal property-enumeration assertions during
import; exit status, extraction integrity and the numerical checks are assessed
separately, and the original logs are retained.

## Independent benchmarks

Run `./pcb-pi test --integration`. Native outputs and numerical summaries are
written to `results/validation/`. The suite includes:

- A uniform trace against `rho * length / area` and a plated annulus against
  cylindrical resistance; potential reconstruction and electrical power balance.
- Unequal parallel conductors with a 2:1 expected current split, a prescribed
  uniform-temperature resistance, zero electrical load and zero heat input.
- A fixed-conductivity heated bar, explicit spatial auxiliary heat, and an
  insulated-side conductor with fixed end temperatures and a closed-form
  temperature-dependent-resistivity solution.
- Analytic complex RLC impedance, capacitor antiresonance, a deliberately violated
  target envelope, a multi-terminal circuit, and independent transfer impedance.
- RC and RLC pulse responses against independent ODE integrations, plus a
  finite-edge sustained step against the closed-form RC ramp response.
- Invalid and missing parameters, missing/DNP pads, incorrect connectivity,
  shorts, disconnected loads, missing fills, unsupported vias, singular circuits,
  nonconvergence, insufficient settling and cache corruption.

The release suite includes 45 tests. Representative independent results:

| Benchmark | Result | Relative error |
|---|---:|---:|
| Uniform trace | 1.724 mOhm | < 0.00000001% |
| Plated barrel | 0.516855 mOhm | 0.07143% |
| Parallel current sharing | 2:1 | < 0.00000001% |
| Fixed-conductivity coupled bar | 32.5 C | < 0.00000001% in temperature rise |
| Temperature-dependent heated conductor | 33.9462 C | 0.02318% in temperature rise |
| AC RLC and antiresonance | Complex impedance | < 0.000000001% |
| RC finite-edge pulse/step | Voltage response | 0.000706% of excursion |
| RLC finite-edge pulse | Voltage response | 0.001641% of excursion |

Initial resistance tolerance is 1%, current/electrical energy imbalance 0.5%,
thermal energy imbalance 1%, and analytic circuit-response error 0.1%.
Transient comparisons normalize by the voltage excursion rather than the supply
voltage; absolute tolerances handle near-zero values. Reported DC current density
at a sharp corner is mesh-sensitive and is not an automatic current-rating limit.

## Integrated synthetic PCB

```bash
./pcb-pi run examples/synthetic-pdn/case.json --analysis all
./pcb-pi verify examples/synthetic-pdn/case.json
python3 scripts/verify-convergence.py --scales 1 .25 .1875
```

The case includes supply and return paths, plated through vias, two loads and
baseline/alternative capacitor populations. All analyses keep the board file
unchanged. Native Elmer results, ngspice decks, CSV metrics, plots and a ParaView
state remain in the ignored result directory.

At least three meshes refine both planar size and through-thickness resolution. The final
refinement must change load voltage drop, copper loss and maximum temperature
rise by at most 2%. Polygon chord error is held fixed, so this study does not
establish convergence of every possible curved-geometry approximation.
Without `--scales`, the script tries relative scales 1, 0.75, 0.5, 0.375, 0.25, 0.1875 and 0.125, stopping after at
least three meshes once the final comparison passes. Failure at the last scale
remains an explicit failed verification. It runs all four analyses on the accepted
mesh. `results/refinement/summary.json` records its output directory; use that
directory's `report.md` for the converged example.
The release check uses the explicit scales shown above, which bracket the final
converged refinement found during the adaptive study.
The quick-start result under `results/examples/synthetic-pdn/` is the initial
coarse run, and is not by itself evidence of mesh convergence.

Each transient population is also solved with half the maximum electrical
timestep. Droop/overshoot must agree within 1% or 1 nV; settling must agree within
1% or 1 ps. Band crossing is interpolated between samples. A waveform that has not settled cannot
claim a converged settling time.

`verify` validates input, implementation, runtime and artifact identities; it
does not silently substitute hash checks for mesh convergence. The reports keep
solver completion, numerical verification and design-limit compliance separate.
A target-envelope failure is tested as a valid solve with a failed design limit.

Compact results from the completed release validation are recorded in
[evidence.json](validation/evidence.json). Full native results and logs are
reproducible through the commands above and are intentionally not committed.

The accepted scale-0.1875 mesh has 198,322 nodes and 1,066,824 tetrahedra. All four
analyses passed their numerical checks and the synthetic case's design limits.
At 25 C ambient and 2 A total load, the observed results are:

| Quantity | Cold DC | Converged electrothermal |
|---|---:|---:|
| Load 1 differential voltage | 3.194837 V | 3.174478 V |
| Load 2 differential voltage | 3.203354 V | 3.184081 V |
| Supply and return copper loss | 0.163513 W | 0.203361 W |
| Maximum temperature | 25 C imposed | 106.234 C solved |

The electrothermal solve converged in 15 iterations, with a final temperature
residual of 0.058915 C. Thermal energy imbalance was 0.0000126%; electrical
current and power imbalances were below 0.000001%.
The final mesh refinement changed temperature rise by 1.5866%, copper loss by
0.3500%, and load voltage drops by at most 0.3565%. Halving the electrical timestep
changed the transient droop, overshoot and settling metrics by at most 0.00794%
across both capacitor populations and both loads.

A separate temporary checkout passed 22 unit tests and a fresh synthetic DC run
with only its own source, inputs and outputs mounted. Its input and implementation
identity matched the main checkout. Artifact verification rejected a deliberately
altered report and passed again after the original bytes were restored.

## Board presentation update

The example board was subsequently centered on landscape A4, with its auxiliary
and grid origins translated together and its reference labels separated. Hidden
fabrication values retain their original component values. The generator now
reproduces this presentation. KiCad 10.0.8 reports zero DRC violations and zero
unconnected items.

[board-presentation.json](validation/board-presentation.json) records the original
and updated board hashes and the geometry comparison. The full extracted model
before and after the edit is identical after reversing the 138.5 mm X / 99 mm Y
translation. The original validated model also matches; the record identifies an
earlier KiCad save's change to an unused dielectric-loss parameter.
Silkscreen is not part of the modeled physics, and the case uses pad-mapped ports
and uniform exterior cooling, so this rigid translation preserves its physical
problem. No numerical simulations were rerun for these presentation edits.

The numerical results above and `evidence.json` retain their original input and
generator hashes. Saved result plots remain in the original coordinate system.
The CLI hashes the whole board file: `verify` therefore reports changed inputs,
and a subsequent `run` recomputes the stages. This separate equivalence record
does not change the cache rules or relabel the historical results as a fresh run.

## Scope of the evidence

The supported physics is DC copper conduction, declared linear circuit models
and steady electrothermal feedback. AC/transient analyses hold the converged
temperature fixed. The evidence does not validate skin effect, layout-extracted
inductance, time-resolved thermal coupling, CUDA, radiation, fluid cooling,
package junction temperature or measured hardware behavior. Qucs-S is optional
for inspecting the generated decks and was not part of the numerical checks.
