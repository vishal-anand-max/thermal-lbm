# ZSC 3D Fixed - Build & Run Instructions

## Quick Start

```bash
# Compile
gfortran -O2 -o zsc3d zsc3d_fixed.f90

# Run Laplace test (quick, ~5 min)
cd laplace_test
cp ../zsc3d .
./zsc3d

# Run straight channel (longer, hours)
cd straight_channel
cp ../zsc3d .
./zsc3d
```

## What Was Fixed (vs Arup's original)

1. **Wall BCs** — Original had sequential overwrite bug (`g(1)=g(2); g(2)=g(1)` makes both equal).
   Fixed with proper bounce-back using `opp_g()` array.

2. **Inlet BC** — Zou-He with parabolic velocity profile (from Sudhakar's paper).
   For buoyancy-only rise, set `u_inlet_max = 0.0`.

3. **Outlet BC** — Linear extrapolation from interior (`2*f(N-1) - f(N-2)`).

4. **Geometry** — Parameterized via `is_fluid` mask instead of hardcoded ±50 offsets.

5. **Neighbor array** — Properly clamped at boundaries (was commented out).

## Validation Targets

### Laplace Test (mode 3)
- Static drop in periodic box, no gravity
- Measure P_in - P_out from stats/VTK
- Should match: ΔP = 2σ/R (3D sphere)
- With σ=0.001, R=15: expect ΔP ≈ 0.000133

### Straight Channel (mode 1) — Sudhakar Fig 1f
- Kerosene drop (ρ=787) rising in water (ρ=1000)
- Track drop shape, rise velocity, film thickness
- Compare against Sudhakar & Das, I&EC Research 2020, Figure 1f

## Known Issues / Things to Watch

1. **Gravity direction**: For straight channel rising drop, gravity should oppose
   the drop motion. With `sinthita=1, costhita=0`, gravity acts in +z on the
   light phase (drop). If drop should rise in +z, gravity helps it — this may
   need sign adjustment depending on your physical setup.

2. **The code is NOT yet tested** — it was written from analysis of the original
   buggy code + Sudhakar's published methodology. The physics core (collision)
   is identical to Arup's original. Only BCs and geometry handling changed.

3. **Periodic BC for Laplace test** may need debugging — the array slice syntax
   for periodic copy may need explicit loops if your compiler complains.

4. **properties.in format** has 2 new lines vs Arup's original:
   - `geom_mode` (1/2/3)
   - `u_inlet_max` (Zou-He inlet velocity)

5. **Performance**: This is single-threaded serial Fortran. For the 52×52×302
   channel at 100K steps, expect several hours on a modern CPU. The Laplace
   64³ test should finish in minutes.

## File Structure

```
zsc3d_fixed.f90          # Main code (all in one file)
straight_channel/
  properties.in          # Sudhakar validation case
  discrete.in            # Single drop
laplace_test/
  properties.in          # Quick Laplace test
  discrete.in            # Static drop in periodic box
```

## Output Files

- `stats.out` — Per-timestep: step, ux, uy, uz, volume_ratio, effective_radius
- `DATA_XXXXXXXX.vtk` — Paraview-compatible VTK with phi, pressure, velocity
