# Physics conventions

MDMTrace models the active MDM spectrometer deck. It is useful independently for field tuning, legacy tracing, map generation, standalone map tracing, comparison, and ion-optics fitting.

## Verified numerical conventions

- nominal dipole sector angle: **100 degrees**;
- dipole reference radius: **160 cm**;
- RAYTRACE speed-of-light convention: `c = 3.0e10 cm/s`;
- field-map positions: cm;
- field-map magnetic fields: Tesla;
- multipole RK4 integration step: 0.1 cm;
- dipole RK4 integration step: 0.1 cm;
- plane-crossing tolerance: `1.0e-8 cm`;
- maximum magnetic integration steps per region: 200000;
- plane-crossing refinement steps: 60;
- relative map-probe comparison tolerance: `2.0e-6` times `max(1, |a|, |b|)`.

The 160 value is a radius in centimetres, not an angle. The sector angle is 100 degrees.

## Active beamline

The standalone field-map transport follows the active sequence:

```text
63.5 cm drift
entrance collimator
18.075 cm drift
entrance multipole
dipole entrance fringe, sector, and exit fringe
36.7 cm zero-field second-multipole drift
exit collimator
96.13 cm final drift
```

The current aperture values and zero-field drifts are copied from `dat/rayin.dat` into the transport source. The second multipole remains disabled: it contributes no magnetic field and has no generated map.

## Magnet settings

Probe mode uses:

```text
B_dipole [T] = mdmDipoleProbe [Gauss] * 1.034e-4
```

Field mode derives the map settings as:

```text
dipoleProbe = mdmDipoleField / 1.034
multipoleProbe = dipoleProbe * 0.71
```

The multipole strengths use the existing calibration ratios in `MdmTrace`. Loaded maps are never rescaled silently. All four maps must contain mutually consistent probe metadata, and the requested probes must match within the source tolerance above.

## Transport

The legacy reference is MIT RAYTRACE through `MdmTrace`. The map generator links that Fortran-backed library and samples its magnetic-field routines directly.

The standalone map reader and map transport are C++ and do not call the Fortran transport. They use straight propagation in zero-field drifts, aperture stops at the copied collimators, and fixed-step RK4 integration through trilinearly interpolated map fields.

For kinetic energy `T` and transported rest mass `m`:

```text
p = sqrt((2*m + T) * T)
v/c = p / (m + T)
```

The transported mass is the AME2020 neutral atomic mass minus `Q` electron masses:

```text
ionMassMeV = neutralAtomicMassU * 931.49410242 - Q * 0.510998950
```

The transport is magnetic only; it has no energy-loss model.

## Map coordinates

The multipole map is centred on the magnet. Local `+z` follows the beam, while `x` and `y` are transverse.

The dipole is split into entrance-fringe, sector, and exit-fringe maps. The entrance-centred physical convention is:

```text
+z: incoming beam
+x: left in top view
+y: upward
central trajectory: bends toward -x
```

The sector map stores `(dr, y, s)` in cm, where `dr = r - 160 cm` and `s = 160 cm * theta`. Stored vector components remain in the dipole-local Cartesian frame.

## Ion-optics coordinates

The phase-space vector is exactly:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

Momentum deviation is a percentage:

```text
deltaP/P0 [%] = 100 * (p - p0) / p0
```

The longitudinal output coordinate is calculated from time of flight:

```text
L [mm] = -v0 * (tof - tof0) / gamma0^2 * 10
```

Later arrival therefore gives negative `L`. Input `L` is a formal coordinate, and the static transport keeps the `L -> L` identity term equal to one.
