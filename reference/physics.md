# Physics Conventions

This project models the current MDM spectrometer deck, not a generic RAYTRACE beamline. The goal is to reproduce the active MDM transport and export magnetic field maps for Geant4 use.

## Beamline Sequence

The implemented validation sequence is:

```text
DRIF 63.5
entrance COLL
DRIF 18.075
entrance multipole
dipole
exit COLL
DRIF 96.13
```

The numerical drift lengths, apertures, fringe parameters, and magnet geometry are read from or copied from the active `dat/rayin.dat` deck. The second multipole is not treated as an active field element.

## Magnet Settings

Hall-probe mode:

```text
B_dipole[T] = mdmDipoleProbe * 1.034e-4
```

Field mode:

```text
dipoleProbe = mdmDipoleField / 1.034
multipoleProbe = dipoleProbe * 0.71
```

The multipole strengths are then set from the same calibration ratios used by `MdmTrace`.

All tools follow the same `usingProbe` convention:
- `usingProbe=true`: use `mdmDipoleProbe` and `mdmMultipoleProbe`
- `usingProbe=false`: use `mdmDipoleField` and derive the probe values

Field-map tracers reject maps if the requested magnetic setting does not match the metadata stored in the map headers.

## Transport Model

The legacy reference is MIT RAYTRACE through `MdmTrace`.

The field-map tracer uses the same visible beamline logic but does not call Fortran:
- straight-line transport in drifts
- aperture stops at collimators
- RK4 magnetic integration through generated field maps
- stopped rays receive the same sentinel-style output used by the legacy examples

Particle speed is computed relativistically from kinetic energy:

```text
m = ionMassMeV
v/c = sqrt((2*m + T) * T) / (m + T)
```

For `scatteredIon`, `ionMassMeV` is the AME2020 neutral atomic mass minus `chargeState` electron masses:

```text
ionMassMeV = neutralAtomicMassU * 931.49410242 - chargeState * 0.510998950
```

The charge state is also the magnetic charge used in the Lorentz force.

The transport is magnetic only. There is no energy loss model.

## Coordinates

### Multipole Map

The multipole map uses a centered local frame:

```text
origin: magnet center
z: beam axis
x,y: transverse axes
beam direction: -z to +z
```

This follows the RAYTRACE-style local magnet convention.

### Dipole Maps

The dipole maps use an entrance-centered frame chosen for Geant4 convenience:

```text
origin: dipole entrance center
+z: incoming beam direction
+x: left in the top view
+y: upward
central trajectory: bends toward -x
```

The dipole is split into:
- entrance fringe
- sector body
- exit fringe

The split maps make the fringe boundaries explicit and avoid forcing one large rectangular map to describe regions with different natural coordinates.

## Second Multipole

The second multipole is ignored as a magnetic field source in this project. It may still appear in historical deck geometry, but no second-multipole map is generated and no second-multipole field is applied by the field-map tracer.

## Longitudinal Coordinate `L`

Ion-optics output uses a LISE-style 6D vector:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

`L` is derived from time-of-flight relative to the reference ray:

```text
L[mm] = -v0_ref * (tof - tof_ref) / gamma0_ref^2 * 10
```

The sign convention means later arrival gives negative `L`. Input `L` is a formal coordinate in the fitted transfer map; the static magnetic transport does not change an initial time offset, so `L -> L` is fixed to 1.
