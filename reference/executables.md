# Executables

Run the tools from `build/` unless you explicitly set paths. The generated field maps are normally written into the current working directory, and the field-map trace tools load maps from the current working directory by default.

## `FindMdmField`

Tune the scalar MDM dipole field for a requested ion and input angle:

```bash
./FindMdmField ../config/MDMFindField.json
```

Run this first when setting up a new ion. The output field or equivalent probe values can then be copied into `config/MDM.json` before generating field maps.

The tool uses legacy RAYTRACE because it can vary the magnetic field directly. It computes an initial rigidity field estimate, then minimizes:

```text
chi2 = (X1/positionScaleCm)^2
     + (Y1/positionScaleCm)^2
     + (AngX1/angleScaleDeg)^2
     + (AngY1/angleScaleDeg)^2
```

Output includes:
- initial rigidity field
- tuned `mdmDipoleField` in Gauss
- equivalent `mdmDipoleProbe`
- equivalent `mdmMultipoleProbe`
- final ray coordinates and angles
- final weighted residual

## `MdmFieldMapGenerator`

Generate direct RAYTRACE-sampled magnetic field maps:

```bash
./MdmFieldMapGenerator ../config/MDM.json
```

Outputs:
- `Multipole.bin`
- `DipoleEntrance.bin`
- `DipoleSector.bin`
- `DipoleExit.bin`

The generator uses RAYTRACE field routines only for stored node values. It does not use interpolation-derived adaptive refinement. Grid spacing is set by `fieldMapSpacingMm`.

The second multipole is not exported.

## `MdmTraceExample`

Run the legacy RAYTRACE transport through `MdmTrace`:

```bash
./MdmTraceExample ../config/MDM.json
./MdmTraceExample ../config/MDMScan.json
```

Each ray prints one final line:

```text
Scattered Angle X: ...deg  Scattered Angle Y: ...deg  Scattered Energy: ...MeV  X1: ...cm  Y1: ...cm  AngX1: ...deg  AngY1: ...deg
```

This is the reference result for validating the generated maps.

## `MdmFieldMapTraceExample`

Run transport through the generated field maps:

```bash
./MdmFieldMapTraceExample ../config/MDM.json
./MdmFieldMapTraceExample ../config/MDMScan.json
```

The stdout format matches `MdmTraceExample`. Diagnostics and map-setting mismatch messages are printed separately.

The field-map tracer is Fortran-free. It uses:
- generated maps for the entrance multipole and dipole
- straight drifts where no field exists
- collimator apertures copied from `rayin.dat`
- no second multipole field

## `Compare`

Run the legacy and field-map transports over the scan config and write ROOT comparison plots:

```bash
./Compare
```

Default config:

```text
config/MDMScan.json
```

Output:

```text
Compare.root
```

The ROOT file contains four canvases:
- `c_X1`
- `c_Y1`
- `c_AngX1`
- `c_AngY1`

Each canvas shows `Legacy` vs `FieldMap`, a fitted linear function, the displayed `R^2`, and the residual `Legacy - FieldMap`. `compareProcesses` in the config controls multiprocessing. If it is missing or zero, the tool uses the available CPU count.

## `GenerateIonOptics`

Fit a LISE-style 6D ion-optical transfer map from field-map transport:

```bash
./GenerateIonOptics
./GenerateIonOptics ../config/MDM.json
```

Default config:

```text
config/MDM.json
```

Output:

```text
IonOpticsMatrix.txt
```

Phase-space vector:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

The fit supports `ionOptics.order = 1..5`. The output starts with a human-readable first-order matrix and determinant, then writes COSY-style coefficient rows.

`L` is derived from time of flight relative to the reference ray. Input `L` is a formal coordinate, so `L -> L` is fixed to 1.
