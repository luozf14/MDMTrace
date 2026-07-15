# Validation

## Build

```bash
cmake -S . -B build
cmake --build build -j4
```

ROOT controls only `Compare` and `GenerateIonOptics`. A machine without ROOT should still configure and build the other four executables. The build also copies `rayin.dat` and `mass_1.mas20.txt` into `build/`.

## Normal transport workflow

```bash
cd build
./FindMdmField ../config/FindField.json
./MdmFieldMapGenerator ../config/MDM.json
./MdmTraceExample ../config/MDM.json
./MdmFieldMapTraceExample ../config/MDM.json
```

The generator must retain the four filenames `Multipole.bin`, `DipoleEntrance.bin`, `DipoleSector.bin`, and `DipoleExit.bin`. Both trace examples must retain the same final-line fields and units so their output can be compared directly:

```text
X1 cm, Y1 cm, AngX1 deg, AngY1 deg
```

The map tracer loads all four required files and fails clearly for a wrong map role, inconsistent probe metadata, or a mismatch between map probes and the requested configuration.

## Scan comparison

```bash
./MdmTraceExample ../config/Compare.json
./MdmFieldMapTraceExample ../config/Compare.json
./Compare ../config/Compare.json
```

Before tracing, each scan path reports:

- energy values;
- horizontal-angle values;
- vertical-angle values;
- total rays.

`Compare` requires ROOT and writes `Compare.root`. It reports rays accepted by both transports, stopped by both, and stopped by only one transport. Its `Legacy - FieldMap` fits, canvases, RMS values, and maximum residuals exclude stopped-ray sentinel positions and use only rays accepted by both transports.

## Ion optics

```bash
./GenerateIonOptics ../config/MDM.json
```

This requires ROOT. Its text output must identify the six-dimensional vector as:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

The fit-grid `delta` values and output coefficients use percent momentum deviation.

## Configuration rejection checks

The three checked-in commented files must parse. Boundary checks include:

- `usingProbe: true` is accepted, while numeric `usingProbe: 1` is rejected;
- missing or malformed ion fields are rejected;
- `A > 0`, `Z >= 0`, `A >= Z`, and `0 <= Q <= Z`;
- kinetic energy is not negative;
- `fieldMapSpacingMm` is positive;
- scan and fit-grid steps are positive and ranges are ordered;
- ion-optics order is an integer from 1 through 5;
- ray, iteration, thread, and process counts are nonnegative integers.

Configuration files support `//` comments, and harmless keys unused by a particular executable remain allowed.

## Preserved physics checks

Source and map headers must remain consistent with:

- 160 cm dipole reference radius;
- 100-degree nominal sector angle;
- `3.0e10 cm/s` RAYTRACE convention;
- 0.1 cm multipole and dipole integration steps;
- cm positions and Tesla fields;
- component-major, x-fastest map payloads;
- trilinear interpolation;
- copied apertures and zero-field drifts;
- inactive second multipole;
- relative probe tolerance `2.0e-6`.

Smaller `fieldMapSpacingMm` values produce larger maps and may improve interpolation agreement. They do not change the underlying RAYTRACE field formulas.
