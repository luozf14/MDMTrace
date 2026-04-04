# MDMTrace

MDMTrace is a C++ interface around the MIT RAYTRACE code plus MDM-specific tools for field-map generation and validation. The active optics deck is [dat/rayin.dat](dat/rayin.dat), and the original RAYTRACE manuals are shipped in the repository root as [raytrace1.pdf](raytrace1.pdf), [raytrace2.pdf](raytrace2.pdf), and [raytrace3.pdf](raytrace3.pdf).

This repository currently models one specific MDM configuration:

- entrance multipole on,
- dipole on,
- second multipole present in the beamline geometry but field-off.

## Overview

The project has three main user-facing executables:

- `MDMTraceExample`: runs the original Fortran RAYTRACE transport through the deck.
- `MDMFieldMapGenerator`: samples the RAYTRACE magnet field formulas and writes `Multipole.bin` and `Dipole.bin`.
- `MDMFieldMapTraceExample`: transports ions through the generated field maps and compares against the original tracer output format.

At the library level, the repo exposes three main C++ interfaces:

- `MDMTrace`: thin wrapper around the original Fortran tracer and common blocks.
- `MDMFieldMap`: binary field-map loader, saver, and trilinear interpolator.
- `MDMFieldMapTrace`: field-map-based transport validator for the current MDM beamline.

## Repository Structure

- `src/`: Fortran RAYTRACE source and the C++ implementations.
- `include/`: public C++ headers.
- `config/`: example JSON config files for the example apps and generator.
- `dat/`: the active RAYTRACE optics deck, `rayin.dat`.
- `MDMTraceExample.C`: example executable for the original Fortran transport.
- `MDMFieldMapTraceExample.cpp`: example executable for the field-map validator.
- `raytrace1.pdf`, `raytrace2.pdf`, `raytrace3.pdf`: original RAYTRACE manuals.

## Build

Use the standard CMake flow:

```bash
cmake -S . -B build
cmake --build build -j4
```

The build produces:

- `build/MDMTraceExample`
- `build/MDMFieldMapGenerator`
- `build/MDMFieldMapTraceExample`

During configuration, CMake copies `dat/rayin.dat` into `build/`. Running the executables from `build/` is the intended default workflow.

## Usage

### `MDMTraceExample`

Purpose: run the original Fortran RAYTRACE transport with the current MDM deck.

Syntax:

```bash
./MDMTraceExample <config-file>
```

Example:

```bash
cd build
./MDMTraceExample ../config/config-MDMTraceExample.json
```

Input:

- shared transport JSON keys described below.

Output:

- confirmation lines about the selected settings,
- one final line per requested scattering angle in the form
  `Scattered Angle: ... X1: ... Y1: ... AngX1: ... AngY1: ...`

Notes:

- This executable uses the original Fortran tracer directly.
- In the example app, `scatteredAngles` is interpreted as a list of horizontal scattering angles; the vertical input angle is fixed to `0`.

### `MDMFieldMapGenerator`

Purpose: sample the entrance multipole and dipole fields from the RAYTRACE formulas and write binary field maps.

Syntax:

```bash
./MDMFieldMapGenerator <config-file>
```

Example:

```bash
cd build
./MDMFieldMapGenerator ../config/config-MDMFieldMapGenerator.json
```

Input:

- generator JSON keys described below.

Output:

- refinement progress on stdout,
- `Multipole.bin`,
- `Dipole.bin`

Notes:

- The generator ignores the second multipole because its field is off in the current deck.
- The dipole generator may emit a best-effort map with a warning if the refinement target is not met before the configured limits.
- Relative output paths are resolved against `outputDirectory`.

### `MDMFieldMapTraceExample`

Purpose: transport ions through `Multipole.bin` and `Dipole.bin` while keeping drifts and collimators from the current deck, then print the same final result format as `MDMTraceExample`.

Syntax:

```bash
./MDMFieldMapTraceExample <config-file>
```

Example with default map paths:

```bash
cd build
./MDMFieldMapTraceExample ../config/config-MDMTraceExample.json
```

Example with explicit map paths:

```bash
cd build
./MDMFieldMapTraceExample ../config/config-MDMFieldMapTraceExample.json
```

Input:

- shared transport JSON keys,
- optional `multipoleMapPath` and `dipoleMapPath`

Output:

- one final line per requested scattering angle in the same form used by `MDMTraceExample`

Notes:

- If `multipoleMapPath` and `dipoleMapPath` are omitted, the executable looks for `Multipole.bin` and `Dipole.bin` in the current working directory.
- The validator checks that the requested magnet settings match the metadata stored in the loaded maps. It rejects mismatches instead of silently rescaling the fields.

## JSON Configuration

### Shared Transport Keys

These keys are used by `MDMTraceExample` and `MDMFieldMapTraceExample`.

| Key | Meaning |
| --- | --- |
| `usingProbe` | If `true`, use `mdmDipoleProbe` and `mdmMultipoleProbe`. If `false`, derive equivalent probe settings from `mdmDipoleField`. |
| `mdmAngle` | MDM spectrometer angle in degrees. |
| `mdmDipoleField` | Legacy dipole field setting used by the existing code path. Ignored when `usingProbe=true`. |
| `mdmDipoleProbe` | Dipole hall-probe value. |
| `mdmMultipoleProbe` | Entrance multipole hall-probe value. |
| `scatteredMass` | Ion mass in AMU. |
| `scatteredCharge` | Ion charge state in units of `e`. |
| `scatteredEnergy` | Ion kinetic energy in MeV. |
| `scatteredAngles` | List of horizontal scattering angles in degrees. The example apps keep the vertical input angle fixed at `0`. |

### Generator Keys

These keys are used by `MDMFieldMapGenerator`.

| Key | Meaning |
| --- | --- |
| `mdmDipoleProbe` | Dipole hall-probe setting used to build the maps. |
| `mdmMultipoleProbe` | Entrance multipole hall-probe setting used to build the maps. |
| `relativeTolerance` | Target relative interpolation error for non-negligible fields. |
| `absoluteTolerance` | Target absolute interpolation error for weak fields. |
| `maxRefinementSteps` | Maximum refinement passes. |
| `maxNodesPerAxis` | Maximum allowed grid size along any axis. |
| `outputDirectory` | Base directory for relative output file names. |
| `multipoleOutput` | Output filename or path for the entrance multipole map. |
| `dipoleOutput` | Output filename or path for the dipole map. |

### Validator-Only Keys

These keys are optional and are used only by `MDMFieldMapTraceExample`.

| Key | Meaning |
| --- | --- |
| `multipoleMapPath` | Path to the multipole map. Defaults to `Multipole.bin` in the current working directory. |
| `dipoleMapPath` | Path to the dipole map. Defaults to `Dipole.bin` in the current working directory. |

## Physics And Modeling Conventions

### Beamline Sequence

For the current deck, the transport sequence is:

1. drift
2. entrance collimator
3. drift
4. first multipole
5. dipole
6. second multipole with field off
7. exit collimator
8. final drift

The field-map validator uses the same sequence. It transports through the two active magnetic elements with the saved maps and handles the zero-field sections separately in C++.

### Magnet Setting Conventions

The implemented hall-probe conversion rules are:

- dipole central field in Tesla:

```text
B_dipole = mdmDipoleProbe * 1.034 * 1e-4
```

- if `usingProbe=false`, the code derives equivalent probe settings as:

```text
dipoleProbe = mdmDipoleField / 1.034
multipoleProbe = dipoleProbe * 0.71
```

- the entrance multipole field strengths are derived from `mdmMultipoleProbe` using the embedded Jeffs calibration ratios already used by `MDMTrace.cpp`

### Transport Model

`MDMTrace` and `MDMFieldMapTrace` both follow the magnetic RAYTRACE particle model:

- particle state is advanced in `(x, y, z, vx, vy, vz)`
- only magnetic fields are used in the current MDM workflow
- ion speed is computed from the same relativistic expression used in the Fortran code:

```text
m = A * 931.48 MeV
v = sqrt((2m + E)E) / (m + E) * c
```

- the field-map validator uses the same four-stage Runge-Kutta structure as RAYTRACE for the magnetic elements
- drifts are transported as straight-line segments
- collimators are enforced from `rayin.dat`

### Coordinate Conventions

#### Multipole Map

- origin at magnet center
- `+z` along the beam direction
- `x` and `y` are transverse
- beam travels from `-z` to `+z`

This follows the centered straight-through convention used for the generated entrance multipole map.

#### Dipole Map

- origin at the dipole entrance center
- `+x` points left in the top view
- `+y` points upward
- `+z` points along the incoming beam
- the reference bend is toward `-x`

This dipole map frame was chosen for Geant4 convenience and is not the same presentation used in the original RAYTRACE 2 manual.

### Second Multipole

The second multipole remains in the beamline geometry, but its field coefficients are zero in the current deck. The generator does not write a second multipole map, and the field-map validator treats that section as zero-field transmission.

## Binary Field Map Format

Both `Multipole.bin` and `Dipole.bin` use the same container format:

1. ASCII header with one `key=value` entry per line
2. header terminator line:

```text
END_HEADER
```

3. binary payload in little-endian `float32`

The payload order is:

1. all `Bx`
2. all `By`
3. all `Bz`

The payload layout is x-fastest, component-major. For grid indices `(ix, iy, iz)`:

```text
linear_index = ix + nx * (iy + ny * iz)
```

The physical coordinate for a grid node is:

```text
x = origin_x + ix * dx
y = origin_y + iy * dy
z = origin_z + iz * dz
```

where `origin_cm = origin_x origin_y origin_z` and `spacing_cm = dx dy dz`.

### Common Header Keys

These keys are written for both map types:

| Key | Meaning |
| --- | --- |
| `version` | File format version. Current value is `1`. |
| `magnet` | Magnet name, currently `Multipole` or `Dipole`. |
| `units_length` | Length unit. Current value is `cm`. |
| `units_field` | Field unit. Current value is `Tesla`. |
| `nx`, `ny`, `nz` | Grid node counts along each axis. |
| `origin_cm` | Grid origin in centimeters. |
| `spacing_cm` | Grid spacing in centimeters. |
| `axis_definition` | Human-readable axis convention for the map. |
| `payload_layout` | Current value is `component_major_x_fastest_float32`. |
| `masked_zero_region` | `true` if regions outside the physical magnet are stored in the grid but evaluate to zero. |
| `relative_tolerance` | Requested relative refinement target used during generation. |
| `absolute_tolerance` | Requested absolute refinement target used during generation. |
| `mdm_dipole_probe` | Dipole hall-probe setting used to generate the map. |
| `mdm_multipole_probe` | Entrance multipole hall-probe setting used to generate the map. |

### Multipole-Specific Header Keys

| Key | Meaning |
| --- | --- |
| `multipole_aperture_radius_cm` | Circular aperture radius used to mask the map outside the active bore. |
| `multipole_transition_planes_cm` | Transition plane locations used during map generation and validation. |

### Dipole-Specific Header Keys

| Key | Meaning |
| --- | --- |
| `dipole_gap_cm` | Dipole gap used for the vertical acceptance mask. |
| `dipole_center_x_cm` | Bend-center x coordinate in the dipole local frame. |
| `dipole_center_z_cm` | Bend-center z coordinate in the dipole local frame. |
| `dipole_inner_radius_cm` | Inner radial boundary of the stored bend strip. |
| `dipole_outer_radius_cm` | Outer radial boundary of the stored bend strip. |
| `dipole_sector_angle_deg` | Central bend angle in degrees. |
| `dipole_alpha_deg` | Entrance face angle. |
| `dipole_beta_deg` | Exit face angle. |
| `dipole_z11_cm` | Entrance fringe extent parameter from the deck. |
| `dipole_z12_cm` | Entrance fringe extent parameter from the deck. |
| `dipole_z21_cm` | Exit fringe extent parameter from the deck. |
| `dipole_z22_cm` | Exit fringe extent parameter from the deck. |
| `dipole_strip_half_width_cm` | Half-width of the stored radial strip around the reference bend radius. |

### `masked_zero_region`

The stored arrays are rectangular, but not every grid point is physically inside a magnet. Points outside the valid region are treated as zero field:

- multipole: outside the circular aperture
- dipole: outside the valid sector or fringe region defined by the deck geometry

This allows downstream code to keep a regular grid while still honoring the physical magnet shape.

## Validation Workflow

The intended workflow is:

```bash
cmake -S . -B build
cmake --build build -j4
cd build
./MDMFieldMapGenerator ../config/config-MDMFieldMapGenerator.json
./MDMTraceExample ../config/config-MDMTraceExample.json
./MDMFieldMapTraceExample ../config/config-MDMTraceExample.json
```

Compare the final values:

- `X1`
- `Y1`
- `AngX1`
- `AngY1`

`MDMFieldMapTraceExample` is designed to make that comparison straightforward by printing the same final result format as `MDMTraceExample`.

## Limitations

- This is not a generic RAYTRACE deck runner. The current field-map workflow is specific to the MDM deck shipped in this repository.
- The field maps contain magnetic fields only. Drifts, collimators, and the inactive second multipole are handled separately by the validator or downstream transport code.
- The second multipole is currently treated as zero-field.
- The example apps interpret `scatteredAngles` as a horizontal-angle sweep only.
- The dipole map generator may write a best-effort map with a warning if refinement limits are reached before the requested tolerance.
