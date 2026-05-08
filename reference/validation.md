# Validation

The main validation is a direct comparison between the legacy RAYTRACE transport and the generated field-map transport.

## Basic Workflow

Build and generate maps:

```bash
cmake -S . -B build
cmake --build build -j4
cd build
./MdmFieldMapGenerator ../config/MDM.json
```

Run one normal ray through both transport paths:

```bash
./MdmTraceExample ../config/MDM.json
./MdmFieldMapTraceExample ../config/MDM.json
```

The two programs print the same line format. Compare:
- `X1`
- `Y1`
- `AngX1`
- `AngY1`

## Scan Validation

Use the scan config:

```bash
./MdmTraceExample ../config/MDMScan.json
./MdmFieldMapTraceExample ../config/MDMScan.json
```

`config/MDMScan.json` can include:
- `scatteredAngleGrid`
- `scatteredEnergyGrid`
- explicit `scatteredAnglePairs`

The output order is deterministic, so lines can be compared directly.

## ROOT Comparison

Run:

```bash
./Compare
```

This reads `config/MDMScan.json`, runs both transports, and writes:

```text
Compare.root
```

The ROOT file contains four canvases:
- `c_X1`
- `c_Y1`
- `c_AngX1`
- `c_AngY1`

Each canvas contains:
- top pad: field-map result versus legacy result
- fitted linear function
- displayed `R^2`
- bottom pad: residual `Legacy - FieldMap`

The compact stdout summary reports RMS and maximum residuals.

## Map Compatibility Checks

The field-map tracer checks that the requested magnetic setting matches the metadata stored in the map headers.

If `usingProbe=true`, the config probes must match:

```text
mdmDipoleProbe
mdmMultipoleProbe
```

If `usingProbe=false`, the equivalent probes are derived from `mdmDipoleField` and compared against the map metadata.

Maps are not silently rescaled. Regenerate maps when the requested magnetic setting changes.

## Known Limitations

The field-map validator is specific to the current MDM layout. It is not a general RAYTRACE deck interpreter.

Generated maps contain magnetic fields only. Drifts, apertures, and the inactive second multipole are handled by the validator code, not by the field-map files.

The field maps are sampled on regular grids. Consumers interpolate between stored RAYTRACE samples. Reducing `fieldMapSpacingMm` can improve agreement with legacy transport at the cost of larger files and slower generation.

The ion-optics fitter is a least-squares map over the chosen fit grid. High-order fits need enough accepted rays and can produce large text output.
