# MDMTrace

MDMTrace is a standalone physics program for the current MDM spectrometer. It combines the MIT RAYTRACE transport with small C++ tools for:

- tuning the MDM dipole field;
- tracing rays with the legacy RAYTRACE model;
- generating magnetic-field maps;
- tracing rays with the standalone field-map transport;
- comparing legacy and field-map transport;
- fitting ion-optical transfer maps.

The active RAYTRACE deck is `dat/rayin.dat`. The active model includes the entrance multipole and dipole; the second multipole remains disabled. The four generated maps may be used directly by MDMTrace or by another consumer such as MdmSim.

## Build

```bash
cmake -S . -B build
cmake --build build -j4
```

C++17 and the legacy Fortran compiler flags are retained. ROOT is needed only for `Compare` and `GenerateIonOptics`. If ROOT is absent, CMake reports that those two targets are skipped and still builds `FindMdmField`, `MdmTraceExample`, `MdmFieldMapGenerator`, and `MdmFieldMapTraceExample`.

`MdmFieldMapGenerator` links the legacy Fortran-backed `MdmTrace` library because it samples the RAYTRACE field routines. `MdmFieldMapTraceExample` links only the standalone C++ map reader and transport; it does not link the legacy Fortran transport.

## Quick start

Run normal commands from `build/`:

```bash
cd build
./FindMdmField ../config/MDMFindField.json
./MdmFieldMapGenerator ../config/MDM.json
./MdmTraceExample ../config/MDM.json
./MdmFieldMapTraceExample ../config/MDM.json
./Compare
./GenerateIonOptics ../config/MDM.json
```

The generator preserves the established outputs and metadata headers:

```text
Multipole.bin
DipoleEntrance.bin
DipoleSector.bin
DipoleExit.bin
```

Relative map paths and output paths are resolved from the process working directory. Running from `build/` therefore creates and reads the four maps there.

## Configuration

The three supported configuration files are:

- `config/MDM.json`: normal tracing, generation, and ion-optics fitting;
- `config/MDMScan.json`: self-contained energy/angle scans and comparison;
- `config/MDMFindField.json`: field tuning.

They are JSON-with-comments files and retain the `.json` extension. Comments beside dimensional values state their units and physics meaning. `usingProbe` must be the JSON Boolean `true` or `false`; numeric `0` and `1` are rejected.

The ion-optics phase-space vector is exactly:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

Here `deltaP/P0` is percent momentum deviation, not a unitless fraction.

## Documentation

- [Configuration](reference/configuration.md)
- [Executables](reference/executables.md)
- [Physics conventions](reference/physics.md)
- [Field-map format](reference/field-map-format.md)
- [Validation](reference/validation.md)
- [Field-map visualization](reference/visualization.md)
