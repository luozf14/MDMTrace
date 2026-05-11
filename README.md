# MdmTrace

MdmTrace is a C++ interface around the MIT RAYTRACE code plus MDM-specific tools for generating magnetic field maps, validating them against the legacy transport, and fitting ion-optical transfer maps.

The normal goal of this project is to export field maps that can be used in a Geant4 model of the MDM spectrometer. The active RAYTRACE deck is `dat/rayin.dat`.

Current model:
- entrance multipole on
- dipole on
- second multipole ignored as a field element
- field-map transport uses generated maps plus hard-coded zero-field drifts and apertures from the deck

## Build

```bash
cmake -S . -B build
cmake --build build -j4
```

ROOT is required for `Compare` and `GenerateIonOptics`. The field-map loader and field-map transport code do not depend on ROOT or Fortran.

## Quickstart

Run the normal field-map and validation workflow from `build/`:

```bash
cd build
./FindMdmField ../config/MDMFindField.json
./MdmFieldMapGenerator ../config/MDM.json
./MdmTraceExample ../config/MDM.json
./MdmFieldMapTraceExample ../config/MDM.json
./Compare
```

Use `FindMdmField` first when you need a dipole setting for a new ion. Put the tuned field or equivalent probe values into `config/MDM.json`, then generate maps and validate transport.

The generator writes:
- `Multipole.bin`
- `DipoleEntrance.bin`
- `DipoleSector.bin`
- `DipoleExit.bin`

The trace examples print the same final-result line format so the legacy RAYTRACE result and field-map result can be compared directly.

## Config Files

There are three self-contained JSON configs:

- `config/MDM.json`: normal setup. Used by `MdmFieldMapGenerator`, `MdmTraceExample`, `MdmFieldMapTraceExample`, and `GenerateIonOptics`.
- `config/MDMScan.json`: larger angle/energy scan. Used by `Compare` by default and by optional scan runs of the two trace examples.
- `config/MDMFindField.json`: field-tuning setup for `FindMdmField`.

Important conventions:
- JSON keys keep the original `mdm...` names, for example `mdmDipoleProbe` and `mdmDipoleField`.
- `scatteredIon` defines the isotope and transported charge state; the tools read AME2020 atomic masses from `dat/mass_1.mas20.txt`.
- If `usingProbe` is true, tools use `mdmDipoleProbe` and `mdmMultipoleProbe`.
- If `usingProbe` is false, tools use `mdmDipoleField` and derive the equivalent probes with the project calibration rules.
- Relative map paths are resolved from the current working directory, so the usual workflow is to run tools from `build/`.

## Reference

Detailed documentation is split into:

- [Executables](reference/executables.md): command-line usage and outputs for the generator, trace examples, comparison tool, ion-optics fitter, and field tuner.
- [Configuration](reference/configuration.md): JSON keys, scan-grid rules, map paths, ion-optics settings, and field-finder settings.
- [Physics Conventions](reference/physics.md): beamline sequence, magnet-setting convention, coordinate systems, transport model, second-multipole behavior, and `L` definition.
- [Field Map Format](reference/field-map-format.md): binary file layout, header keys, `MdmFieldMap` usage notes, split-dipole metadata, and masked-zero regions.
- [Validation Workflow](reference/validation.md): legacy-vs-field-map comparison steps, ROOT comparison plots, map compatibility checks, and known limitations.
