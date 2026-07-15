# Executables

Run the normal commands from `build/`. Relative map paths are resolved from the process working directory, so the default names refer to maps in `build/`.

## `FindMdmField`

```bash
./FindMdmField ../config/MDMFindField.json
```

This tool uses the legacy Fortran-backed RAYTRACE transport to tune one scalar dipole field. It starts from the 160 cm rigidity estimate, performs the existing coarse scan and golden-section refinement, and prints:

- convergence status and refinement iteration count;
- tuned `mdmDipoleField` in Gauss;
- derived dipole and first-multipole probes in Gauss;
- final `X1`, `Y1` residuals in cm;
- final `AngX1`, `AngY1` residuals in degrees;
- final weighted objective.

Field and probe values use at least 12 significant digits so they can be copied into `config/MDM.json` before regenerating maps. The tool does not change `MDM.json`.

## `MdmTraceExample`

```bash
./MdmTraceExample ../config/MDM.json
./MdmTraceExample ../config/MDMScan.json
```

This is legacy RAYTRACE transport through the Fortran-linked `MdmTrace` library. Before tracing it reports scan dimensions and total rays. Each final ray uses this stable line format:

```text
Scattered Angle X: ...deg  Scattered Angle Y: ...deg  Scattered Energy: ...MeV  X1: ...cm  Y1: ...cm  AngX1: ...deg  AngY1: ...deg
```

## `MdmFieldMapGenerator`

```bash
./MdmFieldMapGenerator ../config/MDM.json
```

The generator links the Fortran-backed `MdmTrace` library and evaluates the RAYTRACE field routines directly at each grid node. It writes:

```text
Multipole.bin
DipoleEntrance.bin
DipoleSector.bin
DipoleExit.bin
```

The four existing `key=value` metadata headers and binary payload layout are retained. `fieldMapSpacingMm` must be positive. The generator exports the entrance multipole and the three dipole regions; the disabled second multipole is not exported.

## `MdmFieldMapTraceExample`

```bash
./MdmFieldMapTraceExample ../config/MDM.json
./MdmFieldMapTraceExample ../config/MDMScan.json
```

This executable uses the standalone C++ map reader and map transport. It does not link the legacy Fortran transport. It loads all four required map files, checks their roles and probe metadata, reports scan dimensions, and prints the same final-line format as `MdmTraceExample`.

## `Compare`

ROOT is required to build this executable.

```bash
./Compare
./Compare ../config/MDMScan.json
```

It traces the self-contained scan with both transports and writes `Compare.root`. `compareProcesses` is a nonnegative integer; zero selects the automatic CPU-count choice. Before tracing, the program reports the number of energy, horizontal-angle, and vertical-angle values and the total ray count.

The program separately reports rays accepted by both transports, stopped by both, and stopped by only one transport. Stopped-ray sentinel positions are excluded from numerical fits, plots, RMS values, and maxima. The ROOT file contains canvases for `X1`, `Y1`, `AngX1`, and `AngY1` using only rays accepted by both transports, with residual convention `Legacy - FieldMap`.

## `GenerateIonOptics`

ROOT is required to build this executable.

```bash
./GenerateIonOptics ../config/MDM.json
```

It fits polynomial orders 1 through 5 using standalone field-map transport and writes `IonOpticsMatrix.txt`. The exact input/output phase-space vector is:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

`deltaP/P0` is percent momentum deviation. `L` is derived from time of flight relative to the reference ray, while input `L` is a formal coordinate with the identity term fixed to one.
