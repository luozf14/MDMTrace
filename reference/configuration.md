# Configuration

The project uses self-contained JSON files. There is no include or inheritance system. Each executable reads the top-level keys it needs and ignores extra keys.

## Config Files

`config/MDM.json`

Normal setup used by:
- `MdmFieldMapGenerator`
- `MdmTraceExample`
- `MdmFieldMapTraceExample`
- `GenerateIonOptics`

`config/MDMScan.json`

Scan setup used by:
- `Compare`
- optional scan runs of `MdmTraceExample`
- optional scan runs of `MdmFieldMapTraceExample`

`config/MDMFindField.json`

Field-finder setup used by:
- `FindMdmField`

## Shared Transport Keys

Common keys:

| Key | Meaning |
| --- | --- |
| `usingProbe` | If true, use hall-probe values. If false, use `mdmDipoleField` and derive probes. |
| `mdmAngle` | Spectrometer angle in degrees. |
| `mdmDipoleField` | Dipole field in Gauss, used when `usingProbe=false`. |
| `mdmDipoleProbe` | Dipole hall-probe value. |
| `mdmMultipoleProbe` | Entrance multipole hall-probe value. |
| `scatteredIon` | Isotope and transported charge state. |
| `massTablePath` | Optional AME2020 mass table path. Default is `dat/mass_1.mas20.txt` in the source tree. |
| `scatteredEnergy` | Kinetic energy in MeV. |

Ion definition:

```json
"scatteredIon": {
  "massNumber": 12,
  "atomicNumber": 6,
  "chargeState": 5
}
```

The transported rest mass is:

```text
ionMassMeV = neutralAtomicMassU * 931.49410242 - chargeState * 0.510998950
```

The neutral atomic mass comes from AME2020. Electron binding energies are ignored.

Probe convention:

```text
dipole central field [T] = mdmDipoleProbe * 1.034e-4
dipoleProbe = mdmDipoleField / 1.034
multipoleProbe = dipoleProbe * 0.71
```

When `usingProbe=true`, the explicit probe values are used. When `usingProbe=false`, `mdmDipoleField` is used and the probes are derived with the formulas above.

## Ray Scan Keys

Legacy single-axis scan:

```json
"scatteredAngles": [0.0, 1.0, 2.0]
```

This means `(xAngleDeg, yAngleDeg) = (value, 0)`.

Explicit angle pairs:

```json
"scatteredAnglePairs": [
  [0.5, -0.5],
  [1.0, 0.0]
]
```

2D angle grid:

```json
"scatteredAngleGrid": {
  "xMin": 0.5,
  "xMax": 1.5,
  "xStep": 0.5,
  "yMin": -0.5,
  "yMax": 0.5,
  "yStep": 0.5
}
```

Energy grid:

```json
"scatteredEnergyGrid": {
  "eMin": 17.9,
  "eMax": 18.1,
  "eStep": 0.1
}
```

Ordering:
- angle rays are generated from `scatteredAngles`, then `scatteredAnglePairs`, then `scatteredAngleGrid`
- `scatteredAngleGrid` is x-major
- if `scatteredEnergyGrid` exists, each angle ray is repeated for all energies in ascending order
- grid endpoints are included; if a step does not land exactly on the maximum, the maximum is appended

## Field-Map Generator Keys

Used by `MdmFieldMapGenerator`:

| Key | Meaning |
| --- | --- |
| `fieldMapSpacingMm` | Uniform sampling spacing in mm. |
| `outputDirectory` | Directory for generated maps. |
| `multipoleOutput` | Multipole map filename. |
| `dipoleEntranceOutput` | Dipole entrance-fringe map filename. |
| `dipoleSectorOutput` | Dipole sector map filename. |
| `dipoleExitOutput` | Dipole exit-fringe map filename. |

The stored values are direct RAYTRACE field evaluations at each grid node. The generator only changes coordinates and serializes the resulting field values.

## Field-Map Load Path Keys

Used by `MdmFieldMapTraceExample`, `Compare`, and `GenerateIonOptics`:

| Key | Default |
| --- | --- |
| `multipoleMapPath` | `Multipole.bin` |
| `dipoleEntranceMapPath` | `DipoleEntrance.bin` |
| `dipoleSectorMapPath` | `DipoleSector.bin` |
| `dipoleExitMapPath` | `DipoleExit.bin` |

Relative paths are resolved from the current working directory.

## Ion-Optics Keys

The `ionOptics` block is used by `GenerateIonOptics`.

Main keys:

| Key | Meaning |
| --- | --- |
| `order` | Polynomial order, from 1 through 5. |
| `threads` | Worker thread count. Missing or zero uses hardware concurrency. |
| `maxRays` | Optional cap on generated fit rays. |
| `outputPath` | Output text file, default `IonOpticsMatrix.txt`. |

Reference ray:

```json
"reference": {
  "xMm": 0.0,
  "thetaXMrad": 0.0,
  "yMm": 0.0,
  "thetaYMrad": 0.0,
  "energyMeV": 18.0
}
```

Fit grid:

```json
"fitGrid": {
  "xMinMm": -1.0,
  "xMaxMm": 1.0,
  "xStepMm": 1.0,
  "thetaXMinMrad": -1.0,
  "thetaXMaxMrad": 1.0,
  "thetaXStepMrad": 1.0,
  "yMinMm": -1.0,
  "yMaxMm": 1.0,
  "yStepMm": 1.0,
  "thetaYMinMrad": -1.0,
  "thetaYMaxMrad": 1.0,
  "thetaYStepMrad": 1.0,
  "deltaMin": -0.1,
  "deltaMax": 0.1,
  "deltaStep": 0.1
}
```

The fitted variable is:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

Input `L` is not scanned. It is a formal coordinate, so all fitted monomials with powers of input `L` are zero except the identity term `L -> L = 1`.

Polynomial orders:

| Order | Active fitted terms | Printed 6D terms |
| --- | ---: | ---: |
| 1 | 5 | 6 |
| 2 | 20 | 27 |
| 3 | 55 | 83 |
| 4 | 125 | 209 |
| 5 | 251 | 461 |

Coefficients multiply monomials directly. There are no factorial factors.

## Field-Finder Keys

The `fieldFinder` block is used by `FindMdmField`.

| Key | Default | Meaning |
| --- | ---: | --- |
| `inputAngleXDeg` | `0` | Input horizontal angle. |
| `inputAngleYDeg` | `0` | Input vertical angle. |
| `positionScaleCm` | `0.1` | Position scale in chi2. |
| `angleScaleDeg` | `0.1` | Angle scale in chi2. |
| `searchHalfWidthFraction` | `0.25` | Initial search half-width around rigidity estimate. |
| `coarseSamples` | `101` | Number of coarse scan points. |
| `fieldToleranceGauss` | `1e-3` | Golden-section field tolerance. |
| `maxIterations` | `100` | Maximum refinement iterations. |

Rigidity estimate:

```text
m = ionMassMeV
p = sqrt((2*m + T) * T)
BRho[kG cm] = p / (0.299792458 * chargeState)
mdmDipoleField[G] = 1000 * BRho / 160
```
