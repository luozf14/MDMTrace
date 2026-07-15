# Configuration

MDMTrace uses three self-contained JSON-with-comments files. They keep the `.json` extension, and each executable reads only the keys it uses. Harmless extra keys in `config/MDM.json` are therefore allowed.

The checked-in files are the canonical, readable examples:

- `config/MDM.json` for normal tracing, map generation, and ion optics;
- `config/MDMScan.json` for scan tracing and `Compare`;
- `config/MDMFindField.json` for `FindMdmField`.

Comments may use `//`. `usingProbe` must be a JSON Boolean:

```jsonc
{
  // true uses the two explicit probe readings.
  // false uses mdmDipoleField and derives both probes.
  "usingProbe": true,

  "mdmAngle": 0.0,                 // deg
  "mdmDipoleField": 2417.50,       // Gauss
  "mdmDipoleProbe": 2338.00,       // Gauss
  "mdmMultipoleProbe": 1659.98     // Gauss
}
```

Numeric `0` or `1` is not accepted for `usingProbe`.

## Shared ion and transport values

```jsonc
"scatteredIon": {
  "massNumber": 12,                // A
  "atomicNumber": 6,               // Z
  "chargeState": 5                 // Q, in units of e
},
"scatteredEnergy": 15.0            // MeV kinetic energy
```

The checks are `A > 0`, `Z >= 0`, `A >= Z`, and `0 <= Q <= Z`. Kinetic energy must not be negative. The neutral atomic mass is read from `dat/mass_1.mas20.txt`, and the transported rest mass is:

```text
ionMassMeV = neutralAtomicMassU * 931.49410242 - Q * 0.510998950
```

Magnet calibration is:

```text
dipole central field [T] = mdmDipoleProbe * 1.034e-4
dipoleProbe [Gauss] = mdmDipoleField [Gauss] / 1.034
multipoleProbe [Gauss] = dipoleProbe [Gauss] * 0.71
```

## `config/MDM.json`

The normal configuration keeps a flat shared transport section. Map generation uses positive `fieldMapSpacingMm` and these unchanged output names:

```jsonc
"fieldMapSpacingMm": 2.5,           // mm

"outputDirectory": ".",
"multipoleOutput": "Multipole.bin",
"dipoleEntranceOutput": "DipoleEntrance.bin",
"dipoleSectorOutput": "DipoleSector.bin",
"dipoleExitOutput": "DipoleExit.bin",

"multipoleMapPath": "Multipole.bin",
"dipoleEntranceMapPath": "DipoleEntrance.bin",
"dipoleSectorMapPath": "DipoleSector.bin",
"dipoleExitMapPath": "DipoleExit.bin"
```

Relative paths are resolved from the current working directory. Normal commands run from `build/`, so `.` normally means `build/`.

The ion-optics section is readable and unit-labelled:

```jsonc
"ionOptics": {
  "reference": {
    "xMm": 0.0,                    // mm
    "thetaXMrad": 0.0,             // mrad
    "yMm": 0.0,                    // mm
    "thetaYMrad": 0.0,             // mrad
    "energyMeV": 15.0              // MeV kinetic energy
  },
  "method": "fit",
  "order": 1,                      // integer, 1 through 5
  "threads": 0,                    // nonnegative integer; 0 is automatic
  "maxRays": 100000000,            // nonnegative integer ray limit
  "fitGrid": {
    "xMinMm": -5.0,                // mm
    "xMaxMm": 5.0,                 // mm
    "xStepMm": 1.0,                // mm
    "thetaXMinMrad": -34.5,        // mrad
    "thetaXMaxMrad": 34.5,         // mrad
    "thetaXStepMrad": 2.5,         // mrad
    "yMinMm": -5.0,                // mm
    "yMaxMm": 5.0,                 // mm
    "yStepMm": 1.0,                // mm
    "thetaYMinMrad": -34.5,        // mrad
    "thetaYMaxMrad": 34.5,         // mrad
    "thetaYStepMrad": 2.5,         // mrad
    "deltaMin": -5.0,              // % momentum deviation
    "deltaMax": 5.0,               // % momentum deviation
    "deltaStep": 0.25              // % momentum deviation
  },
  "outputPath": "IonOpticsMatrix.txt"
}
```

Every grid step must be positive and every minimum must be no greater than its maximum. The phase-space vector is:

```text
[x mm, thetaX mrad, y mm, thetaY mrad, L mm, deltaP/P0 %]
```

The configured `delta` values mean `100 * (p - p0) / p0`. Input `L` is not scanned; it is the formal longitudinal coordinate.

## `config/MDMScan.json`

The scan file contains its own magnet and ion settings. Its grids are:

```jsonc
"scatteredEnergyGrid": {
  "eMin": 13.0,                    // MeV
  "eMax": 17.0,                    // MeV
  "eStep": 0.5                     // MeV
},
"scatteredAngleGrid": {
  "xMin": -1.9,                    // deg
  "xMax": 1.9,                     // deg
  "xStep": 0.1,                    // deg
  "yMin": -1.9,                    // deg
  "yMax": 1.9,                     // deg
  "yStep": 0.1                     // deg
},
"compareProcesses": 0              // nonnegative integer; 0 is automatic
```

Steps must be positive and ranges must be ordered. Endpoints are included; if repeated stepping does not land exactly on a maximum, the maximum is appended. Before tracing, the trace examples and `Compare` report the energy count, horizontal-angle count, vertical-angle count, and total ray count.

`scatteredAngles: [0.0]` remains the normal one-dimensional horizontal-angle form used by `config/MDM.json`.

## `config/MDMFindField.json`

```jsonc
"fieldFinder": {
  "inputAngleXDeg": 0.0,           // deg
  "inputAngleYDeg": 0.0,           // deg
  "positionScaleCm": 0.1,          // cm
  "angleScaleDeg": 0.1,            // deg
  "searchHalfWidthFraction": 0.25,
  "coarseSamples": 101,            // integer sample count
  "fieldToleranceGauss": 0.001,    // Gauss
  "maxIterations": 1000            // nonnegative integer
}
```

The weighted objective is:

```text
(X1/positionScaleCm)^2 + (Y1/positionScaleCm)^2
+ (AngX1/angleScaleDeg)^2 + (AngY1/angleScaleDeg)^2
```

`FindMdmField` prints the tuned field, both derived probes, convergence status, iteration count, final position and angle residuals, and the final weighted objective. It does not edit another configuration file.
