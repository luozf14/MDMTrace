# Field Map Format

The generated field maps are dense binary files with a small ASCII header followed by float32 field arrays.

Generated files:

```text
Multipole.bin
DipoleEntrance.bin
DipoleSector.bin
DipoleExit.bin
```

## File Layout

Each file has:

```text
key=value
key=value
...
END_HEADER
binary payload
```

The payload is little-endian `float32`, component-major, with x-fastest indexing:

```text
Bx[nx*ny*nz]
By[nx*ny*nz]
Bz[nx*ny*nz]
```

Indexing:

```text
i = ix + nx * (iy + ny * iz)
x = origin_x + ix * dx
y = origin_y + iy * dy
z = origin_z + iz * dz
```

Units in the files:

```text
position: cm
field: Tesla
```

## `MdmFieldMap`

`include/MdmFieldMap.h` and `src/MdmFieldMap.cpp` are intentionally standalone C++ standard-library code. They do not depend on RAYTRACE, ROOT, JsonCpp, or Geant4.

Minimal usage:

```cpp
#include "MdmFieldMap.h"

MdmFieldMap map("Multipole.bin");
Vec3 b = map.FieldTesla(x_cm, y_cm, z_cm);
```

External-use rules:
- convert the Geant4/global position into the local magnet coordinates in cm
- call `FieldTesla`
- rotate the returned local field vector back into the external frame
- `FieldTesla` returns zero outside the stored box or outside the map aperture/mask
- metadata is available in `map.h` and extra untyped fields are in `map.h.extra`

This class is the intended file to copy into a future Geant4 project. A Geant4 field class should mainly handle coordinate transforms and unit conversion around it.

## Typed Header Fields

`MdmFieldMap` reads common fields into:

```cpp
struct MdmFieldMapHeader {
  std::string magnet;
  int nx, ny, nz;
  Vec3 origin_cm;
  Vec3 step_cm;
  std::string axis_definition;
  std::string coordinate_system;
  std::string payload_layout;
  double mdm_dipole_probe;
  double mdm_multipole_probe;
  double aperture_radius_cm;
  std::map<std::string, std::string> extra;
};
```

Common header keys:

| Key | Meaning |
| --- | --- |
| `version` | Format version. |
| `magnet` | Magnet or region name. |
| `units_length` | Usually `cm`. |
| `units_field` | Usually `Tesla`. |
| `nx`, `ny`, `nz` | Grid dimensions. |
| `origin_cm` | Local coordinate of node `(0,0,0)`. |
| `spacing_cm` | Uniform grid spacing. |
| `axis_definition` | Short human-readable axis description. |
| `coordinate_system` | Local coordinate-system note. |
| `payload_layout` | Should state component-major, x-fastest. |
| `sampling_method` | `direct_raytrace`. |
| `masked_zero_region` | Whether zero values encode outside-region masking. |
| `mdm_dipole_probe` | Dipole probe setting used to generate the map. |
| `mdm_multipole_probe` | Multipole probe setting used to generate the map. |
| `requested_spacing_mm` | Requested generator spacing. |

Unknown keys are preserved in `extra`.

## Multipole Headers

Multipole-specific fields include:

| Key | Meaning |
| --- | --- |
| `multipole_aperture_radius_cm` | Circular aperture used by the map. |
| `multipole_transition_planes_cm` | Longitudinal fringe/transition plane metadata. |

The multipole map frame is centered on the magnet with beam direction along local `z`.

## Dipole Headers

Dipole-specific fields include:

| Key | Meaning |
| --- | --- |
| `dipole_region` | `entrance`, `sector`, or `exit`. |
| `field_component_frame` | Frame of stored vector components. |
| `raytrace_sector_frame` | Notes on the RAYTRACE sector coordinate relation. |
| `dipole_gap_cm` | Vertical gap used for map bounds. |
| `dipole_reference_radius_cm` | Reference bend radius. |
| `dipole_sector_angle_deg` | Nominal sector bend angle. |
| `dipole_alpha_deg` | Entrance edge angle. |
| `dipole_beta_deg` | Exit edge angle. |
| `dipole_z11_cm`, `dipole_z12_cm` | Entrance fringe limits. |
| `dipole_z21_cm`, `dipole_z22_cm` | Exit fringe limits. |
| `dipole_strip_half_width_cm` | Radial half-width of sector strip. |

The dipole local frame is entrance-centered:

```text
+z incoming beam
+x left in top view
+y up
bend toward -x
```

## Masked Zero Regions

The maps are stored on rectangular grids. Points outside the physical region are written as zero so the consumer can use dense arrays.

Examples:
- multipole nodes outside the circular aperture return zero
- dipole split maps only represent their own entrance, sector, or exit region
- outside the stored grid, `MdmFieldMap::FieldTesla` returns zero

Zero masking is not a smoothing operation. Stored nonzero nodes are direct RAYTRACE field samples.
