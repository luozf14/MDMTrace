# Field-map format

MDMTrace preserves four dense binary map outputs:

```text
Multipole.bin
DipoleEntrance.bin
DipoleSector.bin
DipoleExit.bin
```

MDMTrace uses these maps itself for standalone map transport. MdmSim is one possible external consumer, but the format is not tied to it.

## File layout

Each file contains an ASCII `key=value` header followed by the established binary payload:

```text
key=value
key=value
...
END_HEADER
binary payload
```

The payload is component-major `float32`:

```text
Bx[nx*ny*nz]
By[nx*ny*nz]
Bz[nx*ny*nz]
```

The x index changes fastest:

```text
i = ix + nx * (iy + ny * iz)
x = origin_x + ix * dx
y = origin_y + iy * dy
z = origin_z + iz * dz
```

Positions are in cm and fields are in Tesla. `MdmFieldMap::FieldTesla` uses trilinear interpolation and returns zero outside the stored grid or multipole aperture mask.

## Common metadata

The existing headers include:

| Key | Meaning |
| --- | --- |
| `version` | Map-format version. |
| `magnet` | Required role: `Multipole`, `DipoleEntrance`, `DipoleSector`, or `DipoleExit`. |
| `units_length` | `cm`. |
| `units_field` | `Tesla`. |
| `nx`, `ny`, `nz` | Grid node counts. |
| `origin_cm` | Position of node `(0,0,0)` in cm. |
| `spacing_cm` | Grid step in cm. |
| `payload_layout` | `component_major_x_fastest_float32`. |
| `axis_definition` | Human-readable coordinate definition. |
| `coordinate_system` | Named local coordinate system. |
| `sampling_method` | `direct_raytrace`. |
| `masked_zero_region` | Whether stored zeros mark a physical mask. |
| `mdm_dipole_probe` | Dipole probe used for generation, in Gauss. |
| `mdm_multipole_probe` | First-multipole probe used for generation, in Gauss. |
| `requested_spacing_mm` | Requested generator spacing in mm. |

The loader preserves additional header fields. The standalone transport checks the four `magnet` roles, checks that probe values agree across all maps, and checks them against the requested configuration.

## Multipole metadata

The `Multipole` map uses centred local Cartesian coordinates with `+z` along the beam. Its metadata includes:

- `multipole_aperture_radius_cm`;
- `multipole_transition_planes_cm`;
- deck fringe and length fields in cm.

Nodes outside the circular aperture are stored as zero. The disabled second multipole is not represented.

## Dipole metadata

The dipole is split into three roles. Their `dipole_region` values are `entrance_fringe`, `sector`, and `exit_fringe`. Shared metadata includes:

- `dipole_gap_cm`;
- `dipole_reference_radius_cm` = 160 cm for the current deck;
- `dipole_sector_angle_deg` = 100 degrees for the current deck;
- `dipole_alpha_deg` and `dipole_beta_deg`;
- entrance and exit fringe limits in cm;
- `dipole_strip_half_width_cm`;
- `field_component_frame = dipole_local_cartesian`.

The sector coordinate is `(dr, y, s)` in cm with `dr = r - RB` and `s = RB * theta`. The entrance and exit maps use their respective local Cartesian frames. Vector components in all three maps remain in the dipole-local Cartesian frame.

## Standalone reader

`include/MdmFieldMap.h` and `src/MdmFieldMap.cpp` use only the C++ standard library. They do not depend on RAYTRACE, ROOT, JsonCpp, or Geant4.

```cpp
#include "MdmFieldMap.h"

MdmFieldMap map("Multipole.bin");
Vec3 fieldTesla = map.FieldTesla(xCm, yCm, zCm);
```

An external consumer is responsible for transforming its global position into the documented map frame and rotating the returned local field vector back to its own frame.
