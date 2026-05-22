# Field-Map Visualization

The `visual/` macros provide ROOT-only diagnostics for generated MDM field
maps. They read the binary field-map format directly and draw interactive ROOT
canvases by default. They do not write PNG, PDF, or ROOT files unless a user
adds that manually from the ROOT session.

## Default Usage

Run from `build/` after generating field maps:

```bash
root -l ../visual/PlotAllFieldMaps.C
```

The default macro opens multipole slice/profile canvases, one combined dipole
top-down `By` canvas, and dipole centerline profiles for:

```text
Multipole.bin
DipoleEntrance.bin
DipoleSector.bin
DipoleExit.bin
```

Individual views can be run directly:

```bash
root -l '../visual/PlotFieldSlices.C("DipoleSector.bin")'
root -l '../visual/PlotFieldProfiles.C("Multipole.bin")'
root -l '../visual/PlotDipoleTopDown.C("DipoleEntrance.bin","DipoleSector.bin","DipoleExit.bin")'
```

For a non-interactive smoke test, add `-b -q`.

## Macros

`FieldMapCommon.C` contains the shared field-map reader and plotting helpers.
It parses the ASCII header, reads the `Bx`, `By`, and `Bz` float arrays, and
uses the stored grid origin and spacing for all axes.

`PlotFieldSlices.C` treats multipole and dipole maps differently. For the
multipole it draws three transverse `x-y` slices: one in the entrance fringe,
one in the uniform region, and one in the exit fringe. The three panels share
the same `|B|` color scale and overlay `Bx,By` arrows for transverse field
direction. For dipole maps it remains a storage-coordinate diagnostic at the
requested `y` value:

- `|B|(x,z)` or `|B|(dr,s)` for total field strength
- `By(x,z)` or `By(dr,s)` for the dominant dipole component

`PlotFieldProfiles.C` draws a longitudinal multipole field-strength summary
instead of using the zero-field symmetry axis. The multipole profile shows
angular-average `|B|` versus `z` at several fixed radii, with dashed lines
marking the transition planes. For dipoles it still draws `Bx`, `By`, `Bz`, and
`|B|` along a longitudinal line, which is useful for checking fringe falloff and
sector flatness.

`PlotDipoleTopDown.C` loads the entrance, sector, and exit dipole maps together
and projects all three into the same physical dipole-local top-down frame. The
horizontal axis is dipole local `x [cm]`, the vertical axis is dipole local
`z [cm]`, and the bend should curve toward `-x`. The color scale shows `By`,
the dominant dipole bending-field component. The plot is sampled in physical
top-down coordinates and uses trilinear interpolation from the underlying field
maps, so valid mapped regions should not show internal grid holes.

`PlotField3D.C` is left as an explicit debug macro, but it is not part of the
default workflow. ROOT's sparse 3D vector display is difficult to interpret for
the curved dipole geometry.

`PlotAllFieldMaps.C` is the default entrypoint. It draws the multipole slice and
profile, the combined dipole top-down view, and the dipole profiles.

## Interpreting Plots

For dipole maps, the storage coordinates are not the best visual coordinates.
Entrance is stored as `(xB,zB)`, sector as `(dr,s)`, and exit as `(xC,zC)`.
The top-down macro transforms all three regions back into physical dipole-local
`x-z` coordinates before plotting, so it shows the actual bend geometry.

`By` should usually dominate in the dipole. The top-down `By` plot is therefore
the easiest first check for the main bending field, fringe regions, and sector
flatness.

For the multipole map, compare the total-field slice with the component
vectors and longitudinal profile to see transverse structure, the circular
aperture mask, and the entrance-uniform-exit field evolution. Zero regions in
the plots can be physical masking or points outside a represented field region;
the maps store those values as zero.

The longitudinal axis label follows the map coordinate system:

- multipole: `z`
- dipole entrance: `zB`
- dipole sector: `s`
- dipole exit: `zC`

The transverse horizontal axis is similarly labeled as `x`, `xB`, `dr`, or
`xC` depending on the map.
