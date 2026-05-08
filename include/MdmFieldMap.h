#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

// Tiny standalone field-map reader/interpolator.
//
// This is the class to copy into an external simulation project.  It has no
// dependency on RAYTRACE, ROOT, JSON, or Geant4.  A Geant4 field class should
// convert the global Geant4 position to the magnet-local coordinates used by
// the map, call FieldTesla(x_cm,y_cm,z_cm), then rotate the returned local
// field vector back to the Geant4 global frame.
//
// Binary map format:
//   ASCII header: key=value lines, terminated by END_HEADER
//   Payload: little-endian float32, component-major and x-fastest:
//     Bx[nx*ny*nz], then By[nx*ny*nz], then Bz[nx*ny*nz]
//   Index: ix + nx*(iy + ny*iz)
//   Position units: cm
//   Field units: Tesla
//   FieldTesla returns zero outside the rectangular map or aperture mask.

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct MdmFieldMapHeader {
  std::string magnet;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  Vec3 origin_cm;
  Vec3 step_cm;
  std::string axis_definition;
  std::string coordinate_system;
  std::string payload_layout = "component_major_x_fastest_float32";
  double mdm_dipole_probe = 0.0;
  double mdm_multipole_probe = 0.0;
  double aperture_radius_cm = 0.0;
  std::map<std::string, std::string> extra;
};

class MdmFieldMap {
 public:
  MdmFieldMap() = default;
  explicit MdmFieldMap(const std::string& file);
  MdmFieldMap(MdmFieldMapHeader header,
              std::vector<float> bx,
              std::vector<float> by,
              std::vector<float> bz);

  void Load(const std::string& file);
  void Save(const std::string& file) const;

  Vec3 FieldTesla(double x_cm, double y_cm, double z_cm) const;
  bool Inside(double x_cm, double y_cm, double z_cm) const;
  std::size_t Index(int ix, int iy, int iz) const;
  std::size_t Size() const;

  MdmFieldMapHeader h;
  std::vector<float> bx;
  std::vector<float> by;
  std::vector<float> bz;
};
