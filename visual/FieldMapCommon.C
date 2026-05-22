#ifndef MDMTRACE_VISUAL_FIELD_MAP_COMMON_C
#define MDMTRACE_VISUAL_FIELD_MAP_COMMON_C

#include <TColor.h>
#include <TStyle.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mdmvis {

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct FieldMapHeader {
  std::string magnet;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  Vec3 origin_cm;
  Vec3 step_cm;
  std::string axis_definition;
  std::string coordinate_system;
  std::map<std::string, std::string> extra;
};

struct FieldMap {
  FieldMapHeader h;
  std::vector<float> bx;
  std::vector<float> by;
  std::vector<float> bz;

  std::size_t Size() const {
    return static_cast<std::size_t>(h.nx) * static_cast<std::size_t>(h.ny) *
           static_cast<std::size_t>(h.nz);
  }

  std::size_t Index(int ix, int iy, int iz) const {
    return static_cast<std::size_t>(ix) +
           static_cast<std::size_t>(h.nx) *
               (static_cast<std::size_t>(iy) +
                static_cast<std::size_t>(h.ny) * static_cast<std::size_t>(iz));
  }

  double X(int ix) const { return h.origin_cm.x + h.step_cm.x * ix; }
  double Y(int iy) const { return h.origin_cm.y + h.step_cm.y * iy; }
  double Z(int iz) const { return h.origin_cm.z + h.step_cm.z * iz; }

  double Bx(int ix, int iy, int iz) const { return bx[Index(ix, iy, iz)]; }
  double By(int ix, int iy, int iz) const { return by[Index(ix, iy, iz)]; }
  double Bz(int ix, int iy, int iz) const { return bz[Index(ix, iy, iz)]; }

  double Bmag(int ix, int iy, int iz) const {
    const double x = Bx(ix, iy, iz);
    const double y = By(ix, iy, iz);
    const double z = Bz(ix, iy, iz);
    return std::sqrt(x * x + y * y + z * z);
  }
};

double Bmag(const Vec3& field) {
  return std::sqrt(field.x * field.x + field.y * field.y +
                   field.z * field.z);
}

std::string Trim(const std::string& value) {
  const char* whitespace = " \t\r\n";
  const std::size_t begin = value.find_first_not_of(whitespace);
  if (begin == std::string::npos) {
    return "";
  }
  const std::size_t end = value.find_last_not_of(whitespace);
  return value.substr(begin, end - begin + 1);
}

Vec3 ParseVec3(const std::string& value) {
  std::istringstream in(value);
  Vec3 v;
  in >> v.x >> v.y >> v.z;
  if (!in) {
    throw std::runtime_error("bad Vec3 header value: " + value);
  }
  return v;
}

FieldMap LoadFieldMap(const char* path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error(std::string("cannot open field map: ") + path);
  }

  FieldMap map;
  std::string line;
  while (std::getline(in, line)) {
    line = Trim(line);
    if (line == "END_HEADER") {
      break;
    }
    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::runtime_error("bad field-map header line: " + line);
    }

    const std::string key = Trim(line.substr(0, eq));
    const std::string value = Trim(line.substr(eq + 1));
    map.h.extra[key] = value;

    if (key == "magnet") {
      map.h.magnet = value;
    } else if (key == "nx") {
      map.h.nx = std::stoi(value);
    } else if (key == "ny") {
      map.h.ny = std::stoi(value);
    } else if (key == "nz") {
      map.h.nz = std::stoi(value);
    } else if (key == "origin_cm") {
      map.h.origin_cm = ParseVec3(value);
    } else if (key == "spacing_cm") {
      map.h.step_cm = ParseVec3(value);
    } else if (key == "axis_definition") {
      map.h.axis_definition = value;
    } else if (key == "coordinate_system") {
      map.h.coordinate_system = value;
    }
  }

  if (!in || line != "END_HEADER") {
    throw std::runtime_error(std::string("missing END_HEADER in: ") + path);
  }
  if (map.h.nx <= 0 || map.h.ny <= 0 || map.h.nz <= 0) {
    throw std::runtime_error(std::string("invalid field-map dimensions in: ") +
                             path);
  }

  const std::size_t size = map.Size();
  map.bx.resize(size);
  map.by.resize(size);
  map.bz.resize(size);

  in.read(reinterpret_cast<char*>(map.bx.data()),
          static_cast<std::streamsize>(map.bx.size() * sizeof(float)));
  in.read(reinterpret_cast<char*>(map.by.data()),
          static_cast<std::streamsize>(map.by.size() * sizeof(float)));
  in.read(reinterpret_cast<char*>(map.bz.data()),
          static_cast<std::streamsize>(map.bz.size() * sizeof(float)));
  if (!in) {
    throw std::runtime_error(std::string("cannot read field-map payload: ") +
                             path);
  }

  return map;
}

int NearestIndex(double q, double origin, double step, int count) {
  int index = static_cast<int>(std::floor((q - origin) / step + 0.5));
  if (index < 0) {
    return 0;
  }
  if (index >= count) {
    return count - 1;
  }
  return index;
}

double AxisLowEdge(double origin, double step) { return origin - 0.5 * step; }

double AxisHighEdge(double origin, double step, int count) {
  return origin + step * static_cast<double>(count - 1) + 0.5 * step;
}

bool IsMultipole(const FieldMap& map) {
  return map.h.coordinate_system == "multipole_local_cartesian" ||
         map.h.magnet == "Multipole";
}

bool InsideMap(const FieldMap& map, double x, double y, double z) {
  const auto inAxis = [](double q, double q0, double dq, int n) {
    const double q1 = q0 + dq * static_cast<double>(n - 1);
    return q >= q0 && q <= q1;
  };
  return inAxis(x, map.h.origin_cm.x, map.h.step_cm.x, map.h.nx) &&
         inAxis(y, map.h.origin_cm.y, map.h.step_cm.y, map.h.ny) &&
         inAxis(z, map.h.origin_cm.z, map.h.step_cm.z, map.h.nz);
}

Vec3 FieldAt(const FieldMap& map, double x, double y, double z) {
  if (!InsideMap(map, x, y, z)) {
    return {};
  }

  const auto index1 = [](double q, double q0, double dq, int n, int& i0,
                         int& i1, double& t) {
    const double u = (q - q0) / dq;
    const double uc = std::min(std::max(u, 0.0),
                               static_cast<double>(n - 1) - 1.0e-12);
    i0 = static_cast<int>(std::floor(uc));
    i1 = std::min(i0 + 1, n - 1);
    t = uc - static_cast<double>(i0);
  };

  const auto lerp = [](double a, double b, double t) {
    return a + (b - a) * t;
  };

  int x0, x1, y0, y1, z0, z1;
  double tx, ty, tz;
  index1(x, map.h.origin_cm.x, map.h.step_cm.x, map.h.nx, x0, x1, tx);
  index1(y, map.h.origin_cm.y, map.h.step_cm.y, map.h.ny, y0, y1, ty);
  index1(z, map.h.origin_cm.z, map.h.step_cm.z, map.h.nz, z0, z1, tz);

  const auto interpolate = [&](const std::vector<float>& component) {
    const auto sample = [&](int ix, int iy, int iz) {
      return static_cast<double>(component[map.Index(ix, iy, iz)]);
    };
    const double c00 = lerp(sample(x0, y0, z0), sample(x1, y0, z0), tx);
    const double c10 = lerp(sample(x0, y1, z0), sample(x1, y1, z0), tx);
    const double c01 = lerp(sample(x0, y0, z1), sample(x1, y0, z1), tx);
    const double c11 = lerp(sample(x0, y1, z1), sample(x1, y1, z1), tx);
    return lerp(lerp(c00, c10, ty), lerp(c01, c11, ty), tz);
  };

  return {interpolate(map.bx), interpolate(map.by), interpolate(map.bz)};
}

int RepresentativeZIndex(const FieldMap& map) {
  double best = -1.0;
  int bestIz = 0;
  for (int iz = 0; iz < map.h.nz; ++iz) {
    double planeMax = 0.0;
    for (int iy = 0; iy < map.h.ny; ++iy) {
      for (int ix = 0; ix < map.h.nx; ++ix) {
        planeMax = std::max(planeMax, map.Bmag(ix, iy, iz));
      }
    }
    if (planeMax > best) {
      best = planeMax;
      bestIz = iz;
    }
  }
  return bestIz;
}

double HeaderDouble(const FieldMap& map, const std::string& key) {
  const auto item = map.h.extra.find(key);
  if (item == map.h.extra.end()) {
    throw std::runtime_error("missing field-map header key: " + key);
  }
  return std::stod(item->second);
}

std::vector<double> HeaderDoubles(const FieldMap& map,
                                  const std::string& key) {
  const auto item = map.h.extra.find(key);
  if (item == map.h.extra.end()) {
    throw std::runtime_error("missing field-map header key: " + key);
  }
  std::istringstream in(item->second);
  std::vector<double> values;
  double value = 0.0;
  while (in >> value) {
    values.push_back(value);
  }
  return values;
}

double DegreesToRadians(double degrees) {
  return degrees * 3.14159265358979323846 / 180.0;
}

Vec3 EntranceFrameToDipoleLocal(double xB, double y, double zB,
                                double alphaDeg) {
  const double alpha = DegreesToRadians(alphaDeg);
  const double cosAlpha = std::cos(alpha);
  const double sinAlpha = std::sin(alpha);
  return {-xB * cosAlpha + zB * sinAlpha, y,
          -xB * sinAlpha - zB * cosAlpha};
}

Vec3 DipoleLocalToEntranceFrame(double x, double y, double z,
                                double alphaDeg) {
  const double alpha = DegreesToRadians(alphaDeg);
  const double cosAlpha = std::cos(alpha);
  const double sinAlpha = std::sin(alpha);
  return {-x * cosAlpha - z * sinAlpha, y, x * sinAlpha - z * cosAlpha};
}

Vec3 SectorFrameToDipoleLocal(double dr, double y, double s, double radius) {
  const double theta = s / radius;
  const double r = radius + dr;
  return {-radius + r * std::cos(theta), y, r * std::sin(theta)};
}

Vec3 DipoleLocalToSectorFrame(double x, double y, double z, double radius) {
  const double r = std::sqrt((x + radius) * (x + radius) + z * z);
  const double theta = std::atan2(z, x + radius);
  return {r - radius, y, radius * theta};
}

Vec3 ExitFrameToDipoleLocal(double xC,
                            double y,
                            double zC,
                            double radius,
                            double phiDeg,
                            double alphaDeg,
                            double betaDeg) {
  const double phi = DegreesToRadians(phiDeg);
  const double rotation = DegreesToRadians(phiDeg - alphaDeg - betaDeg);
  const double cosRot = std::cos(rotation);
  const double sinRot = std::sin(rotation);
  const double cosPb = std::cos(DegreesToRadians(phiDeg / 2.0 - betaDeg));
  const double sinPb = std::sin(DegreesToRadians(phiDeg / 2.0 - betaDeg));
  const double sinHalfPhi = std::sin(phi / 2.0);
  const double tx = 2.0 * radius * sinHalfPhi * sinPb;
  const double tz = 2.0 * radius * sinHalfPhi * cosPb;

  const double xB = -cosRot * (xC + tx) + sinRot * (zC + tz);
  const double zB = -sinRot * (xC + tx) - cosRot * (zC + tz);
  return EntranceFrameToDipoleLocal(xB, y, zB, alphaDeg);
}

Vec3 DipoleLocalToExitFrame(double x,
                            double y,
                            double z,
                            double radius,
                            double phiDeg,
                            double alphaDeg,
                            double betaDeg) {
  const Vec3 entrance = DipoleLocalToEntranceFrame(x, y, z, alphaDeg);
  const double phi = DegreesToRadians(phiDeg);
  const double rotation = DegreesToRadians(phiDeg - alphaDeg - betaDeg);
  const double cosRot = std::cos(rotation);
  const double sinRot = std::sin(rotation);
  const double cosPb = std::cos(DegreesToRadians(phiDeg / 2.0 - betaDeg));
  const double sinPb = std::sin(DegreesToRadians(phiDeg / 2.0 - betaDeg));
  const double sinHalfPhi = std::sin(phi / 2.0);
  const double tx = 2.0 * radius * sinHalfPhi * sinPb;
  const double tz = 2.0 * radius * sinHalfPhi * cosPb;

  const double xC = -entrance.z * sinRot - entrance.x * cosRot - tx;
  const double zC = -entrance.z * cosRot + entrance.x * sinRot - tz;
  return {xC, y, zC};
}

std::string SanitizeName(const std::string& text) {
  std::string result;
  for (char c : text) {
    result += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
  }
  return result.empty() ? "field_map" : result;
}

std::string UniqueName(const std::string& prefix, const FieldMap& map) {
  static int counter = 0;
  std::ostringstream out;
  out << prefix << "_" << SanitizeName(map.h.magnet) << "_" << ++counter;
  return out.str();
}

std::string XAxisTitle(const FieldMap& map) {
  if (map.h.coordinate_system == "dipole_sector_frame") {
    return "dr [cm]";
  }
  if (map.h.coordinate_system == "dipole_entrance_frame") {
    return "xB [cm]";
  }
  if (map.h.coordinate_system == "dipole_exit_frame") {
    return "xC [cm]";
  }
  return "x [cm]";
}

std::string LongitudinalAxisTitle(const FieldMap& map) {
  if (map.h.coordinate_system == "dipole_sector_frame") {
    return "s [cm]";
  }
  if (map.h.coordinate_system == "dipole_entrance_frame") {
    return "zB [cm]";
  }
  if (map.h.coordinate_system == "dipole_exit_frame") {
    return "zC [cm]";
  }
  return "z [cm]";
}

void UseFieldMapStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
}

int PaletteColor(double value, double maxValue) {
  if (maxValue <= 0.0) {
    return kBlack;
  }
  double fraction = value / maxValue;
  if (fraction < 0.0) {
    fraction = 0.0;
  }
  if (fraction > 1.0) {
    fraction = 1.0;
  }
  const int colors = TColor::GetNumberOfColors();
  const int index = static_cast<int>(fraction * static_cast<double>(colors - 1));
  return TColor::GetColorPalette(index);
}

}  // namespace mdmvis

#endif
