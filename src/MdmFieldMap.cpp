#include "MdmFieldMap.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

std::string Trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) {
    return "";
  }
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

Vec3 Triple(const std::string& s) {
  Vec3 v;
  std::istringstream in(s);
  in >> v.x >> v.y >> v.z;
  return v;
}

void WriteIf(std::ostream& out, const char* key, const std::string& value) {
  if (!value.empty()) {
    out << key << "=" << value << "\n";
  }
}

bool WrittenByHeader(const std::string& key) {
  return key == "version" || key == "magnet" || key == "units_length" ||
         key == "units_field" || key == "nx" || key == "ny" || key == "nz" ||
         key == "origin_cm" || key == "spacing_cm" ||
         key == "payload_layout" || key == "axis_definition" ||
         key == "coordinate_system" || key == "mdm_dipole_probe" ||
         key == "mdm_multipole_probe" ||
         key == "multipole_aperture_radius_cm";
}

double Lerp(double a, double b, double t) { return a + (b - a) * t; }

}  // namespace

MdmFieldMap::MdmFieldMap(const std::string& file) { Load(file); }

MdmFieldMap::MdmFieldMap(MdmFieldMapHeader header,
                         std::vector<float> bxIn,
                         std::vector<float> byIn,
                         std::vector<float> bzIn)
    : h(std::move(header)),
      bx(std::move(bxIn)),
      by(std::move(byIn)),
      bz(std::move(bzIn)) {
  if (bx.size() != Size() || by.size() != Size() || bz.size() != Size()) {
    throw std::runtime_error("field-map payload size does not match header");
  }
}

void MdmFieldMap::Load(const std::string& file) {
  std::ifstream in(file, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot open field map: " + file);
  }

  std::map<std::string, std::string> kv;
  std::string line;
  while (std::getline(in, line)) {
    line = Trim(line);
    if (line == "END_HEADER") {
      break;
    }
    if (line.empty()) {
      continue;
    }
    const auto eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::runtime_error("bad field-map header line: " + line);
    }
    kv[Trim(line.substr(0, eq))] = Trim(line.substr(eq + 1));
  }

  h = MdmFieldMapHeader{};
  h.magnet = kv.at("magnet");
  h.nx = std::stoi(kv.at("nx"));
  h.ny = std::stoi(kv.at("ny"));
  h.nz = std::stoi(kv.at("nz"));
  h.origin_cm = Triple(kv.at("origin_cm"));
  h.step_cm = Triple(kv.at("spacing_cm"));
  if (kv.count("payload_layout")) {
    h.payload_layout = kv["payload_layout"];
  }
  if (kv.count("axis_definition")) {
    h.axis_definition = kv["axis_definition"];
  }
  if (kv.count("coordinate_system")) {
    h.coordinate_system = kv["coordinate_system"];
  }
  if (kv.count("mdm_dipole_probe")) {
    h.mdm_dipole_probe = std::stod(kv["mdm_dipole_probe"]);
  }
  if (kv.count("mdm_multipole_probe")) {
    h.mdm_multipole_probe = std::stod(kv["mdm_multipole_probe"]);
  }
  if (kv.count("multipole_aperture_radius_cm")) {
    h.aperture_radius_cm = std::stod(kv["multipole_aperture_radius_cm"]);
  }
  for (const auto& item : kv) {
    if (!WrittenByHeader(item.first)) {
      h.extra[item.first] = item.second;
    }
  }

  bx.assign(Size(), 0.0f);
  by.assign(Size(), 0.0f);
  bz.assign(Size(), 0.0f);
  in.read(reinterpret_cast<char*>(bx.data()),
          static_cast<std::streamsize>(bx.size() * sizeof(float)));
  in.read(reinterpret_cast<char*>(by.data()),
          static_cast<std::streamsize>(by.size() * sizeof(float)));
  in.read(reinterpret_cast<char*>(bz.data()),
          static_cast<std::streamsize>(bz.size() * sizeof(float)));
  if (!in) {
    throw std::runtime_error("cannot read field-map payload: " + file);
  }
}

void MdmFieldMap::Save(const std::string& file) const {
  std::ofstream out(file, std::ios::binary);
  if (!out) {
    throw std::runtime_error("cannot create field map: " + file);
  }
  out << std::setprecision(17);

  const auto version = h.extra.find("version");
  out << "version=" << (version == h.extra.end() ? "1" : version->second)
      << "\n";
  out << "magnet=" << h.magnet << "\n";
  out << "units_length=cm\n";
  out << "units_field=Tesla\n";
  out << "nx=" << h.nx << "\n";
  out << "ny=" << h.ny << "\n";
  out << "nz=" << h.nz << "\n";
  out << "origin_cm=" << h.origin_cm.x << " " << h.origin_cm.y << " "
      << h.origin_cm.z << "\n";
  out << "spacing_cm=" << h.step_cm.x << " " << h.step_cm.y << " "
      << h.step_cm.z << "\n";
  WriteIf(out, "payload_layout", h.payload_layout);
  WriteIf(out, "axis_definition", h.axis_definition);
  WriteIf(out, "coordinate_system", h.coordinate_system);
  out << "mdm_dipole_probe=" << h.mdm_dipole_probe << "\n";
  out << "mdm_multipole_probe=" << h.mdm_multipole_probe << "\n";
  if (h.aperture_radius_cm > 0.0) {
    out << "multipole_aperture_radius_cm=" << h.aperture_radius_cm << "\n";
  }
  for (const auto& item : h.extra) {
    if (!WrittenByHeader(item.first)) {
      out << item.first << "=" << item.second << "\n";
    }
  }

  out << "END_HEADER\n";
  out.write(reinterpret_cast<const char*>(bx.data()),
            static_cast<std::streamsize>(bx.size() * sizeof(float)));
  out.write(reinterpret_cast<const char*>(by.data()),
            static_cast<std::streamsize>(by.size() * sizeof(float)));
  out.write(reinterpret_cast<const char*>(bz.data()),
            static_cast<std::streamsize>(bz.size() * sizeof(float)));
}

bool MdmFieldMap::Inside(double x, double y, double z) const {
  const auto inAxis = [](double q, double q0, double dq, int n) {
    const double q1 = q0 + dq * static_cast<double>(n - 1);
    return q >= q0 && q <= q1;
  };
  if (!inAxis(x, h.origin_cm.x, h.step_cm.x, h.nx) ||
      !inAxis(y, h.origin_cm.y, h.step_cm.y, h.ny) ||
      !inAxis(z, h.origin_cm.z, h.step_cm.z, h.nz)) {
    return false;
  }
  return h.aperture_radius_cm <= 0.0 ||
         x * x + y * y <= h.aperture_radius_cm * h.aperture_radius_cm;
}

Vec3 MdmFieldMap::FieldTesla(double x, double y, double z) const {
  if (!Inside(x, y, z)) {
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

  int x0, x1, y0, y1, z0, z1;
  double tx, ty, tz;
  index1(x, h.origin_cm.x, h.step_cm.x, h.nx, x0, x1, tx);
  index1(y, h.origin_cm.y, h.step_cm.y, h.ny, y0, y1, ty);
  index1(z, h.origin_cm.z, h.step_cm.z, h.nz, z0, z1, tz);

  const auto interpolate = [&](const std::vector<float>& c) {
    const auto s = [&](int ix, int iy, int iz) {
      return static_cast<double>(c[Index(ix, iy, iz)]);
    };
    const double c00 = Lerp(s(x0, y0, z0), s(x1, y0, z0), tx);
    const double c10 = Lerp(s(x0, y1, z0), s(x1, y1, z0), tx);
    const double c01 = Lerp(s(x0, y0, z1), s(x1, y0, z1), tx);
    const double c11 = Lerp(s(x0, y1, z1), s(x1, y1, z1), tx);
    return Lerp(Lerp(c00, c10, ty), Lerp(c01, c11, ty), tz);
  };

  return {interpolate(bx), interpolate(by), interpolate(bz)};
}

std::size_t MdmFieldMap::Index(int ix, int iy, int iz) const {
  return static_cast<std::size_t>(ix) +
         static_cast<std::size_t>(h.nx) *
             (static_cast<std::size_t>(iy) +
              static_cast<std::size_t>(h.ny) * static_cast<std::size_t>(iz));
}

std::size_t MdmFieldMap::Size() const {
  return static_cast<std::size_t>(h.nx) * static_cast<std::size_t>(h.ny) *
         static_cast<std::size_t>(h.nz);
}
