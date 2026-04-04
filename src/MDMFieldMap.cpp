#include "MDMFieldMap.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::array<double, 3> ParseTriple(const std::string& value) {
  std::istringstream stream(value);
  std::array<double, 3> triple{};
  if (!(stream >> triple[0] >> triple[1] >> triple[2])) {
    throw std::runtime_error("Failed to parse triple value: " + value);
  }
  return triple;
}

std::size_t ParseSize(const std::string& value) {
  return static_cast<std::size_t>(std::stoull(value));
}

double ParseDouble(const std::string& value) { return std::stod(value); }

std::string Trim(const std::string& value) {
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return "";
  }
  const auto end = value.find_last_not_of(" \t\r\n");
  return value.substr(begin, end - begin + 1);
}

double Lerp(double a, double b, double t) {
  return a + (b - a) * t;
}

double DegreesToRadians(double degrees) {
  return degrees * 3.14159265358979323846 / 180.0;
}

bool InsideDipoleRegion(const std::map<std::string, std::string>& fields,
                        double xCm,
                        double yCm,
                        double zCm) {
  const auto find = [&](const std::string& key) -> const std::string* {
    const auto it = fields.find(key);
    return (it == fields.end()) ? nullptr : &it->second;
  };

  const auto* gapValue = find("dipole_gap_cm");
  if (gapValue == nullptr) {
    return true;
  }

  const double gap = ParseDouble(*gapValue);
  const double radius = ParseDouble(*find("dipole_outer_radius_cm")) -
                        ParseDouble(*find("dipole_strip_half_width_cm"));
  const double phiDeg = ParseDouble(*find("dipole_sector_angle_deg"));
  const double alphaDeg = ParseDouble(*find("dipole_alpha_deg"));
  const double betaDeg = ParseDouble(*find("dipole_beta_deg"));
  const double z11 = ParseDouble(*find("dipole_z11_cm"));
  const double z12 = ParseDouble(*find("dipole_z12_cm"));
  const double z21 = ParseDouble(*find("dipole_z21_cm"));
  const double z22 = ParseDouble(*find("dipole_z22_cm"));
  const double stripHalfWidth = ParseDouble(*find("dipole_strip_half_width_cm"));

  if (std::abs(yCm) > gap / 2.0) {
    return false;
  }

  const double alpha = DegreesToRadians(alphaDeg);
  const double cosAlpha = std::cos(alpha);
  const double sinAlpha = std::sin(alpha);
  const double xB = -xCm * cosAlpha - zCm * sinAlpha;
  const double zB = xCm * sinAlpha - zCm * cosAlpha;

  const double phi = DegreesToRadians(phiDeg);
  const double rotation = DegreesToRadians(phiDeg - alphaDeg - betaDeg);
  const double cosRot = std::cos(rotation);
  const double sinRot = std::sin(rotation);
  const double cosPb = std::cos(DegreesToRadians(phiDeg / 2.0 - betaDeg));
  const double sinPb = std::sin(DegreesToRadians(phiDeg / 2.0 - betaDeg));
  const double sinHalfPhi = std::sin(phi / 2.0);
  const double tx = 2.0 * radius * sinHalfPhi * sinPb;
  const double tz = 2.0 * radius * sinHalfPhi * cosPb;
  const double xC = -zB * sinRot - xB * cosRot - tx;
  const double zC = -zB * cosRot + xB * sinRot - tz;

  if (std::abs(xB) <= stripHalfWidth && zB >= z12 && zB <= z11) {
    return true;
  }

  if (std::abs(xC) <= stripHalfWidth && zC >= z21 && zC <= z22) {
    return true;
  }

  const double radialDistance =
      std::sqrt((xCm + radius) * (xCm + radius) + zCm * zCm);
  const double dr = radialDistance - radius;
  const double theta = std::atan2(zCm, xCm + radius);
  return dr >= -stripHalfWidth && dr <= stripHalfWidth && theta >= 0.0 &&
         theta <= phi;
}

}  // namespace

MDMFieldMap::MDMFieldMap(MDMFieldMapMetadata metadata,
                         std::vector<float> bx,
                         std::vector<float> by,
                         std::vector<float> bz)
    : metadata_(std::move(metadata)),
      bx_(std::move(bx)),
      by_(std::move(by)),
      bz_(std::move(bz)) {
  const std::size_t expectedSize = metadata_.nx * metadata_.ny * metadata_.nz;
  if (bx_.size() != expectedSize || by_.size() != expectedSize ||
      bz_.size() != expectedSize) {
    throw std::runtime_error("Field map payload size does not match metadata");
  }
}

MDMFieldMap MDMFieldMap::Load(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) {
    throw std::runtime_error("Unable to open field map: " + path);
  }

  MDMFieldMapMetadata metadata;
  std::string line;
  while (std::getline(stream, line)) {
    line = Trim(line);
    if (line == "END_HEADER") {
      break;
    }
    if (line.empty()) {
      continue;
    }

    const auto separator = line.find('=');
    if (separator == std::string::npos) {
      throw std::runtime_error("Malformed header line: " + line);
    }

    const std::string key = Trim(line.substr(0, separator));
    const std::string value = Trim(line.substr(separator + 1));
    metadata.fields[key] = value;
  }

  metadata.magnetName = metadata.fields.at("magnet");
  metadata.nx = ParseSize(metadata.fields.at("nx"));
  metadata.ny = ParseSize(metadata.fields.at("ny"));
  metadata.nz = ParseSize(metadata.fields.at("nz"));
  metadata.originCm = ParseTriple(metadata.fields.at("origin_cm"));
  metadata.spacingCm = ParseTriple(metadata.fields.at("spacing_cm"));

  const std::size_t payloadSize = metadata.nx * metadata.ny * metadata.nz;
  std::vector<float> bx(payloadSize);
  std::vector<float> by(payloadSize);
  std::vector<float> bz(payloadSize);

  stream.read(reinterpret_cast<char*>(bx.data()),
              static_cast<std::streamsize>(payloadSize * sizeof(float)));
  stream.read(reinterpret_cast<char*>(by.data()),
              static_cast<std::streamsize>(payloadSize * sizeof(float)));
  stream.read(reinterpret_cast<char*>(bz.data()),
              static_cast<std::streamsize>(payloadSize * sizeof(float)));

  if (!stream) {
    throw std::runtime_error("Failed to read field map payload: " + path);
  }

  return MDMFieldMap(std::move(metadata), std::move(bx), std::move(by),
                     std::move(bz));
}

void MDMFieldMap::Save(const std::string& path) const {
  std::ofstream stream(path, std::ios::binary);
  if (!stream) {
    throw std::runtime_error("Unable to create field map: " + path);
  }

  stream << "version=" << metadata_.fields.at("version") << "\n";
  stream << "magnet=" << metadata_.magnetName << "\n";
  stream << "units_length=cm\n";
  stream << "units_field=Tesla\n";
  stream << "nx=" << metadata_.nx << "\n";
  stream << "ny=" << metadata_.ny << "\n";
  stream << "nz=" << metadata_.nz << "\n";
  stream << "origin_cm=" << metadata_.originCm[0] << " " << metadata_.originCm[1]
         << " " << metadata_.originCm[2] << "\n";
  stream << "spacing_cm=" << metadata_.spacingCm[0] << " "
         << metadata_.spacingCm[1] << " " << metadata_.spacingCm[2] << "\n";

  for (const auto& [key, value] : metadata_.fields) {
    if (key == "version" || key == "magnet" || key == "nx" || key == "ny" ||
        key == "nz" || key == "origin_cm" || key == "spacing_cm") {
      continue;
    }
    stream << key << "=" << value << "\n";
  }

  stream << "END_HEADER\n";
  stream.write(reinterpret_cast<const char*>(bx_.data()),
               static_cast<std::streamsize>(bx_.size() * sizeof(float)));
  stream.write(reinterpret_cast<const char*>(by_.data()),
               static_cast<std::streamsize>(by_.size() * sizeof(float)));
  stream.write(reinterpret_cast<const char*>(bz_.data()),
               static_cast<std::streamsize>(bz_.size() * sizeof(float)));

  if (!stream) {
    throw std::runtime_error("Failed to write field map payload: " + path);
  }
}

std::array<double, 3> MDMFieldMap::Evaluate(double xCm,
                                            double yCm,
                                            double zCm) const {
  const auto within = [&](double coord, double origin, double spacing,
                          std::size_t count) {
    const double maxCoord = origin + spacing * static_cast<double>(count - 1);
    return coord >= origin && coord <= maxCoord;
  };

  if (!within(xCm, metadata_.originCm[0], metadata_.spacingCm[0], metadata_.nx) ||
      !within(yCm, metadata_.originCm[1], metadata_.spacingCm[1], metadata_.ny) ||
      !within(zCm, metadata_.originCm[2], metadata_.spacingCm[2], metadata_.nz)) {
    return {0.0, 0.0, 0.0};
  }

  const auto aperture = metadata_.fields.find("multipole_aperture_radius_cm");
  if (aperture != metadata_.fields.end()) {
    const double radiusCm = std::stod(aperture->second);
    if (xCm * xCm + yCm * yCm > radiusCm * radiusCm) {
      return {0.0, 0.0, 0.0};
    }
  }

  if (!InsideDipoleRegion(metadata_.fields, xCm, yCm, zCm)) {
    return {0.0, 0.0, 0.0};
  }

  const auto positionToIndex = [](double coord, double origin, double spacing,
                                  std::size_t count, std::size_t* i0,
                                  std::size_t* i1, double* fraction) {
    const double scaled = (coord - origin) / spacing;
    const double clamped = std::min(
        std::max(scaled, 0.0), static_cast<double>(count - 1) - 1.0e-12);
    *i0 = static_cast<std::size_t>(std::floor(clamped));
    *i1 = std::min(*i0 + 1, count - 1);
    *fraction = clamped - static_cast<double>(*i0);
  };

  std::size_t x0 = 0;
  std::size_t x1 = 0;
  std::size_t y0 = 0;
  std::size_t y1 = 0;
  std::size_t z0 = 0;
  std::size_t z1 = 0;
  double tx = 0.0;
  double ty = 0.0;
  double tz = 0.0;

  positionToIndex(xCm, metadata_.originCm[0], metadata_.spacingCm[0],
                  metadata_.nx, &x0, &x1, &tx);
  positionToIndex(yCm, metadata_.originCm[1], metadata_.spacingCm[1],
                  metadata_.ny, &y0, &y1, &ty);
  positionToIndex(zCm, metadata_.originCm[2], metadata_.spacingCm[2],
                  metadata_.nz, &z0, &z1, &tz);

  const auto interpolateComponent = [&](const std::vector<float>& component) {
    const auto sample = [&](std::size_t ix, std::size_t iy, std::size_t iz) {
      return static_cast<double>(component[Index(ix, iy, iz)]);
    };

    const double c000 = sample(x0, y0, z0);
    const double c100 = sample(x1, y0, z0);
    const double c010 = sample(x0, y1, z0);
    const double c110 = sample(x1, y1, z0);
    const double c001 = sample(x0, y0, z1);
    const double c101 = sample(x1, y0, z1);
    const double c011 = sample(x0, y1, z1);
    const double c111 = sample(x1, y1, z1);

    const double c00 = Lerp(c000, c100, tx);
    const double c10 = Lerp(c010, c110, tx);
    const double c01 = Lerp(c001, c101, tx);
    const double c11 = Lerp(c011, c111, tx);
    const double c0 = Lerp(c00, c10, ty);
    const double c1 = Lerp(c01, c11, ty);
    return Lerp(c0, c1, tz);
  };

  return {interpolateComponent(bx_), interpolateComponent(by_),
          interpolateComponent(bz_)};
}

const MDMFieldMapMetadata& MDMFieldMap::GetMetadata() const { return metadata_; }

const std::vector<float>& MDMFieldMap::GetBx() const { return bx_; }

const std::vector<float>& MDMFieldMap::GetBy() const { return by_; }

const std::vector<float>& MDMFieldMap::GetBz() const { return bz_; }

std::size_t MDMFieldMap::Index(std::size_t ix,
                               std::size_t iy,
                               std::size_t iz) const {
  return ix + metadata_.nx * (iy + metadata_.ny * iz);
}
