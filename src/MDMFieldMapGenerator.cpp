#include "MDMFieldMap.h"
#include "MDMFieldMapInterop.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "json.h"

namespace {

constexpr std::size_t kMultipoleElement = 4;
constexpr std::size_t kDipoleElement = 5;
constexpr double kRelativeToleranceDefault = 5.0e-3;
constexpr double kAbsoluteToleranceDefault = 1.0e-4;
constexpr std::size_t kMaxRefinementStepsDefault = 8;
constexpr std::size_t kMaxNodesPerAxisDefault = 257;
constexpr double kDipoleStripHalfWidthCm = 30.0;
constexpr double kPi = 3.14159265358979323846;

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct GridSpec {
  std::array<double, 3> min{};
  std::array<double, 3> max{};
  std::array<std::size_t, 3> counts{};
};

struct GeneratorConfig {
  double mdmDipoleProbe = 0.0;
  double mdmMultipoleProbe = 0.0;
  double relativeTolerance = kRelativeToleranceDefault;
  double absoluteTolerance = kAbsoluteToleranceDefault;
  std::size_t maxRefinementSteps = kMaxRefinementStepsDefault;
  std::size_t maxNodesPerAxis = kMaxNodesPerAxisDefault;
  std::filesystem::path outputDirectory = ".";
  std::filesystem::path multipoleOutput = "Multipole.bin";
  std::filesystem::path dipoleOutput = "Dipole.bin";
};

struct ValidationReport {
  bool withinTolerance = false;
  double maxRelativeError = 0.0;
  double maxAbsoluteError = 0.0;
  std::array<double, 3> axisScores{0.0, 0.0, 0.0};
};

using Evaluator = std::function<Vec3(double, double, double)>;
using ValidationRegion =
    std::function<bool(double, double, double, const std::array<double, 3>&)>;

double DegreesToRadians(double degrees) { return degrees * kPi / 180.0; }

double DeckValue(std::size_t elementOneBased, std::size_t fieldOneBased) {
  return blck0_.DATA[elementOneBased - 1][fieldOneBased - 1];
}

Vec3 EvaluateMultipoleDirect(double x, double y, double z) {
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  mdmfm_eval_entry_multipole_field(x, y, z, &bx, &by, &bz);
  return {bx, by, bz};
}

Vec3 EvaluateDipoleDirect(double x, double y, double z) {
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  mdmfm_eval_dipole_field(x, y, z, &bx, &by, &bz);
  return {bx, by, bz};
}

double Norm(const Vec3& value) {
  return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}

double MetricScore(const Vec3& direct,
                   const Vec3& interpolated,
                   double absoluteTolerance,
                   double relativeTolerance,
                   double* relativeError,
                   double* absoluteError) {
  const Vec3 delta{interpolated.x - direct.x, interpolated.y - direct.y,
                   interpolated.z - direct.z};
  *absoluteError = Norm(delta);
  const double directNorm = Norm(direct);
  if (directNorm >= 0.02) {
    *relativeError = *absoluteError / directNorm;
    return *relativeError / relativeTolerance;
  }
  *relativeError = 0.0;
  return *absoluteError / absoluteTolerance;
}

std::size_t FlattenIndex(const GridSpec& grid,
                         std::size_t ix,
                         std::size_t iy,
                         std::size_t iz) {
  return ix + grid.counts[0] * (iy + grid.counts[1] * iz);
}

std::array<double, 3> GridSpacing(const GridSpec& grid) {
  return {
      (grid.max[0] - grid.min[0]) / static_cast<double>(grid.counts[0] - 1),
      (grid.max[1] - grid.min[1]) / static_cast<double>(grid.counts[1] - 1),
      (grid.max[2] - grid.min[2]) / static_cast<double>(grid.counts[2] - 1)};
}

double GridCoordinate(const GridSpec& grid, std::size_t axis, std::size_t index) {
  const auto spacing = GridSpacing(grid);
  return grid.min[axis] + spacing[axis] * static_cast<double>(index);
}

std::size_t RoundUpEven(std::size_t value) {
  return (value % 2 == 0) ? value : value + 1;
}

std::size_t InitialNodeCount(double range, double characteristicScale) {
  const double targetStep =
      std::clamp(characteristicScale / 8.0, 0.25, 5.0);
  std::size_t segments =
      std::max<std::size_t>(8, RoundUpEven(static_cast<std::size_t>(
                                 std::ceil(range / targetStep))));
  return segments + 1;
}

GridSpec BuildMultipoleBounds() {
  const double l = DeckValue(kMultipoleElement, 12);
  const double rad = DeckValue(kMultipoleElement, 13);
  const double z11 = DeckValue(kMultipoleElement, 19);
  const double z12 = DeckValue(kMultipoleElement, 20);
  const double z21 = DeckValue(kMultipoleElement, 21);
  const double z22 = DeckValue(kMultipoleElement, 22);
  const double zMin = -(l / 2.0 + std::max(z11, -z12));
  const double zMax = +(l / 2.0 + std::max(z22, -z21));

  GridSpec grid;
  grid.min = {-rad, -rad, zMin};
  grid.max = {+rad, +rad, zMax};
  grid.counts = {
      InitialNodeCount(grid.max[0] - grid.min[0], rad),
      InitialNodeCount(grid.max[1] - grid.min[1], rad),
      InitialNodeCount(grid.max[2] - grid.min[2], rad)};
  return grid;
}

std::array<double, 2> EntranceFrameToLocal(double xB, double zB, double alphaDeg) {
  const double alpha = DegreesToRadians(alphaDeg);
  const double cosAlpha = std::cos(alpha);
  const double sinAlpha = std::sin(alpha);
  return {-xB * cosAlpha + zB * sinAlpha, -xB * sinAlpha - zB * cosAlpha};
}

std::array<double, 2> ExitFrameToLocal(double xC,
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
  return EntranceFrameToLocal(xB, zB, alphaDeg);
}

bool InsideDipoleRegion(double xLocal,
                        double yLocal,
                        double zLocal,
                        double gap,
                        double radius,
                        double phiDeg,
                        double alphaDeg,
                        double betaDeg,
                        double z11,
                        double z12,
                        double z21,
                        double z22,
                        double stripHalfWidth,
                        double radialMargin,
                        double verticalMargin,
                        double axialMargin) {
  if (std::abs(yLocal) > gap / 2.0 - verticalMargin) {
    return false;
  }

  const double alpha = DegreesToRadians(alphaDeg);
  const double cosAlpha = std::cos(alpha);
  const double sinAlpha = std::sin(alpha);
  const double xB = -xLocal * cosAlpha - zLocal * sinAlpha;
  const double zB = xLocal * sinAlpha - zLocal * cosAlpha;

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

  const double effectiveStrip = stripHalfWidth - radialMargin;
  if (effectiveStrip <= 0.0) {
    return false;
  }

  if (std::abs(xB) <= effectiveStrip && zB >= z12 + axialMargin &&
      zB <= z11 - axialMargin) {
    return true;
  }

  if (std::abs(xC) <= effectiveStrip && zC >= z21 + axialMargin &&
      zC <= z22 - axialMargin) {
    return true;
  }

  const double radialDistance =
      std::sqrt((xLocal + radius) * (xLocal + radius) + zLocal * zLocal);
  const double dr = radialDistance - radius;
  const double theta = std::atan2(zLocal, xLocal + radius);
  const double angleMargin = (radius > 0.0) ? radialMargin / radius : 0.0;
  if (zB > z12 - axialMargin || zC > z21 - axialMargin) {
    return false;
  }
  return dr >= -effectiveStrip && dr <= effectiveStrip &&
         theta >= angleMargin && theta <= phi - angleMargin;
}

GridSpec BuildDipoleBounds() {
  const double d = DeckValue(kDipoleElement, 13);
  const double rb = DeckValue(kDipoleElement, 14);
  const double phi = DeckValue(kDipoleElement, 16);
  const double alpha = DeckValue(kDipoleElement, 17);
  const double beta = DeckValue(kDipoleElement, 18);
  const double z11 = DeckValue(kDipoleElement, 25);
  const double z12 = DeckValue(kDipoleElement, 26);
  const double z21 = DeckValue(kDipoleElement, 27);
  const double z22 = DeckValue(kDipoleElement, 28);

  std::vector<std::array<double, 2>> points;
  const double innerRadius = rb - kDipoleStripHalfWidthCm;
  const double outerRadius = rb + kDipoleStripHalfWidthCm;
  const double phiRad = DegreesToRadians(phi);

  std::vector<double> angles{0.0, phiRad};
  if (phiRad > kPi / 2.0) {
    angles.push_back(kPi / 2.0);
  }

  for (const double radius : {innerRadius, outerRadius}) {
    for (const double angle : angles) {
      points.push_back(
          {-rb + radius * std::cos(angle), radius * std::sin(angle)});
    }
  }

  for (const double xB : {-kDipoleStripHalfWidthCm, kDipoleStripHalfWidthCm}) {
    for (const double zB : {z12, z11}) {
      points.push_back(EntranceFrameToLocal(xB, zB, alpha));
    }
  }

  for (const double xC : {-kDipoleStripHalfWidthCm, kDipoleStripHalfWidthCm}) {
    for (const double zC : {z21, z22}) {
      points.push_back(ExitFrameToLocal(xC, zC, rb, phi, alpha, beta));
    }
  }

  double xMin = std::numeric_limits<double>::max();
  double xMax = std::numeric_limits<double>::lowest();
  double zMin = std::numeric_limits<double>::max();
  double zMax = std::numeric_limits<double>::lowest();
  for (const auto& point : points) {
    xMin = std::min(xMin, point[0]);
    xMax = std::max(xMax, point[0]);
    zMin = std::min(zMin, point[1]);
    zMax = std::max(zMax, point[1]);
  }

  GridSpec grid;
  grid.min = {xMin, -d / 2.0, zMin};
  grid.max = {xMax, d / 2.0, zMax};
  grid.counts = {
      InitialNodeCount(grid.max[0] - grid.min[0], kDipoleStripHalfWidthCm),
      InitialNodeCount(grid.max[1] - grid.min[1], d),
      InitialNodeCount(grid.max[2] - grid.min[2], kDipoleStripHalfWidthCm)};
  return grid;
}

GeneratorConfig ReadConfig(const std::string& path) {
  std::ifstream stream(path.c_str());
  if (!stream) {
    throw std::runtime_error("Unable to open config file: " + path);
  }

  Json::Value config;
  stream >> config;

  GeneratorConfig result;
  result.mdmDipoleProbe = config["mdmDipoleProbe"].asDouble();
  result.mdmMultipoleProbe = config["mdmMultipoleProbe"].asDouble();

  if (config.isMember("relativeTolerance")) {
    result.relativeTolerance = config["relativeTolerance"].asDouble();
  }
  if (config.isMember("absoluteTolerance")) {
    result.absoluteTolerance = config["absoluteTolerance"].asDouble();
  }
  if (config.isMember("maxRefinementSteps")) {
    result.maxRefinementSteps =
        static_cast<std::size_t>(config["maxRefinementSteps"].asUInt64());
  }
  if (config.isMember("maxNodesPerAxis")) {
    result.maxNodesPerAxis =
        static_cast<std::size_t>(config["maxNodesPerAxis"].asUInt64());
  }
  if (config.isMember("outputDirectory")) {
    result.outputDirectory = config["outputDirectory"].asString();
  }
  if (config.isMember("multipoleOutput")) {
    result.multipoleOutput = config["multipoleOutput"].asString();
  }
  if (config.isMember("dipoleOutput")) {
    result.dipoleOutput = config["dipoleOutput"].asString();
  }

  return result;
}

MDMFieldMap BuildMap(const std::string& magnetName,
                     const GridSpec& grid,
                     const Evaluator& evaluator,
                     const GeneratorConfig& config,
                     const std::map<std::string, std::string>& extraFields) {
  const auto spacing = GridSpacing(grid);
  const std::size_t payloadSize =
      grid.counts[0] * grid.counts[1] * grid.counts[2];
  std::vector<float> bx(payloadSize, 0.0f);
  std::vector<float> by(payloadSize, 0.0f);
  std::vector<float> bz(payloadSize, 0.0f);

  for (std::size_t iz = 0; iz < grid.counts[2]; ++iz) {
    const double z = GridCoordinate(grid, 2, iz);
    for (std::size_t iy = 0; iy < grid.counts[1]; ++iy) {
      const double y = GridCoordinate(grid, 1, iy);
      for (std::size_t ix = 0; ix < grid.counts[0]; ++ix) {
        const double x = GridCoordinate(grid, 0, ix);
        const Vec3 value = evaluator(x, y, z);
        const std::size_t index = FlattenIndex(grid, ix, iy, iz);
        bx[index] = static_cast<float>(value.x);
        by[index] = static_cast<float>(value.y);
        bz[index] = static_cast<float>(value.z);
      }
    }
  }

  MDMFieldMapMetadata metadata;
  metadata.magnetName = magnetName;
  metadata.nx = grid.counts[0];
  metadata.ny = grid.counts[1];
  metadata.nz = grid.counts[2];
  metadata.originCm = {grid.min[0], grid.min[1], grid.min[2]};
  metadata.spacingCm = spacing;
  metadata.fields["version"] = "1";
  metadata.fields["axis_definition"] = extraFields.at("axis_definition");
  metadata.fields["payload_layout"] = "component_major_x_fastest_float32";
  metadata.fields["masked_zero_region"] = extraFields.at("masked_zero_region");
  metadata.fields["relative_tolerance"] = extraFields.at("relative_tolerance");
  metadata.fields["absolute_tolerance"] = extraFields.at("absolute_tolerance");
  metadata.fields["mdm_dipole_probe"] = extraFields.at("mdm_dipole_probe");
  metadata.fields["mdm_multipole_probe"] =
      extraFields.at("mdm_multipole_probe");
  for (const auto& [key, value] : extraFields) {
    if (metadata.fields.find(key) == metadata.fields.end()) {
      metadata.fields[key] = value;
    }
  }

  return MDMFieldMap(std::move(metadata), std::move(bx), std::move(by),
                     std::move(bz));
}

ValidationReport ValidateMap(const MDMFieldMap& map,
                             const Evaluator& evaluator,
                             const ValidationRegion& validationRegion,
                             double absoluteTolerance,
                             double relativeTolerance) {
  ValidationReport report;
  const auto& metadata = map.GetMetadata();
  const GridSpec grid{{metadata.originCm[0], metadata.originCm[1], metadata.originCm[2]},
                      {metadata.originCm[0] + metadata.spacingCm[0] * static_cast<double>(metadata.nx - 1),
                       metadata.originCm[1] + metadata.spacingCm[1] * static_cast<double>(metadata.ny - 1),
                       metadata.originCm[2] + metadata.spacingCm[2] * static_cast<double>(metadata.nz - 1)},
                      {metadata.nx, metadata.ny, metadata.nz}};
  const auto spacing = GridSpacing(grid);

  const std::array<std::array<double, 3>, 4> samples{{
      {0.5, 0.5, 0.5},
      {0.5, 0.25, 0.25},
      {0.25, 0.5, 0.25},
      {0.25, 0.25, 0.5},
  }};

  for (std::size_t iz = 0; iz + 1 < grid.counts[2]; ++iz) {
    const double z0 = GridCoordinate(grid, 2, iz);
    for (std::size_t iy = 0; iy + 1 < grid.counts[1]; ++iy) {
      const double y0 = GridCoordinate(grid, 1, iy);
      for (std::size_t ix = 0; ix + 1 < grid.counts[0]; ++ix) {
        const double x0 = GridCoordinate(grid, 0, ix);
        for (std::size_t sampleIndex = 0; sampleIndex < samples.size();
             ++sampleIndex) {
          const auto& sample = samples[sampleIndex];
          const double x = x0 + sample[0] * spacing[0];
          const double y = y0 + sample[1] * spacing[1];
          const double z = z0 + sample[2] * spacing[2];
          if (validationRegion && !validationRegion(x, y, z, spacing)) {
            continue;
          }
          const Vec3 direct = evaluator(x, y, z);
          const auto interpolated = map.Evaluate(x, y, z);
          const Vec3 interpolatedVec{interpolated[0], interpolated[1],
                                     interpolated[2]};

          double relativeError = 0.0;
          double absoluteError = 0.0;
          const double score =
              MetricScore(direct, interpolatedVec, absoluteTolerance,
                          relativeTolerance, &relativeError, &absoluteError);
          if (Norm(direct) >= 0.02) {
            report.maxRelativeError =
                std::max(report.maxRelativeError, relativeError);
          } else {
            report.maxAbsoluteError =
                std::max(report.maxAbsoluteError, absoluteError);
          }
          if (sampleIndex > 0) {
            report.axisScores[sampleIndex - 1] =
                std::max(report.axisScores[sampleIndex - 1], score);
          }
        }
      }
    }
  }

  report.withinTolerance = report.maxAbsoluteError <= absoluteTolerance &&
                           report.maxRelativeError <= relativeTolerance;
  return report;
}

GridSpec RefineGrid(const GridSpec& current,
                    const ValidationReport& report,
                    std::size_t maxNodesPerAxis) {
  GridSpec next = current;
  const auto axisScores = report.axisScores;
  std::array<std::size_t, 3> order{0, 1, 2};
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return axisScores[lhs] > axisScores[rhs];
  });

  if (axisScores[order[0]] > 1.25 * axisScores[order[1]]) {
    next.counts[order[0]] =
        std::min(maxNodesPerAxis, 2 * (current.counts[order[0]] - 1) + 1);
  } else {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      next.counts[axis] =
          std::min(maxNodesPerAxis, 2 * (current.counts[axis] - 1) + 1);
    }
  }
  return next;
}

std::string FormatDouble(double value) {
  std::ostringstream stream;
  stream << std::setprecision(10) << value;
  return stream.str();
}

std::string FormatGridCounts(const GridSpec& grid) {
  std::ostringstream stream;
  stream << grid.counts[0] << " x " << grid.counts[1] << " x " << grid.counts[2];
  return stream.str();
}

void VerifyReferenceSamples(const MDMFieldMap& map,
                            const Evaluator& evaluator,
                            const std::string& label) {
  const auto& metadata = map.GetMetadata();
  const std::array<Vec3, 2> samplePoints{{
      {metadata.originCm[0], metadata.originCm[1], metadata.originCm[2]},
      {metadata.originCm[0] + metadata.spacingCm[0] * 0.5 * static_cast<double>(metadata.nx - 1),
       metadata.originCm[1] + metadata.spacingCm[1] * 0.5 * static_cast<double>(metadata.ny - 1),
       metadata.originCm[2] + metadata.spacingCm[2] * 0.5 * static_cast<double>(metadata.nz - 1)},
  }};

  for (const auto& point : samplePoints) {
    const Vec3 direct = evaluator(point.x, point.y, point.z);
    const auto loaded = map.Evaluate(point.x, point.y, point.z);
    const Vec3 loadedVec{loaded[0], loaded[1], loaded[2]};
    const Vec3 delta{loadedVec.x - direct.x, loadedVec.y - direct.y,
                     loadedVec.z - direct.z};
    if (Norm(delta) > 2.0e-4) {
      throw std::runtime_error(label + " round-trip verification failed");
    }
  }
}

void EnsureParentDirectory(const std::filesystem::path& path) {
  const auto parent = path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Usage: ./MDMFieldMapGenerator <config-file>" << std::endl;
    return 0;
  }

  try {
    const GeneratorConfig config = ReadConfig(argv[1]);
    mdmfm_init();
    mdmfm_set_probes(config.mdmDipoleProbe, config.mdmMultipoleProbe);

    std::filesystem::path multipolePath = config.multipoleOutput;
    if (multipolePath.is_relative()) {
      multipolePath = config.outputDirectory / multipolePath;
    }
    std::filesystem::path dipolePath = config.dipoleOutput;
    if (dipolePath.is_relative()) {
      dipolePath = config.outputDirectory / dipolePath;
    }

    const auto generateMagnet =
        [&](const std::string& magnetName,
            const GridSpec& initialGrid,
            const Evaluator& evaluator,
            const ValidationRegion& validationRegion,
            bool forceFullRefinement,
            bool allowBestEffort,
            const std::map<std::string, std::string>& staticFields,
            const std::filesystem::path& outputPath) {
          GridSpec grid = initialGrid;
          MDMFieldMap fieldMap;
          ValidationReport report;
          MDMFieldMap bestFieldMap;
          ValidationReport bestReport;
          bool haveBestFieldMap = false;

          for (std::size_t step = 0; step < config.maxRefinementSteps; ++step) {
            std::map<std::string, std::string> metadata = staticFields;
            metadata["relative_tolerance"] = FormatDouble(config.relativeTolerance);
            metadata["absolute_tolerance"] = FormatDouble(config.absoluteTolerance);
            metadata["mdm_dipole_probe"] = FormatDouble(config.mdmDipoleProbe);
            metadata["mdm_multipole_probe"] =
                FormatDouble(config.mdmMultipoleProbe);
            fieldMap = BuildMap(magnetName, grid, evaluator, config, metadata);
            report =
                ValidateMap(fieldMap, evaluator, validationRegion,
                            config.absoluteTolerance, config.relativeTolerance);
            const bool betterRelative =
                !haveBestFieldMap ||
                report.maxRelativeError <
                    bestReport.maxRelativeError - 1.0e-12;
            const bool comparableRelative =
                !haveBestFieldMap ||
                std::abs(report.maxRelativeError - bestReport.maxRelativeError) <=
                    1.0e-12;
            const bool betterAbsolute =
                !haveBestFieldMap ||
                report.maxAbsoluteError <
                    bestReport.maxAbsoluteError - 1.0e-12;
            if (betterRelative || (comparableRelative && betterAbsolute)) {
              bestFieldMap = fieldMap;
              bestReport = report;
              haveBestFieldMap = true;
            }

            std::cout << magnetName << " refinement step " << step + 1 << ": "
                      << FormatGridCounts(grid) << ", max rel err "
                      << report.maxRelativeError << ", max abs err "
                      << report.maxAbsoluteError << std::endl;

            if (report.withinTolerance) {
              break;
            }

            if (allowBestEffort && haveBestFieldMap &&
                report.maxRelativeError > bestReport.maxRelativeError * 1.10 &&
                step > 0) {
              break;
            }

            const GridSpec refined = forceFullRefinement
                                         ? GridSpec{
                                               grid.min,
                                               grid.max,
                                               {std::min(config.maxNodesPerAxis,
                                                         2 * (grid.counts[0] - 1) + 1),
                                                std::min(config.maxNodesPerAxis,
                                                         2 * (grid.counts[1] - 1) + 1),
                                                std::min(config.maxNodesPerAxis,
                                                         2 * (grid.counts[2] - 1) + 1)}}
                                         : RefineGrid(grid, report, config.maxNodesPerAxis);
            if (refined.counts == grid.counts) {
              throw std::runtime_error(
                  "Reached node-count limit before satisfying tolerance for " +
                  magnetName);
            }
            grid = refined;
          }

          if (!report.withinTolerance) {
            if (!allowBestEffort || !haveBestFieldMap) {
              throw std::runtime_error("Failed to converge for " + magnetName);
            }
            fieldMap = bestFieldMap;
            report = bestReport;
            std::cerr << "WARNING: Using best-effort " << magnetName
                      << " map after validation shortfall. max rel err "
                      << report.maxRelativeError << ", max abs err "
                      << report.maxAbsoluteError << std::endl;
          }

          EnsureParentDirectory(outputPath);
          fieldMap.Save(outputPath.string());
          const MDMFieldMap reloaded = MDMFieldMap::Load(outputPath.string());
          VerifyReferenceSamples(reloaded, evaluator, magnetName);
          std::cout << "Wrote " << outputPath << std::endl;
        };

    const GridSpec multipoleGrid = BuildMultipoleBounds();
    const double multipoleRadius = DeckValue(kMultipoleElement, 13);
    const double multipoleLength = DeckValue(kMultipoleElement, 12);
    const double multipoleZ12 = DeckValue(kMultipoleElement, 20);
    const double multipoleZ21 = DeckValue(kMultipoleElement, 21);
    const double multipoleEnterPlane = -0.5 * multipoleLength - multipoleZ12;
    const double multipoleExitPlane = 0.5 * multipoleLength + multipoleZ21;
    const ValidationRegion multipoleValidationRegion =
        [multipoleRadius, multipoleEnterPlane, multipoleExitPlane](
            double x, double y, double z,
                          const std::array<double, 3>& spacing) {
          const double margin = std::hypot(spacing[0], spacing[1]);
          const double effectiveRadius = multipoleRadius - margin;
          return effectiveRadius > 0.0 &&
                 x * x + y * y <= effectiveRadius * effectiveRadius &&
                 std::abs(z - multipoleEnterPlane) >= spacing[2] &&
                 std::abs(z - multipoleExitPlane) >= spacing[2];
        };
    std::map<std::string, std::string> multipoleFields{
        {"axis_definition",
         "origin=center; +z=beam; +x,+y transverse; beam from -z to +z"},
        {"masked_zero_region", "true"},
        {"multipole_aperture_radius_cm", FormatDouble(multipoleRadius)},
        {"multipole_transition_planes_cm",
         FormatDouble(multipoleEnterPlane) + " " +
             FormatDouble(multipoleExitPlane)},
    };

    generateMagnet("Multipole", multipoleGrid, EvaluateMultipoleDirect,
                   multipoleValidationRegion,
                   false, false,
                   multipoleFields, multipolePath);

    const double dipoleGap = DeckValue(kDipoleElement, 13);
    const double dipoleRadius = DeckValue(kDipoleElement, 14);
    const double dipolePhi = DeckValue(kDipoleElement, 16);
    const double dipoleAlpha = DeckValue(kDipoleElement, 17);
    const double dipoleBeta = DeckValue(kDipoleElement, 18);
    const double dipoleZ11 = DeckValue(kDipoleElement, 25);
    const double dipoleZ12 = DeckValue(kDipoleElement, 26);
    const double dipoleZ21 = DeckValue(kDipoleElement, 27);
    const double dipoleZ22 = DeckValue(kDipoleElement, 28);
    const GridSpec dipoleGrid = BuildDipoleBounds();
    const ValidationRegion dipoleValidationRegion =
        [dipoleGap, dipoleRadius, dipolePhi, dipoleAlpha, dipoleBeta,
         dipoleZ11, dipoleZ12, dipoleZ21, dipoleZ22](
            double x, double y, double z,
            const std::array<double, 3>& spacing) {
          const double radialMargin = 2.0 * std::hypot(spacing[0], spacing[2]);
          return InsideDipoleRegion(
              x, y, z, dipoleGap, dipoleRadius, dipolePhi, dipoleAlpha,
              dipoleBeta, dipoleZ11, dipoleZ12, dipoleZ21, dipoleZ22,
              kDipoleStripHalfWidthCm, radialMargin, 2.0 * spacing[1],
              2.0 * spacing[2]);
        };
    std::map<std::string, std::string> dipoleFields{
        {"axis_definition",
         "origin=entrance_center; +x=left; +y=up; +z=incoming_beam; bend_toward=-x"},
        {"masked_zero_region", "true"},
        {"dipole_gap_cm", FormatDouble(dipoleGap)},
        {"dipole_center_x_cm", FormatDouble(-dipoleRadius)},
        {"dipole_center_z_cm", "0"},
        {"dipole_inner_radius_cm",
         FormatDouble(dipoleRadius - kDipoleStripHalfWidthCm)},
        {"dipole_outer_radius_cm",
         FormatDouble(dipoleRadius + kDipoleStripHalfWidthCm)},
        {"dipole_sector_angle_deg", FormatDouble(dipolePhi)},
        {"dipole_alpha_deg", FormatDouble(dipoleAlpha)},
        {"dipole_beta_deg", FormatDouble(dipoleBeta)},
        {"dipole_z11_cm", FormatDouble(dipoleZ11)},
        {"dipole_z12_cm", FormatDouble(dipoleZ12)},
        {"dipole_z21_cm", FormatDouble(dipoleZ21)},
        {"dipole_z22_cm", FormatDouble(dipoleZ22)},
        {"dipole_strip_half_width_cm", FormatDouble(kDipoleStripHalfWidthCm)},
    };

    generateMagnet("Dipole", dipoleGrid, EvaluateDipoleDirect,
                   dipoleValidationRegion,
                   false, true,
                   dipoleFields, dipolePath);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
}
