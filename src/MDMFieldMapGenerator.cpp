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
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "json.h"

namespace {

constexpr std::size_t kMultipoleElement = 4;
constexpr std::size_t kDipoleElement = 5;
constexpr double kDipoleStripHalfWidthCm = 30.0;
constexpr double kPi = 3.14159265358979323846;
constexpr double kSampleRoundTripToleranceTesla = 2.0e-5;

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
  std::filesystem::path outputDirectory = ".";
  std::filesystem::path multipoleOutput = "Multipole.bin";
  std::filesystem::path dipoleEntranceOutput = "DipoleEntrance.bin";
  std::filesystem::path dipoleSectorOutput = "DipoleSector.bin";
  std::filesystem::path dipoleExitOutput = "DipoleExit.bin";
};

using Evaluator = std::function<Vec3(double, double, double)>;

double DegreesToRadians(double degrees) { return degrees * kPi / 180.0; }

double DeckValue(std::size_t elementOneBased, std::size_t fieldOneBased) {
  return blck0_.DATA[elementOneBased - 1][fieldOneBased - 1];
}

double AbsDeckStep(const std::string& name, double value) {
  const double step = std::abs(value);
  if (step <= 0.0) {
    throw std::runtime_error("Invalid deck step '" + name +
                             "': expected non-zero positive magnitude");
  }
  return step;
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

std::size_t FlattenIndex(const GridSpec& grid,
                         std::size_t ix,
                         std::size_t iy,
                         std::size_t iz) {
  return ix + grid.counts[0] * (iy + grid.counts[1] * iz);
}

std::array<double, 3> GridSpacing(const GridSpec& grid) {
  return {(grid.max[0] - grid.min[0]) / static_cast<double>(grid.counts[0] - 1),
          (grid.max[1] - grid.min[1]) / static_cast<double>(grid.counts[1] - 1),
          (grid.max[2] - grid.min[2]) / static_cast<double>(grid.counts[2] - 1)};
}

double GridCoordinate(const GridSpec& grid, std::size_t axis, std::size_t index) {
  const auto spacing = GridSpacing(grid);
  return grid.min[axis] + spacing[axis] * static_cast<double>(index);
}

std::size_t CountFromRangeAndStep(double range, double targetStep) {
  if (range < 0.0) {
    throw std::runtime_error("Invalid grid range: max must be >= min");
  }
  if (targetStep <= 0.0) {
    throw std::runtime_error("Grid target step must be positive");
  }
  const std::size_t segments =
      static_cast<std::size_t>(std::ceil(range / targetStep));
  return std::max<std::size_t>(2, segments + 1);
}

GridSpec BuildGridFromStepTargets(const std::array<double, 3>& min,
                                  const std::array<double, 3>& max,
                                  const std::array<double, 3>& targetStep) {
  GridSpec grid;
  grid.min = min;
  grid.max = max;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    grid.counts[axis] =
        CountFromRangeAndStep(max[axis] - min[axis], targetStep[axis]);
  }
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

void EnsureParentDirectory(const std::filesystem::path& path) {
  const auto parent = path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
}

std::filesystem::path WithSuffixBeforeExtension(const std::filesystem::path& path,
                                                const std::string& suffix) {
  const std::string stem = path.stem().string();
  const std::string extension = path.extension().string();
  if (extension.empty()) {
    return path.parent_path() / (path.filename().string() + suffix);
  }
  return path.parent_path() / (stem + suffix + extension);
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

  if (config.isMember("outputDirectory")) {
    result.outputDirectory = config["outputDirectory"].asString();
  }
  if (config.isMember("multipoleOutput")) {
    result.multipoleOutput = config["multipoleOutput"].asString();
  }

  const bool hasSplitOutputs =
      config.isMember("dipoleEntranceOutput") ||
      config.isMember("dipoleSectorOutput") ||
      config.isMember("dipoleExitOutput");

  if (hasSplitOutputs) {
    if (config.isMember("dipoleEntranceOutput")) {
      result.dipoleEntranceOutput = config["dipoleEntranceOutput"].asString();
    }
    if (config.isMember("dipoleSectorOutput")) {
      result.dipoleSectorOutput = config["dipoleSectorOutput"].asString();
    }
    if (config.isMember("dipoleExitOutput")) {
      result.dipoleExitOutput = config["dipoleExitOutput"].asString();
    }
  } else if (config.isMember("dipoleOutput")) {
    const std::filesystem::path base = config["dipoleOutput"].asString();
    result.dipoleEntranceOutput = WithSuffixBeforeExtension(base, "Entrance");
    result.dipoleSectorOutput = WithSuffixBeforeExtension(base, "Sector");
    result.dipoleExitOutput = WithSuffixBeforeExtension(base, "Exit");
  }

  return result;
}

MDMFieldMap BuildMap(const std::string& magnetName,
                     const GridSpec& grid,
                     const Evaluator& evaluator,
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
  metadata.fields["payload_layout"] = "component_major_x_fastest_float32";
  metadata.fields["sampling_method"] = "direct_raytrace";
  metadata.fields["axis_definition"] = extraFields.at("axis_definition");
  metadata.fields["masked_zero_region"] = extraFields.at("masked_zero_region");
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

void VerifyDirectSamples(const MDMFieldMap& map,
                         const GridSpec& grid,
                         const Evaluator& evaluator,
                         const std::string& label) {
  const auto& bx = map.GetBx();
  const auto& by = map.GetBy();
  const auto& bz = map.GetBz();

  std::vector<std::array<std::size_t, 3>> samples{
      {0, 0, 0},
      {grid.counts[0] - 1, grid.counts[1] - 1, grid.counts[2] - 1},
      {grid.counts[0] / 2, grid.counts[1] / 2, grid.counts[2] / 2},
  };

  std::mt19937 rng(42);
  std::uniform_int_distribution<std::size_t> xDist(0, grid.counts[0] - 1);
  std::uniform_int_distribution<std::size_t> yDist(0, grid.counts[1] - 1);
  std::uniform_int_distribution<std::size_t> zDist(0, grid.counts[2] - 1);
  for (std::size_t i = 0; i < 8; ++i) {
    samples.push_back({xDist(rng), yDist(rng), zDist(rng)});
  }

  for (const auto& idx : samples) {
    const double x = GridCoordinate(grid, 0, idx[0]);
    const double y = GridCoordinate(grid, 1, idx[1]);
    const double z = GridCoordinate(grid, 2, idx[2]);
    const Vec3 expected = evaluator(x, y, z);
    const std::size_t flat = FlattenIndex(grid, idx[0], idx[1], idx[2]);
    const double dx = static_cast<double>(bx[flat]) - expected.x;
    const double dy = static_cast<double>(by[flat]) - expected.y;
    const double dz = static_cast<double>(bz[flat]) - expected.z;
    const double diff = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (diff > kSampleRoundTripToleranceTesla) {
      throw std::runtime_error(label + " direct-sample check failed");
    }
  }
}

void VerifyRoundTripPayload(const MDMFieldMap& original,
                            const MDMFieldMap& reloaded,
                            const std::string& label) {
  if (original.GetBx() != reloaded.GetBx() ||
      original.GetBy() != reloaded.GetBy() ||
      original.GetBz() != reloaded.GetBz()) {
    throw std::runtime_error(label + " payload round-trip mismatch");
  }
}

std::filesystem::path ResolveOutputPath(const std::filesystem::path& base,
                                        const std::filesystem::path& configured) {
  if (configured.is_relative()) {
    return base / configured;
  }
  return configured;
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

    const std::filesystem::path multipolePath =
        ResolveOutputPath(config.outputDirectory, config.multipoleOutput);
    const std::filesystem::path dipoleEntrancePath =
        ResolveOutputPath(config.outputDirectory, config.dipoleEntranceOutput);
    const std::filesystem::path dipoleSectorPath =
        ResolveOutputPath(config.outputDirectory, config.dipoleSectorOutput);
    const std::filesystem::path dipoleExitPath =
        ResolveOutputPath(config.outputDirectory, config.dipoleExitOutput);

    const auto writeMap = [&](const std::string& label,
                              const GridSpec& grid,
                              const Evaluator& evaluator,
                              const std::map<std::string, std::string>& metadata,
                              const std::filesystem::path& path) {
      MDMFieldMap map = BuildMap(label, grid, evaluator, metadata);
      VerifyDirectSamples(map, grid, evaluator, label);

      EnsureParentDirectory(path);
      map.Save(path.string());
      const MDMFieldMap reloaded = MDMFieldMap::Load(path.string());
      VerifyRoundTripPayload(map, reloaded, label);

      const auto spacing = GridSpacing(grid);
      std::cout << label << " grid " << FormatGridCounts(grid)
                << " spacing(cm)=[" << spacing[0] << ", " << spacing[1] << ", "
                << spacing[2] << "]" << std::endl;
      std::cout << "Wrote " << path << std::endl;
    };

    const double multipoleLF1 =
        AbsDeckStep("multipole LF1", DeckValue(kMultipoleElement, 1));
    const double multipoleLU1 =
        AbsDeckStep("multipole LU1", DeckValue(kMultipoleElement, 2));
    const double multipoleLF2 =
        AbsDeckStep("multipole LF2", DeckValue(kMultipoleElement, 3));
    const double multipoleStep =
        std::min({multipoleLF1, multipoleLU1, multipoleLF2});
    const double multipoleLength = DeckValue(kMultipoleElement, 12);
    const double multipoleRadius = DeckValue(kMultipoleElement, 13);
    const double multipoleZ11 = DeckValue(kMultipoleElement, 19);
    const double multipoleZ12 = DeckValue(kMultipoleElement, 20);
    const double multipoleZ21 = DeckValue(kMultipoleElement, 21);
    const double multipoleZ22 = DeckValue(kMultipoleElement, 22);
    const double multipoleZMin =
        -(multipoleLength / 2.0 + std::max(multipoleZ11, -multipoleZ12));
    const double multipoleZMax =
        +(multipoleLength / 2.0 + std::max(multipoleZ22, -multipoleZ21));
    const double multipoleEnterPlane = -0.5 * multipoleLength - multipoleZ12;
    const double multipoleExitPlane = 0.5 * multipoleLength + multipoleZ21;

    const GridSpec multipoleGrid = BuildGridFromStepTargets(
        {-multipoleRadius, -multipoleRadius, multipoleZMin},
        {multipoleRadius, multipoleRadius, multipoleZMax},
        {multipoleStep, multipoleStep, multipoleStep});

    std::map<std::string, std::string> multipoleFields{
        {"axis_definition",
         "origin=center; +z=beam; +x,+y transverse; beam from -z to +z"},
        {"masked_zero_region", "true"},
        {"multipole_aperture_radius_cm", FormatDouble(multipoleRadius)},
        {"multipole_transition_planes_cm",
         FormatDouble(multipoleEnterPlane) + " " +
             FormatDouble(multipoleExitPlane)},
        {"deck_lf1_cm", FormatDouble(multipoleLF1)},
        {"deck_lu1_cm", FormatDouble(multipoleLU1)},
        {"deck_lf2_cm", FormatDouble(multipoleLF2)},
        {"mdm_dipole_probe", FormatDouble(config.mdmDipoleProbe)},
        {"mdm_multipole_probe", FormatDouble(config.mdmMultipoleProbe)},
    };

    writeMap("Multipole", multipoleGrid, EvaluateMultipoleDirect, multipoleFields,
             multipolePath);

    const double dipoleLF1 =
        AbsDeckStep("dipole LF1", DeckValue(kDipoleElement, 1));
    const double dipoleLU1 =
        AbsDeckStep("dipole LU1", DeckValue(kDipoleElement, 2));
    const double dipoleLF2 =
        AbsDeckStep("dipole LF2", DeckValue(kDipoleElement, 3));
    const double dipoleDG = AbsDeckStep("dipole DG", DeckValue(kDipoleElement, 4));
    const double dipoleGap = DeckValue(kDipoleElement, 13);
    const double dipoleRadius = DeckValue(kDipoleElement, 14);
    const double dipolePhiDeg = DeckValue(kDipoleElement, 16);
    const double dipolePhiRad = DegreesToRadians(dipolePhiDeg);
    const double dipoleAlphaDeg = DeckValue(kDipoleElement, 17);
    const double dipoleBetaDeg = DeckValue(kDipoleElement, 18);
    const double dipoleZ11 = DeckValue(kDipoleElement, 25);
    const double dipoleZ12 = DeckValue(kDipoleElement, 26);
    const double dipoleZ21 = DeckValue(kDipoleElement, 27);
    const double dipoleZ22 = DeckValue(kDipoleElement, 28);

    const GridSpec dipoleEntranceGrid = BuildGridFromStepTargets(
        {-kDipoleStripHalfWidthCm, -dipoleGap / 2.0, dipoleZ12},
        {kDipoleStripHalfWidthCm, dipoleGap / 2.0, dipoleZ11},
        {dipoleDG, dipoleDG, dipoleLF1});
    const GridSpec dipoleSectorGrid = BuildGridFromStepTargets(
        {-kDipoleStripHalfWidthCm, -dipoleGap / 2.0, 0.0},
        {kDipoleStripHalfWidthCm, dipoleGap / 2.0, dipoleRadius * dipolePhiRad},
        {dipoleDG, dipoleDG, dipoleLU1});
    const GridSpec dipoleExitGrid = BuildGridFromStepTargets(
        {-kDipoleStripHalfWidthCm, -dipoleGap / 2.0, dipoleZ21},
        {kDipoleStripHalfWidthCm, dipoleGap / 2.0, dipoleZ22},
        {dipoleDG, dipoleDG, dipoleLF2});

    const Evaluator dipoleEntranceEvaluator =
        [dipoleAlphaDeg](double xB, double y, double zB) {
          const auto local = EntranceFrameToLocal(xB, zB, dipoleAlphaDeg);
          return EvaluateDipoleDirect(local[0], y, local[1]);
        };
    const Evaluator dipoleSectorEvaluator =
        [dipoleRadius](double dr, double y, double s) {
          const double theta = s / dipoleRadius;
          const double radius = dipoleRadius + dr;
          const double xLocal = -dipoleRadius + radius * std::cos(theta);
          const double zLocal = radius * std::sin(theta);
          return EvaluateDipoleDirect(xLocal, y, zLocal);
        };
    const Evaluator dipoleExitEvaluator =
        [dipoleRadius, dipolePhiDeg, dipoleAlphaDeg, dipoleBetaDeg](
            double xC, double y, double zC) {
          const auto local = ExitFrameToLocal(
              xC, zC, dipoleRadius, dipolePhiDeg, dipoleAlphaDeg, dipoleBetaDeg);
          return EvaluateDipoleDirect(local[0], y, local[1]);
        };

    const std::map<std::string, std::string> sharedDipoleFields{
        {"masked_zero_region", "false"},
        {"field_component_frame", "dipole_local_cartesian"},
        {"dipole_gap_cm", FormatDouble(dipoleGap)},
        {"dipole_reference_radius_cm", FormatDouble(dipoleRadius)},
        {"dipole_sector_angle_deg", FormatDouble(dipolePhiDeg)},
        {"dipole_alpha_deg", FormatDouble(dipoleAlphaDeg)},
        {"dipole_beta_deg", FormatDouble(dipoleBetaDeg)},
        {"dipole_z11_cm", FormatDouble(dipoleZ11)},
        {"dipole_z12_cm", FormatDouble(dipoleZ12)},
        {"dipole_z21_cm", FormatDouble(dipoleZ21)},
        {"dipole_z22_cm", FormatDouble(dipoleZ22)},
        {"dipole_strip_half_width_cm", FormatDouble(kDipoleStripHalfWidthCm)},
        {"deck_lf1_cm", FormatDouble(dipoleLF1)},
        {"deck_lu1_cm", FormatDouble(dipoleLU1)},
        {"deck_lf2_cm", FormatDouble(dipoleLF2)},
        {"deck_dg_cm", FormatDouble(dipoleDG)},
        {"mdm_dipole_probe", FormatDouble(config.mdmDipoleProbe)},
        {"mdm_multipole_probe", FormatDouble(config.mdmMultipoleProbe)},
    };

    {
      auto fields = sharedDipoleFields;
      fields["dipole_region"] = "entrance_fringe";
      fields["axis_definition"] =
          "coords=(xB,y,zB)_cm; fields=(Bx,By,Bz)_dipole_local";
      writeMap("DipoleEntrance", dipoleEntranceGrid, dipoleEntranceEvaluator,
               fields, dipoleEntrancePath);
    }
    {
      auto fields = sharedDipoleFields;
      fields["dipole_region"] = "sector";
      fields["axis_definition"] =
          "coords=(dr,y,s)_cm; dr=r-RB; s=RB*theta; fields=(Bx,By,Bz)_dipole_local";
      writeMap("DipoleSector", dipoleSectorGrid, dipoleSectorEvaluator, fields,
               dipoleSectorPath);
    }
    {
      auto fields = sharedDipoleFields;
      fields["dipole_region"] = "exit_fringe";
      fields["axis_definition"] =
          "coords=(xC,y,zC)_cm; fields=(Bx,By,Bz)_dipole_local";
      writeMap("DipoleExit", dipoleExitGrid, dipoleExitEvaluator, fields,
               dipoleExitPath);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
}
