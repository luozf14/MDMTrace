#pragma once

#include "json.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mdm {

struct RayInput {
  double xDeg = 0.0;
  double yDeg = 0.0;
  double energyMeV = 0.0;
};

inline double RequireNumber(const Json::Value& value,
                            const std::string& context) {
  if (!value.isNumeric()) {
    throw std::runtime_error(context + " must be numeric");
  }
  const double result = value.asDouble();
  if (!std::isfinite(result)) {
    throw std::runtime_error(context + " must be finite");
  }
  return result;
}

inline std::vector<double> BuildGridAxis(double minValue,
                                         double maxValue,
                                         double step,
                                         const std::string& gridName,
                                         const std::string& axisName) {
  if (step <= 0.0 || !std::isfinite(step)) {
    throw std::runtime_error(gridName + " " + axisName +
                             "Step must be positive");
  }
  if (!std::isfinite(minValue) || !std::isfinite(maxValue)) {
    throw std::runtime_error(gridName + " " + axisName +
                             " range must be finite");
  }
  if (maxValue < minValue) {
    throw std::runtime_error(gridName + " " + axisName + "Max must be >= " +
                             axisName + "Min");
  }

  std::vector<double> values;
  values.push_back(minValue);

  const double range = maxValue - minValue;
  const double tolerance = 1.0e-12 * std::max(1.0, std::fabs(range));
  if (range == 0.0) {
    values[0] = maxValue;
    return values;
  }

  double current = minValue + step;
  while (current < maxValue - tolerance) {
    values.push_back(current);
    current += step;
    if (values.size() > 1000000u) {
      throw std::runtime_error(gridName + " " + axisName +
                               " creates too many scan points");
    }
  }

  if (std::fabs(values.back() - maxValue) > tolerance) {
    values.push_back(maxValue);
  } else {
    values.back() = maxValue;
  }
  return values;
}

inline void AppendLegacyAngles(const Json::Value& value,
                               std::vector<RayInput>& rays) {
  if (!value.isArray()) {
    throw std::runtime_error("scatteredAngles must be an array");
  }
  for (Json::ArrayIndex index = 0; index < value.size(); ++index) {
    RayInput ray;
    ray.xDeg = RequireNumber(value[index], "scatteredAngles entry");
    ray.yDeg = 0.0;
    rays.push_back(ray);
  }
}

inline void AppendAnglePairs(const Json::Value& value,
                             std::vector<RayInput>& rays) {
  if (!value.isArray()) {
    throw std::runtime_error("scatteredAnglePairs must be an array");
  }
  for (Json::ArrayIndex index = 0; index < value.size(); ++index) {
    const Json::Value& pair = value[index];
    if (!pair.isArray() || pair.size() != 2u) {
      std::ostringstream message;
      message << "scatteredAnglePairs[" << index
              << "] must be [xDeg, yDeg]";
      throw std::runtime_error(message.str());
    }
    RayInput ray;
    ray.xDeg = RequireNumber(pair[0u], "scatteredAnglePairs x angle");
    ray.yDeg = RequireNumber(pair[1u], "scatteredAnglePairs y angle");
    rays.push_back(ray);
  }
}

inline void AppendAngleGrid(const Json::Value& value,
                            std::vector<RayInput>& rays) {
  if (!value.isObject()) {
    throw std::runtime_error("scatteredAngleGrid must be an object");
  }

  const double xMin = RequireNumber(value["xMin"], "scatteredAngleGrid xMin");
  const double xMax = RequireNumber(value["xMax"], "scatteredAngleGrid xMax");
  const double xStep =
      RequireNumber(value["xStep"], "scatteredAngleGrid xStep");
  const double yMin = RequireNumber(value["yMin"], "scatteredAngleGrid yMin");
  const double yMax = RequireNumber(value["yMax"], "scatteredAngleGrid yMax");
  const double yStep =
      RequireNumber(value["yStep"], "scatteredAngleGrid yStep");

  const std::vector<double> xValues =
      BuildGridAxis(xMin, xMax, xStep, "scatteredAngleGrid", "x");
  const std::vector<double> yValues =
      BuildGridAxis(yMin, yMax, yStep, "scatteredAngleGrid", "y");

  for (double xDeg : xValues) {
    for (double yDeg : yValues) {
      RayInput ray;
      ray.xDeg = xDeg;
      ray.yDeg = yDeg;
      rays.push_back(ray);
    }
  }
}

inline std::vector<double> ParseEnergies(const Json::Value& config) {
  if (!config.isMember("scatteredEnergyGrid")) {
    return {RequireNumber(config["scatteredEnergy"], "scatteredEnergy")};
  }

  const Json::Value& grid = config["scatteredEnergyGrid"];
  if (!grid.isObject()) {
    throw std::runtime_error("scatteredEnergyGrid must be an object");
  }

  const double eMin = RequireNumber(grid["eMin"], "scatteredEnergyGrid eMin");
  const double eMax = RequireNumber(grid["eMax"], "scatteredEnergyGrid eMax");
  const double eStep =
      RequireNumber(grid["eStep"], "scatteredEnergyGrid eStep");
  return BuildGridAxis(eMin, eMax, eStep, "scatteredEnergyGrid", "e");
}

inline std::vector<RayInput> ParseRayInputs(const Json::Value& config) {
  std::vector<RayInput> angleRays;
  if (config.isMember("scatteredAngles")) {
    AppendLegacyAngles(config["scatteredAngles"], angleRays);
  }
  if (config.isMember("scatteredAnglePairs")) {
    AppendAnglePairs(config["scatteredAnglePairs"], angleRays);
  }
  if (config.isMember("scatteredAngleGrid")) {
    AppendAngleGrid(config["scatteredAngleGrid"], angleRays);
  }
  if (angleRays.empty()) {
    throw std::runtime_error(
        "No input rays configured. Provide scatteredAngles, "
        "scatteredAnglePairs, or scatteredAngleGrid.");
  }

  const std::vector<double> energies = ParseEnergies(config);
  std::vector<RayInput> rays;
  rays.reserve(angleRays.size() * energies.size());
  for (const RayInput& angleRay : angleRays) {
    for (double energyMeV : energies) {
      RayInput ray = angleRay;
      ray.energyMeV = energyMeV;
      rays.push_back(ray);
    }
  }
  return rays;
}

}  // namespace mdm
