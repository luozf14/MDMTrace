#ifndef MDM_RAY_SCAN_H
#define MDM_RAY_SCAN_H

#include "json.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mdm {

struct RayAngle {
  double xDeg = 0.0;
  double yDeg = 0.0;
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
                                         const std::string& axisName) {
  if (step <= 0.0 || !std::isfinite(step)) {
    throw std::runtime_error("scatteredAngleGrid " + axisName +
                             "Step must be positive");
  }
  if (!std::isfinite(minValue) || !std::isfinite(maxValue)) {
    throw std::runtime_error("scatteredAngleGrid " + axisName +
                             " range must be finite");
  }
  if (maxValue < minValue) {
    throw std::runtime_error("scatteredAngleGrid " + axisName +
                             "Max must be >= " + axisName + "Min");
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
      throw std::runtime_error("scatteredAngleGrid " + axisName +
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
                               std::vector<RayAngle>& rays) {
  if (!value.isArray()) {
    throw std::runtime_error("scatteredAngles must be an array");
  }
  for (Json::ArrayIndex index = 0; index < value.size(); ++index) {
    RayAngle ray;
    ray.xDeg = RequireNumber(value[index], "scatteredAngles entry");
    ray.yDeg = 0.0;
    rays.push_back(ray);
  }
}

inline void AppendAnglePairs(const Json::Value& value,
                             std::vector<RayAngle>& rays) {
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
    RayAngle ray;
    ray.xDeg = RequireNumber(pair[0u], "scatteredAnglePairs x angle");
    ray.yDeg = RequireNumber(pair[1u], "scatteredAnglePairs y angle");
    rays.push_back(ray);
  }
}

inline void AppendAngleGrid(const Json::Value& value,
                            std::vector<RayAngle>& rays) {
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

  const std::vector<double> xValues = BuildGridAxis(xMin, xMax, xStep, "x");
  const std::vector<double> yValues = BuildGridAxis(yMin, yMax, yStep, "y");

  for (double xDeg : xValues) {
    for (double yDeg : yValues) {
      RayAngle ray;
      ray.xDeg = xDeg;
      ray.yDeg = yDeg;
      rays.push_back(ray);
    }
  }
}

inline std::vector<RayAngle> ParseRayAngles(const Json::Value& config) {
  std::vector<RayAngle> rays;
  if (config.isMember("scatteredAngles")) {
    AppendLegacyAngles(config["scatteredAngles"], rays);
  }
  if (config.isMember("scatteredAnglePairs")) {
    AppendAnglePairs(config["scatteredAnglePairs"], rays);
  }
  if (config.isMember("scatteredAngleGrid")) {
    AppendAngleGrid(config["scatteredAngleGrid"], rays);
  }
  if (rays.empty()) {
    throw std::runtime_error(
        "No input rays configured. Provide scatteredAngles, "
        "scatteredAnglePairs, or scatteredAngleGrid.");
  }
  return rays;
}

}  // namespace mdm

#endif
