#include "MdmConfig.h"
#include "MdmIonConfig.h"
#include "MdmTrace.h"
#include "json.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr double kRigidityConstant = 0.299792458;
constexpr double kDipoleRadiusCm = 160.0;
constexpr double kStoppedPositionCm = 1.0e9;
constexpr double kStoppedPenalty = 1.0e60;
constexpr int kMaxBoundaryExpansions = 8;

struct Config {
  double mdmAngleDeg = 0.0;
  MdmIon ion;
  double energyMeV = 0.0;
  double inputAngleXDeg = 0.0;
  double inputAngleYDeg = 0.0;
  double positionScaleCm = 0.1;
  double angleScaleDeg = 0.1;
  double searchHalfWidthFraction = 0.25;
  int coarseSamples = 101;
  double fieldToleranceGauss = 1.0e-3;
  int maxIterations = 100;
};

struct RayResult {
  double fieldGauss = 0.0;
  double xCm = 0.0;
  double yCm = 0.0;
  double angXDeg = 0.0;
  double angYDeg = 0.0;
  double chi2 = kStoppedPenalty;
  bool stopped = true;
};

double GetDouble(const Json::Value& value,
                 const char* key,
                 double fallback = 0.0) {
  return value.isObject() && value.isMember(key) ? value[key].asDouble()
                                                 : fallback;
}

Config ParseConfig(const Json::Value& json) {
  Config cfg;
  cfg.mdmAngleDeg = GetDouble(json, "mdmAngle");
  cfg.ion = mdm::ParseScatteredIon(json);
  cfg.energyMeV = GetDouble(json, "scatteredEnergy");

  const Json::Value& finder = json["fieldFinder"];
  if (!finder.isObject()) {
    throw std::runtime_error("fieldFinder must be an object");
  }
  cfg.inputAngleXDeg = GetDouble(finder, "inputAngleXDeg", 0.0);
  cfg.inputAngleYDeg = GetDouble(finder, "inputAngleYDeg", 0.0);
  cfg.positionScaleCm = GetDouble(finder, "positionScaleCm", 0.1);
  cfg.angleScaleDeg = GetDouble(finder, "angleScaleDeg", 0.1);
  cfg.searchHalfWidthFraction =
      GetDouble(finder, "searchHalfWidthFraction", 0.25);
  cfg.coarseSamples = static_cast<int>(
      mdm::GetNonnegativeInteger(finder, "coarseSamples", 101));
  cfg.fieldToleranceGauss =
      GetDouble(finder, "fieldToleranceGauss", 1.0e-3);
  cfg.maxIterations = static_cast<int>(
      mdm::GetNonnegativeInteger(finder, "maxIterations", 100));
  return cfg;
}

void CheckConfig(const Config& cfg) {
  if (cfg.ion.ionMassMeV <= 0.0 || cfg.ion.chargeState <= 0 ||
      cfg.energyMeV <= 0.0) {
    throw std::runtime_error("scatteredIon and scatteredEnergy must be "
                             "physical");
  }
  if (cfg.positionScaleCm <= 0.0 || cfg.angleScaleDeg <= 0.0) {
    throw std::runtime_error("fieldFinder scales must be positive");
  }
  if (cfg.searchHalfWidthFraction <= 0.0) {
    throw std::runtime_error("fieldFinder.searchHalfWidthFraction must be "
                             "positive");
  }
  if (cfg.coarseSamples < 3) {
    throw std::runtime_error("fieldFinder.coarseSamples must be at least 3");
  }
  if (cfg.fieldToleranceGauss <= 0.0) {
    throw std::runtime_error(
        "fieldFinder.fieldToleranceGauss must be positive");
  }
}

double InitialFieldGauss(const Config& cfg) {
  const double momentumMeVPerC =
      std::sqrt((2.0 * cfg.ion.ionMassMeV + cfg.energyMeV) * cfg.energyMeV);
  const double bRhoKgCm =
      momentumMeVPerC / (kRigidityConstant * cfg.ion.chargeState);
  return 1000.0 * bRhoKgCm / kDipoleRadiusCm;
}

RayResult TraceAtField(const Config& cfg, double fieldGauss) {
  MdmTrace trace;
  trace.SetMdmAngle(cfg.mdmAngleDeg);
  trace.SetMdmDipoleField(fieldGauss);
  trace.SetScatteredIon(cfg.ion);
  trace.SetScatteredEnergy(cfg.energyMeV);
  trace.SetScatteredAngle(cfg.inputAngleXDeg, cfg.inputAngleYDeg);
  trace.SendRay();

  RayResult result;
  result.fieldGauss = fieldGauss;
  trace.GetPositionAngleFirstWire(result.xCm, result.yCm, result.angXDeg,
                                  result.angYDeg);
  result.stopped = !std::isfinite(result.xCm) || !std::isfinite(result.yCm) ||
                   !std::isfinite(result.angXDeg) ||
                   !std::isfinite(result.angYDeg) ||
                   std::abs(result.xCm) > kStoppedPositionCm ||
                   std::abs(result.yCm) > kStoppedPositionCm;
  if (!result.stopped) {
    const double x = result.xCm / cfg.positionScaleCm;
    const double y = result.yCm / cfg.positionScaleCm;
    const double ax = result.angXDeg / cfg.angleScaleDeg;
    const double ay = result.angYDeg / cfg.angleScaleDeg;
    result.chi2 = x * x + y * y + ax * ax + ay * ay;
  }
  return result;
}

RayResult CoarseScan(const Config& cfg,
                     double lowField,
                     double highField,
                     int& bestIndex) {
  RayResult best;
  bestIndex = 0;
  for (int i = 0; i < cfg.coarseSamples; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(cfg.coarseSamples - 1);
    const double field = lowField + fraction * (highField - lowField);
    const RayResult result = TraceAtField(cfg, field);
    if (result.chi2 < best.chi2) {
      best = result;
      bestIndex = i;
    }
  }
  return best;
}

RayResult GoldenSectionSearch(const Config& cfg,
                              double lowField,
                              double highField,
                              int& iterations,
                              bool& converged) {
  const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
  double a = lowField;
  double b = highField;
  double c = b - gr * (b - a);
  double d = a + gr * (b - a);
  RayResult fc = TraceAtField(cfg, c);
  RayResult fd = TraceAtField(cfg, d);

  iterations = 0;
  while (std::abs(b - a) > cfg.fieldToleranceGauss &&
         iterations < cfg.maxIterations) {
    if (fc.chi2 < fd.chi2) {
      b = d;
      d = c;
      fd = fc;
      c = b - gr * (b - a);
      fc = TraceAtField(cfg, c);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + gr * (b - a);
      fd = TraceAtField(cfg, d);
    }
    ++iterations;
  }
  converged = std::abs(b - a) <= cfg.fieldToleranceGauss;
  return fc.chi2 < fd.chi2 ? fc : fd;
}

RayResult TuneField(const Config& cfg,
                    double initialField,
                    int& coarseEvaluations,
                    int& refinementIterations,
                    int& expansions,
                    bool& converged) {
  double halfWidth = cfg.searchHalfWidthFraction * initialField;
  double low = std::max(1.0, initialField - halfWidth);
  double high = initialField + halfWidth;

  int bestIndex = 0;
  RayResult coarseBest = CoarseScan(cfg, low, high, bestIndex);
  coarseEvaluations = cfg.coarseSamples;
  expansions = 0;
  while ((bestIndex == 0 || bestIndex == cfg.coarseSamples - 1) &&
         expansions < kMaxBoundaryExpansions) {
    halfWidth *= 2.0;
    low = std::max(1.0, initialField - halfWidth);
    high = initialField + halfWidth;
    coarseBest = CoarseScan(cfg, low, high, bestIndex);
    coarseEvaluations += cfg.coarseSamples;
    ++expansions;
  }

  const double coarseStep = (high - low) / (cfg.coarseSamples - 1);
  const double bracketLow =
      std::max(1.0, coarseBest.fieldGauss - coarseStep);
  const double bracketHigh = coarseBest.fieldGauss + coarseStep;
  RayResult refined = GoldenSectionSearch(cfg, bracketLow, bracketHigh,
                                          refinementIterations, converged);
  return refined.chi2 < coarseBest.chi2 ? refined : coarseBest;
}

void PrintResult(const Config& cfg,
                 double initialField,
                 const RayResult& initial,
                 const RayResult& best,
                 int coarseEvaluations,
                 int refinementIterations,
                 int expansions,
                 bool converged) {
  const double dipoleProbe = best.fieldGauss / 1.034;
  const double multipoleProbe = dipoleProbe * 0.71;

  std::printf("Initial rigidity field: %.12f G\n", initialField);
  std::printf("Ion: A=%d Z=%d q=%d neutralMass=%.12f u ionMass=%.8f MeV\n",
              cfg.ion.massNumber, cfg.ion.atomicNumber, cfg.ion.chargeState,
              cfg.ion.neutralMassU, cfg.ion.ionMassMeV);
  std::printf("Initial residual: %.8e\n", initial.chi2);
  std::printf("Best mdmDipoleField: %.12f G\n", best.fieldGauss);
  std::printf("Best mdmDipoleProbe: %.12f G\n", dipoleProbe);
  std::printf("Best mdmMultipoleProbe: %.12f G\n", multipoleProbe);
  std::printf("Final X1 residual: %.12g cm\n", best.xCm);
  std::printf("Final Y1 residual: %.12g cm\n", best.yCm);
  std::printf("Final AngX1 residual: %.12g deg\n", best.angXDeg);
  std::printf("Final AngY1 residual: %.12g deg\n", best.angYDeg);
  std::printf("Final weighted objective: %.12g\n", best.chi2);
  std::printf("Convergence status: %s\n",
              converged ? "converged" : "maximum iterations reached");
  std::printf("Iteration count: %d\n", refinementIterations);
  std::printf("Improvement factor: %.8e\n", initial.chi2 / best.chi2);
  std::printf("Input ray: AngX %.8f deg, AngY %.8f deg, Energy %.8f MeV\n",
              cfg.inputAngleXDeg, cfg.inputAngleYDeg, cfg.energyMeV);
  std::printf("Search summary: coarse evaluations %d, boundary expansions %d, "
              "refinement iterations %d\n",
              coarseEvaluations, expansions, refinementIterations);
}

}  // namespace

int Run(int argc, char* argv[]) {
  const std::string configPath =
      argc > 1 ? argv[1] : "../config/FindField.json";
  const Config cfg = ParseConfig(mdm::ReadConfig(configPath));
  CheckConfig(cfg);

  const double initialField = InitialFieldGauss(cfg);
  const RayResult initial = TraceAtField(cfg, initialField);
  int coarseEvaluations = 0;
  int refinementIterations = 0;
  int expansions = 0;
  bool converged = false;
  const RayResult best = TuneField(cfg, initialField, coarseEvaluations,
                                  refinementIterations, expansions, converged);
  if (best.stopped) {
    std::cerr << "ERROR: best ray is stopped; no valid field found\n";
    return 1;
  }
  PrintResult(cfg, initialField, initial, best, coarseEvaluations,
              refinementIterations, expansions, converged);
  return 0;
}

int main(int argc, char* argv[]) {
  try {
    return Run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
}
