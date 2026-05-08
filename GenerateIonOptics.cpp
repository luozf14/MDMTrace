#include "MDMFieldMapTrace.h"
#include "json.h"

#include <TDecompSVD.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifndef MDMTRACE_SOURCE_DIR
#define MDMTRACE_SOURCE_DIR "."
#endif

namespace {

constexpr double kAmuMeV = 931.48;
constexpr double kMradPerDegree = 1000.0 * 3.14159265358979323846 / 180.0;
constexpr double kTraceMradPerDegree = 17.453;
constexpr double kStoppedCm = 1.0e9;
constexpr double kPercentToFraction = 0.01;
constexpr std::size_t kInputDim = 5;
constexpr std::size_t kOutputDim = 5;
constexpr std::size_t kPhysicsOutputDim = 4;

struct Config {
  bool usingProbe = false;
  double mdmAngleDeg = 0.0;
  double dipoleField = 0.0;
  double dipoleProbe = 0.0;
  double multipoleProbe = 0.0;
  double massAmu = 0.0;
  double charge = 0.0;
  std::string multipoleMap = "Multipole.bin";
  std::string dipoleEntranceMap = "DipoleEntrance.bin";
  std::string dipoleSectorMap = "DipoleSector.bin";
  std::string dipoleExitMap = "DipoleExit.bin";
};

struct PhaseSpace {
  double xMm = 0.0;
  double thetaXMrad = 0.0;
  double yMm = 0.0;
  double thetaYMrad = 0.0;
  double deltaPercent = 0.0;
};

struct Axis {
  double min = 0.0;
  double max = 0.0;
  double step = 1.0;
  std::vector<double> values;
};

struct FitGrid {
  Axis xMm;
  Axis thetaXMrad;
  Axis yMm;
  Axis thetaYMrad;
  Axis deltaPercent;
};

struct IonOpticsConfig {
  std::string method = "fit";
  unsigned int order = 2;
  PhaseSpace reference;
  double referenceEnergyMeV = 0.0;
  FitGrid grid;
  unsigned int requestedThreads = 0;
  unsigned int maxRays = 1000000;
  std::string outputPath = "IonOpticsMatrix.txt";
};

struct TraceResult {
  bool stopped = false;
  std::array<double, kOutputDim> output{};
};

struct FitStats {
  std::size_t total = 0;
  std::size_t accepted = 0;
  std::size_t stopped = 0;
  std::size_t basisTerms = 0;
  unsigned int threads = 1;
  std::string firstSkippedReason;
  std::array<double, 4> rmsResidual{};
  std::array<double, 4> maxAbsResidual{};
};

struct FitSample {
  std::array<double, kInputDim> inputDelta{};
  std::array<double, kOutputDim> outputDelta{};
};

struct SampleSlot {
  bool accepted = false;
  bool stopped = false;
  FitSample sample;
};

using Matrix5 = std::array<std::array<double, kInputDim>, kOutputDim>;
using BasisVector = std::vector<double>;

struct QuadraticTerm {
  std::size_t first = 0;
  std::size_t second = 0;
};

struct TransferMap {
  unsigned int order = 1;
  Matrix5 first{};
  std::array<std::array<double, 15>, kOutputDim> second{};
};

Json::Value ReadJson(const std::string& path) {
  std::ifstream stream(path.c_str());
  if (!stream) {
    throw std::runtime_error("Cannot open config: " + path);
  }
  Json::Value config;
  stream >> config;
  return config;
}

double GetDouble(const Json::Value& value,
                 const char* key,
                 double fallback = 0.0) {
  return value.isObject() && value.isMember(key) ? value[key].asDouble()
                                                 : fallback;
}

std::string GetString(const Json::Value& value,
                      const char* key,
                      const char* fallback) {
  return value.isObject() && value.isMember(key) ? value[key].asString()
                                                 : fallback;
}

unsigned int GetUnsigned(const Json::Value& value,
                         const char* key,
                         unsigned int fallback = 0) {
  if (!value.isObject() || !value.isMember(key)) {
    return fallback;
  }
  const int raw = value[key].asInt();
  if (raw < 0) {
    throw std::runtime_error(std::string(key) + " must be non-negative");
  }
  return static_cast<unsigned int>(raw);
}

std::vector<double> BuildAxis(double minValue,
                              double maxValue,
                              double step,
                              const std::string& name) {
  if (step <= 0.0 || maxValue < minValue) {
    throw std::runtime_error("invalid ionOptics.fitGrid axis: " + name);
  }
  std::vector<double> values;
  values.push_back(minValue);
  const double range = maxValue - minValue;
  const double tolerance = 1.0e-12 * std::max(1.0, std::abs(range));
  double current = minValue + step;
  while (current < maxValue - tolerance) {
    values.push_back(current);
    current += step;
  }
  if (std::abs(values.back() - maxValue) > tolerance) {
    values.push_back(maxValue);
  } else {
    values.back() = maxValue;
  }
  return values;
}

Axis ParseAxis(const Json::Value& grid,
               const char* minKey,
               const char* maxKey,
               const char* stepKey,
               double defaultMin,
               double defaultMax,
               double defaultStep,
               const std::string& name) {
  Axis axis;
  axis.min = GetDouble(grid, minKey, defaultMin);
  axis.max = GetDouble(grid, maxKey, defaultMax);
  axis.step = GetDouble(grid, stepKey, defaultStep);
  axis.values = BuildAxis(axis.min, axis.max, axis.step, name);
  return axis;
}

Config ParseConfig(const Json::Value& json) {
  Config cfg;
  cfg.usingProbe = json.isMember("usingProbe") && json["usingProbe"].asBool();
  cfg.mdmAngleDeg = GetDouble(json, "mdmAngle");
  cfg.dipoleField = GetDouble(json, "mdmDipoleField");
  cfg.dipoleProbe = GetDouble(json, "mdmDipoleProbe");
  cfg.multipoleProbe = GetDouble(json, "mdmMultipoleProbe");
  cfg.massAmu = GetDouble(json, "scatteredMass");
  cfg.charge = GetDouble(json, "scatteredCharge");
  cfg.multipoleMap = GetString(json, "multipoleMapPath", "Multipole.bin");
  cfg.dipoleEntranceMap =
      GetString(json, "dipoleEntranceMapPath", "DipoleEntrance.bin");
  cfg.dipoleSectorMap =
      GetString(json, "dipoleSectorMapPath", "DipoleSector.bin");
  cfg.dipoleExitMap = GetString(json, "dipoleExitMapPath", "DipoleExit.bin");
  return cfg;
}

IonOpticsConfig ParseIonOpticsConfig(const Json::Value& json) {
  IonOpticsConfig cfg;
  cfg.referenceEnergyMeV = GetDouble(json, "scatteredEnergy");

  const Json::Value& optics = json["ionOptics"];
  cfg.method = GetString(optics, "method", "fit");
  cfg.order = GetUnsigned(optics, "order", 2);
  const Json::Value& reference = optics["reference"];
  if (reference.isObject()) {
    cfg.reference.xMm = GetDouble(reference, "xMm");
    cfg.reference.thetaXMrad = GetDouble(reference, "thetaXMrad");
    cfg.reference.yMm = GetDouble(reference, "yMm");
    cfg.reference.thetaYMrad = GetDouble(reference, "thetaYMrad");
    cfg.referenceEnergyMeV =
        GetDouble(reference, "energyMeV", cfg.referenceEnergyMeV);
  }

  const Json::Value& grid = optics["fitGrid"];
  cfg.grid.xMm =
      ParseAxis(grid, "xMinMm", "xMaxMm", "xStepMm", -2.0, 2.0, 1.0, "x");
  cfg.grid.thetaXMrad =
      ParseAxis(grid, "thetaXMinMrad", "thetaXMaxMrad", "thetaXStepMrad",
                -5.0, 5.0, 2.5, "thetaX");
  cfg.grid.yMm =
      ParseAxis(grid, "yMinMm", "yMaxMm", "yStepMm", -2.0, 2.0, 1.0, "y");
  cfg.grid.thetaYMrad =
      ParseAxis(grid, "thetaYMinMrad", "thetaYMaxMrad", "thetaYStepMrad",
                -5.0, 5.0, 2.5, "thetaY");
  cfg.grid.deltaPercent =
      ParseAxis(grid, "deltaMin", "deltaMax", "deltaStep", -0.2, 0.2, 0.2,
                "deltaP/P0 percent");
  cfg.requestedThreads = GetUnsigned(optics, "threads", 0);
  cfg.maxRays = GetUnsigned(optics, "maxRays", cfg.maxRays);
  cfg.outputPath = GetString(optics, "outputPath", cfg.outputPath.c_str());
  return cfg;
}

std::size_t GridRayCount(const FitGrid& grid);

void CheckConfig(const Config& cfg, const IonOpticsConfig& optics) {
  if (optics.method != "fit") {
    throw std::runtime_error("GenerateIonOptics only supports method=fit");
  }
  if (optics.order < 1 || optics.order > 2) {
    throw std::runtime_error("ionOptics.order must be 1 or 2");
  }
  if (cfg.massAmu <= 0.0) {
    throw std::runtime_error("scatteredMass must be positive");
  }
  if (cfg.charge == 0.0) {
    throw std::runtime_error("scatteredCharge must be non-zero");
  }
  if (optics.referenceEnergyMeV <= 0.0) {
    throw std::runtime_error("reference energy must be positive");
  }
  const std::size_t totalRays = GridRayCount(optics.grid);
  if (totalRays > optics.maxRays) {
    throw std::runtime_error("ionOptics.fitGrid creates " +
                             std::to_string(totalRays) +
                             " rays, exceeding ionOptics.maxRays=" +
                             std::to_string(optics.maxRays));
  }
}

void ConfigureTrace(MDMFieldMapTrace& trace, const Config& cfg) {
  trace.SetMDMAngle(cfg.mdmAngleDeg);
  cfg.usingProbe ? trace.SetMDMProbe(cfg.dipoleProbe, cfg.multipoleProbe)
                 : trace.SetMDMDipoleField(cfg.dipoleField);
  trace.SetScatteredMass(cfg.massAmu);
  trace.SetScatteredCharge(cfg.charge);
  trace.LoadFieldMaps(cfg.multipoleMap, cfg.dipoleEntranceMap,
                      cfg.dipoleSectorMap, cfg.dipoleExitMap);
}

double EnergyFromDeltaPercent(double referenceEnergyMeV,
                              double massAmu,
                              double deltaPercent) {
  const double massMeV = massAmu * kAmuMeV;
  const double p0 =
      std::sqrt((2.0 * massMeV + referenceEnergyMeV) * referenceEnergyMeV);
  const double p = p0 * (1.0 + deltaPercent * kPercentToFraction);
  return std::sqrt(p * p + massMeV * massMeV) - massMeV;
}

TraceResult TraceRay(MDMFieldMapTrace& trace,
                     const Config& cfg,
                     const IonOpticsConfig& optics,
                     const PhaseSpace& input) {
  trace.SetInitialPosition(input.xMm / 10.0, input.yMm / 10.0);
  trace.SetScatteredAngle(cfg.mdmAngleDeg + input.thetaXMrad /
                                             kTraceMradPerDegree,
                          input.thetaYMrad / kTraceMradPerDegree);
  trace.SetScatteredEnergy(
      EnergyFromDeltaPercent(optics.referenceEnergyMeV, cfg.massAmu,
                             input.deltaPercent));
  trace.SendRay();

  double xCm = 0.0;
  double yCm = 0.0;
  double angXDeg = 0.0;
  double angYDeg = 0.0;
  trace.GetPositionAngleFirstWire(xCm, yCm, angXDeg, angYDeg);

  TraceResult result;
  result.stopped = !std::isfinite(xCm) || !std::isfinite(yCm) ||
                   !std::isfinite(angXDeg) || !std::isfinite(angYDeg) ||
                   std::abs(xCm) > kStoppedCm || std::abs(yCm) > kStoppedCm;
  result.output = {10.0 * xCm, angXDeg * kMradPerDegree, 10.0 * yCm,
                   angYDeg * kMradPerDegree, input.deltaPercent};
  return result;
}

std::array<double, kInputDim> InputDelta(const PhaseSpace& input,
                                         const PhaseSpace& reference) {
  return {input.xMm - reference.xMm,
          input.thetaXMrad - reference.thetaXMrad,
          input.yMm - reference.yMm,
          input.thetaYMrad - reference.thetaYMrad,
          input.deltaPercent - reference.deltaPercent};
}

std::vector<QuadraticTerm> BuildQuadraticTerms() {
  std::vector<QuadraticTerm> terms;
  terms.reserve(15);
  for (std::size_t first = 0; first < kInputDim; ++first) {
    for (std::size_t second = first; second < kInputDim; ++second) {
      terms.push_back({first, second});
    }
  }
  return terms;
}

std::vector<const char*> InputLabels() {
  return {"x[mm]", "thetaX[mrad]", "y[mm]", "thetaY[mrad]",
          "deltaP/P0[%]"};
}

std::vector<std::string> QuadraticLabels() {
  const std::vector<const char*> labels = InputLabels();
  std::vector<std::string> result;
  for (const QuadraticTerm& term : BuildQuadraticTerms()) {
    result.push_back(std::string(labels[term.first]) + "*" +
                     labels[term.second]);
  }
  return result;
}

BasisVector BuildBasis(const std::array<double, kInputDim>& inputDelta,
                       unsigned int order) {
  BasisVector basis;
  basis.reserve(order >= 2 ? 20 : 5);
  for (double value : inputDelta) {
    basis.push_back(value);
  }
  if (order >= 2) {
    for (const QuadraticTerm& term : BuildQuadraticTerms()) {
      basis.push_back(inputDelta[term.first] * inputDelta[term.second]);
    }
  }
  return basis;
}

std::array<double, kOutputDim> EvaluateMap(
    const TransferMap& map,
    const std::array<double, kInputDim>& inputDelta) {
  std::array<double, kOutputDim> result{};
  for (std::size_t row = 0; row < kOutputDim; ++row) {
    for (std::size_t col = 0; col < kInputDim; ++col) {
      result[row] += map.first[row][col] * inputDelta[col];
    }
  }
  if (map.order >= 2) {
    const std::vector<QuadraticTerm> terms = BuildQuadraticTerms();
    for (std::size_t row = 0; row < kOutputDim; ++row) {
      for (std::size_t termIndex = 0; termIndex < terms.size(); ++termIndex) {
        const QuadraticTerm& term = terms[termIndex];
        result[row] += map.second[row][termIndex] * inputDelta[term.first] *
                       inputDelta[term.second];
      }
    }
  }
  return result;
}

std::size_t GridRayCount(const FitGrid& grid) {
  return grid.xMm.values.size() * grid.thetaXMrad.values.size() *
         grid.yMm.values.size() * grid.thetaYMrad.values.size() *
         grid.deltaPercent.values.size();
}

std::vector<PhaseSpace> BuildInputGrid(const FitGrid& grid) {
  std::vector<PhaseSpace> inputs;
  inputs.reserve(GridRayCount(grid));
  for (double xMm : grid.xMm.values) {
    for (double thetaXMrad : grid.thetaXMrad.values) {
      for (double yMm : grid.yMm.values) {
        for (double thetaYMrad : grid.thetaYMrad.values) {
          for (double deltaPercent : grid.deltaPercent.values) {
            inputs.push_back(
                {xMm, thetaXMrad, yMm, thetaYMrad, deltaPercent});
          }
        }
      }
    }
  }
  return inputs;
}

unsigned int ResolveThreadCount(unsigned int requestedThreads,
                                std::size_t totalRays) {
  if (totalRays == 0) {
    return 1;
  }
  unsigned int threads = requestedThreads;
  if (threads == 0) {
    threads = std::thread::hardware_concurrency();
    if (threads == 0) {
      threads = 1;
    }
  }
  return std::min<std::size_t>(threads, totalRays);
}

std::vector<FitSample> TraceGrid(const std::vector<PhaseSpace>& inputs,
                                 const Config& cfg,
                                 const IonOpticsConfig& optics,
                                 const std::array<double, kOutputDim>& referenceOutput,
                                 FitStats& stats) {
  stats.total = inputs.size();
  stats.threads = ResolveThreadCount(optics.requestedThreads, inputs.size());

  std::vector<std::unique_ptr<MDMFieldMapTrace>> traces;
  traces.reserve(stats.threads);
  for (unsigned int worker = 0; worker < stats.threads; ++worker) {
    auto trace = std::make_unique<MDMFieldMapTrace>();
    ConfigureTrace(*trace, cfg);
    const TraceResult warmup = TraceRay(*trace, cfg, optics, optics.reference);
    if (warmup.stopped) {
      throw std::runtime_error("worker reference ray stopped before fit tracing");
    }
    traces.push_back(std::move(trace));
  }

  std::vector<SampleSlot> slots(inputs.size());
  std::mutex skippedReasonMutex;
  const auto recordSkippedReason = [&](const std::string& reason) {
    std::lock_guard<std::mutex> lock(skippedReasonMutex);
    if (stats.firstSkippedReason.empty()) {
      stats.firstSkippedReason = reason;
    }
  };
  const auto runRange = [&](MDMFieldMapTrace& trace, std::size_t begin,
                            std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      TraceResult traceResult;
      try {
        traceResult = TraceRay(trace, cfg, optics, inputs[i]);
      } catch (const std::exception& error) {
        recordSkippedReason(error.what());
        slots[i].stopped = true;
        continue;
      }
      if (traceResult.stopped) {
        recordSkippedReason("sentinel or non-finite output");
        slots[i].stopped = true;
        continue;
      }

      const std::array<double, kInputDim> inputDelta =
          InputDelta(inputs[i], optics.reference);
      std::array<double, kOutputDim> outputDelta{};
      for (std::size_t row = 0; row < kOutputDim; ++row) {
        outputDelta[row] = traceResult.output[row] - referenceOutput[row];
      }
      slots[i].accepted = true;
      slots[i].sample = {inputDelta, outputDelta};
    }
  };

  if (stats.threads == 1) {
    runRange(*traces.front(), 0, inputs.size());
  } else {
    std::exception_ptr threadException;
    std::mutex exceptionMutex;
    std::vector<std::thread> workers;
    workers.reserve(stats.threads);

    for (unsigned int worker = 0; worker < stats.threads; ++worker) {
      const std::size_t begin = inputs.size() * worker / stats.threads;
      const std::size_t end = inputs.size() * (worker + 1) / stats.threads;
      workers.emplace_back([&, worker, begin, end]() {
        try {
          runRange(*traces[worker], begin, end);
        } catch (...) {
          std::lock_guard<std::mutex> lock(exceptionMutex);
          if (!threadException) {
            threadException = std::current_exception();
          }
        }
      });
    }

    for (std::thread& worker : workers) {
      worker.join();
    }
    if (threadException) {
      std::rethrow_exception(threadException);
    }
  }

  std::vector<FitSample> samples;
  samples.reserve(inputs.size());
  for (const SampleSlot& slot : slots) {
    if (slot.accepted) {
      samples.push_back(slot.sample);
    } else if (slot.stopped) {
      ++stats.stopped;
    }
  }
  stats.accepted = samples.size();
  return samples;
}

TransferMap SolveMapWithRootSvd(const std::vector<FitSample>& samples,
                                unsigned int order) {
  const std::size_t basisTerms = BuildBasis(samples.front().inputDelta, order).size();
  TMatrixD normal(static_cast<Int_t>(basisTerms), static_cast<Int_t>(basisTerms));
  std::array<std::vector<double>, kOutputDim> rhs;
  for (std::vector<double>& row : rhs) {
    row.assign(basisTerms, 0.0);
  }

  for (const FitSample& sample : samples) {
    const BasisVector basis = BuildBasis(sample.inputDelta, order);
    for (Int_t i = 0; i < static_cast<Int_t>(basisTerms); ++i) {
      for (Int_t j = 0; j < static_cast<Int_t>(basisTerms); ++j) {
        normal(i, j) += basis[i] * basis[j];
      }
      for (Int_t row = 0; row < static_cast<Int_t>(kOutputDim); ++row) {
        rhs[row][i] += sample.outputDelta[row] * basis[i];
      }
    }
  }

  TransferMap map;
  map.order = order;
  TDecompSVD svd(normal);
  for (Int_t row = 0; row < static_cast<Int_t>(kPhysicsOutputDim); ++row) {
    TVectorD output(static_cast<Int_t>(basisTerms));
    for (Int_t col = 0; col < static_cast<Int_t>(basisTerms); ++col) {
      output[col] = rhs[row][col];
    }
    Bool_t ok = false;
    const TVectorD coeff = svd.Solve(output, ok);
    if (!ok) {
      throw std::runtime_error("ROOT TDecompSVD failed to solve fit matrix");
    }
    for (Int_t col = 0; col < static_cast<Int_t>(kInputDim); ++col) {
      map.first[row][col] = coeff[col];
    }
    if (order >= 2) {
      for (Int_t col = 0; col < 15; ++col) {
        map.second[row][col] = coeff[static_cast<Int_t>(kInputDim) + col];
      }
    }
  }
  map.first[4] = {0.0, 0.0, 0.0, 0.0, 1.0};
  return map;
}

TransferMap FitTransferMap(const std::vector<FitSample>& samples,
                           unsigned int order,
                           FitStats& stats) {
  if (samples.empty()) {
    throw std::runtime_error("not enough accepted rays for ion-optics fit");
  }
  stats.basisTerms = BuildBasis(samples.front().inputDelta, order).size();
  if (stats.accepted < stats.basisTerms) {
    throw std::runtime_error("not enough accepted rays for ion-optics fit");
  }
  const TransferMap map = SolveMapWithRootSvd(samples, order);
  for (const FitSample& sample : samples) {
    const std::array<double, kOutputDim> predictedDelta =
        EvaluateMap(map, sample.inputDelta);
    for (std::size_t row = 0; row < kPhysicsOutputDim; ++row) {
      const double residual = predictedDelta[row] - sample.outputDelta[row];
      stats.rmsResidual[row] += residual * residual;
      stats.maxAbsResidual[row] =
          std::max(stats.maxAbsResidual[row], std::abs(residual));
    }
  }
  for (double& rms : stats.rmsResidual) {
    rms = std::sqrt(rms / static_cast<double>(stats.accepted));
  }
  return map;
}

std::string DefaultConfigPath() {
  return std::string(MDMTRACE_SOURCE_DIR) + "/config/MDM.json";
}

void WriteAxis(std::ostream& out, const char* label, const Axis& axis) {
  out << "# " << label << ": min=" << axis.min << " max=" << axis.max
      << " step=" << axis.step << " n=" << axis.values.size() << "\n";
}

void WriteReport(std::ostream& out,
                 const std::string& configPath,
                 const Config& cfg,
                 const IonOpticsConfig& optics,
                 const std::array<double, kOutputDim>& referenceOutput,
                 const TransferMap& map,
                 const FitStats& stats) {
  const std::vector<const char*> columns = InputLabels();
  const std::array<const char*, kOutputDim> rows{
      "X1[mm]", "AngX1[mrad]", "Y1[mm]", "AngY1[mrad]", "deltaP/P0[%]"};
  const std::vector<std::string> quadraticLabels = QuadraticLabels();
  const std::array<const char*, 4> residualLabels{"X1[mm]", "AngX1[mrad]",
                                                  "Y1[mm]", "AngY1[mrad]"};

  out << "# GenerateIonOptics fit-based ion-optics map\n";
  out << "# method: " << optics.method << "\n";
  out << "# order: " << map.order << "\n";
  out << "# config: " << configPath << "\n";
  out << "# maps: " << cfg.multipoleMap << ", " << cfg.dipoleEntranceMap
      << ", " << cfg.dipoleSectorMap << ", " << cfg.dipoleExitMap << "\n";
  out << "# vector columns: x[mm] thetaX[mrad] y[mm] thetaY[mrad] "
         "deltaP/P0[%]\n";
  out << "# momentum coordinate: deltaP/P0 is percent; p = p0*(1 + "
         "deltaP/P0[%]/100)\n";
  out << "# reference: x=" << optics.reference.xMm
      << " mm thetaX=" << optics.reference.thetaXMrad
      << " mrad y=" << optics.reference.yMm
      << " mm thetaY=" << optics.reference.thetaYMrad
      << " mrad energy=" << optics.referenceEnergyMeV << " MeV\n";
  out << "# reference output: X1=" << referenceOutput[0]
      << " mm AngX1=" << referenceOutput[1]
      << " mrad Y1=" << referenceOutput[2]
      << " mm AngY1=" << referenceOutput[3] << " mrad\n";
  WriteAxis(out, "fitGrid x[mm]", optics.grid.xMm);
  WriteAxis(out, "fitGrid thetaX[mrad]", optics.grid.thetaXMrad);
  WriteAxis(out, "fitGrid y[mm]", optics.grid.yMm);
  WriteAxis(out, "fitGrid thetaY[mrad]", optics.grid.thetaYMrad);
  WriteAxis(out, "fitGrid deltaP/P0[%]", optics.grid.deltaPercent);
  out << "# solver: ROOT TDecompSVD normal-matrix solve\n";
  out << "# threads: " << stats.threads << "\n";
  out << "# maxRays: " << optics.maxRays << "\n";
  out << "# basisTerms: " << stats.basisTerms << "\n";
  out << "# fit rays: total=" << stats.total << " accepted=" << stats.accepted
      << " stopped=" << stats.stopped << "\n";
  if (!stats.firstSkippedReason.empty()) {
    out << "# first skipped-ray reason: " << stats.firstSkippedReason << "\n";
  }
  out << "# residual convention: fitted prediction - traced field-map output\n";
  for (std::size_t i = 0; i < 4; ++i) {
    out << "# residual " << residualLabels[i] << ": RMS="
        << stats.rmsResidual[i] << " MaxAbs=" << stats.maxAbsResidual[i]
        << "\n";
  }
  out << "# expansion: output-referenceOutput = R*q";
  if (map.order >= 2) {
    out << " + T*(q_j*q_k)";
  }
  out << "\n";
  out << "# q = [x[mm], thetaX[mrad], y[mm], thetaY[mrad], "
         "deltaP/P0[%]]\n";
  out << "# second-order coefficients multiply q_j*q_k directly; no 1/2 "
         "factor is applied\n";
  out << "# first-order units are row unit divided by column unit\n";
  out << std::scientific << std::setprecision(10);
  out << "# First-order R matrix\n";
  out << std::setw(16) << "";
  for (const char* column : columns) {
    out << std::setw(18) << column;
  }
  out << "\n";
  for (std::size_t row = 0; row < 5; ++row) {
    out << std::setw(16) << rows[row];
    for (std::size_t column = 0; column < 5; ++column) {
      out << std::setw(18) << map.first[row][column];
    }
    out << "\n";
  }

  if (map.order >= 2) {
    out << "# Second-order T coefficients\n";
    out << "# second-order units are row unit divided by product of the two "
           "input units\n";
    out << std::setw(16) << "";
    for (const std::string& label : quadraticLabels) {
      out << std::setw(32) << label;
    }
    out << "\n";
    for (std::size_t row = 0; row < 5; ++row) {
      out << std::setw(16) << rows[row];
      for (std::size_t term = 0; term < quadraticLabels.size(); ++term) {
        out << std::setw(32) << map.second[row][term];
      }
      out << "\n";
    }
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  const std::string configPath = argc > 1 ? argv[1] : DefaultConfigPath();
  const Json::Value json = ReadJson(configPath);
  const Config cfg = ParseConfig(json);
  const IonOpticsConfig optics = ParseIonOpticsConfig(json);
  CheckConfig(cfg, optics);

  MDMFieldMapTrace trace;
  ConfigureTrace(trace, cfg);
  const TraceResult referenceResult = TraceRay(trace, cfg, optics, optics.reference);
  if (referenceResult.stopped) {
    throw std::runtime_error("reference ray stopped before the first wire");
  }

  FitStats stats;
  const std::vector<PhaseSpace> inputs = BuildInputGrid(optics.grid);
  std::cerr << "Tracing " << inputs.size() << " ion-optics rays\n";
  const std::vector<FitSample> samples =
      TraceGrid(inputs, cfg, optics, referenceResult.output, stats);
  std::cerr << "Fitting " << stats.accepted << " accepted rays with "
            << stats.threads << " thread(s)\n";
  if (!stats.firstSkippedReason.empty()) {
    std::cerr << "First skipped-ray reason: " << stats.firstSkippedReason
              << "\n";
  }
  const TransferMap map = FitTransferMap(samples, optics.order, stats);
  std::cerr << "Writing " << optics.outputPath << "\n";

  std::ofstream output(optics.outputPath.c_str());
  if (!output) {
    throw std::runtime_error("Cannot write output: " + optics.outputPath);
  }
  WriteReport(output, configPath, cfg, optics, referenceResult.output, map,
              stats);
  WriteReport(std::cout, configPath, cfg, optics, referenceResult.output,
              map, stats);
  return 0;
}
