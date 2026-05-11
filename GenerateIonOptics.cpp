#include "MdmFieldMapTrace.h"
#include "MdmIonConfig.h"
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

#ifndef MdmTrace_SOURCE_DIR
#define MdmTrace_SOURCE_DIR "."
#endif

namespace {

constexpr double kMradPerDegree = 1000.0 * 3.14159265358979323846 / 180.0;
constexpr double kTraceMradPerDegree = 17.453;
constexpr double kSpeedOfLightCmPerSecond = 3.0e10;
constexpr double kStoppedCm = 1.0e9;
constexpr double kPercentToFraction = 0.01;
constexpr double kFirstOrderPrintZero = 1.0e-10;
constexpr std::size_t kInputDim = 6;
constexpr std::size_t kOutputDim = 6;
constexpr std::size_t kPhysicsOutputDim = 5;
constexpr std::size_t kLIndex = 4;
constexpr std::size_t kDeltaIndex = 5;

struct Config {
  bool usingProbe = false;
  double mdmAngleDeg = 0.0;
  double dipoleField = 0.0;
  double dipoleProbe = 0.0;
  double multipoleProbe = 0.0;
  MdmIon ion;
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
  double lMm = 0.0;
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
  double tofSeconds = 0.0;
  std::array<double, kOutputDim> output{};
};

struct FitStats {
  std::size_t total = 0;
  std::size_t accepted = 0;
  std::size_t stopped = 0;
  std::size_t basisTerms = 0;
  std::size_t activeFitTerms = 0;
  unsigned int threads = 1;
  std::string firstSkippedReason;
  std::array<double, kPhysicsOutputDim> rmsResidual{};
  std::array<double, kPhysicsOutputDim> maxAbsResidual{};
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

using BasisVector = std::vector<double>;

struct Monomial {
  std::array<unsigned int, kInputDim> powers{};
};

struct TransferMap {
  unsigned int order = 1;
  std::vector<Monomial> terms;
  std::array<std::vector<double>, kOutputDim> coeff;
};

struct LongitudinalReference {
  double tofSeconds = 0.0;
  double speedCmPerSecond = 0.0;
  double gamma = 1.0;
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
  cfg.ion = mdm::ParseScatteredIon(json);
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
    cfg.reference.lMm = GetDouble(reference, "lMm");
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
  if (optics.order < 1 || optics.order > 5) {
    throw std::runtime_error("ionOptics.order must be 1, 2, 3, 4, or 5");
  }
  if (cfg.ion.ionMassMeV <= 0.0 || cfg.ion.chargeState <= 0) {
    throw std::runtime_error("scatteredIon must be physical");
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

void ConfigureTrace(MdmFieldMapTrace& trace, const Config& cfg) {
  trace.SetMdmAngle(cfg.mdmAngleDeg);
  cfg.usingProbe ? trace.SetMdmProbe(cfg.dipoleProbe, cfg.multipoleProbe)
                 : trace.SetMdmDipoleField(cfg.dipoleField);
  trace.SetScatteredIon(cfg.ion);
  trace.LoadFieldMaps(cfg.multipoleMap, cfg.dipoleEntranceMap,
                      cfg.dipoleSectorMap, cfg.dipoleExitMap);
}

double EnergyFromDeltaPercent(double referenceEnergyMeV,
                              double massMeV,
                              double deltaPercent) {
  const double p0 =
      std::sqrt((2.0 * massMeV + referenceEnergyMeV) * referenceEnergyMeV);
  const double p = p0 * (1.0 + deltaPercent * kPercentToFraction);
  return std::sqrt(p * p + massMeV * massMeV) - massMeV;
}

double GammaFromEnergy(double energyMeV, double massMeV) {
  return (massMeV + energyMeV) / massMeV;
}

double SpeedFromEnergy(double energyMeV, double massMeV) {
  const double p = std::sqrt((2.0 * massMeV + energyMeV) * energyMeV);
  return p / (massMeV + energyMeV) * kSpeedOfLightCmPerSecond;
}

TraceResult TraceRay(MdmFieldMapTrace& trace,
                     const Config& cfg,
                     const IonOpticsConfig& optics,
                     const PhaseSpace& input,
                     const LongitudinalReference* longitudinal = nullptr) {
  trace.SetInitialPosition(input.xMm / 10.0, input.yMm / 10.0);
  trace.SetScatteredAngle(cfg.mdmAngleDeg + input.thetaXMrad /
                                             kTraceMradPerDegree,
                          input.thetaYMrad / kTraceMradPerDegree);
  trace.SetScatteredEnergy(
      EnergyFromDeltaPercent(optics.referenceEnergyMeV, cfg.ion.ionMassMeV,
                             input.deltaPercent));
  trace.SendRay();

  double xCm = 0.0;
  double yCm = 0.0;
  double angXDeg = 0.0;
  double angYDeg = 0.0;
  trace.GetPositionAngleFirstWire(xCm, yCm, angXDeg, angYDeg);
  const double tofSeconds = trace.GetTimeOfFlightSeconds();

  TraceResult result;
  result.tofSeconds = tofSeconds;
  result.stopped = !std::isfinite(xCm) || !std::isfinite(yCm) ||
                   !std::isfinite(angXDeg) || !std::isfinite(angYDeg) ||
                   !std::isfinite(tofSeconds) ||
                   std::abs(xCm) > kStoppedCm || std::abs(yCm) > kStoppedCm;
  double lMm = input.lMm;
  if (longitudinal) {
    lMm += -longitudinal->speedCmPerSecond *
           (tofSeconds - longitudinal->tofSeconds) /
           (longitudinal->gamma * longitudinal->gamma) * 10.0;
  }
  result.output = {10.0 * xCm, angXDeg * kMradPerDegree, 10.0 * yCm,
                   angYDeg * kMradPerDegree, lMm, input.deltaPercent};
  return result;
}

std::array<double, kInputDim> InputDelta(const PhaseSpace& input,
                                         const PhaseSpace& reference) {
  return {input.xMm - reference.xMm,
          input.thetaXMrad - reference.thetaXMrad,
          input.yMm - reference.yMm,
          input.thetaYMrad - reference.thetaYMrad,
          input.lMm - reference.lMm,
          input.deltaPercent - reference.deltaPercent};
}

unsigned int Degree(const Monomial& term) {
  unsigned int degree = 0;
  for (unsigned int power : term.powers) {
    degree += power;
  }
  return degree;
}

void AddMonomialsOfDegree(unsigned int remaining,
                          std::size_t firstVariable,
                          Monomial& term,
                          std::vector<Monomial>& terms) {
  if (remaining == 0) {
    terms.push_back(term);
    return;
  }
  for (std::size_t variable = firstVariable; variable < kInputDim; ++variable) {
    ++term.powers[variable];
    AddMonomialsOfDegree(remaining - 1, variable, term, terms);
    --term.powers[variable];
  }
}

std::vector<Monomial> BuildMonomials(unsigned int order) {
  std::vector<Monomial> terms;
  for (unsigned int degree = 1; degree <= order; ++degree) {
    Monomial term;
    AddMonomialsOfDegree(degree, 0, term, terms);
  }
  return terms;
}

std::vector<std::size_t> ActiveFitTermIndices(
    const std::vector<Monomial>& terms) {
  std::vector<std::size_t> indices;
  for (std::size_t i = 0; i < terms.size(); ++i) {
    if (terms[i].powers[kLIndex] == 0) {
      indices.push_back(i);
    }
  }
  return indices;
}

double MonomialValue(const std::array<double, kInputDim>& inputDelta,
                     const Monomial& term) {
  double value = 1.0;
  for (std::size_t variable = 0; variable < kInputDim; ++variable) {
    for (unsigned int power = 0; power < term.powers[variable]; ++power) {
      value *= inputDelta[variable];
    }
  }
  return value;
}

std::size_t LinearTermIndex(const std::vector<Monomial>& terms,
                            std::size_t variable) {
  for (std::size_t index = 0; index < terms.size(); ++index) {
    if (Degree(terms[index]) == 1 && terms[index].powers[variable] == 1) {
      return index;
    }
  }
  return terms.size();
}

std::string ExponentCode(const Monomial& term) {
  std::string code;
  code.reserve(kInputDim);
  for (unsigned int power : term.powers) {
    code += static_cast<char>('0' + power);
  }
  return code;
}

BasisVector BuildBasis(const std::array<double, kInputDim>& inputDelta,
                       const std::vector<Monomial>& terms) {
  BasisVector basis;
  basis.reserve(terms.size());
  for (const Monomial& term : terms) {
    basis.push_back(MonomialValue(inputDelta, term));
  }
  return basis;
}

BasisVector BuildBasis(const std::array<double, kInputDim>& inputDelta,
                       const std::vector<Monomial>& terms,
                       const std::vector<std::size_t>& termIndices) {
  BasisVector basis;
  basis.reserve(termIndices.size());
  for (std::size_t term : termIndices) {
    basis.push_back(MonomialValue(inputDelta, terms[term]));
  }
  return basis;
}

std::array<double, kOutputDim> EvaluateMap(
    const TransferMap& map,
    const std::array<double, kInputDim>& inputDelta) {
  std::array<double, kOutputDim> result{};
  const BasisVector basis = BuildBasis(inputDelta, map.terms);
  for (std::size_t row = 0; row < kOutputDim; ++row) {
    for (std::size_t term = 0; term < map.terms.size(); ++term) {
      result[row] += map.coeff[row][term] * basis[term];
    }
  }
  return result;
}

std::size_t GridRayCount(const FitGrid& grid) {
  return grid.xMm.values.size() * grid.thetaXMrad.values.size() *
         grid.yMm.values.size() * grid.thetaYMrad.values.size() *
         grid.deltaPercent.values.size();
}

std::vector<PhaseSpace> BuildInputGrid(const FitGrid& grid,
                                       const PhaseSpace& reference) {
  std::vector<PhaseSpace> inputs;
  inputs.reserve(GridRayCount(grid));
  for (double xMm : grid.xMm.values) {
    for (double thetaXMrad : grid.thetaXMrad.values) {
      for (double yMm : grid.yMm.values) {
        for (double thetaYMrad : grid.thetaYMrad.values) {
          for (double deltaPercent : grid.deltaPercent.values) {
            inputs.push_back({xMm, thetaXMrad, yMm, thetaYMrad,
                              reference.lMm, deltaPercent});
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
                                 const LongitudinalReference& longitudinal,
                                 FitStats& stats) {
  stats.total = inputs.size();
  stats.threads = ResolveThreadCount(optics.requestedThreads, inputs.size());

  std::vector<std::unique_ptr<MdmFieldMapTrace>> traces;
  traces.reserve(stats.threads);
  for (unsigned int worker = 0; worker < stats.threads; ++worker) {
    auto trace = std::make_unique<MdmFieldMapTrace>();
    ConfigureTrace(*trace, cfg);
    const TraceResult warmup =
        TraceRay(*trace, cfg, optics, optics.reference, &longitudinal);
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
  const auto runRange = [&](MdmFieldMapTrace& trace, std::size_t begin,
                            std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      TraceResult traceResult;
      try {
        traceResult = TraceRay(trace, cfg, optics, inputs[i], &longitudinal);
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

std::vector<double> BuildColumnScales(const std::vector<FitSample>& samples,
                                      const std::vector<Monomial>& terms,
                                      const std::vector<std::size_t>& activeTerms,
                                      unsigned int threads) {
  std::vector<std::vector<double>> localMax(
      threads, std::vector<double>(activeTerms.size(), 0.0));

  const auto runRange = [&](unsigned int worker, std::size_t begin,
                            std::size_t end) {
    for (std::size_t sampleIndex = begin; sampleIndex < end; ++sampleIndex) {
      const BasisVector basis =
          BuildBasis(samples[sampleIndex].inputDelta, terms, activeTerms);
      for (std::size_t term = 0; term < activeTerms.size(); ++term) {
        localMax[worker][term] =
            std::max(localMax[worker][term], std::abs(basis[term]));
      }
    }
  };

  if (threads == 1) {
    runRange(0, 0, samples.size());
  } else {
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned int worker = 0; worker < threads; ++worker) {
      const std::size_t begin = samples.size() * worker / threads;
      const std::size_t end = samples.size() * (worker + 1) / threads;
      workers.emplace_back(runRange, worker, begin, end);
    }
    for (std::thread& worker : workers) {
      worker.join();
    }
  }

  std::vector<double> scales(activeTerms.size(), 1.0);
  for (std::size_t term = 0; term < activeTerms.size(); ++term) {
    double maximum = 0.0;
    for (unsigned int worker = 0; worker < threads; ++worker) {
      maximum = std::max(maximum, localMax[worker][term]);
    }
    if (maximum > 0.0) {
      scales[term] = 1.0 / maximum;
    }
  }
  return scales;
}

TransferMap SolveMapWithRootSvd(const std::vector<FitSample>& samples,
                                unsigned int order,
                                unsigned int threads) {
  const std::vector<Monomial> terms = BuildMonomials(order);
  const std::vector<std::size_t> activeTerms = ActiveFitTermIndices(terms);
  const std::size_t basisTerms = activeTerms.size();
  const std::vector<double> columnScale =
      BuildColumnScales(samples, terms, activeTerms, threads);

  std::vector<std::vector<double>> localNormal(
      threads, std::vector<double>(basisTerms * basisTerms, 0.0));
  std::vector<std::array<std::vector<double>, kPhysicsOutputDim>> localRhs(
      threads);
  for (auto& rhs : localRhs) {
    for (std::vector<double>& row : rhs) {
      row.assign(basisTerms, 0.0);
    }
  }

  const auto runRange = [&](unsigned int worker, std::size_t begin,
                            std::size_t end) {
    std::vector<double>& normal = localNormal[worker];
    auto& rhs = localRhs[worker];
    for (std::size_t sampleIndex = begin; sampleIndex < end; ++sampleIndex) {
      BasisVector basis =
          BuildBasis(samples[sampleIndex].inputDelta, terms, activeTerms);
      for (std::size_t term = 0; term < basisTerms; ++term) {
        basis[term] *= columnScale[term];
      }
      for (std::size_t i = 0; i < basisTerms; ++i) {
        const double bi = basis[i];
        for (std::size_t j = i; j < basisTerms; ++j) {
          normal[i * basisTerms + j] += bi * basis[j];
        }
        for (std::size_t row = 0; row < kPhysicsOutputDim; ++row) {
          rhs[row][i] += samples[sampleIndex].outputDelta[row] * bi;
        }
      }
    }
  };

  if (threads == 1) {
    runRange(0, 0, samples.size());
  } else {
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned int worker = 0; worker < threads; ++worker) {
      const std::size_t begin = samples.size() * worker / threads;
      const std::size_t end = samples.size() * (worker + 1) / threads;
      workers.emplace_back(runRange, worker, begin, end);
    }
    for (std::thread& worker : workers) {
      worker.join();
    }
  }

  TMatrixD normal(static_cast<Int_t>(basisTerms),
                  static_cast<Int_t>(basisTerms));
  std::array<std::vector<double>, kPhysicsOutputDim> rhs;
  for (std::vector<double>& row : rhs) {
    row.assign(basisTerms, 0.0);
  }
  for (std::size_t i = 0; i < basisTerms; ++i) {
    for (std::size_t j = i; j < basisTerms; ++j) {
      double value = 0.0;
      for (unsigned int worker = 0; worker < threads; ++worker) {
        value += localNormal[worker][i * basisTerms + j];
      }
      normal(static_cast<Int_t>(i), static_cast<Int_t>(j)) = value;
      normal(static_cast<Int_t>(j), static_cast<Int_t>(i)) = value;
    }
    for (std::size_t row = 0; row < kPhysicsOutputDim; ++row) {
      for (unsigned int worker = 0; worker < threads; ++worker) {
        rhs[row][i] += localRhs[worker][row][i];
      }
    }
  }

  TransferMap map;
  map.order = order;
  map.terms = terms;
  for (std::vector<double>& row : map.coeff) {
    row.assign(terms.size(), 0.0);
  }

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
    for (Int_t col = 0; col < static_cast<Int_t>(basisTerms); ++col) {
      map.coeff[row][activeTerms[col]] = coeff[col] * columnScale[col];
    }
  }
  const std::size_t lTerm = LinearTermIndex(map.terms, kLIndex);
  if (lTerm < map.terms.size()) {
    map.coeff[kLIndex][lTerm] = 1.0;
  }
  const std::size_t deltaTerm = LinearTermIndex(map.terms, kDeltaIndex);
  if (deltaTerm < map.terms.size()) {
    map.coeff[kDeltaIndex][deltaTerm] = 1.0;
  }
  return map;
}

TransferMap FitTransferMap(const std::vector<FitSample>& samples,
                           unsigned int order,
                           FitStats& stats) {
  if (samples.empty()) {
    throw std::runtime_error("not enough accepted rays for ion-optics fit");
  }
  const std::vector<Monomial> terms = BuildMonomials(order);
  stats.basisTerms = terms.size();
  stats.activeFitTerms = ActiveFitTermIndices(terms).size();
  if (stats.accepted < stats.activeFitTerms) {
    throw std::runtime_error("not enough accepted rays for ion-optics fit");
  }
  const TransferMap map = SolveMapWithRootSvd(samples, order, stats.threads);
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
  return std::string(MdmTrace_SOURCE_DIR) + "/config/MDM.json";
}

void WriteAxis(std::ostream& out, const char* label, const Axis& axis) {
  out << "# " << label << ": min=" << axis.min << " max=" << axis.max
      << " step=" << axis.step << " n=" << axis.values.size() << "\n";
}

double FirstOrderPrintValue(double value) {
  return std::abs(value) < kFirstOrderPrintZero ? 0.0 : value;
}

double FirstOrderElement(const TransferMap& map,
                         std::size_t row,
                         std::size_t column) {
  const std::size_t term = LinearTermIndex(map.terms, column);
  return term < map.terms.size() ? map.coeff[row][term] : 0.0;
}

double FirstOrderDeterminant(const TransferMap& map) {
  std::array<std::array<double, kInputDim>, kInputDim> matrix{};
  for (std::size_t row = 0; row < kInputDim; ++row) {
    for (std::size_t column = 0; column < kInputDim; ++column) {
      matrix[row][column] = FirstOrderElement(map, row, column);
    }
  }

  double determinant = 1.0;
  for (std::size_t pivot = 0; pivot < kInputDim; ++pivot) {
    std::size_t best = pivot;
    for (std::size_t row = pivot + 1; row < kInputDim; ++row) {
      if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot])) {
        best = row;
      }
    }
    if (std::abs(matrix[best][pivot]) < 1.0e-30) {
      return 0.0;
    }
    if (best != pivot) {
      std::swap(matrix[best], matrix[pivot]);
      determinant = -determinant;
    }
    const double pivotValue = matrix[pivot][pivot];
    determinant *= pivotValue;
    for (std::size_t row = pivot + 1; row < kInputDim; ++row) {
      const double factor = matrix[row][pivot] / pivotValue;
      for (std::size_t column = pivot + 1; column < kInputDim; ++column) {
        matrix[row][column] -= factor * matrix[pivot][column];
      }
    }
  }
  return determinant;
}

void WriteReport(std::ostream& out,
                 const std::string& configPath,
                 const Config& cfg,
                 const IonOpticsConfig& optics,
                 const std::array<double, kOutputDim>& referenceOutput,
                 const TransferMap& map,
                 const FitStats& stats) {
  const std::array<const char*, kOutputDim> columns{
      "X1[mm]", "AngX1[mrad]", "Y1[mm]", "AngY1[mrad]",
      "L[mm]", "deltaP/P0[%]"};
  const std::array<const char*, kPhysicsOutputDim> residualLabels{
      "X1[mm]", "AngX1[mrad]", "Y1[mm]", "AngY1[mrad]", "L[mm]"};

  out << "# GenerateIonOptics fit-based ion-optics map\n";
  out << "# method: " << optics.method << "\n";
  out << "# order: " << map.order << "\n";
  out << "# config: " << configPath << "\n";
  out << "# maps: " << cfg.multipoleMap << ", " << cfg.dipoleEntranceMap
      << ", " << cfg.dipoleSectorMap << ", " << cfg.dipoleExitMap << "\n";
  out << "# ion: A=" << cfg.ion.massNumber << " Z=" << cfg.ion.atomicNumber
      << " q=" << cfg.ion.chargeState << " neutralMassU="
      << cfg.ion.neutralMassU << " ionMassMeV=" << cfg.ion.ionMassMeV
      << "\n";
  out << "# vector columns: x[mm] thetaX[mrad] y[mm] thetaY[mrad] "
         "L[mm] deltaP/P0[%]\n";
  out << "# momentum coordinate: deltaP/P0 is percent; p = p0*(1 + "
         "deltaP/P0[%]/100)\n";
  out << "# L convention: L = -v0_ref*(tof-tof_ref)/gamma0_ref^2; unit mm\n";
  out << "# input L is a formal identity coordinate; no L scan grid is traced\n";
  out << "# reference: x=" << optics.reference.xMm
      << " mm thetaX=" << optics.reference.thetaXMrad
      << " mrad y=" << optics.reference.yMm
      << " mm thetaY=" << optics.reference.thetaYMrad
      << " mrad L=" << optics.reference.lMm
      << " mm energy=" << optics.referenceEnergyMeV << " MeV\n";
  out << "# reference output: X1=" << referenceOutput[0]
      << " mm AngX1=" << referenceOutput[1]
      << " mrad Y1=" << referenceOutput[2]
      << " mm AngY1=" << referenceOutput[3]
      << " mrad L=" << referenceOutput[4]
      << " mm deltaP/P0=" << referenceOutput[5] << " %\n";
  WriteAxis(out, "fitGrid x[mm]", optics.grid.xMm);
  WriteAxis(out, "fitGrid thetaX[mrad]", optics.grid.thetaXMrad);
  WriteAxis(out, "fitGrid y[mm]", optics.grid.yMm);
  WriteAxis(out, "fitGrid thetaY[mrad]", optics.grid.thetaYMrad);
  WriteAxis(out, "fitGrid deltaP/P0[%]", optics.grid.deltaPercent);
  out << "# solver: ROOT TDecompSVD normal-matrix solve\n";
  out << "# threads: " << stats.threads << "\n";
  out << "# maxRays: " << optics.maxRays << "\n";
  out << "# basisTerms: " << stats.basisTerms << "\n";
  out << "# activeFitTerms: " << stats.activeFitTerms << "\n";
  out << "# fit rays: total=" << stats.total << " accepted=" << stats.accepted
      << " stopped=" << stats.stopped << "\n";
  if (!stats.firstSkippedReason.empty()) {
    out << "# first skipped-ray reason: " << stats.firstSkippedReason << "\n";
  }
  out << "# residual convention: fitted prediction - traced field-map output\n";
  for (std::size_t i = 0; i < kPhysicsOutputDim; ++i) {
    out << "# residual " << residualLabels[i] << ": RMS="
        << stats.rmsResidual[i] << " MaxAbs=" << stats.maxAbsResidual[i]
        << "\n";
  }
  out << "# expansion: output-referenceOutput = sum C_m * monomial_m(q)\n";
  out << "# q = [x[mm], thetaX[mrad], y[mm], thetaY[mrad], "
         "L[mm], deltaP/P0[%]]\n";
  out << "# fitted terms have L exponent zero; all L-input terms are zero "
         "except L->L=1\n";
  out << "# coefficients multiply monomials directly; no factorial factors "
         "are applied\n";
  out << "# human-readable first-order coefficients with abs(value)<1e-10 "
         "are printed as zero\n";
  out << std::uppercase << std::scientific << std::setprecision(10);
  out << "# Human-readable first-order R matrix\n";
  out << std::setw(16) << "";
  for (const char* column : columns) {
    out << std::setw(18) << column;
  }
  out << "\n";
  for (std::size_t row = 0; row < kOutputDim; ++row) {
    out << std::setw(16) << columns[row];
    for (std::size_t column = 0; column < kInputDim; ++column) {
      const double value = FirstOrderElement(map, row, column);
      out << std::setw(18) << FirstOrderPrintValue(value);
    }
    out << "\n";
  }
  out << "# det(first-order R) = " << FirstOrderDeterminant(map) << "\n";

  out << "# COSY-style map columns: X1[mm] AngX1[mrad] Y1[mm] "
         "AngY1[mrad] L[mm] deltaP/P0[%] exponents\n";
  out << "# exponent order: x thetaX y thetaY L deltaP/P0[%]\n";
  out << "#";
  for (const char* column : columns) {
    out << std::setw(18) << column;
  }
  out << std::setw(18) << "exponents" << "\n";
  for (std::size_t term = 0; term < map.terms.size(); ++term) {
    for (std::size_t row = 0; row < kOutputDim; ++row) {
      out << std::setw(18) << map.coeff[row][term];
    }
    out << "  " << ExponentCode(map.terms[term]);
    out << "\n";
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  const std::string configPath = argc > 1 ? argv[1] : DefaultConfigPath();
  const Json::Value json = ReadJson(configPath);
  const Config cfg = ParseConfig(json);
  const IonOpticsConfig optics = ParseIonOpticsConfig(json);
  CheckConfig(cfg, optics);

  MdmFieldMapTrace trace;
  ConfigureTrace(trace, cfg);
  const TraceResult referenceResult = TraceRay(trace, cfg, optics, optics.reference);
  if (referenceResult.stopped) {
    throw std::runtime_error("reference ray stopped before the first wire");
  }
  const LongitudinalReference longitudinal{
      referenceResult.tofSeconds,
      SpeedFromEnergy(optics.referenceEnergyMeV, cfg.ion.ionMassMeV),
      GammaFromEnergy(optics.referenceEnergyMeV, cfg.ion.ionMassMeV)};

  FitStats stats;
  const std::vector<PhaseSpace> inputs =
      BuildInputGrid(optics.grid, optics.reference);
  std::cerr << "Tracing " << inputs.size() << " ion-optics rays\n";
  const std::vector<FitSample> samples =
      TraceGrid(inputs, cfg, optics, referenceResult.output, longitudinal,
                stats);
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
