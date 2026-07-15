#include "MdmConfig.h"
#include "MdmFieldMapTrace.h"
#include "MdmIonConfig.h"
#include "MdmRayScan.h"
#include "MdmTrace.h"

#include "json.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

#ifndef MdmTrace_SOURCE_DIR
#define MdmTrace_SOURCE_DIR "."
#endif

namespace {

struct Config {
  bool usingProbe = false;
  double mdmAngle = 0.0;
  double dipoleField = 0.0;
  double dipoleProbe = 0.0;
  double multipoleProbe = 0.0;
  MdmIon ion;
  std::string multipoleMap = "Multipole.bin";
  std::string dipoleEntranceMap = "DipoleEntrance.bin";
  std::string dipoleSectorMap = "DipoleSector.bin";
  std::string dipoleExitMap = "DipoleExit.bin";
};

struct Result {
  double x = 0.0;
  double y = 0.0;
  double ax = 0.0;
  double ay = 0.0;
  bool stopped = false;
};

struct Row {
  mdm::RayInput ray;
  Result legacy;
  Result fieldMap;
};

struct Quantity {
  const char* key;
  const char* title;
  const char* unit;
  std::function<double(const Result&)> value;
};

struct LinearFit {
  double intercept = 0.0;
  double slope = 0.0;
  double r2 = 0.0;
};

constexpr double kStoppedPositionCm = 1.0e9;

bool IsStopped(const Result& result) {
  return !std::isfinite(result.x) || !std::isfinite(result.y) ||
         !std::isfinite(result.ax) || !std::isfinite(result.ay) ||
         std::abs(result.x) >= kStoppedPositionCm ||
         std::abs(result.y) >= kStoppedPositionCm;
}

double GetDouble(const Json::Value& c, const char* key, double fallback = 0.0) {
  return c.isMember(key) ? c[key].asDouble() : fallback;
}

std::string GetString(const Json::Value& c,
                      const char* key,
                      const char* fallback) {
  return c.isMember(key) ? c[key].asString() : fallback;
}

Config ParseConfig(const Json::Value& c) {
  Config cfg;
  cfg.usingProbe = mdm::RequireBoolean(c, "usingProbe");
  cfg.mdmAngle = GetDouble(c, "mdmAngle");
  cfg.dipoleField = GetDouble(c, "mdmDipoleField");
  cfg.dipoleProbe = GetDouble(c, "mdmDipoleProbe");
  cfg.multipoleProbe = GetDouble(c, "mdmMultipoleProbe");
  cfg.ion = mdm::ParseScatteredIon(c);
  cfg.multipoleMap = GetString(c, "multipoleMapPath", "Multipole.bin");
  cfg.dipoleEntranceMap =
      GetString(c, "dipoleEntranceMapPath", "DipoleEntrance.bin");
  cfg.dipoleSectorMap = GetString(c, "dipoleSectorMapPath", "DipoleSector.bin");
  cfg.dipoleExitMap = GetString(c, "dipoleExitMapPath", "DipoleExit.bin");
  return cfg;
}

void Configure(MdmTrace& trace, const Config& cfg) {
  trace.SetMdmAngle(cfg.mdmAngle);
  cfg.usingProbe ? trace.SetMdmProbe(cfg.dipoleProbe, cfg.multipoleProbe)
                 : trace.SetMdmDipoleField(cfg.dipoleField);
  trace.SetScatteredIon(cfg.ion);
}

void Configure(MdmFieldMapTrace& trace, const Config& cfg) {
  trace.SetMdmAngle(cfg.mdmAngle);
  cfg.usingProbe ? trace.SetMdmProbe(cfg.dipoleProbe, cfg.multipoleProbe)
                 : trace.SetMdmDipoleField(cfg.dipoleField);
  trace.SetScatteredIon(cfg.ion);
  trace.LoadFieldMaps(cfg.multipoleMap, cfg.dipoleEntranceMap,
                      cfg.dipoleSectorMap, cfg.dipoleExitMap);
}

Result Run(MdmTrace& trace, const mdm::RayInput& ray) {
  Result r;
  trace.SetScatteredAngle(ray.xDeg, ray.yDeg);
  trace.SetScatteredEnergy(ray.energyMeV);
  trace.SendRay();
  trace.GetPositionAngleFirstWire(r.x, r.y, r.ax, r.ay);
  r.stopped = IsStopped(r);
  return r;
}

Result Run(MdmFieldMapTrace& trace, const mdm::RayInput& ray) {
  Result r;
  trace.SetScatteredAngle(ray.xDeg, ray.yDeg);
  trace.SetScatteredEnergy(ray.energyMeV);
  trace.SendRay();
  trace.GetPositionAngleFirstWire(r.x, r.y, r.ax, r.ay);
  r.stopped = IsStopped(r);
  return r;
}

unsigned int ResolveProcessCount(const Json::Value& json,
                                 std::size_t totalRays) {
  if (totalRays == 0) {
    return 1;
  }
  unsigned int processes = 0;
  processes =
      mdm::GetNonnegativeInteger(json, "compareProcesses", processes);
  if (processes == 0) {
    processes = std::thread::hardware_concurrency();
    if (processes == 0) {
      processes = 1;
    }
  }
  return std::min<std::size_t>(processes, totalRays);
}

void WriteWorkerRows(const std::string& configPath,
                     std::size_t begin,
                     std::size_t end,
                     const std::string& outputPath) {
  const Json::Value json = mdm::ReadConfig(configPath);
  const Config cfg = ParseConfig(json);
  const std::vector<mdm::RayInput> rays = mdm::ParseRayInputs(json);

  MdmTrace legacy;
  MdmFieldMapTrace fieldMap;
  Configure(legacy, cfg);
  Configure(fieldMap, cfg);

  std::ofstream out(outputPath.c_str());
  out << std::setprecision(17);
  for (std::size_t i = begin; i < end; ++i) {
    const Result legacyResult = Run(legacy, rays[i]);
    const Result fieldMapResult = Run(fieldMap, rays[i]);
    out << i << " " << legacyResult.x << " " << legacyResult.y << " "
        << legacyResult.ax << " " << legacyResult.ay << " "
        << fieldMapResult.x << " " << fieldMapResult.y << " "
        << fieldMapResult.ax << " " << fieldMapResult.ay << "\n";
  }
}

void ReadWorkerRows(const std::string& path, std::vector<Row>& rows) {
  std::ifstream in(path.c_str());
  std::size_t index = 0;
  while (in >> index) {
    in >> rows[index].legacy.x >> rows[index].legacy.y >>
        rows[index].legacy.ax >> rows[index].legacy.ay >>
        rows[index].fieldMap.x >> rows[index].fieldMap.y >>
        rows[index].fieldMap.ax >> rows[index].fieldMap.ay;
    rows[index].legacy.stopped = IsStopped(rows[index].legacy);
    rows[index].fieldMap.stopped = IsStopped(rows[index].fieldMap);
  }
}

std::vector<Row> TraceRowsSerial(const std::vector<mdm::RayInput>& rays,
                                 const Config& cfg) {
  std::vector<Row> rows(rays.size());
  MdmTrace legacy;
  MdmFieldMapTrace fieldMap;
  Configure(legacy, cfg);
  Configure(fieldMap, cfg);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    rows[i].ray = rays[i];
    rows[i].legacy = Run(legacy, rays[i]);
    rows[i].fieldMap = Run(fieldMap, rays[i]);
  }
  return rows;
}

std::vector<Row> TraceRowsParallel(const std::vector<mdm::RayInput>& rays,
                                   const Config& cfg,
                                   const std::string& executable,
                                   const std::string& configPath,
                                   unsigned int processes) {
  if (processes <= 1) {
    return TraceRowsSerial(rays, cfg);
  }

  std::vector<Row> rows(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    rows[i].ray = rays[i];
  }

  std::vector<pid_t> pids;
  std::vector<std::string> files;
  pids.reserve(processes);
  files.reserve(processes);
  for (unsigned int worker = 0; worker < processes; ++worker) {
    const std::size_t begin = rays.size() * worker / processes;
    const std::size_t end = rays.size() * (worker + 1) / processes;
    if (begin == end) {
      continue;
    }

    std::ostringstream fileName;
    fileName << "/private/tmp/Compare_" << getpid() << "_" << worker << ".txt";
    files.push_back(fileName.str());

    const std::string beginText = std::to_string(begin);
    const std::string endText = std::to_string(end);
    const pid_t pid = fork();
    if (pid == 0) {
      execlp(executable.c_str(), executable.c_str(), "--worker",
             configPath.c_str(), beginText.c_str(), endText.c_str(),
             files.back().c_str(), static_cast<char*>(nullptr));
      _exit(127);
    }
    if (pid < 0) {
      throw std::runtime_error("fork failed");
    }
    pids.push_back(pid);
  }

  for (pid_t pid : pids) {
    int status = 0;
    waitpid(pid, &status, 0);
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
      throw std::runtime_error("Compare worker process failed");
    }
  }
  for (const std::string& file : files) {
    ReadWorkerRows(file, rows);
    std::remove(file.c_str());
  }
  return rows;
}

std::pair<double, double> Range(const std::vector<double>& values,
                               bool includeZero = false) {
  double lo = *std::min_element(values.begin(), values.end());
  double hi = *std::max_element(values.begin(), values.end());
  if (includeZero) {
    lo = std::min(lo, 0.0);
    hi = std::max(hi, 0.0);
  }
  const double span = hi - lo;
  const double pad = span > 0.0 ? 0.08 * span : 0.08 * std::max(1.0, std::abs(lo));
  return {lo - pad, hi + pad};
}

LinearFit FitLine(const std::vector<double>& x, const std::vector<double>& y) {
  LinearFit fit;
  const double n = static_cast<double>(x.size());
  double meanX = 0.0;
  double meanY = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    meanX += x[i];
    meanY += y[i];
  }
  meanX /= n;
  meanY /= n;

  double sxx = 0.0;
  double sxy = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double dx = x[i] - meanX;
    sxx += dx * dx;
    sxy += dx * (y[i] - meanY);
  }
  fit.slope = sxx > 0.0 ? sxy / sxx : 0.0;
  fit.intercept = meanY - fit.slope * meanX;

  double ssResidual = 0.0;
  double ssTotal = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double fitted = fit.intercept + fit.slope * x[i];
    ssResidual += (y[i] - fitted) * (y[i] - fitted);
    ssTotal += (y[i] - meanY) * (y[i] - meanY);
  }
  fit.r2 = ssTotal > 0.0 ? 1.0 - ssResidual / ssTotal : 1.0;
  return fit;
}

void MakePlot(const std::vector<Row>& rows, const Quantity& q) {
  std::vector<double> legacy;
  std::vector<double> fieldMap;
  std::vector<double> residual;
  legacy.reserve(rows.size());
  fieldMap.reserve(rows.size());
  residual.reserve(rows.size());
  for (const Row& row : rows) {
    legacy.push_back(q.value(row.legacy));
    fieldMap.push_back(q.value(row.fieldMap));
    residual.push_back(legacy.back() - fieldMap.back());
  }

  std::vector<double> both = legacy;
  both.insert(both.end(), fieldMap.begin(), fieldMap.end());
  const auto xyRange = Range(both);
  const auto residualRange = Range(residual, true);
  const LinearFit fit = FitLine(legacy, fieldMap);

  TCanvas canvas((std::string("c_") + q.key).c_str(),
                 (std::string(q.title) +
                  " legacy vs field-map (both accepted)")
                     .c_str(),
                 900, 900);
  TPad top("top", "", 0.0, 0.32, 1.0, 1.0);
  TPad bottom("bottom", "", 0.0, 0.0, 1.0, 0.32);
  top.SetBottomMargin(0.03);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  bottom.SetTopMargin(0.04);
  bottom.SetBottomMargin(0.30);
  bottom.SetLeftMargin(0.12);
  bottom.SetRightMargin(0.04);
  top.Draw();
  bottom.Draw();

  top.cd();
  TH2D topFrame("topFrame",
                (std::string(q.title) + " (both accepted);Legacy " + q.title +
                 " [" + q.unit + "];FieldMap " + q.title + " [" + q.unit +
                 "]")
                    .c_str(),
                1, xyRange.first, xyRange.second, 1, xyRange.first,
                xyRange.second);
  topFrame.SetStats(false);
  topFrame.GetXaxis()->SetLabelSize(0.0);
  topFrame.GetXaxis()->SetTitleSize(0.0);
  topFrame.Draw();
  TGraph scatter(static_cast<int>(rows.size()), legacy.data(), fieldMap.data());
  scatter.SetMarkerStyle(20);
  scatter.Draw("P SAME");
  TLine fitLine(xyRange.first, fit.intercept + fit.slope * xyRange.first,
                xyRange.second, fit.intercept + fit.slope * xyRange.second);
  fitLine.SetLineColor(kRed);
  fitLine.SetLineWidth(2);
  fitLine.Draw("SAME");
  TPaveText fitText(0.15, 0.74, 0.62, 0.90, "NDC");
  fitText.SetFillColor(0);
  fitText.SetBorderSize(0);
  fitText.SetTextAlign(12);
  fitText.SetTextSize(0.045);
  std::ostringstream equation;
  equation << "Fit: y = " << std::setprecision(5) << fit.intercept << " + "
           << fit.slope << " x";
  std::ostringstream r2Text;
  r2Text << "R^{2} = " << std::setprecision(6) << fit.r2;
  fitText.AddText(equation.str().c_str());
  fitText.AddText(r2Text.str().c_str());
  fitText.Draw("SAME");

  bottom.cd();
  TH2D bottomFrame(
      "bottomFrame",
      (std::string(";Legacy ") + q.title + " [" + q.unit +
       "];Legacy - FieldMap [" + q.unit + "]")
          .c_str(),
      1, xyRange.first, xyRange.second, 1, residualRange.first,
      residualRange.second);
  bottomFrame.SetStats(false);
  bottomFrame.GetXaxis()->SetTitleSize(0.11);
  bottomFrame.GetXaxis()->SetLabelSize(0.10);
  bottomFrame.GetYaxis()->SetTitleSize(0.10);
  bottomFrame.GetYaxis()->SetLabelSize(0.09);
  bottomFrame.GetYaxis()->SetTitleOffset(0.55);
  bottomFrame.Draw();
  TGraph residualGraph(static_cast<int>(rows.size()), legacy.data(),
                       residual.data());
  residualGraph.SetMarkerStyle(20);
  residualGraph.Draw("P SAME");
  TLine zero(xyRange.first, 0.0, xyRange.second, 0.0);
  zero.SetLineColor(kGray + 2);
  zero.SetLineStyle(2);
  zero.Draw("SAME");

  canvas.Write();
}

std::string DefaultConfigPath() {
  return std::string(MdmTrace_SOURCE_DIR) + "/config/MDMScan.json";
}

}  // namespace

int Run(int argc, char* argv[]) {
  if (argc == 6 && std::string(argv[1]) == "--worker") {
    WriteWorkerRows(argv[2],
                    static_cast<std::size_t>(std::stoull(argv[3])),
                    static_cast<std::size_t>(std::stoull(argv[4])),
                    argv[5]);
    return 0;
  }

  const std::string configPath = argc > 1 ? argv[1] : DefaultConfigPath();
  const Json::Value json = mdm::ReadConfig(configPath);
  const Config cfg = ParseConfig(json);
  const std::vector<mdm::RayInput> rays = mdm::ParseRayInputs(json);
  const mdm::RayScanCounts counts = mdm::CountRayScan(json, rays.size());
  const unsigned int processes = ResolveProcessCount(json, rays.size());
  std::cerr << "Energy values: " << counts.energies << "\n"
            << "Horizontal-angle values: " << counts.horizontalAngles << "\n"
            << "Vertical-angle values: " << counts.verticalAngles << "\n"
            << "Total rays: " << counts.rays << "\n"
            << "Tracing " << rays.size() << " rays with " << processes
            << " process(es)\n";
  const std::vector<Row> rows =
      TraceRowsParallel(rays, cfg, argv[0], configPath, processes);

  std::vector<Row> comparableRows;
  comparableRows.reserve(rows.size());
  std::size_t bothStopped = 0;
  std::size_t legacyOnlyStopped = 0;
  std::size_t fieldMapOnlyStopped = 0;
  const Row* firstLegacyOnlyStopped = nullptr;
  const Row* firstFieldMapOnlyStopped = nullptr;
  for (const Row& row : rows) {
    if (!row.legacy.stopped && !row.fieldMap.stopped) {
      comparableRows.push_back(row);
    } else if (row.legacy.stopped && row.fieldMap.stopped) {
      ++bothStopped;
    } else if (row.legacy.stopped) {
      ++legacyOnlyStopped;
      if (firstLegacyOnlyStopped == nullptr) {
        firstLegacyOnlyStopped = &row;
      }
    } else {
      ++fieldMapOnlyStopped;
      if (firstFieldMapOnlyStopped == nullptr) {
        firstFieldMapOnlyStopped = &row;
      }
    }
  }
  if (comparableRows.empty()) {
    throw std::runtime_error("No rays reached the first wire in both transports");
  }

  std::cout << "Stop classification:\n"
            << "  Both accepted: " << comparableRows.size() << "\n"
            << "  Both stopped: " << bothStopped << "\n"
            << "  Legacy-only stopped: " << legacyOnlyStopped << "\n"
            << "  FieldMap-only stopped: " << fieldMapOnlyStopped << "\n";
  const auto printFirstMismatch = [](const char* label, const Row* row) {
    if (row != nullptr) {
      std::cout << label << ": AngleX=" << row->ray.xDeg
                << " deg, AngleY=" << row->ray.yDeg
                << " deg, Energy=" << row->ray.energyMeV << " MeV\n";
    }
  };
  printFirstMismatch("First legacy-only stop", firstLegacyOnlyStopped);
  printFirstMismatch("First FieldMap-only stop", firstFieldMapOnlyStopped);

  const std::vector<Quantity> quantities{
      {"X1", "X", "cm", [](const Result& r) { return r.x; }},
      {"Y1", "Y", "cm", [](const Result& r) { return r.y; }},
      {"AngX1", "AngX", "deg", [](const Result& r) { return r.ax; }},
      {"AngY1", "AngY", "deg", [](const Result& r) { return r.ay; }},
  };

  TFile output("Compare.root", "RECREATE");
  gStyle->SetOptStat(0);
  for (const Quantity& q : quantities) {
    MakePlot(comparableRows, q);
  }
  output.Close();

  std::cout << "Wrote Compare.root with 4 canvases and "
            << comparableRows.size()
            << " rays accepted by both transports using " << processes
            << " process(es)\n";
  std::cout << "Residual convention: Legacy - FieldMap\n"
            << "Residual statistics include only rays accepted by both "
               "transports.\n";
  std::cout << std::fixed << std::setprecision(6);
  for (const Quantity& q : quantities) {
    double sumSq = 0.0;
    double maxAbs = 0.0;
    for (const Row& row : comparableRows) {
      const double residual = q.value(row.legacy) - q.value(row.fieldMap);
      sumSq += residual * residual;
      maxAbs = std::max(maxAbs, std::abs(residual));
    }
    std::cout << std::setw(8) << q.key << " RMS="
              << std::sqrt(sumSq /
                           static_cast<double>(comparableRows.size()))
              << " MaxAbs=" << maxAbs << " " << q.unit << "\n";
  }
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
