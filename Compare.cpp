#include "MDMFieldMapTrace.h"
#include "MDMRayScan.h"
#include "MDMTrace.h"

#include "json.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef MDMTRACE_SOURCE_DIR
#define MDMTRACE_SOURCE_DIR "."
#endif

namespace {

struct Config {
  bool usingProbe = false;
  double mdmAngle = 0.0;
  double dipoleField = 0.0;
  double dipoleProbe = 0.0;
  double multipoleProbe = 0.0;
  double mass = 0.0;
  double charge = 0.0;
  double energy = 0.0;
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
};

struct Row {
  mdm::RayAngle ray;
  Result legacy;
  Result fieldMap;
};

struct Quantity {
  const char* key;
  const char* title;
  const char* unit;
  std::function<double(const Result&)> value;
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
  cfg.usingProbe = c.isMember("usingProbe") && c["usingProbe"].asBool();
  cfg.mdmAngle = GetDouble(c, "mdmAngle");
  cfg.dipoleField = GetDouble(c, "mdmDipoleField");
  cfg.dipoleProbe = GetDouble(c, "mdmDipoleProbe");
  cfg.multipoleProbe = GetDouble(c, "mdmMultipoleProbe");
  cfg.mass = GetDouble(c, "scatteredMass");
  cfg.charge = GetDouble(c, "scatteredCharge");
  cfg.energy = GetDouble(c, "scatteredEnergy");
  cfg.multipoleMap = GetString(c, "multipoleMapPath", "Multipole.bin");
  cfg.dipoleEntranceMap =
      GetString(c, "dipoleEntranceMapPath", "DipoleEntrance.bin");
  cfg.dipoleSectorMap = GetString(c, "dipoleSectorMapPath", "DipoleSector.bin");
  cfg.dipoleExitMap = GetString(c, "dipoleExitMapPath", "DipoleExit.bin");
  return cfg;
}

void Configure(MDMTrace& trace, const Config& cfg) {
  trace.SetMDMAngle(cfg.mdmAngle);
  cfg.usingProbe ? trace.SetMDMProbe(cfg.dipoleProbe, cfg.multipoleProbe)
                 : trace.SetMDMDipoleField(cfg.dipoleField);
  trace.SetScatteredMass(cfg.mass);
  trace.SetScatteredCharge(cfg.charge);
  trace.SetScatteredEnergy(cfg.energy);
}

void Configure(MDMFieldMapTrace& trace, const Config& cfg) {
  trace.SetMDMAngle(cfg.mdmAngle);
  cfg.usingProbe ? trace.SetMDMProbe(cfg.dipoleProbe, cfg.multipoleProbe)
                 : trace.SetMDMDipoleField(cfg.dipoleField);
  trace.SetScatteredMass(cfg.mass);
  trace.SetScatteredCharge(cfg.charge);
  trace.SetScatteredEnergy(cfg.energy);
  trace.LoadFieldMaps(cfg.multipoleMap, cfg.dipoleEntranceMap,
                      cfg.dipoleSectorMap, cfg.dipoleExitMap);
}

Result Run(MDMTrace& trace, const mdm::RayAngle& ray) {
  Result r;
  trace.SetScatteredAngle(ray.xDeg, ray.yDeg);
  trace.SendRay();
  trace.GetPositionAngleFirstWire(r.x, r.y, r.ax, r.ay);
  return r;
}

Result Run(MDMFieldMapTrace& trace, const mdm::RayAngle& ray) {
  Result r;
  trace.SetScatteredAngle(ray.xDeg, ray.yDeg);
  trace.SendRay();
  trace.GetPositionAngleFirstWire(r.x, r.y, r.ax, r.ay);
  return r;
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

  TCanvas canvas((std::string("c_") + q.key).c_str(),
                 (std::string(q.title) + " legacy vs field-map").c_str(), 900,
                 900);
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
                (std::string(q.title) + ";Legacy " + q.title + " [" + q.unit +
                 "];FieldMap " + q.title + " [" + q.unit + "]")
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
  TLine diagonal(xyRange.first, xyRange.first, xyRange.second, xyRange.second);
  diagonal.SetLineColor(kRed);
  diagonal.SetLineWidth(2);
  diagonal.Draw("SAME");

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
  return std::string(MDMTRACE_SOURCE_DIR) +
         "/config/config-MDMTraceAngleGrid.json";
}

}  // namespace

int main(int argc, char* argv[]) {
  const Json::Value json = ReadJson(argc > 1 ? argv[1] : DefaultConfigPath());
  const Config cfg = ParseConfig(json);
  const std::vector<mdm::RayAngle> rays = mdm::ParseRayAngles(json);
  std::vector<Row> rows;
  rows.reserve(rays.size());

  {
    MDMTrace trace;
    Configure(trace, cfg);
    for (const mdm::RayAngle& ray : rays) {
      rows.push_back({ray, Run(trace, ray), {}});
    }
  }
  {
    MDMFieldMapTrace trace;
    Configure(trace, cfg);
    for (Row& row : rows) {
      row.fieldMap = Run(trace, row.ray);
    }
  }

  const std::vector<Quantity> quantities{
      {"X1", "X", "cm", [](const Result& r) { return r.x; }},
      {"Y1", "Y", "cm", [](const Result& r) { return r.y; }},
      {"AngX1", "AngX", "deg", [](const Result& r) { return r.ax; }},
      {"AngY1", "AngY", "deg", [](const Result& r) { return r.ay; }},
  };

  TFile output("Compare.root", "RECREATE");
  gStyle->SetOptStat(0);
  for (const Quantity& q : quantities) {
    MakePlot(rows, q);
  }
  output.Close();

  std::cout << "Wrote Compare.root with 4 canvases and " << rows.size()
            << " rays\n";
  std::cout << "Residual convention: Legacy - FieldMap\n";
  std::cout << std::fixed << std::setprecision(6);
  for (const Quantity& q : quantities) {
    double sumSq = 0.0;
    double maxAbs = 0.0;
    for (const Row& row : rows) {
      const double residual = q.value(row.legacy) - q.value(row.fieldMap);
      sumSq += residual * residual;
      maxAbs = std::max(maxAbs, std::abs(residual));
    }
    std::cout << std::setw(8) << q.key << " RMS="
              << std::sqrt(sumSq / static_cast<double>(rows.size()))
              << " MaxAbs=" << maxAbs << " " << q.unit << "\n";
  }
  return 0;
}
