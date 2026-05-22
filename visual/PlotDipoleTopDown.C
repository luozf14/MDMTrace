#include "FieldMapCommon.C"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

using mdmvis::FieldMap;
using mdmvis::Vec3;

struct DipoleProjection {
  std::string label;
  int color = kBlack;
  FieldMap map;
};

Vec3 TopDownPoint(const FieldMap& map, int ix, int iy, int iz) {
  const double radius = mdmvis::HeaderDouble(map, "dipole_reference_radius_cm");
  const double phiDeg = mdmvis::HeaderDouble(map, "dipole_sector_angle_deg");
  const double alphaDeg = mdmvis::HeaderDouble(map, "dipole_alpha_deg");
  const double betaDeg = mdmvis::HeaderDouble(map, "dipole_beta_deg");

  if (map.h.coordinate_system == "dipole_entrance_frame") {
    return mdmvis::EntranceFrameToDipoleLocal(map.X(ix), map.Y(iy), map.Z(iz),
                                             alphaDeg);
  }
  if (map.h.coordinate_system == "dipole_sector_frame") {
    return mdmvis::SectorFrameToDipoleLocal(map.X(ix), map.Y(iy), map.Z(iz),
                                           radius);
  }
  if (map.h.coordinate_system == "dipole_exit_frame") {
    return mdmvis::ExitFrameToDipoleLocal(map.X(ix), map.Y(iy), map.Z(iz),
                                          radius, phiDeg, alphaDeg, betaDeg);
  }
  throw std::runtime_error("expected a dipole field map: " + map.h.magnet);
}

Vec3 StoragePointFromTopDown(const FieldMap& map,
                             double x,
                             double y,
                             double z) {
  const double radius = mdmvis::HeaderDouble(map, "dipole_reference_radius_cm");
  const double phiDeg = mdmvis::HeaderDouble(map, "dipole_sector_angle_deg");
  const double alphaDeg = mdmvis::HeaderDouble(map, "dipole_alpha_deg");
  const double betaDeg = mdmvis::HeaderDouble(map, "dipole_beta_deg");

  if (map.h.coordinate_system == "dipole_entrance_frame") {
    return mdmvis::DipoleLocalToEntranceFrame(x, y, z, alphaDeg);
  }
  if (map.h.coordinate_system == "dipole_sector_frame") {
    return mdmvis::DipoleLocalToSectorFrame(x, y, z, radius);
  }
  if (map.h.coordinate_system == "dipole_exit_frame") {
    return mdmvis::DipoleLocalToExitFrame(x, y, z, radius, phiDeg, alphaDeg,
                                          betaDeg);
  }
  throw std::runtime_error("expected a dipole field map: " + map.h.magnet);
}

void ExpandBounds(const Vec3& point,
                  double& minX,
                  double& maxX,
                  double& minZ,
                  double& maxZ) {
  minX = std::min(minX, point.x);
  maxX = std::max(maxX, point.x);
  minZ = std::min(minZ, point.z);
  maxZ = std::max(maxZ, point.z);
}

TGraph* MakeOutline(const FieldMap& map, int color) {
  std::vector<double> xs;
  std::vector<double> zs;
  const int iy = mdmvis::NearestIndex(0.0, map.h.origin_cm.y, map.h.step_cm.y,
                                      map.h.ny);

  const auto addPoint = [&](int ix, int iz) {
    const Vec3 point = TopDownPoint(map, ix, iy, iz);
    xs.push_back(point.x);
    zs.push_back(point.z);
  };

  if (map.h.coordinate_system == "dipole_sector_frame") {
    for (int iz = 0; iz < map.h.nz; ++iz) addPoint(0, iz);
    for (int ix = 1; ix < map.h.nx; ++ix) addPoint(ix, map.h.nz - 1);
    for (int iz = map.h.nz - 2; iz >= 0; --iz) addPoint(map.h.nx - 1, iz);
    for (int ix = map.h.nx - 2; ix >= 0; --ix) addPoint(ix, 0);
    addPoint(0, 0);
  } else {
    addPoint(0, 0);
    addPoint(map.h.nx - 1, 0);
    addPoint(map.h.nx - 1, map.h.nz - 1);
    addPoint(0, map.h.nz - 1);
    addPoint(0, 0);
  }

  TGraph* graph = new TGraph(static_cast<int>(xs.size()), xs.data(), zs.data());
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  graph->SetFillStyle(0);
  return graph;
}

}  // namespace

void PlotDipoleTopDown(
    const char* entrancePath = "DipoleEntrance.bin",
    const char* sectorPath = "DipoleSector.bin",
    const char* exitPath = "DipoleExit.bin",
    double yCm = 0.0) {
  using namespace mdmvis;
  UseFieldMapStyle();

  std::vector<DipoleProjection> regions{
      {"entrance", kBlue + 1, LoadFieldMap(entrancePath)},
      {"sector", kBlack, LoadFieldMap(sectorPath)},
      {"exit", kRed + 1, LoadFieldMap(exitPath)},
  };

  double minX = 1.0e99;
  double maxX = -1.0e99;
  double minZ = 1.0e99;
  double maxZ = -1.0e99;
  double minStep = 1.0e99;

  for (const DipoleProjection& region : regions) {
    const FieldMap& map = region.map;
    minStep = std::min(minStep, std::min(map.h.step_cm.x, map.h.step_cm.z));
    const int iy = NearestIndex(yCm, map.h.origin_cm.y, map.h.step_cm.y,
                                map.h.ny);
    for (int iz = 0; iz < map.h.nz; ++iz) {
      for (int ix = 0; ix < map.h.nx; ++ix) {
        ExpandBounds(TopDownPoint(map, ix, iy, iz), minX, maxX, minZ, maxZ);
      }
    }
  }

  const double margin = 0.05 * std::max(maxX - minX, maxZ - minZ);
  minX -= margin;
  maxX += margin;
  minZ -= margin;
  maxZ += margin;

  const int nX = std::max(100, static_cast<int>((maxX - minX) / minStep) + 1);
  const int nZ = std::max(100, static_cast<int>((maxZ - minZ) / minStep) + 1);

  TH2D* hBy = new TH2D("h_dipole_top_down_by",
                       "Dipole top-down By;dipole local x [cm];dipole local z [cm];By [T]",
                       nX, minX, maxX, nZ, minZ, maxZ);
  hBy->SetContour(100);

  std::vector<const DipoleProjection*> sampleOrder{&regions[0], &regions[2],
                                                   &regions[1]};
  for (int ix = 1; ix <= hBy->GetNbinsX(); ++ix) {
    const double x = hBy->GetXaxis()->GetBinCenter(ix);
    for (int iz = 1; iz <= hBy->GetNbinsY(); ++iz) {
      const double z = hBy->GetYaxis()->GetBinCenter(iz);
      for (const DipoleProjection* region : sampleOrder) {
        const FieldMap& map = region->map;
        const Vec3 storage = StoragePointFromTopDown(map, x, yCm, z);
        if (!InsideMap(map, storage.x, storage.y, storage.z)) {
          continue;
        }
        hBy->SetBinContent(ix, iz, FieldAt(map, storage.x, storage.y,
                                           storage.z)
                                       .y);
        break;
      }
    }
  }

  TCanvas* canvas =
      new TCanvas("c_dipole_top_down",
                  "Dipole top-down field in physical coordinates", 900, 800);
  hBy->Draw("COLZ");

  TLegend* legend = new TLegend(0.72, 0.75, 0.90, 0.90);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  for (const DipoleProjection& region : regions) {
    TGraph* outline = MakeOutline(region.map, region.color);
    outline->Draw("L same");
    legend->AddEntry(outline, region.label.c_str(), "l");
  }
  legend->Draw();

  canvas->Update();
}
