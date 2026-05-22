#include "FieldMapCommon.C"

#include <TArrow.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TPad.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {

int PlaneMaxIndexInRange(const mdmvis::FieldMap& map,
                         double minZ,
                         double maxZ) {
  int bestIndex = mdmvis::NearestIndex(0.5 * (minZ + maxZ), map.h.origin_cm.z,
                                       map.h.step_cm.z, map.h.nz);
  double best = -1.0;
  for (int iz = 0; iz < map.h.nz; ++iz) {
    const double z = map.Z(iz);
    if (z < minZ || z > maxZ) {
      continue;
    }
    double planeMax = 0.0;
    for (int iy = 0; iy < map.h.ny; ++iy) {
      for (int ix = 0; ix < map.h.nx; ++ix) {
        planeMax = std::max(planeMax, map.Bmag(ix, iy, iz));
      }
    }
    if (planeMax > best) {
      best = planeMax;
      bestIndex = iz;
    }
  }
  return bestIndex;
}

std::vector<int> MultipoleSlicePlanes(const mdmvis::FieldMap& map) {
  const std::vector<double> transitions =
      mdmvis::HeaderDoubles(map, "multipole_transition_planes_cm");
  if (transitions.size() < 2) {
    return {mdmvis::RepresentativeZIndex(map)};
  }

  const double zMin = map.Z(0);
  const double zMax = map.Z(map.h.nz - 1);
  return {PlaneMaxIndexInRange(map, zMin, transitions[0]),
          mdmvis::NearestIndex(0.5 * (transitions[0] + transitions[1]),
                               map.h.origin_cm.z, map.h.step_cm.z, map.h.nz),
          PlaneMaxIndexInRange(map, transitions[1], zMax)};
}

TH2D* BuildMultipoleSlice(const mdmvis::FieldMap& map,
                          int iz,
                          const std::string& label) {
  const double z = map.Z(iz);
  const int displayScale = 4;
  const int nx = displayScale * (map.h.nx - 1) + 1;
  const int ny = displayScale * (map.h.ny - 1) + 1;
  std::ostringstream title;
  title << std::fixed << std::setprecision(2) << map.h.magnet << " " << label
        << " |B| at z = " << z << " cm;x [cm];y [cm];|B| [T]";

  TH2D* hist =
      new TH2D(mdmvis::UniqueName("h_multipole_xy_bmag", map).c_str(),
               title.str().c_str(), nx,
               mdmvis::AxisLowEdge(map.h.origin_cm.x, map.h.step_cm.x),
               mdmvis::AxisHighEdge(map.h.origin_cm.x, map.h.step_cm.x,
                                     map.h.nx),
               ny,
               mdmvis::AxisLowEdge(map.h.origin_cm.y, map.h.step_cm.y),
               mdmvis::AxisHighEdge(map.h.origin_cm.y, map.h.step_cm.y,
                                     map.h.ny));
  hist->SetMinimum(0.0);
  for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
    const double y = hist->GetYaxis()->GetBinCenter(iy);
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
      const double x = hist->GetXaxis()->GetBinCenter(ix);
      const mdmvis::Vec3 field =
          mdmvis::InsideMap(map, x, y, z) ? mdmvis::FieldAt(map, x, y, z)
                                          : mdmvis::Vec3{};
      hist->SetBinContent(ix, iy, mdmvis::Bmag(field));
    }
  }
  return hist;
}

double MaxTransverseField(const mdmvis::FieldMap& map,
                          const std::vector<int>& planes) {
  double maxTransverse = 0.0;
  for (int iz : planes) {
    const double z = map.Z(iz);
    for (int iy = 0; iy < map.h.ny; ++iy) {
      for (int ix = 0; ix < map.h.nx; ++ix) {
        const mdmvis::Vec3 field =
            mdmvis::FieldAt(map, map.X(ix), map.Y(iy), z);
        maxTransverse =
            std::max(maxTransverse, std::sqrt(field.x * field.x +
                                              field.y * field.y));
      }
    }
  }
  return maxTransverse;
}

void DrawMultipoleArrows(const mdmvis::FieldMap& map,
                         int iz,
                         double maxTransverse) {
  const double z = map.Z(iz);
  const int nArrow = 17;
  const double xMin = map.h.origin_cm.x;
  const double xMax = map.X(map.h.nx - 1);
  const double yMin = map.h.origin_cm.y;
  const double yMax = map.Y(map.h.ny - 1);
  const double dxGrid = (xMax - xMin) / static_cast<double>(nArrow - 1);
  const double dyGrid = (yMax - yMin) / static_cast<double>(nArrow - 1);
  const double arrowBase = 0.45 * std::min(dxGrid, dyGrid);
  for (int iy = 0; iy < nArrow; ++iy) {
    const double y = yMin + dyGrid * iy;
    for (int ix = 0; ix < nArrow; ++ix) {
      const double x = xMin + dxGrid * ix;
      if (!mdmvis::InsideMap(map, x, y, z)) {
        continue;
      }
      const mdmvis::Vec3 field = mdmvis::FieldAt(map, x, y, z);
      const double transverse =
          std::sqrt(field.x * field.x + field.y * field.y);
      if (transverse <= 0.0 || maxTransverse <= 0.0) {
        continue;
      }
      const double length = arrowBase * transverse / maxTransverse;
      const double dx = length * field.x / transverse;
      const double dy = length * field.y / transverse;
      TArrow* arrow = new TArrow(x - 0.5 * dx, y - 0.5 * dy, x + 0.5 * dx,
                                 y + 0.5 * dy, 0.0035, "|>");
      arrow->SetLineColor(kWhite);
      arrow->SetFillColor(kWhite);
      arrow->SetLineWidth(2);
      arrow->Draw();
    }
  }
}

}  // namespace

void PlotFieldSlices(const char* mapPath = "DipoleSector.bin",
                     double yCm = 0.0) {
  using namespace mdmvis;
  UseFieldMapStyle();

  const FieldMap map = LoadFieldMap(mapPath);
  if (IsMultipole(map)) {
    const std::vector<int> planes = MultipoleSlicePlanes(map);
    const std::vector<std::string> labels{"entrance fringe", "uniform",
                                          "exit fringe"};
    std::vector<TH2D*> slices;
    double maxField = 0.0;
    for (std::size_t i = 0; i < planes.size(); ++i) {
      TH2D* hist = BuildMultipoleSlice(
          map, planes[i], i < labels.size() ? labels[i] : "slice");
      slices.push_back(hist);
      maxField = std::max(maxField, hist->GetMaximum());
    }
    const double maxTransverse = MaxTransverseField(map, planes);
    for (TH2D* hist : slices) {
      hist->SetMinimum(0.0);
      hist->SetMaximum(maxField);
    }

    TCanvas* canvas =
        new TCanvas(UniqueName("c_multipole_xy_slices", map).c_str(),
                    (map.h.magnet + " transverse field slices").c_str(), 1800,
                    650);
    canvas->Divide(static_cast<int>(slices.size()), 1);
    for (std::size_t i = 0; i < slices.size(); ++i) {
      canvas->cd(static_cast<int>(i + 1));
      gPad->SetRightMargin(0.16);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.10);
      gPad->SetFixedAspectRatio();
      slices[i]->Draw("COLZ");
      DrawMultipoleArrows(map, planes[i], maxTransverse);
    }

    canvas->Update();
    return;
  }

  const int iy = NearestIndex(yCm, map.h.origin_cm.y, map.h.step_cm.y, map.h.ny);
  const double y = map.Y(iy);

  std::ostringstream suffix;
  suffix << std::fixed << std::setprecision(2) << y;

  const std::string magName = UniqueName("h_bmag_slice", map);
  TH2D* hMag =
      new TH2D(magName.c_str(),
               (map.h.magnet + " |B| at y = " + suffix.str() +
                " cm;" + XAxisTitle(map) + ";" + LongitudinalAxisTitle(map) +
                ";|B| [T]")
                   .c_str(),
               map.h.nx, AxisLowEdge(map.h.origin_cm.x, map.h.step_cm.x),
               AxisHighEdge(map.h.origin_cm.x, map.h.step_cm.x, map.h.nx),
               map.h.nz, AxisLowEdge(map.h.origin_cm.z, map.h.step_cm.z),
               AxisHighEdge(map.h.origin_cm.z, map.h.step_cm.z, map.h.nz));

  const std::string byName = UniqueName("h_by_slice", map);
  TH2D* hBy =
      new TH2D(byName.c_str(),
               (map.h.magnet + " By at y = " + suffix.str() +
                " cm;" + XAxisTitle(map) + ";" + LongitudinalAxisTitle(map) +
                ";By [T]")
                   .c_str(),
               map.h.nx, AxisLowEdge(map.h.origin_cm.x, map.h.step_cm.x),
               AxisHighEdge(map.h.origin_cm.x, map.h.step_cm.x, map.h.nx),
               map.h.nz, AxisLowEdge(map.h.origin_cm.z, map.h.step_cm.z),
               AxisHighEdge(map.h.origin_cm.z, map.h.step_cm.z, map.h.nz));

  for (int iz = 0; iz < map.h.nz; ++iz) {
    for (int ix = 0; ix < map.h.nx; ++ix) {
      const Vec3 field = FieldAt(map, map.X(ix), y, map.Z(iz));
      hMag->SetBinContent(ix + 1, iz + 1, Bmag(field));
      hBy->SetBinContent(ix + 1, iz + 1, field.y);
    }
  }

  const std::string canvasName = UniqueName("c_field_slices", map);
  TCanvas* canvas =
      new TCanvas(canvasName.c_str(),
                  (map.h.magnet + " field slices").c_str(), 1300, 600);
  canvas->Divide(2, 1);
  canvas->cd(1);
  hMag->Draw("COLZ");
  canvas->cd(2);
  hBy->Draw("COLZ");
  canvas->Update();
}
