#include "FieldMapCommon.C"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMultiGraph.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

void PlotFieldProfiles(const char* mapPath = "DipoleSector.bin",
                       double xCm = 0.0,
                       double yCm = 0.0) {
  using namespace mdmvis;
  UseFieldMapStyle();

  const FieldMap map = LoadFieldMap(mapPath);
  if (IsMultipole(map)) {
    const double rMax =
        map.h.extra.count("multipole_aperture_radius_cm") > 0
            ? HeaderDouble(map, "multipole_aperture_radius_cm")
            : std::min(std::abs(map.h.origin_cm.x),
                       std::abs(map.h.origin_cm.y));
    const int zDisplayScale = 1;
    const int nZ = zDisplayScale * (map.h.nz - 1) + 1;
    std::vector<double> z(nZ);
    const std::vector<double> radiusFractions{0.25, 0.50, 0.75, 0.92};
    std::vector<double> radii;
    for (double fraction : radiusFractions) {
      radii.push_back(fraction * rMax);
    }
    const std::vector<int> colors{kBlue + 1, kGreen + 2, kOrange + 7,
                                  kRed + 1};
    std::vector<std::vector<double>> samples(radii.size(),
                                             std::vector<double>(nZ));

    const int nPhi = 360;
    for (int iz = 0; iz < nZ; ++iz) {
      const double fraction =
          static_cast<double>(iz) / static_cast<double>(nZ - 1);
      z[iz] = map.Z(0) + fraction * (map.Z(map.h.nz - 1) - map.Z(0));
      for (std::size_t ir = 0; ir < radii.size(); ++ir) {
        double sum = 0.0;
        int count = 0;
        for (int iphi = 0; iphi < nPhi; ++iphi) {
          const double phi =
              2.0 * 3.14159265358979323846 * static_cast<double>(iphi) /
              static_cast<double>(nPhi);
          const double x = radii[ir] * std::cos(phi);
          const double y = radii[ir] * std::sin(phi);
          if (!InsideMap(map, x, y, z[iz])) {
            continue;
          }
          sum += Bmag(FieldAt(map, x, y, z[iz]));
          ++count;
        }
        samples[ir][iz] = count > 0 ? sum / static_cast<double>(count) : 0.0;
      }
    }

    TMultiGraph* graphs = new TMultiGraph();
    graphs->SetName(UniqueName("mg_multipole_z_profile", map).c_str());
    graphs->SetTitle(
        (map.h.magnet +
         " fixed-radius field strength vs z;z [cm];angular-average |B| [T]")
            .c_str());

    std::vector<TGraph*> radiusGraphs;
    for (std::size_t ir = 0; ir < radii.size(); ++ir) {
      TGraph* graph = new TGraph(nZ, z.data(), samples[ir].data());
      graph->SetName(UniqueName("g_multipole_radius_profile", map).c_str());
      graph->SetLineColor(colors[ir]);
      graph->SetLineWidth(2);
      graphs->Add(graph, "L");
      radiusGraphs.push_back(graph);
    }

    TCanvas* canvas =
        new TCanvas(UniqueName("c_multipole_z_profile", map).c_str(),
                    (map.h.magnet + " field strength vs z").c_str(), 900,
                    650);
    graphs->Draw("A");

    TLegend* legend = new TLegend(0.62, 0.68, 0.90, 0.90);
    for (std::size_t ir = 0; ir < radii.size(); ++ir) {
      std::ostringstream label;
      label << std::fixed << std::setprecision(2) << "r = "
            << radiusFractions[ir] << "R";
      legend->AddEntry(radiusGraphs[ir], label.str().c_str(), "l");
    }

    const std::vector<double> transitions =
        map.h.extra.count("multipole_transition_planes_cm") > 0
            ? HeaderDoubles(map, "multipole_transition_planes_cm")
            : std::vector<double>{};
    double yMax = 0.0;
    for (const std::vector<double>& values : samples) {
      yMax = std::max(yMax, *std::max_element(values.begin(), values.end()));
    }
    for (double transition : transitions) {
      if (transition < z.front() || transition > z.back()) {
          continue;
        }
      TLine* line = new TLine(transition, 0.0, transition, yMax);
      line->SetLineColor(kGray + 2);
      line->SetLineStyle(2);
      line->Draw();
    }
    legend->Draw();
    canvas->Update();
    return;
  }

  const int ix = NearestIndex(xCm, map.h.origin_cm.x, map.h.step_cm.x, map.h.nx);
  const int iy = NearestIndex(yCm, map.h.origin_cm.y, map.h.step_cm.y, map.h.ny);

  std::vector<double> z(map.h.nz);
  std::vector<double> bx(map.h.nz);
  std::vector<double> by(map.h.nz);
  std::vector<double> bz(map.h.nz);
  std::vector<double> bmag(map.h.nz);

  for (int iz = 0; iz < map.h.nz; ++iz) {
    z[iz] = map.Z(iz);
    const Vec3 field = FieldAt(map, map.X(ix), map.Y(iy), z[iz]);
    bx[iz] = field.x;
    by[iz] = field.y;
    bz[iz] = field.z;
    bmag[iz] = Bmag(field);
  }

  TGraph* gBx = new TGraph(map.h.nz, z.data(), bx.data());
  TGraph* gBy = new TGraph(map.h.nz, z.data(), by.data());
  TGraph* gBz = new TGraph(map.h.nz, z.data(), bz.data());
  TGraph* gMag = new TGraph(map.h.nz, z.data(), bmag.data());

  gBx->SetName(UniqueName("g_bx_profile", map).c_str());
  gBy->SetName(UniqueName("g_by_profile", map).c_str());
  gBz->SetName(UniqueName("g_bz_profile", map).c_str());
  gMag->SetName(UniqueName("g_bmag_profile", map).c_str());

  gBx->SetLineColor(kRed + 1);
  gBy->SetLineColor(kBlue + 1);
  gBz->SetLineColor(kGreen + 2);
  gMag->SetLineColor(kBlack);
  gBx->SetLineWidth(2);
  gBy->SetLineWidth(2);
  gBz->SetLineWidth(2);
  gMag->SetLineWidth(3);

  std::ostringstream point;
  point << std::fixed << std::setprecision(2) << " at " << XAxisTitle(map)
        << " = " << map.X(ix) << ", y = " << map.Y(iy) << " cm";

  TMultiGraph* graphs = new TMultiGraph();
  graphs->SetName(UniqueName("mg_field_profile", map).c_str());
  graphs->SetTitle((map.h.magnet + " field profile" + point.str() + ";" +
                    LongitudinalAxisTitle(map) + ";B [T]")
                       .c_str());
  graphs->Add(gBx, "L");
  graphs->Add(gBy, "L");
  graphs->Add(gBz, "L");
  graphs->Add(gMag, "L");

  TCanvas* canvas =
      new TCanvas(UniqueName("c_field_profile", map).c_str(),
                  (map.h.magnet + " field profile").c_str(), 900, 650);
  graphs->Draw("A");

  const std::string coord = map.h.coordinate_system;
  if (coord == "dipole_sector_frame") {
    const double radius = HeaderDouble(map, "dipole_reference_radius_cm");
    const double phiDeg = HeaderDouble(map, "dipole_sector_angle_deg");
    const double z12 = HeaderDouble(map, "dipole_z12_cm");
    const double z21 = HeaderDouble(map, "dipole_z21_cm");
    const double phi = DegreesToRadians(phiDeg);
    const double profileRadius = radius + map.X(ix);
    const double entranceToSector =
        radius * std::asin(-z12 / profileRadius);
    const double sectorToExit = radius * (phi + z21 / profileRadius);
    const double yMax =
        std::max(*std::max_element(bmag.begin(), bmag.end()),
                 std::max(*std::max_element(bx.begin(), bx.end()),
                          *std::max_element(by.begin(), by.end())));
    const double yMin =
        std::min(*std::min_element(bmag.begin(), bmag.end()),
                 std::min(*std::min_element(bx.begin(), bx.end()),
                          *std::min_element(by.begin(), by.end())));
    const std::vector<double> transitions{entranceToSector, sectorToExit};
    const std::vector<std::string> labels{"entrance -> sector",
                                          "sector -> exit"};
    for (std::size_t i = 0; i < transitions.size(); ++i) {
      const double transition = transitions[i];
      if (transition < z.front() || transition > z.back()) {
        continue;
      }
      TLine* line = new TLine(transition, yMin, transition, yMax);
      line->SetLineColor(kGray + 2);
      line->SetLineStyle(2);
      line->Draw();

      TLatex* label = new TLatex(
          transition + 0.01 * (z.back() - z.front()),
          yMin + 0.50 * (yMax - yMin), labels[i].c_str());
      label->SetTextAngle(0.0);
      label->SetTextSize(0.025);
      label->SetTextColor(kBlack);
      label->Draw();
    }
  }

  TLegend* legend = new TLegend(0.78, 0.72, 0.92, 0.90);
  legend->AddEntry(gBx, "Bx", "l");
  legend->AddEntry(gBy, "By", "l");
  legend->AddEntry(gBz, "Bz", "l");
  legend->AddEntry(gMag, "|B|", "l");
  legend->Draw();

  canvas->Update();
}
