#include "FieldMapCommon.C"

#include <TCanvas.h>
#include <TH3D.h>
#include <TPolyLine3D.h>

#include <algorithm>
#include <cmath>
#include <string>

void PlotField3D(const char* mapPath = "DipoleSector.bin",
                 int maxVectors = 1200) {
  using namespace mdmvis;
  UseFieldMapStyle();

  const FieldMap map = LoadFieldMap(mapPath);
  const double totalNodes = static_cast<double>(map.Size());
  int stride =
      std::max(1, static_cast<int>(std::ceil(std::cbrt(totalNodes /
                                                       std::max(1, maxVectors)))));

  double maxB = 0.0;
  for (int iz = 0; iz < map.h.nz; iz += stride) {
    for (int iy = 0; iy < map.h.ny; iy += stride) {
      for (int ix = 0; ix < map.h.nx; ix += stride) {
        maxB = std::max(maxB, map.Bmag(ix, iy, iz));
      }
    }
  }

  const std::string canvasName = UniqueName("c_field_3d", map);
  TCanvas* canvas =
      new TCanvas(canvasName.c_str(),
                  (map.h.magnet + " sparse 3D field").c_str(), 900, 750);

  TH3D* frame =
      new TH3D(UniqueName("h3_field_frame", map).c_str(),
               (map.h.magnet + " sparse 3D field;" + XAxisTitle(map) +
                ";y [cm];" + LongitudinalAxisTitle(map))
                   .c_str(),
               2, AxisLowEdge(map.h.origin_cm.x, map.h.step_cm.x),
               AxisHighEdge(map.h.origin_cm.x, map.h.step_cm.x, map.h.nx), 2,
               AxisLowEdge(map.h.origin_cm.y, map.h.step_cm.y),
               AxisHighEdge(map.h.origin_cm.y, map.h.step_cm.y, map.h.ny), 2,
               AxisLowEdge(map.h.origin_cm.z, map.h.step_cm.z),
               AxisHighEdge(map.h.origin_cm.z, map.h.step_cm.z, map.h.nz));
  frame->SetStats(0);
  frame->Draw("BOX");

  const double minStep =
      std::min(map.h.step_cm.x, std::min(map.h.step_cm.y, map.h.step_cm.z));
  const double baseLength = 0.75 * minStep * stride;

  for (int iz = 0; iz < map.h.nz; iz += stride) {
    for (int iy = 0; iy < map.h.ny; iy += stride) {
      for (int ix = 0; ix < map.h.nx; ix += stride) {
        const double mag = map.Bmag(ix, iy, iz);
        if (mag <= 0.0 || !std::isfinite(mag)) {
          continue;
        }

        const double x0 = map.X(ix);
        const double y0 = map.Y(iy);
        const double z0 = map.Z(iz);
        const double length = baseLength * (0.25 + 0.75 * mag / maxB);
        const double scale = 0.5 * length / mag;
        const double x1 = x0 + scale * map.Bx(ix, iy, iz);
        const double y1 = y0 + scale * map.By(ix, iy, iz);
        const double z1 = z0 + scale * map.Bz(ix, iy, iz);
        const double x2 = x0 - scale * map.Bx(ix, iy, iz);
        const double y2 = y0 - scale * map.By(ix, iy, iz);
        const double z2 = z0 - scale * map.Bz(ix, iy, iz);

        TPolyLine3D* line = new TPolyLine3D(2);
        line->SetPoint(0, x2, y2, z2);
        line->SetPoint(1, x1, y1, z1);
        line->SetLineColor(PaletteColor(mag, maxB));
        line->SetLineWidth(2);
        line->Draw("same");
      }
    }
  }

  canvas->Update();
}
