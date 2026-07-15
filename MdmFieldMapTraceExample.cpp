#include "MdmConfig.h"
#include "MdmFieldMapTrace.h"
#include "MdmIonConfig.h"
#include "MdmRayScan.h"
#include "json.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

double GetDouble(const Json::Value& config, const char* key) {
  return config.isMember(key) ? config[key].asDouble() : 0.0;
}

std::string GetString(const Json::Value& config,
                      const char* key,
                      const char* fallback) {
  return config.isMember(key) ? config[key].asString() : fallback;
}

}  // namespace

int Run(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: ./MdmFieldMapTraceExample <config-file>" << std::endl;
    return 1;
  }

  const std::string configFileName = argv[1];
  const Json::Value config = mdm::ReadConfig(configFileName);
  const bool usingProbe = mdm::RequireBoolean(config, "usingProbe");
  const MdmIon ion = mdm::ParseScatteredIon(config);
  const std::vector<mdm::RayInput> rays = mdm::ParseRayInputs(config);
  const mdm::RayScanCounts counts = mdm::CountRayScan(config, rays.size());
  std::cerr << "Energy values: " << counts.energies << "\n"
            << "Horizontal-angle values: " << counts.horizontalAngles << "\n"
            << "Vertical-angle values: " << counts.verticalAngles << "\n"
            << "Total rays: " << counts.rays << "\n";

  MdmFieldMapTrace trace;
  trace.SetMdmAngle(GetDouble(config, "mdmAngle"));
  if (usingProbe) {
    trace.SetMdmProbe(GetDouble(config, "mdmDipoleProbe"),
                      GetDouble(config, "mdmMultipoleProbe"));
  } else {
    trace.SetMdmDipoleField(GetDouble(config, "mdmDipoleField"));
  }
  trace.SetScatteredIon(ion);

  const std::string multipolePath =
      GetString(config, "multipoleMapPath", "Multipole.bin");
  trace.LoadFieldMaps(multipolePath,
                      GetString(config, "dipoleEntranceMapPath",
                                "DipoleEntrance.bin"),
                      GetString(config, "dipoleSectorMapPath",
                                "DipoleSector.bin"),
                      GetString(config, "dipoleExitMapPath", "DipoleExit.bin"));

  for (const mdm::RayInput& ray : rays) {
    trace.SetScatteredAngle(ray.xDeg, ray.yDeg);
    trace.SetScatteredEnergy(ray.energyMeV);
    trace.SendRay();

    double x1 = 0.0;
    double y1 = 0.0;
    double angX1 = 0.0;
    double angY1 = 0.0;
    trace.GetPositionAngleFirstWire(x1, y1, angX1, angY1);
    std::printf("Scattered Angle X: %.4fdeg  Scattered Angle Y: %.4fdeg  "
                "Scattered Energy: %.4fMeV  X1: %.4fcm  Y1: %.4fcm  "
                "AngX1: %.4fdeg  AngY1: %.4fdeg\n",
                ray.xDeg, ray.yDeg, ray.energyMeV, x1, y1, angX1, angY1);
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
