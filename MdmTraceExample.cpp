#include "MdmRayScan.h"
#include "MdmTrace.h"
#include "json.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

double GetDouble(const Json::Value& config, const char* key) {
  return config.isMember(key) ? config[key].asDouble() : 0.0;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Please specify the config file.\n"
              << "Usage: ./MdmTraceExample <config-file>" << std::endl;
    return 0;
  }

  const std::string configFileName = argv[1];
  std::ifstream configStream(configFileName.c_str());
  if (!configStream) {
    std::cerr << "ERROR: Unable to open config file: " << configFileName
              << std::endl;
    return 1;
  }

  Json::Value config;
  configStream >> config;
  const bool usingProbe =
      config.isMember("usingProbe") && config["usingProbe"].asBool();
  const std::vector<mdm::RayInput> rays = mdm::ParseRayInputs(config);

  MdmTrace trace;
  trace.SetMdmAngle(GetDouble(config, "mdmAngle"));

  if (usingProbe) {
    trace.SetMdmProbe(GetDouble(config, "mdmDipoleProbe"),
                      GetDouble(config, "mdmMultipoleProbe"));
  } else {
    trace.SetMdmDipoleField(GetDouble(config, "mdmDipoleField"));
  }

  trace.SetScatteredMass(GetDouble(config, "scatteredMass"));
  trace.SetScatteredCharge(GetDouble(config, "scatteredCharge"));

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
