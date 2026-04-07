#include "MDMFieldMapTrace.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "json.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: ./MDMFieldMapTraceExample <config-file>" << std::endl;
    return 1;
  }

  try {
    Json::Value config;
    const std::string configFileName = argv[1];
    std::ifstream configStream(configFileName.c_str());
    if (!configStream) {
      std::cerr << "ERROR: Unable to open config file: " << configFileName
                << std::endl;
      return 1;
    }
    configStream >> config;

    bool usingProbe = false;
    double mdmAngle = 0.0;
    double mdmDipoleField = 0.0;
    double mdmDipoleProbe = 0.0;
    double mdmMultipoleProbe = 0.0;
    double scatteredMass = 0.0;
    double scatteredCharge = 0.0;
    double scatteredEnergy = 0.0;
    std::string multipoleMapPath = "Multipole.bin";
    std::string dipoleMapPath = "Dipole.bin";
    std::string dipoleEntranceMapPath = "DipoleEntrance.bin";
    std::string dipoleSectorMapPath = "DipoleSector.bin";
    std::string dipoleExitMapPath = "DipoleExit.bin";
    bool hasSplitDipolePaths = false;
    bool hasLegacyDipolePath = false;
    std::vector<double> scatteredAngles;

    for (Json::Value::iterator it = config.begin(); it != config.end(); ++it) {
      if (false) {
      } else if (it.key().asString() == "usingProbe") {
        usingProbe = it->asBool();
      } else if (it.key().asString() == "mdmAngle") {
        mdmAngle = it->asDouble();
      } else if (it.key().asString() == "mdmDipoleField") {
        mdmDipoleField = it->asDouble();
      } else if (it.key().asString() == "mdmDipoleProbe") {
        mdmDipoleProbe = it->asDouble();
      } else if (it.key().asString() == "mdmMultipoleProbe") {
        mdmMultipoleProbe = it->asDouble();
      } else if (it.key().asString() == "scatteredMass") {
        scatteredMass = it->asDouble();
      } else if (it.key().asString() == "scatteredCharge") {
        scatteredCharge = it->asDouble();
      } else if (it.key().asString() == "scatteredEnergy") {
        scatteredEnergy = it->asDouble();
      } else if (it.key().asString() == "scatteredAngles") {
        for (unsigned int index = 0; index < it->size(); ++index) {
          scatteredAngles.push_back((*it)[index].asDouble());
        }
      } else if (it.key().asString() == "multipoleMapPath") {
        multipoleMapPath = it->asString();
      } else if (it.key().asString() == "dipoleMapPath") {
        dipoleMapPath = it->asString();
        hasLegacyDipolePath = true;
      } else if (it.key().asString() == "dipoleEntranceMapPath") {
        dipoleEntranceMapPath = it->asString();
        hasSplitDipolePaths = true;
      } else if (it.key().asString() == "dipoleSectorMapPath") {
        dipoleSectorMapPath = it->asString();
        hasSplitDipolePaths = true;
      } else if (it.key().asString() == "dipoleExitMapPath") {
        dipoleExitMapPath = it->asString();
        hasSplitDipolePaths = true;
      }
    }

    MDMFieldMapTrace trace;
    trace.SetMDMAngle(mdmAngle);
    if (usingProbe) {
      trace.SetMDMProbe(mdmDipoleProbe, mdmMultipoleProbe);
    } else {
      trace.SetMDMDipoleField(mdmDipoleField);
    }
    trace.SetScatteredMass(scatteredMass);
    trace.SetScatteredCharge(scatteredCharge);
    trace.SetScatteredEnergy(scatteredEnergy);
    if (hasSplitDipolePaths || !hasLegacyDipolePath) {
      trace.LoadFieldMaps(multipoleMapPath, dipoleEntranceMapPath,
                          dipoleSectorMapPath, dipoleExitMapPath);
    } else {
      trace.LoadFieldMaps(multipoleMapPath, dipoleMapPath);
    }

    for (const double angle : scatteredAngles) {
      trace.SetScatteredAngle(angle);
      trace.SendRay();

      double x1 = 0.0;
      double y1 = 0.0;
      double angX1 = 0.0;
      double angY1 = 0.0;
      trace.GetPositionAngleFirstWire(x1, y1, angX1, angY1);
      std::printf(
          "Scattered Angle: %.2fdeg  X1: %.2fcm  Y1: %.2fcm  AngX1: %.2fdeg  "
          "AngY1: %.2fdeg\n",
          trace.GetScatteredAngle(), x1, y1, angX1, angY1);
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
}
