#pragma once

#include "MdmFieldMap.h"
#include "MdmIon.h"

#include <string>

class MdmFieldMapTrace {
 public:
  MdmFieldMapTrace();

  void LoadFieldMaps(const std::string& multipolePath,
                     const std::string& dipoleEntrancePath,
                     const std::string& dipoleSectorPath,
                     const std::string& dipoleExitPath);

  void SetMdmAngle(double angleDeg);
  double GetMdmAngle() const;

  void SetMdmProbe(double dipoleProbe, double multipoleProbe);
  void SetMdmDipoleField(double dipoleField);

  void SetScatteredIon(const MdmIon& ion);

  void SetScatteredEnergy(double energyMeV);
  double GetScatteredEnergy() const;

  void SetInitialPosition(double xCm, double yCm);

  void SetScatteredAngle(double xAngleDeg);
  void SetScatteredAngle(double xAngleDeg, double yAngleDeg);
  double GetScatteredAngle() const;

  void SendRay();

  void GetPositionAngleFirstWire(double& posX,
                                 double& posY,
                                 double& angX,
                                 double& angY) const;
  double GetTimeOfFlightSeconds() const;

 private:
  void ValidateLoadedMaps() const;

  double mdmAngleDeg_ = 0.0;
  double ionMassMeV_ = 0.0;
  double ionChargeState_ = 0.0;
  double scatteredEnergyMeV_ = 0.0;
  double initialXcm_ = 0.0;
  double initialYcm_ = 0.0;
  double scatteredAnglesDeg_[2]{0.0, 0.0};

  double requestedDipoleProbe_ = 0.0;
  double requestedMultipoleProbe_ = 0.0;
  bool requestedProbesSet_ = false;

  bool mapsLoaded_ = false;
  MdmFieldMap multipoleMap_;
  MdmFieldMap dipoleEntranceMap_;
  MdmFieldMap dipoleSectorMap_;
  MdmFieldMap dipoleExitMap_;

  double firstWireX_ = 0.0;
  double firstWireY_ = 0.0;
  double firstWireAngXDeg_ = 0.0;
  double firstWireAngYDeg_ = 0.0;
  double timeOfFlightSeconds_ = 0.0;
};
