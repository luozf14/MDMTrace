#ifndef MDMFIELDMAPTRACE_H
#define MDMFIELDMAPTRACE_H

#include "MDMFieldMap.h"

#include <string>

class MDMFieldMapTrace {
 public:
  MDMFieldMapTrace();

  void LoadFieldMaps(const std::string& multipolePath = "Multipole.bin",
                     const std::string& dipolePath = "Dipole.bin");
  void LoadFieldMaps(const std::string& multipolePath,
                     const std::string& dipoleEntrancePath,
                     const std::string& dipoleSectorPath,
                     const std::string& dipoleExitPath);

  void SetMDMAngle(double angleDeg);
  double GetMDMAngle() const;

  void SetMDMProbe(double dipoleProbe, double multipoleProbe);
  void SetMDMDipoleField(double dipoleField);

  void SetScatteredMass(double massAmu);
  double GetScatteredMass() const;

  void SetScatteredCharge(double charge);
  double GetScatteredCharge() const;

  void SetScatteredEnergy(double energyMeV);
  double GetScatteredEnergy() const;

  void SetScatteredAngle(double xAngleDeg);
  void SetScatteredAngle(double xAngleDeg, double yAngleDeg);
  double GetScatteredAngle() const;

  void SendRay();

  void GetPositionAngleFirstWire(double& posX,
                                 double& posY,
                                 double& angX,
                                 double& angY) const;

 private:
  void ValidateLoadedMaps() const;

  double mdmAngleDeg_ = 0.0;
  double scatteredMassAmu_ = 0.0;
  double scatteredCharge_ = 0.0;
  double scatteredEnergyMeV_ = 0.0;
  double scatteredAnglesDeg_[2]{0.0, 0.0};

  double requestedDipoleProbe_ = 0.0;
  double requestedMultipoleProbe_ = 0.0;
  bool requestedProbesSet_ = false;

  bool mapsLoaded_ = false;
  bool usingSplitDipoleMaps_ = false;
  MDMFieldMap multipoleMap_;
  MDMFieldMap dipoleMap_;
  MDMFieldMap dipoleEntranceMap_;
  MDMFieldMap dipoleSectorMap_;
  MDMFieldMap dipoleExitMap_;

  double firstWireX_ = 0.0;
  double firstWireY_ = 0.0;
  double firstWireAngXDeg_ = 0.0;
  double firstWireAngYDeg_ = 0.0;
};

#endif
