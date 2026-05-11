#pragma once

#include "MdmIon.h"

class MdmTrace {
 public:
  MdmTrace();

  void SetBeamEnergy(double);
  double GetBeamEnergy() const;

  void SetMdmAngle(double);
  double GetMdmAngle() const;
  void SetMdmBRho(double);
  void SetMdmDipoleField(double);
  void SetMdmField(double, double);
  void SetMdmProbe(double, double);
  double GetMdmDipoleField() const;

  void SetScatteredAngle(double);
  void SetScatteredAngle(double, double);
  double GetScatteredAngle() const;

  void SetScatteredEnergy(double);
  double GetScatteredEnergy() const;

  void SetQValue(double);
  double GetQValue() const;
  void SetResidualEnergy(double);
  double GetResidualEnergy() const;

  void SetTargetMass(double);
  double GetTargetMass() const;
  void SetProjectileMass(double);
  double GetProjectileMass() const;

  void SetScatteredIon(const MdmIon&);

  double GetEnergyAfterKinematics() const;
  void SendRayWithKinematics();
  void SendRay();

  void GetPositionAngleFirstWire(double&, double&) const;
  void GetPositionAngleFirstWire(double&, double&, double&, double&) const;
  void GetOxfordWirePositions(double&, double&, double&, double&);

  double GetFirstWireX() const;
  double GetFirstWireY() const;
  double GetFirstWireXAngle() const;
  double GetFirstWireYAngle() const;

 private:
  static void SetMultipoleProbe(double multipoleProbe);

  static double jeffParams_[6];
  static double oxfordWireSpacing_[3];
  double beamEnergy_;
  double scatteredEnergy_;
  double scatteredAngles_[2];
};
