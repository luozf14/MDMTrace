#include "MdmTrace.h"

#include <cmath>

extern "C" {
  void raytrace_(int*);

  extern struct {
    double DATA[200][75];
    double ITITLE[200];
  } blck0_;

  extern struct {
    double XI[1000];
    double YI[1000];
    double ZI[1000];
    double VXI[1000];
    double VYI[1000];
    double VZI[1000];
    double DELP[1000];
  } blck1_;

  extern struct {
    double XO[1000];
    double YO[1000];
    double ZO[1000];
    double VXO[1000];
    double VYO[1000];
    double VZO[1000];
    double RTL[1000];
    double RLL[1000];
  } blck2_;

  extern struct {
    double ENERGY;
    double VEL;
    double PMASS;
    double Q0;
  } blck4_;

  extern struct {
    double THTSPEC;
    double TRGT1;
    double AM[4];
    double QVALUE;
    double EEXC;
    double THETACAL[10];
    double EKINE;
  } kineblck_;
}

double MdmTrace::jeffParams_[6] = {
    -0.51927, 0.038638, 0.028404, -0.022797, -0.019275, 0.755583};
double MdmTrace::oxfordWireSpacing_[3] = {15.1, 16.3, 16.3};

MdmTrace::MdmTrace() {
  int flag = 0;
  raytrace_(&flag);
  kineblck_.TRGT1 = 0.;
  beamEnergy_ = 0;
  scatteredEnergy_ = 0;
}

void MdmTrace::SetBeamEnergy(double energy) {
  beamEnergy_ = energy;
}

double MdmTrace::GetBeamEnergy() const {
  return beamEnergy_;
}

void MdmTrace::SetMdmAngle(double angle) {
  kineblck_.THTSPEC = angle;
}

double MdmTrace::GetMdmAngle() const {
  return kineblck_.THTSPEC;
}

void MdmTrace::SetMdmBRho(double bRho) {
  double field = bRho / 160. * 1000;
  SetMdmDipoleField(field);
}

void MdmTrace::SetMdmDipoleField(double field) {
  blck0_.DATA[4][14] = field * 1.e-4;
  SetMultipoleProbe(field / 1.034 * 0.71);
}

void MdmTrace::SetMdmField(double dipoleField, double multipoleField) {
  blck0_.DATA[4][14] = dipoleField * 1.e-4;
  SetMultipoleProbe(multipoleField / 1.034);
}

void MdmTrace::SetMdmProbe(double dipoleProbe, double multipoleProbe) {
  blck0_.DATA[4][14] = dipoleProbe * 1.034 * 1.e-4;
  SetMultipoleProbe(multipoleProbe);
}

void MdmTrace::SetMultipoleProbe(double multipoleProbe) {
  double BQR = -multipoleProbe * 1e-4 * jeffParams_[5];
  double BHR = BQR * jeffParams_[1] / jeffParams_[0];
  double BOR = BQR * jeffParams_[2] / jeffParams_[0];
  double BDR = BQR * jeffParams_[3] / jeffParams_[0];
  double BDDR = BQR * jeffParams_[4] / jeffParams_[0];
  blck0_.DATA[3][13] = BQR;
  blck0_.DATA[3][14] = BHR;
  blck0_.DATA[3][15] = BOR;
  blck0_.DATA[3][16] = BDR;
  blck0_.DATA[3][17] = BDDR;
}

double MdmTrace::GetMdmDipoleField() const {
  return blck0_.DATA[4][14] * 1.e4;
}

void MdmTrace::SetScatteredAngle(double angle) {
  kineblck_.THETACAL[0] = angle;
  scatteredAngles_[0] = angle;
  scatteredAngles_[1] = 0.;
}

void MdmTrace::SetScatteredAngle(double xAngle, double yAngle) {
  kineblck_.THETACAL[0] = xAngle;
  scatteredAngles_[0] = xAngle;
  scatteredAngles_[1] = yAngle;
}

double MdmTrace::GetScatteredAngle() const {
  return kineblck_.THETACAL[0];
}

void MdmTrace::SetQValue(double qValue) {
  kineblck_.QVALUE =qValue;
}

double MdmTrace::GetQValue() const {
  return kineblck_.QVALUE;
}

void MdmTrace::SetResidualEnergy(double energy) {
  kineblck_.EEXC = energy;
}

double MdmTrace::GetResidualEnergy() const {
  return kineblck_.EEXC;
}

void MdmTrace::SetScatteredEnergy(double energy) {
  scatteredEnergy_ = energy;
}

double MdmTrace::GetScatteredEnergy() const {
  return scatteredEnergy_;
}

void MdmTrace::SetTargetMass(double mass) {
    kineblck_.AM[1] = mass;
}

double MdmTrace::GetTargetMass() const {
  return kineblck_.AM[1];
}

void MdmTrace::SetProjectileMass(double mass) {
    kineblck_.AM[0] = mass;
}

double MdmTrace::GetProjectileMass() const {
  return kineblck_.AM[0];
}

void MdmTrace::SetScatteredMass(double mass) {
    blck4_.PMASS = mass;
}

double MdmTrace::GetScatteredMass() const {
  return blck4_.PMASS;
}

void MdmTrace::SetScatteredCharge(double charge) {
    blck4_.Q0 = charge;
}

double MdmTrace::GetScatteredCharge() const {
  return blck4_.Q0;
}

double MdmTrace::GetEnergyAfterKinematics() const {
  return kineblck_.EKINE*(1.+blck1_.DELP[0]/100.);
}

void MdmTrace::SendRayWithKinematics() {
    int flag = 1;
    blck4_.ENERGY = beamEnergy_;
    raytrace_(&flag);
}

void MdmTrace::SendRay() {
  int flag = 2;
  blck4_.ENERGY = scatteredEnergy_;
  blck1_.XI[0]=0.;
  blck1_.YI[0]=0.;
  blck1_.ZI[0]=0.;
  blck1_.VXI[0]=17.453*(scatteredAngles_[0]-kineblck_.THTSPEC);
  blck1_.VYI[0]=17.453*(scatteredAngles_[1]);
  blck1_.VZI[0]=0.;
  blck1_.DELP[0]=0.;
  raytrace_(&flag);
}

void MdmTrace::GetPositionAngleFirstWire(double& pos, double& ang) const {
  pos  = blck2_.XO[0];
  ang  = blck2_.VXO[0]/1000.*180./3.14159;
}

void MdmTrace::GetPositionAngleFirstWire(double& posX, double& posY, double& angX, double& angY) const {
  posX  = blck2_.XO[0];
  posY  = blck2_.YO[0];
  angX  = blck2_.VXO[0]/1000.*180./3.14159;
  angY  = blck2_.VYO[0]/1000.*180./3.14159;
}

void MdmTrace::GetOxfordWirePositions(double& x1,double& x2,double& x3,double& x4) {
  double oxfordWire1Pos = blck2_.XO[0];
  double oxfordWire1Ang = blck2_.VXO[0];
  double tanAngle = std::tan(1e-3*oxfordWire1Ang);

  x1 = oxfordWire1Pos;
  x2 = oxfordWire1Pos+tanAngle*oxfordWireSpacing_[0];
  x3 = oxfordWire1Pos+tanAngle*(oxfordWireSpacing_[0]+oxfordWireSpacing_[1]);
  x4 = oxfordWire1Pos+tanAngle*(oxfordWireSpacing_[0]+oxfordWireSpacing_[1]+oxfordWireSpacing_[2]);
}

double MdmTrace::GetFirstWireX() const{return blck2_.XO[0];}
double MdmTrace::GetFirstWireY() const{return blck2_.YO[0];}
double MdmTrace::GetFirstWireXAngle() const{return blck2_.VXO[0]/1000.*180./3.14159;}
double MdmTrace::GetFirstWireYAngle() const{return blck2_.VYO[0]/1000.*180./3.14159;}
