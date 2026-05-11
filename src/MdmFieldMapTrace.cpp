#include "MdmFieldMapTrace.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDegreesPerRadian = 180.0 / kPi;
constexpr double kRadiansPerDegree = kPi / 180.0;
constexpr double kMilliradiansPerDegree = 17.453;
constexpr double kSpeedOfLightCmPerSecond = 3.0e10;
constexpr double kMultipoleStepCm = 0.1;
constexpr double kDipoleStepCm = 0.1;
constexpr double kPlaneToleranceCm = 1.0e-8;
constexpr double kProbeTolerance = 2.0e-6;
constexpr double kStoppedPosition = 1.0e10;
constexpr std::size_t kMaxIntegrationSteps = 200000;
constexpr std::size_t kPlaneRefinementSteps = 60;

struct State {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
};

struct CollimatorGeometry {
  double type = 0.0;
  double xCenter = 0.0;
  double yCenter = 0.0;
  double xMax = 0.0;
  double yMax = 0.0;
};

struct MultipoleGeometry {
  double a = 0.0;
  double b = 0.0;
  double length = 0.0;
  double z11 = 0.0;
};

struct DipoleGeometry {
  double a = 0.0;
  double b = 0.0;
  double radius = 0.0;
  double phiDeg = 0.0;
  double alphaDeg = 0.0;
  double betaDeg = 0.0;
  double z11 = 0.0;
  double xcr1 = 0.0;
  double xcr2 = 0.0;
};

struct BeamlineGeometry {
  double entranceDrift = 0.0;
  CollimatorGeometry entranceCollimator;
  double intermediateDrift = 0.0;
  MultipoleGeometry multipole;
  DipoleGeometry dipole;
  double secondMultipoleDrift = 0.0;
  CollimatorGeometry exitCollimator;
  double finalDrift = 0.0;
};

struct PlaneCrossingResult {
  State state;
  double elapsedSeconds = 0.0;
};

struct IntegrationResult {
  State state;
  double elapsedSeconds = 0.0;
};

const BeamlineGeometry kMdmGeometry{
    63.5,
    {0.0, 0.0, 0.0, 2.27965, 2.2},
    18.075,
    {1.925, 113.2, 26.0, 20.0},
    {26.0, 32.55, 160.0, 100.0, 0.0, 0.0, 46.0, 0.0, 0.0},
    36.7,
    {0.0, -0.08, 0.0, 29.61, 20.21},
    96.13};

double DegreesToRadians(double degrees) {
  return degrees * kRadiansPerDegree;
}

bool NearlyEqual(double lhs, double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= kProbeTolerance * scale;
}

State AddScaled(const State& lhs, const State& rhs, double scale) {
  return {lhs.x + scale * rhs.x,
          lhs.y + scale * rhs.y,
          lhs.z + scale * rhs.z,
          lhs.vx + scale * rhs.vx,
          lhs.vy + scale * rhs.vy,
          lhs.vz + scale * rhs.vz};
}

State CombineForRK4(const State& k1,
                    const State& k2,
                    const State& k3,
                    const State& k4) {
  return {(k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) / 6.0,
          (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y) / 6.0,
          (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z) / 6.0,
          (k1.vx + 2.0 * k2.vx + 2.0 * k3.vx + k4.vx) / 6.0,
          (k1.vy + 2.0 * k2.vy + 2.0 * k3.vy + k4.vy) / 6.0,
          (k1.vz + 2.0 * k2.vz + 2.0 * k3.vz + k4.vz) / 6.0};
}

using FieldFunction = std::function<Vec3(const State&)>;
using PlaneFunction = std::function<double(const State&)>;

State Derivative(const State& state,
                 const Vec3& fieldTesla,
                 double kFactor) {
  return {state.vx,
          state.vy,
          state.vz,
          kFactor * (state.vy * fieldTesla.z - state.vz * fieldTesla.y),
          kFactor * (state.vz * fieldTesla.x - state.vx * fieldTesla.z),
          kFactor * (state.vx * fieldTesla.y - state.vy * fieldTesla.x)};
}

State RK4Step(const State& state,
              double dt,
              const FieldFunction& fieldFunction,
              double kFactor) {
  const State k1 = Derivative(state, fieldFunction(state), kFactor);
  const State k2 =
      Derivative(AddScaled(state, k1, 0.5 * dt),
                 fieldFunction(AddScaled(state, k1, 0.5 * dt)), kFactor);
  const State k3 =
      Derivative(AddScaled(state, k2, 0.5 * dt),
                 fieldFunction(AddScaled(state, k2, 0.5 * dt)), kFactor);
  const State k4 =
      Derivative(AddScaled(state, k3, dt),
                 fieldFunction(AddScaled(state, k3, dt)), kFactor);
  return AddScaled(state, CombineForRK4(k1, k2, k3, k4), dt);
}

bool PlaneCrossed(double left, double right) {
  return (left <= 0.0 && right >= 0.0) || (left >= 0.0 && right <= 0.0);
}

PlaneCrossingResult RefinePlaneCrossing(const State& leftInitial,
                                        const State& rightInitial,
                                        double stepDt,
                                        const FieldFunction& fieldFunction,
                                        const PlaneFunction& planeFunction,
                                        double kFactor) {
  State left = leftInitial;
  State right = rightInitial;
  double leftValue = planeFunction(left);
  double rightValue = planeFunction(right);
  double leftDt = 0.0;
  double rightDt = stepDt;
  double intervalDt = stepDt;

  for (std::size_t iteration = 0; iteration < kPlaneRefinementSteps;
       ++iteration) {
    const double midDt = 0.5 * intervalDt;
    const State mid = RK4Step(left, midDt, fieldFunction, kFactor);
    const double midValue = planeFunction(mid);
    if (std::abs(midValue) <= kPlaneToleranceCm) {
      return {mid, leftDt + midDt};
    }
    if (PlaneCrossed(leftValue, midValue)) {
      right = mid;
      rightValue = midValue;
      rightDt = leftDt + midDt;
    } else {
      left = mid;
      leftValue = midValue;
      leftDt += midDt;
    }
    intervalDt = midDt;
  }

  return (std::abs(leftValue) < std::abs(rightValue))
             ? PlaneCrossingResult{left, leftDt}
             : PlaneCrossingResult{right, rightDt};
}

IntegrationResult IntegrateToPlane(const State& initialState,
                                   const FieldFunction& fieldFunction,
                                   const PlaneFunction& planeFunction,
                                   double stepCm,
                                   double speedCmPerSecond,
                                   double kFactor,
                                   const std::string& label) {
  State current = initialState;
  double currentPlaneValue = planeFunction(current);
  if (std::abs(currentPlaneValue) <= kPlaneToleranceCm) {
    return {current, 0.0};
  }

  const double dt = stepCm / speedCmPerSecond;
  double elapsed = 0.0;
  for (std::size_t step = 0; step < kMaxIntegrationSteps; ++step) {
    const State next = RK4Step(current, dt, fieldFunction, kFactor);
    const double nextPlaneValue = planeFunction(next);
    if (std::abs(nextPlaneValue) <= kPlaneToleranceCm) {
      return {next, elapsed + dt};
    }
    if (PlaneCrossed(currentPlaneValue, nextPlaneValue)) {
      const PlaneCrossingResult crossing = RefinePlaneCrossing(
          current, next, dt, fieldFunction, planeFunction, kFactor);
      return {crossing.state, elapsed + crossing.elapsedSeconds};
    }
    current = next;
    currentPlaneValue = nextPlaneValue;
    elapsed += dt;
  }

  throw std::runtime_error("Failed to reach " + label + " exit plane");
}

double OutputAngleXDegrees(const State& state) {
  return std::atan2(state.vx, state.vz) * kDegreesPerRadian;
}

double OutputAngleYDegrees(const State& state, double speedCmPerSecond) {
  return std::asin(state.vy / speedCmPerSecond) * kDegreesPerRadian;
}

bool InMapBounds(const MdmFieldMap& map, double x, double y, double z) {
  const auto within = [](double coord, double origin, double spacing,
                         int count) {
    const double maxCoord = origin + spacing * static_cast<double>(count - 1);
    const double tolerance = 1.0e-9;
    return coord >= origin - tolerance && coord <= maxCoord + tolerance;
  };
  return within(x, map.h.origin_cm.x, map.h.step_cm.x, map.h.nx) &&
         within(y, map.h.origin_cm.y, map.h.step_cm.y, map.h.ny) &&
         within(z, map.h.origin_cm.z, map.h.step_cm.z, map.h.nz);
}

}  // namespace

MdmFieldMapTrace::MdmFieldMapTrace() = default;

void MdmFieldMapTrace::LoadFieldMaps(const std::string& multipolePath,
                                     const std::string& dipoleEntrancePath,
                                     const std::string& dipoleSectorPath,
                                     const std::string& dipoleExitPath) {
  multipoleMap_.Load(multipolePath);
  dipoleEntranceMap_.Load(dipoleEntrancePath);
  dipoleSectorMap_.Load(dipoleSectorPath);
  dipoleExitMap_.Load(dipoleExitPath);

  mapsLoaded_ = true;
  ValidateLoadedMaps();
}

void MdmFieldMapTrace::SetMdmAngle(double angleDeg) { mdmAngleDeg_ = angleDeg; }

double MdmFieldMapTrace::GetMdmAngle() const { return mdmAngleDeg_; }

void MdmFieldMapTrace::SetMdmProbe(double dipoleProbe, double multipoleProbe) {
  requestedDipoleProbe_ = dipoleProbe;
  requestedMultipoleProbe_ = multipoleProbe;
  requestedProbesSet_ = true;
  if (mapsLoaded_) {
    ValidateLoadedMaps();
  }
}

void MdmFieldMapTrace::SetMdmDipoleField(double dipoleField) {
  requestedDipoleProbe_ = dipoleField / 1.034;
  requestedMultipoleProbe_ = requestedDipoleProbe_ * 0.71;
  requestedProbesSet_ = true;
  if (mapsLoaded_) {
    ValidateLoadedMaps();
  }
}

void MdmFieldMapTrace::SetScatteredIon(const MdmIon& ion) {
  ionMassMeV_ = ion.ionMassMeV;
  ionChargeState_ = ion.chargeState;
}

void MdmFieldMapTrace::SetScatteredEnergy(double energyMeV) {
  scatteredEnergyMeV_ = energyMeV;
}

double MdmFieldMapTrace::GetScatteredEnergy() const {
  return scatteredEnergyMeV_;
}

void MdmFieldMapTrace::SetInitialPosition(double xCm, double yCm) {
  initialXcm_ = xCm;
  initialYcm_ = yCm;
}

void MdmFieldMapTrace::SetScatteredAngle(double xAngleDeg) {
  SetScatteredAngle(xAngleDeg, 0.0);
}

void MdmFieldMapTrace::SetScatteredAngle(double xAngleDeg, double yAngleDeg) {
  scatteredAnglesDeg_[0] = xAngleDeg;
  scatteredAnglesDeg_[1] = yAngleDeg;
}

double MdmFieldMapTrace::GetScatteredAngle() const {
  return scatteredAnglesDeg_[0];
}

void MdmFieldMapTrace::ValidateLoadedMaps() const {
  if (!mapsLoaded_) {
    return;
  }

  std::vector<const MdmFieldMap*> mapsToCheck;
  mapsToCheck.push_back(&multipoleMap_);
  mapsToCheck.push_back(&dipoleEntranceMap_);
  mapsToCheck.push_back(&dipoleSectorMap_);
  mapsToCheck.push_back(&dipoleExitMap_);

  const double referenceDipoleProbe = mapsToCheck.front()->h.mdm_dipole_probe;
  const double referenceMultipoleProbe =
      mapsToCheck.front()->h.mdm_multipole_probe;

  for (const MdmFieldMap* map : mapsToCheck) {
    const double dipoleProbe = map->h.mdm_dipole_probe;
    const double multipoleProbe = map->h.mdm_multipole_probe;
    if (!NearlyEqual(referenceDipoleProbe, dipoleProbe) ||
        !NearlyEqual(referenceMultipoleProbe, multipoleProbe)) {
      throw std::runtime_error("Loaded field maps were generated with different "
                               "magnet settings");
    }
  }

  if (!requestedProbesSet_) {
    return;
  }
  if (!NearlyEqual(requestedDipoleProbe_, referenceDipoleProbe) ||
      !NearlyEqual(requestedMultipoleProbe_, referenceMultipoleProbe)) {
    throw std::runtime_error(
        "Requested probe settings do not match field-map metadata");
  }
}

void MdmFieldMapTrace::SendRay() {
  timeOfFlightSeconds_ = 0.0;
  if (!mapsLoaded_) {
    throw std::runtime_error("LoadFieldMaps must be called before SendRay");
  }
  if (!requestedProbesSet_) {
    throw std::runtime_error(
        "Magnet settings must be configured before SendRay");
  }
  ValidateLoadedMaps();
  if (ionMassMeV_ <= 0.0) {
    throw std::runtime_error("Ion mass must be positive");
  }
  if (ionChargeState_ == 0.0) {
    throw std::runtime_error("Ion charge state must be non-zero");
  }
  if (scatteredEnergyMeV_ <= 0.0) {
    throw std::runtime_error("Scattered energy must be positive");
  }

  const double massMeV = ionMassMeV_;
  const double totalEnergyMeV = massMeV + scatteredEnergyMeV_;
  const double speedCmPerSecond =
      std::sqrt((2.0 * massMeV + scatteredEnergyMeV_) * scatteredEnergyMeV_) /
      totalEnergyMeV * kSpeedOfLightCmPerSecond;
  const double kFactor = (ionChargeState_ / totalEnergyMeV) * 9.0e10;

  const double thetaXMrad =
      kMilliradiansPerDegree * (scatteredAnglesDeg_[0] - mdmAngleDeg_);
  const double thetaYMrad = kMilliradiansPerDegree * scatteredAnglesDeg_[1];
  const double thetaXRad = thetaXMrad / 1000.0;
  const double thetaYRad = thetaYMrad / 1000.0;

  State state;
  state.x = initialXcm_;
  state.y = initialYcm_;
  state.z = 0.0;
  state.vx = speedCmPerSecond * std::sin(thetaXRad) * std::cos(thetaYRad);
  state.vy = speedCmPerSecond * std::sin(thetaYRad);
  state.vz = speedCmPerSecond * std::cos(thetaXRad) * std::cos(thetaYRad);

  const auto stopRay = [&]() {
    firstWireX_ = kStoppedPosition;
    firstWireY_ = kStoppedPosition;
    firstWireAngXDeg_ = 0.0;
    firstWireAngYDeg_ = 0.0;
    timeOfFlightSeconds_ = std::numeric_limits<double>::quiet_NaN();
  };

  const auto propagateDrift = [&](double distanceCm) {
    const double dt = distanceCm / std::abs(state.vz);
    timeOfFlightSeconds_ += dt;
    state.x += dt * state.vx;
    state.y += dt * state.vy;
    state.z = 0.0;
  };

  const auto applyCollimator = [&](const CollimatorGeometry& collimator) {
    if (collimator.type == 0.0) {
      return std::abs(state.x - collimator.xCenter) <= collimator.xMax &&
             std::abs(state.y - collimator.yCenter) <= collimator.yMax;
    }

    const double xc = (state.x - collimator.xCenter) / collimator.xMax;
    const double yc = (state.y - collimator.yCenter) / collimator.yMax;
    return xc * xc + yc * yc <= 1.0;
  };

  propagateDrift(kMdmGeometry.entranceDrift);
  if (!applyCollimator(kMdmGeometry.entranceCollimator)) {
    stopRay();
    return;
  }

  propagateDrift(kMdmGeometry.intermediateDrift);

  {
    const MultipoleGeometry& multipole = kMdmGeometry.multipole;
    const double a = multipole.a;
    const double b = multipole.b;
    const double l = multipole.length;
    const double z11 = multipole.z11;

    State localState;
    localState.x = state.x;
    localState.y = state.y;
    localState.z = state.z - (a + 0.5 * l);
    localState.vx = state.vx;
    localState.vy = state.vy;
    localState.vz = state.vz;

    // RAYTRACE element reference planes can sit inside the fringe model. Match
    // POLES by backing up to Z11 before starting magnetic integration.
    if (std::abs(localState.vz) <= 1.0e-30) {
      throw std::runtime_error("Multipole local vz is zero");
    }
    const double multipoleStartZ = -0.5 * l - z11;
    const double startDt = (multipoleStartZ - localState.z) / localState.vz;
    timeOfFlightSeconds_ += startDt;
    localState.x += startDt * localState.vx;
    localState.y += startDt * localState.vy;
    localState.z = multipoleStartZ;

    const FieldFunction fieldFunction = [&](const State& sample) {
      return multipoleMap_.FieldTesla(sample.x, sample.y, sample.z);
    };
    const PlaneFunction planeFunction = [&](const State& sample) {
      return sample.z - (b + 0.5 * l);
    };

    const IntegrationResult integrated = IntegrateToPlane(
        localState, fieldFunction, planeFunction, kMultipoleStepCm,
        speedCmPerSecond, kFactor, "multipole");
    timeOfFlightSeconds_ += integrated.elapsedSeconds;
    localState = integrated.state;

    state.x = localState.x;
    state.y = localState.y;
    state.z = 0.0;
    state.vx = localState.vx;
    state.vy = localState.vy;
    state.vz = localState.vz;
  }

  {
    const DipoleGeometry& dipole = kMdmGeometry.dipole;
    const double a = dipole.a;
    const double b = dipole.b;
    const double radius = dipole.radius;
    const double phi = dipole.phiDeg;
    const double alphaDeg = dipole.alphaDeg;
    const double betaDeg = dipole.betaDeg;
    const double alpha = DegreesToRadians(alphaDeg);
    const double beta = DegreesToRadians(betaDeg);
    const double z11 = dipole.z11;
    const double xcr1 = dipole.xcr1;
    const double xcr2 = dipole.xcr2;

    const double cosAlpha = std::cos(alpha);
    const double sinAlpha = std::sin(alpha);
    const double cosBeta = std::cos(beta);
    const double sinBeta = std::sin(beta);
    const double rotation = DegreesToRadians(phi - alphaDeg - betaDeg);
    const double cosRot = std::cos(rotation);
    const double sinRot = std::sin(rotation);
    const double cosPb = std::cos(DegreesToRadians(0.5 * phi - betaDeg));
    const double sinPb = std::sin(DegreesToRadians(0.5 * phi - betaDeg));
    const double sinHalfPhi = std::sin(DegreesToRadians(0.5 * phi));
    const double tx = 2.0 * radius * sinHalfPhi * sinPb;
    const double tz = 2.0 * radius * sinHalfPhi * cosPb;

    const auto localToOutputD = [&](const State& localState) {
      const double xB = -localState.x * cosAlpha - localState.z * sinAlpha;
      const double zB = localState.x * sinAlpha - localState.z * cosAlpha;
      const double vxB = -localState.vx * cosAlpha - localState.vz * sinAlpha;
      const double vzB = localState.vx * sinAlpha - localState.vz * cosAlpha;

      const double xC = -zB * sinRot - xB * cosRot - tx;
      const double zC = -zB * cosRot + xB * sinRot - tz;
      const double vxC = -vzB * sinRot - vxB * cosRot;
      const double vzC = -vzB * cosRot + vxB * sinRot;

      State output;
      output.x = zC * sinBeta + xC * cosBeta - xcr2;
      output.y = localState.y;
      output.z = zC * cosBeta - xC * sinBeta - b;
      output.vx = vzC * sinBeta + vxC * cosBeta;
      output.vy = localState.vy;
      output.vz = vzC * cosBeta - vxC * sinBeta;
      return output;
    };

    State localState;
    const double xB = (a - state.z) * sinAlpha - (state.x + xcr1) * cosAlpha;
    const double zB = (a - state.z) * cosAlpha + (state.x + xcr1) * sinAlpha;
    const double vxB = -state.vz * sinAlpha - state.vx * cosAlpha;
    const double vzB = -state.vz * cosAlpha + state.vx * sinAlpha;

    localState.x = -xB * cosAlpha + zB * sinAlpha;
    localState.y = state.y;
    localState.z = -xB * sinAlpha - zB * cosAlpha;
    localState.vx = -vxB * cosAlpha + vzB * sinAlpha;
    localState.vy = state.vy;
    localState.vz = -vxB * sinAlpha - vzB * cosAlpha;

    // Match DIPOLE: transform to the entrance B frame, then back up to Z11
    // before following the split maps through fringe, sector, and exit.
    const double entranceZB =
        localState.x * sinAlpha - localState.z * cosAlpha;
    const double entranceVzB =
        localState.vx * sinAlpha - localState.vz * cosAlpha;
    if (std::abs(entranceVzB) <= 1.0e-30) {
      throw std::runtime_error("Dipole entrance-frame vz is zero");
    }
    const double startDt = (z11 - entranceZB) / entranceVzB;
    timeOfFlightSeconds_ += startDt;
    localState.x += startDt * localState.vx;
    localState.y += startDt * localState.vy;
    localState.z += startDt * localState.vz;

    const FieldFunction fieldFunction = [&](const State& sample) {
      const double xB = -sample.x * cosAlpha - sample.z * sinAlpha;
      const double zB = sample.x * sinAlpha - sample.z * cosAlpha;
      if (InMapBounds(dipoleEntranceMap_, xB, sample.y, zB)) {
        return dipoleEntranceMap_.FieldTesla(xB, sample.y, zB);
      }

      const double xC = -zB * sinRot - xB * cosRot - tx;
      const double zC = -zB * cosRot + xB * sinRot - tz;
      if (InMapBounds(dipoleExitMap_, xC, sample.y, zC)) {
        return dipoleExitMap_.FieldTesla(xC, sample.y, zC);
      }

      const double radialDistance =
          std::sqrt((sample.x + radius) * (sample.x + radius) +
                    sample.z * sample.z);
      const double dr = radialDistance - radius;
      const double theta = std::atan2(sample.z, sample.x + radius);
      const double s = radius * theta;
      if (InMapBounds(dipoleSectorMap_, dr, sample.y, s)) {
        return dipoleSectorMap_.FieldTesla(dr, sample.y, s);
      }
      return Vec3{};
    };
    const PlaneFunction planeFunction = [&](const State& sample) {
      return localToOutputD(sample).z;
    };

    const IntegrationResult integrated =
        IntegrateToPlane(localState, fieldFunction, planeFunction,
                         kDipoleStepCm, speedCmPerSecond, kFactor, "dipole");
    timeOfFlightSeconds_ += integrated.elapsedSeconds;
    localState = integrated.state;

    state = localToOutputD(localState);
    state.z = 0.0;
  }

  propagateDrift(kMdmGeometry.secondMultipoleDrift);

  if (!applyCollimator(kMdmGeometry.exitCollimator)) {
    stopRay();
    return;
  }

  propagateDrift(kMdmGeometry.finalDrift);

  firstWireX_ = state.x;
  firstWireY_ = state.y;
  firstWireAngXDeg_ = OutputAngleXDegrees(state);
  firstWireAngYDeg_ = OutputAngleYDegrees(state, speedCmPerSecond);
}

void MdmFieldMapTrace::GetPositionAngleFirstWire(double& posX,
                                                 double& posY,
                                                 double& angX,
                                                 double& angY) const {
  posX = firstWireX_;
  posY = firstWireY_;
  angX = firstWireAngXDeg_;
  angY = firstWireAngYDeg_;
}

double MdmFieldMapTrace::GetTimeOfFlightSeconds() const {
  return timeOfFlightSeconds_;
}
