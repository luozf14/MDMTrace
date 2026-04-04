#include "MDMFieldMapTrace.h"

#include "MDMFieldMapInterop.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr std::size_t kEntranceDriftElement = 1;
constexpr std::size_t kEntranceCollimatorElement = 2;
constexpr std::size_t kIntermediateDriftElement = 3;
constexpr std::size_t kMultipoleElement = 4;
constexpr std::size_t kDipoleElement = 5;
constexpr std::size_t kSecondMultipoleElement = 6;
constexpr std::size_t kExitCollimatorElement = 7;
constexpr std::size_t kFinalDriftElement = 8;

constexpr double kPi = 3.14159265358979323846;
constexpr double kDegreesPerRadian = 180.0 / kPi;
constexpr double kRadiansPerDegree = kPi / 180.0;
constexpr double kMilliradiansPerDegree = 17.453;
constexpr double kSpeedOfLightCmPerSecond = 3.0e10;
constexpr double kMultipoleStepCm = 0.1;
constexpr double kDipoleStepCm = 0.1;
constexpr double kPlaneToleranceCm = 1.0e-8;
constexpr double kProbeTolerance = 1.0e-6;
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

double DegreesToRadians(double degrees) {
  return degrees * kRadiansPerDegree;
}

double DeckValue(std::size_t elementOneBased, std::size_t fieldOneBased) {
  return blck0_.DATA[elementOneBased - 1][fieldOneBased - 1];
}

double MetadataDouble(const MDMFieldMap& map, const std::string& key) {
  const auto it = map.GetMetadata().fields.find(key);
  if (it == map.GetMetadata().fields.end()) {
    throw std::runtime_error("Missing field-map metadata entry: " + key);
  }
  return std::stod(it->second);
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

State Scale(const State& state, double scale) {
  return {scale * state.x, scale * state.y, scale * state.z,
          scale * state.vx, scale * state.vy, scale * state.vz};
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

using FieldFunction = std::function<std::array<double, 3>(const State&)>;
using PlaneFunction = std::function<double(const State&)>;

State Derivative(const State& state,
                 const std::array<double, 3>& fieldTesla,
                 double kFactor) {
  return {state.vx,
          state.vy,
          state.vz,
          kFactor * (state.vy * fieldTesla[2] - state.vz * fieldTesla[1]),
          kFactor * (state.vz * fieldTesla[0] - state.vx * fieldTesla[2]),
          kFactor * (state.vx * fieldTesla[1] - state.vy * fieldTesla[0])};
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

State RefinePlaneCrossing(const State& leftInitial,
                          const State& rightInitial,
                          double stepDt,
                          const FieldFunction& fieldFunction,
                          const PlaneFunction& planeFunction,
                          double kFactor) {
  State left = leftInitial;
  State right = rightInitial;
  double leftValue = planeFunction(left);
  double rightValue = planeFunction(right);
  double intervalDt = stepDt;

  for (std::size_t iteration = 0; iteration < kPlaneRefinementSteps;
       ++iteration) {
    const double midDt = 0.5 * intervalDt;
    const State mid = RK4Step(left, midDt, fieldFunction, kFactor);
    const double midValue = planeFunction(mid);
    if (std::abs(midValue) <= kPlaneToleranceCm) {
      return mid;
    }
    if (PlaneCrossed(leftValue, midValue)) {
      right = mid;
      rightValue = midValue;
    } else {
      left = mid;
      leftValue = midValue;
    }
    intervalDt = midDt;
  }

  return (std::abs(leftValue) < std::abs(rightValue)) ? left : right;
}

State IntegrateToPlane(const State& initialState,
                       const FieldFunction& fieldFunction,
                       const PlaneFunction& planeFunction,
                       double stepCm,
                       double speedCmPerSecond,
                       double kFactor,
                       const std::string& label) {
  State current = initialState;
  double currentPlaneValue = planeFunction(current);
  if (std::abs(currentPlaneValue) <= kPlaneToleranceCm) {
    return current;
  }

  const double dt = stepCm / speedCmPerSecond;
  for (std::size_t step = 0; step < kMaxIntegrationSteps; ++step) {
    const State next = RK4Step(current, dt, fieldFunction, kFactor);
    const double nextPlaneValue = planeFunction(next);
    if (std::abs(nextPlaneValue) <= kPlaneToleranceCm) {
      return next;
    }
    if (PlaneCrossed(currentPlaneValue, nextPlaneValue)) {
      return RefinePlaneCrossing(current, next, dt, fieldFunction,
                                 planeFunction, kFactor);
    }
    current = next;
    currentPlaneValue = nextPlaneValue;
  }

  throw std::runtime_error("Failed to reach " + label + " exit plane");
}

double OutputAngleXDegrees(const State& state) {
  return std::atan2(state.vx, state.vz) * kDegreesPerRadian;
}

double OutputAngleYDegrees(const State& state, double speedCmPerSecond) {
  return std::asin(state.vy / speedCmPerSecond) * kDegreesPerRadian;
}

}  // namespace

MDMFieldMapTrace::MDMFieldMapTrace() { mdmfm_init(); }

void MDMFieldMapTrace::LoadFieldMaps(const std::string& multipolePath,
                                     const std::string& dipolePath) {
  multipoleMap_ = MDMFieldMap::Load(multipolePath);
  dipoleMap_ = MDMFieldMap::Load(dipolePath);
  if (multipoleMap_.GetMetadata().magnetName != "Multipole") {
    throw std::runtime_error("Unexpected magnet name in multipole map");
  }
  if (dipoleMap_.GetMetadata().magnetName != "Dipole") {
    throw std::runtime_error("Unexpected magnet name in dipole map");
  }
  mapsLoaded_ = true;
  ValidateLoadedMaps();
}

void MDMFieldMapTrace::SetMDMAngle(double angleDeg) { mdmAngleDeg_ = angleDeg; }

double MDMFieldMapTrace::GetMDMAngle() const { return mdmAngleDeg_; }

void MDMFieldMapTrace::SetMDMProbe(double dipoleProbe, double multipoleProbe) {
  requestedDipoleProbe_ = dipoleProbe;
  requestedMultipoleProbe_ = multipoleProbe;
  requestedProbesSet_ = true;
  if (mapsLoaded_) {
    ValidateLoadedMaps();
  }
}

void MDMFieldMapTrace::SetMDMDipoleField(double dipoleField) {
  requestedDipoleProbe_ = dipoleField / 1.034;
  requestedMultipoleProbe_ = requestedDipoleProbe_ * 0.71;
  requestedProbesSet_ = true;
  if (mapsLoaded_) {
    ValidateLoadedMaps();
  }
}

void MDMFieldMapTrace::SetScatteredMass(double massAmu) {
  scatteredMassAmu_ = massAmu;
}

double MDMFieldMapTrace::GetScatteredMass() const { return scatteredMassAmu_; }

void MDMFieldMapTrace::SetScatteredCharge(double charge) {
  scatteredCharge_ = charge;
}

double MDMFieldMapTrace::GetScatteredCharge() const { return scatteredCharge_; }

void MDMFieldMapTrace::SetScatteredEnergy(double energyMeV) {
  scatteredEnergyMeV_ = energyMeV;
}

double MDMFieldMapTrace::GetScatteredEnergy() const {
  return scatteredEnergyMeV_;
}

void MDMFieldMapTrace::SetScatteredAngle(double xAngleDeg) {
  SetScatteredAngle(xAngleDeg, 0.0);
}

void MDMFieldMapTrace::SetScatteredAngle(double xAngleDeg, double yAngleDeg) {
  scatteredAnglesDeg_[0] = xAngleDeg;
  scatteredAnglesDeg_[1] = yAngleDeg;
}

double MDMFieldMapTrace::GetScatteredAngle() const {
  return scatteredAnglesDeg_[0];
}

void MDMFieldMapTrace::ValidateLoadedMaps() const {
  if (!mapsLoaded_ || !requestedProbesSet_) {
    return;
  }

  const double multipoleDipoleProbe =
      MetadataDouble(multipoleMap_, "mdm_dipole_probe");
  const double multipoleMultipoleProbe =
      MetadataDouble(multipoleMap_, "mdm_multipole_probe");
  const double dipoleDipoleProbe =
      MetadataDouble(dipoleMap_, "mdm_dipole_probe");
  const double dipoleMultipoleProbe =
      MetadataDouble(dipoleMap_, "mdm_multipole_probe");

  if (!NearlyEqual(multipoleDipoleProbe, dipoleDipoleProbe) ||
      !NearlyEqual(multipoleMultipoleProbe, dipoleMultipoleProbe)) {
    throw std::runtime_error("Multipole and dipole maps were generated with "
                             "different magnet settings");
  }

  if (!NearlyEqual(requestedDipoleProbe_, dipoleDipoleProbe) ||
      !NearlyEqual(requestedMultipoleProbe_, dipoleMultipoleProbe)) {
    throw std::runtime_error(
        "Requested probe settings do not match field-map metadata");
  }
}

void MDMFieldMapTrace::SendRay() {
  if (!mapsLoaded_) {
    throw std::runtime_error("LoadFieldMaps must be called before SendRay");
  }
  if (!requestedProbesSet_) {
    throw std::runtime_error(
        "Magnet settings must be configured before SendRay");
  }
  ValidateLoadedMaps();
  if (scatteredMassAmu_ <= 0.0) {
    throw std::runtime_error("Scattered mass must be positive");
  }
  if (scatteredCharge_ == 0.0) {
    throw std::runtime_error("Scattered charge must be non-zero");
  }
  if (scatteredEnergyMeV_ <= 0.0) {
    throw std::runtime_error("Scattered energy must be positive");
  }

  const double massMeV = scatteredMassAmu_ * 931.48;
  const double totalEnergyMeV = massMeV + scatteredEnergyMeV_;
  const double speedCmPerSecond =
      std::sqrt((2.0 * massMeV + scatteredEnergyMeV_) * scatteredEnergyMeV_) /
      totalEnergyMeV * kSpeedOfLightCmPerSecond;
  const double kFactor = (scatteredCharge_ / totalEnergyMeV) * 9.0e10;

  const double thetaXMrad =
      kMilliradiansPerDegree * (scatteredAnglesDeg_[0] - mdmAngleDeg_);
  const double thetaYMrad = kMilliradiansPerDegree * scatteredAnglesDeg_[1];
  const double thetaXRad = thetaXMrad / 1000.0;
  const double thetaYRad = thetaYMrad / 1000.0;

  State state;
  state.x = 0.0;
  state.y = 0.0;
  state.z = 0.0;
  state.vx = speedCmPerSecond * std::sin(thetaXRad) * std::cos(thetaYRad);
  state.vy = speedCmPerSecond * std::sin(thetaYRad);
  state.vz = speedCmPerSecond * std::cos(thetaXRad) * std::cos(thetaYRad);

  const auto stopRay = [&]() {
    firstWireX_ = kStoppedPosition;
    firstWireY_ = kStoppedPosition;
    firstWireAngXDeg_ = 0.0;
    firstWireAngYDeg_ = 0.0;
  };

  const auto propagateDrift = [&](double distanceCm) {
    const double dt = distanceCm / std::abs(state.vz);
    state.x += dt * state.vx;
    state.y += dt * state.vy;
    state.z = 0.0;
  };

  const auto applyCollimator = [&](std::size_t element) {
    const double type = DeckValue(element, 1);
    const double xCenter = DeckValue(element, 2);
    const double yCenter = DeckValue(element, 3);
    const double xMax = DeckValue(element, 4);
    const double yMax = DeckValue(element, 5);

    if (type == 0.0) {
      return std::abs(state.x - xCenter) <= xMax &&
             std::abs(state.y - yCenter) <= yMax;
    }

    const double xc = (state.x - xCenter) / xMax;
    const double yc = (state.y - yCenter) / yMax;
    return xc * xc + yc * yc <= 1.0;
  };

  propagateDrift(DeckValue(kEntranceDriftElement, 1));
  if (!applyCollimator(kEntranceCollimatorElement)) {
    stopRay();
    return;
  }

  propagateDrift(DeckValue(kIntermediateDriftElement, 1));

  {
    const double a = DeckValue(kMultipoleElement, 10);
    const double b = DeckValue(kMultipoleElement, 11);
    const double l = DeckValue(kMultipoleElement, 12);

    State localState;
    localState.x = state.x;
    localState.y = state.y;
    localState.z = state.z - (a + 0.5 * l);
    localState.vx = state.vx;
    localState.vy = state.vy;
    localState.vz = state.vz;

    const FieldFunction fieldFunction = [&](const State& sample) {
      return multipoleMap_.Evaluate(sample.x, sample.y, sample.z);
    };
    const PlaneFunction planeFunction = [&](const State& sample) {
      return sample.z - (b + 0.5 * l);
    };

    localState = IntegrateToPlane(localState, fieldFunction, planeFunction,
                                  kMultipoleStepCm, speedCmPerSecond, kFactor,
                                  "multipole");

    state.x = localState.x;
    state.y = localState.y;
    state.z = 0.0;
    state.vx = localState.vx;
    state.vy = localState.vy;
    state.vz = localState.vz;
  }

  {
    const double a = DeckValue(kDipoleElement, 11);
    const double b = DeckValue(kDipoleElement, 12);
    const double radius = DeckValue(kDipoleElement, 14);
    const double phi = DeckValue(kDipoleElement, 16);
    const double alpha = DegreesToRadians(DeckValue(kDipoleElement, 17));
    const double beta = DegreesToRadians(DeckValue(kDipoleElement, 18));
    const double xcr1 = DeckValue(kDipoleElement, 43);
    const double xcr2 = DeckValue(kDipoleElement, 44);

    const double cosAlpha = std::cos(alpha);
    const double sinAlpha = std::sin(alpha);
    const double cosBeta = std::cos(beta);
    const double sinBeta = std::sin(beta);
    const double phiMinusAlphaMinusBeta =
        DegreesToRadians(phi - DeckValue(kDipoleElement, 17) -
                         DeckValue(kDipoleElement, 18));
    const double cosRot = std::cos(phiMinusAlphaMinusBeta);
    const double sinRot = std::sin(phiMinusAlphaMinusBeta);
    const double cosPb =
        std::cos(DegreesToRadians(0.5 * phi - DeckValue(kDipoleElement, 18)));
    const double sinPb =
        std::sin(DegreesToRadians(0.5 * phi - DeckValue(kDipoleElement, 18)));
    const double sinHalfPhi = std::sin(DegreesToRadians(0.5 * phi));
    const double tx = 2.0 * radius * sinHalfPhi * sinPb;
    const double tz = 2.0 * radius * sinHalfPhi * cosPb;

    const auto localToOutputD = [&](const State& localState) {
      const double xB =
          -localState.x * cosAlpha - localState.z * sinAlpha;
      const double zB =
          localState.x * sinAlpha - localState.z * cosAlpha;
      const double vxB =
          -localState.vx * cosAlpha - localState.vz * sinAlpha;
      const double vzB =
          localState.vx * sinAlpha - localState.vz * cosAlpha;

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
    localState.x = state.x + xcr1;
    localState.y = state.y;
    localState.z = state.z - a;
    localState.vx = state.vx;
    localState.vy = state.vy;
    localState.vz = state.vz;

    const FieldFunction fieldFunction = [&](const State& sample) {
      return dipoleMap_.Evaluate(sample.x, sample.y, sample.z);
    };
    const PlaneFunction planeFunction = [&](const State& sample) {
      return localToOutputD(sample).z;
    };

    localState = IntegrateToPlane(localState, fieldFunction, planeFunction,
                                  kDipoleStepCm, speedCmPerSecond, kFactor,
                                  "dipole");
    state = localToOutputD(localState);
    state.z = 0.0;
  }

  {
    const double a = DeckValue(kSecondMultipoleElement, 10);
    const double b = DeckValue(kSecondMultipoleElement, 11);
    const double l = DeckValue(kSecondMultipoleElement, 12);
    propagateDrift(a + b + l);
  }

  if (!applyCollimator(kExitCollimatorElement)) {
    stopRay();
    return;
  }

  propagateDrift(DeckValue(kFinalDriftElement, 1));

  firstWireX_ = state.x;
  firstWireY_ = state.y;
  firstWireAngXDeg_ = OutputAngleXDegrees(state);
  firstWireAngYDeg_ = OutputAngleYDegrees(state, speedCmPerSecond);
}

void MDMFieldMapTrace::GetPositionAngleFirstWire(double& posX,
                                                 double& posY,
                                                 double& angX,
                                                 double& angY) const {
  posX = firstWireX_;
  posY = firstWireY_;
  angX = firstWireAngXDeg_;
  angY = firstWireAngYDeg_;
}
