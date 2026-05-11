#pragma once

#include <string>

struct MdmIon {
  int massNumber = 0;
  int atomicNumber = 0;
  int chargeState = 0;
  double neutralMassU = 0.0;
  double ionMassMeV = 0.0;

  double RaytraceMassAmu() const;
};

std::string DefaultMassTablePath();
MdmIon LoadMdmIon(int massNumber,
                  int atomicNumber,
                  int chargeState,
                  const std::string& massTablePath = DefaultMassTablePath());
