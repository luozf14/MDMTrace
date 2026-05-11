#include "MdmIon.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef MdmTrace_SOURCE_DIR
#define MdmTrace_SOURCE_DIR "."
#endif

namespace {

constexpr double kAtomicMassUnitMeV = 931.49410242;
constexpr double kElectronMassMeV = 0.510998950;
constexpr double kRaytraceMassUnitMeV = 931.48;

std::string CleanNumber(std::string text) {
  std::string cleaned;
  cleaned.reserve(text.size());
  for (char c : text) {
    if (c != '#') {
      cleaned.push_back(c);
    }
  }
  return cleaned;
}

std::vector<std::string> Tokens(const std::string& line) {
  std::istringstream in(line);
  std::vector<std::string> tokens;
  std::string token;
  while (in >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

bool ToInt(const std::string& text, int& value) {
  std::istringstream in(CleanNumber(text));
  in >> value;
  return !in.fail();
}

bool ToDouble(const std::string& text, double& value) {
  std::istringstream in(CleanNumber(text));
  in >> value;
  return !in.fail();
}

void CheckIonNumbers(int massNumber, int atomicNumber, int chargeState) {
  if (massNumber <= 0 || atomicNumber <= 0 || chargeState <= 0) {
    throw std::runtime_error("scatteredIon requires positive A, Z, and q");
  }
  if (chargeState > atomicNumber) {
    throw std::runtime_error("scatteredIon chargeState cannot exceed Z");
  }
}

}  // namespace

double MdmIon::RaytraceMassAmu() const {
  return ionMassMeV / kRaytraceMassUnitMeV;
}

std::string DefaultMassTablePath() {
  return std::string(MdmTrace_SOURCE_DIR) + "/dat/mass_1.mas20.txt";
}

MdmIon LoadMdmIon(int massNumber,
                  int atomicNumber,
                  int chargeState,
                  const std::string& massTablePath) {
  CheckIonNumbers(massNumber, atomicNumber, chargeState);

  std::ifstream in(massTablePath.c_str());
  if (!in) {
    throw std::runtime_error("Cannot open AME2020 mass table: " +
                             massTablePath);
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.size() < 2) {
      continue;
    }

    // AME2020 is a fixed-width Fortran table with a leading carriage-control
    // character. Dropping that first character makes whitespace token parsing
    // stable for the fields we need: N-Z, N, Z, A, element, ..., mass.
    const std::vector<std::string> t = Tokens(line.substr(1));
    if (t.size() < 8) {
      continue;
    }

    int a = 0;
    int z = 0;
    if (!ToInt(t[2], z) || !ToInt(t[3], a)) {
      continue;
    }
    if (a != massNumber || z != atomicNumber) {
      continue;
    }

    const std::size_t n = t.size();
    double massInteger = 0.0;
    double massMicroU = 0.0;
    if (!ToDouble(t[n - 3], massInteger) || !ToDouble(t[n - 2], massMicroU)) {
      throw std::runtime_error("AME2020 mass entry for A=" +
                               std::to_string(massNumber) +
                               " Z=" + std::to_string(atomicNumber) +
                               " has no numeric atomic mass");
    }
    const double neutralMassU = massInteger + 1.0e-6 * massMicroU;

    MdmIon ion;
    ion.massNumber = massNumber;
    ion.atomicNumber = atomicNumber;
    ion.chargeState = chargeState;
    ion.neutralMassU = neutralMassU;
    ion.ionMassMeV =
        neutralMassU * kAtomicMassUnitMeV - chargeState * kElectronMassMeV;
    return ion;
  }

  throw std::runtime_error("Ion A=" + std::to_string(massNumber) +
                           " Z=" + std::to_string(atomicNumber) +
                           " was not found in " + massTablePath);
}
