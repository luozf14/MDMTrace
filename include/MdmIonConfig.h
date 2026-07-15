#pragma once

#include "MdmIon.h"
#include "json.h"

#include <stdexcept>
#include <string>

namespace mdm {

inline MdmIon ParseScatteredIon(const Json::Value& config) {
  if (!config.isMember("scatteredIon") || !config["scatteredIon"].isObject()) {
    throw std::runtime_error("config must define scatteredIon");
  }

  const Json::Value& ion = config["scatteredIon"];
  for (const char* key : {"massNumber", "atomicNumber", "chargeState"}) {
    if (!ion.isMember(key) || !ion[key].isInt()) {
      throw std::runtime_error(std::string("scatteredIon.") + key +
                               " must be an integer");
    }
  }
  const int massNumber = ion["massNumber"].asInt();
  const int atomicNumber = ion["atomicNumber"].asInt();
  const int chargeState = ion["chargeState"].asInt();
  if (massNumber <= 0) {
    throw std::runtime_error("scatteredIon.massNumber must be greater than zero");
  }
  if (atomicNumber < 0 || massNumber < atomicNumber) {
    throw std::runtime_error(
        "scatteredIon must satisfy atomicNumber >= 0 and massNumber >= atomicNumber");
  }
  if (chargeState < 0 || chargeState > atomicNumber) {
    throw std::runtime_error(
        "scatteredIon.chargeState must be between 0 and atomicNumber");
  }
  return LoadMdmIon(massNumber, atomicNumber, chargeState);
}

}  // namespace mdm
