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
  const std::string massTablePath =
      config.isMember("massTablePath") ? config["massTablePath"].asString()
                                       : DefaultMassTablePath();
  return LoadMdmIon(ion["massNumber"].asInt(),
                    ion["atomicNumber"].asInt(),
                    ion["chargeState"].asInt(),
                    massTablePath);
}

}  // namespace mdm
