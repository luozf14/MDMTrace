#pragma once

#include "json.h"

#include <fstream>
#include <stdexcept>
#include <string>

namespace mdm {

inline Json::Value ReadConfig(const std::string& path) {
  std::ifstream stream(path.c_str());
  if (!stream) {
    throw std::runtime_error("Error in " + path + ":\nCannot open file.");
  }

  Json::CharReaderBuilder builder;
  Json::Value config;
  std::string errors;
  if (!Json::parseFromStream(builder, stream, &config, &errors)) {
    throw std::runtime_error("Error in " + path + ":\n" + errors);
  }
  return config;
}

inline bool RequireBoolean(const Json::Value& object, const char* key) {
  if (!object.isMember(key) || !object[key].isBool()) {
    throw std::runtime_error(std::string(key) + " must be true or false");
  }
  return object[key].asBool();
}

inline unsigned int GetNonnegativeInteger(const Json::Value& object,
                                          const char* key,
                                          unsigned int fallback) {
  if (!object.isObject() || !object.isMember(key)) {
    return fallback;
  }
  if (!object[key].isInt() || object[key].asInt() < 0) {
    throw std::runtime_error(std::string(key) +
                             " must be a non-negative integer");
  }
  return static_cast<unsigned int>(object[key].asInt());
}

}  // namespace mdm
