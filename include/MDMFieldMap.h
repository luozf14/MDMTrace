#ifndef MDMFIELDMAP_H
#define MDMFIELDMAP_H

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

struct MDMFieldMapMetadata {
  std::string magnetName;
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  std::array<double, 3> originCm{0.0, 0.0, 0.0};
  std::array<double, 3> spacingCm{0.0, 0.0, 0.0};
  std::map<std::string, std::string> fields;
};

class MDMFieldMap {
 public:
  MDMFieldMap() = default;
  MDMFieldMap(MDMFieldMapMetadata metadata,
              std::vector<float> bx,
              std::vector<float> by,
              std::vector<float> bz);

  static MDMFieldMap Load(const std::string& path);
  void Save(const std::string& path) const;

  std::array<double, 3> Evaluate(double xCm, double yCm, double zCm) const;

  const MDMFieldMapMetadata& GetMetadata() const;
  const std::vector<float>& GetBx() const;
  const std::vector<float>& GetBy() const;
  const std::vector<float>& GetBz() const;

 private:
  std::size_t Index(std::size_t ix, std::size_t iy, std::size_t iz) const;

  MDMFieldMapMetadata metadata_;
  std::vector<float> bx_;
  std::vector<float> by_;
  std::vector<float> bz_;
};

#endif
