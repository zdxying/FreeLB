/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// geometry_build.cpp
// Reads the shared .geo recipe and builds a .flbgeo file.
//
// The recipe is parsed by hand instead of with iniReader because the order of
// [source.*] sections is part of the semantics: the builder applies sources in
// that order, and a later solid may overwrite an earlier domain cell.
// [patch.*] sections are applied after every source: entity/kind/box select
// Domain cells of the surface band and the patch is stored in the record.

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "head.h"
#include "geometry/basic_geometry.hh"
#include "io/geometry_file.hh"
#include "io/stlreader.hh"
#include "io/voxel_geometry_builder.hh"

using T = FLOAT;

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "[geometry_build] " << message << std::endl;
  std::exit(1);
}

std::string trim(const std::string& text) {
  const std::string::size_type first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  const std::string::size_type last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::string location(const std::string& path, int line) {
  return path + ":" + std::to_string(line);
}

T parseReal(const std::string& text, const std::string& where) {
  if (text.empty()) fail(where + ": expected a number, got an empty value");
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE ||
      !std::isfinite(value)) {
    fail(where + ": invalid number '" + text + "'");
  }
  const T result = static_cast<T>(value);
  if (!std::isfinite(result)) fail(where + ": number '" + text + "' is out of range");
  return result;
}

long long parseInteger(const std::string& text, const std::string& where) {
  if (text.empty()) fail(where + ": expected an integer, got an empty value");
  errno = 0;
  char* end = nullptr;
  const long long value = std::strtoll(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE) {
    fail(where + ": invalid integer '" + text + "'");
  }
  return value;
}

std::uint16_t parseU16(const std::string& text, const std::string& what,
                       const std::string& where) {
  const long long value = parseInteger(text, where + " (" + what + ")");
  if (value < 0 || value > 65535) {
    fail(where + ": " + what + " must be in [0, 65535], got " + text);
  }
  return static_cast<std::uint16_t>(value);
}

std::vector<std::string> splitList(const std::string& text) {
  std::vector<std::string> tokens;
  std::string::size_type start = 0;
  while (true) {
    const std::string::size_type comma = text.find(',', start);
    tokens.push_back(trim(text.substr(start, comma == std::string::npos
                                               ? std::string::npos
                                               : comma - start)));
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  return tokens;
}

std::vector<T> parseRealList(const std::string& text, std::size_t expected,
                             const std::string& where) {
  const std::vector<std::string> tokens = splitList(text);
  if (tokens.size() != expected) {
    fail(where + ": expected " + std::to_string(expected) +
         " comma-separated numbers, got " + std::to_string(tokens.size()));
  }
  std::vector<T> values;
  values.reserve(tokens.size());
  for (const std::string& token : tokens) values.push_back(parseReal(token, where));
  return values;
}

std::vector<int> parseIntegerList(const std::string& text, std::size_t expected,
                                  const std::string& where) {
  const std::vector<std::string> tokens = splitList(text);
  if (tokens.size() != expected) {
    fail(where + ": expected " + std::to_string(expected) +
         " comma-separated integers, got " + std::to_string(tokens.size()));
  }
  std::vector<int> values;
  values.reserve(tokens.size());
  for (const std::string& token : tokens) {
    const long long value = parseInteger(token, where);
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
      fail(where + ": integer '" + token + "' is out of range");
    }
    values.push_back(static_cast<int>(value));
  }
  return values;
}

struct Entry {
  std::string Value;
  int Line;
};

using Entries = std::map<std::string, Entry>;

struct RawSource {
  std::string Name;
  int Line;
  Entries Values;
};

struct RawRecipe {
  bool HasGeometry{false};
  int GeometryLine{0};
  Entries Geometry;
  std::vector<RawSource> Sources;
  std::vector<RawSource> Patches;
};

const Entry& requiredEntry(const Entries& entries, const std::string& key,
                           const std::string& where) {
  const auto it = entries.find(key);
  if (it == entries.end()) fail(where + ": missing required key '" + key + "'");
  return it->second;
}

const Entry* findEntry(const Entries& entries, const std::string& key) {
  const auto it = entries.find(key);
  return it == entries.end() ? nullptr : &it->second;
}

RawRecipe parseRecipe(const std::string& path) {
  std::ifstream in(path);
  if (!in) fail("cannot open recipe '" + path + "'");

  RawRecipe recipe;
  bool inGeometry = false;
  std::size_t currentSource = std::string::npos;
  std::size_t currentPatch = std::string::npos;

  std::string line;
  int lineNumber = 0;
  while (std::getline(in, line)) {
    ++lineNumber;
    const std::string text = trim(line);
    if (text.empty() || text.front() == '#') continue;

    if (text.front() == '[') {
      if (text.back() != ']') {
        fail(location(path, lineNumber) + ": malformed section header '" + text +
             "'");
      }
      const std::string name = trim(text.substr(1, text.size() - 2));
      inGeometry = false;
      currentSource = std::string::npos;
      currentPatch = std::string::npos;
      if (name == "geometry") {
        if (recipe.HasGeometry) {
          fail(location(path, lineNumber) + ": duplicate [geometry] section");
        }
        recipe.HasGeometry = true;
        recipe.GeometryLine = lineNumber;
        inGeometry = true;
      } else if (name.rfind("source.", 0) == 0) {
        const std::string sourceName = trim(name.substr(7));
        if (sourceName.empty()) {
          fail(location(path, lineNumber) + ": [source.] needs a non-empty name");
        }
        for (const RawSource& source : recipe.Sources) {
          if (source.Name == sourceName) {
            fail(location(path, lineNumber) + ": duplicate section [source." +
                 sourceName + "] (first defined at line " +
                 std::to_string(source.Line) + ")");
          }
        }
        recipe.Sources.push_back(RawSource{sourceName, lineNumber, {}});
        currentSource = recipe.Sources.size() - 1;
      } else if (name.rfind("patch.", 0) == 0) {
        const std::string patchName = trim(name.substr(6));
        if (patchName.empty()) {
          fail(location(path, lineNumber) + ": [patch.] needs a non-empty name");
        }
        for (const RawSource& patch : recipe.Patches) {
          if (patch.Name == patchName) {
            fail(location(path, lineNumber) + ": duplicate section [patch." +
                 patchName + "] (first defined at line " +
                 std::to_string(patch.Line) + ")");
          }
        }
        recipe.Patches.push_back(RawSource{patchName, lineNumber, {}});
        currentPatch = recipe.Patches.size() - 1;
      } else {
        fail(location(path, lineNumber) + ": unknown section [" + name + "]");
      }
      continue;
    }

    std::string where;
    Entries* entries = nullptr;
    if (inGeometry) {
      where = location(path, recipe.GeometryLine) + ": [geometry]";
      entries = &recipe.Geometry;
    } else if (currentSource != std::string::npos) {
      const RawSource& source = recipe.Sources[currentSource];
      where = location(path, source.Line) + ": [source." + source.Name + "]";
      entries = &recipe.Sources[currentSource].Values;
    } else if (currentPatch != std::string::npos) {
      const RawSource& patch = recipe.Patches[currentPatch];
      where = location(path, patch.Line) + ": [patch." + patch.Name + "]";
      entries = &recipe.Patches[currentPatch].Values;
    } else {
      fail(location(path, lineNumber) + ": key/value outside any section");
    }

    const std::string::size_type equals = text.find('=');
    if (equals == std::string::npos) {
      fail(location(path, lineNumber) + ": expected 'key = value', got '" + text +
           "'");
    }
    const std::string key = trim(text.substr(0, equals));
    const std::string value = trim(text.substr(equals + 1));
    if (key.empty()) fail(location(path, lineNumber) + ": empty key");
    if (value.empty()) {
      fail(location(path, lineNumber) + ": empty value for key '" + key + "'");
    }

    const auto inserted = entries->emplace(key, Entry{value, lineNumber});
    if (!inserted.second) {
      fail(location(path, lineNumber) + ": duplicate key '" + key + "' (first set at " +
           location(path, inserted.first->second.Line) + ")");
    }
  }

  if (!recipe.HasGeometry) fail("recipe '" + path + "' has no [geometry] section");
  return recipe;
}

struct GeometrySpec {
  int Nx;
  int Ny;
  int Nz;
  Vector<T, 3> DomainMin;
  Vector<T, 3> DomainMax;
  T CellLen;
  T BandWidth;
};

void checkAllowedKeys(const Entries& entries, const std::set<std::string>& common,
                      const std::set<std::string>& specific,
                      const std::set<std::string>& allSpecific,
                      const std::string& path, const std::string& section) {
  for (const auto& item : entries) {
    if (common.count(item.first) != 0 || specific.count(item.first) != 0) continue;
    const std::string at =
      location(path, item.second.Line) + ": " + section + ": ";
    if (allSpecific.count(item.first) != 0) {
      fail(at + "key '" + item.first + "' is not valid for this source shape");
    }
    fail(at + "unknown key '" + item.first + "'");
  }
}

GeometrySpec parseGeometry(const RawRecipe& recipe, const std::string& path) {
  const std::string where =
    location(path, recipe.GeometryLine) + ": [geometry]";
  const std::set<std::string> allowed = {"grid", "domain", "cell_len",
                                         "band_width"};
  for (const auto& item : recipe.Geometry) {
    if (allowed.count(item.first) == 0) {
      fail(location(path, item.second.Line) + ": [geometry]: unknown key '" +
           item.first + "'");
    }
  }

  GeometrySpec spec{};
  const Entry& grid = requiredEntry(recipe.Geometry, "grid", where);
  const Entry& domain = requiredEntry(recipe.Geometry, "domain", where);
  const Entry& cellLen = requiredEntry(recipe.Geometry, "cell_len", where);

  const std::vector<int> dims =
    parseIntegerList(grid.Value, 3, location(path, grid.Line) + ": grid");
  for (int i = 0; i < 3; ++i) {
    if (dims[i] <= 0) {
      fail(location(path, grid.Line) + ": grid dimensions must be positive, got " +
           grid.Value);
    }
  }
  spec.Nx = dims[0];
  spec.Ny = dims[1];
  spec.Nz = dims[2];

  const std::vector<T> box =
    parseRealList(domain.Value, 6, location(path, domain.Line) + ": domain");
  spec.DomainMin = Vector<T, 3>{box[0], box[1], box[2]};
  spec.DomainMax = Vector<T, 3>{box[3], box[4], box[5]};

  spec.CellLen = parseReal(cellLen.Value, location(path, cellLen.Line) + ": cell_len");
  if (!(spec.CellLen > T(0))) {
    fail(location(path, cellLen.Line) + ": cell_len must be positive, got " +
         cellLen.Value);
  }

  spec.BandWidth = T(1.5);
  if (const Entry* band = findEntry(recipe.Geometry, "band_width")) {
    spec.BandWidth =
      parseReal(band->Value, location(path, band->Line) + ": band_width");
    if (!(spec.BandWidth > T(0))) {
      fail(location(path, band->Line) + ": band_width must be positive, got " +
           band->Value);
    }
  }

  const char* axes[3] = {"x", "y", "z"};
  const int counts[3] = {spec.Nx, spec.Ny, spec.Nz};
  for (int i = 0; i < 3; ++i) {
    const T extent = spec.DomainMax[i] - spec.DomainMin[i];
    if (!(extent > T(0))) {
      fail(location(path, domain.Line) + ": domain min must be below max on " +
           std::string(axes[i]) + " (min " + std::to_string(spec.DomainMin[i]) +
           ", max " + std::to_string(spec.DomainMax[i]) + ")");
    }
    const T expected = static_cast<T>(counts[i]) * spec.CellLen;
    const T tolerance = T(8) * std::numeric_limits<T>::epsilon() *
                        std::max(T(1), std::max(extent, expected));
    if (std::abs(extent - expected) > tolerance) {
      fail(location(path, domain.Line) + ": domain extent on " +
           std::string(axes[i]) + " is " + std::to_string(extent) +
           ", expected n * cell_len = " + std::to_string(expected) +
           " (check grid, domain and cell_len)");
    }
  }

  return spec;
}

struct SourceSpec {
  enum class Shape { Box, Sphere, Stl };
  std::string Name;
  int Line;
  Shape Kind;
  std::uint16_t Entity;
  std::uint16_t Material;
  flbgeo::VoxelRole Role;
  AABB<T, 3> Box;
  Vector<T, 3> Center;
  T Radius{};
  std::string Path;
  bool Closed{true};
};

SourceSpec parseSource(const RawSource& raw, const std::string& path) {
  const std::string section = "[source." + raw.Name + "]";
  const std::string where = location(path, raw.Line) + ": " + section;

  const std::set<std::string> common = {"shape", "entity", "material", "role"};
  const std::set<std::string> boxKeys = {"box"};
  const std::set<std::string> sphereKeys = {"center", "radius"};
  const std::set<std::string> stlKeys = {"path", "closed"};
  const std::set<std::string> allSpecific = {"box", "center", "radius", "path",
                                             "closed"};

  const Entry& shapeEntry = requiredEntry(raw.Values, "shape", where);
  const std::string shape = shapeEntry.Value;

  SourceSpec spec{};
  spec.Name = raw.Name;
  spec.Line = raw.Line;

  std::set<std::string> specific;
  if (shape == "box") {
    spec.Kind = SourceSpec::Shape::Box;
    specific = boxKeys;
  } else if (shape == "sphere") {
    spec.Kind = SourceSpec::Shape::Sphere;
    specific = sphereKeys;
  } else if (shape == "stl") {
    spec.Kind = SourceSpec::Shape::Stl;
    specific = stlKeys;
  } else {
    fail(location(path, shapeEntry.Line) + ": " + section +
         ": shape must be box, sphere or stl, got '" + shape + "'");
  }

  checkAllowedKeys(raw.Values, common, specific, allSpecific, path, section);

  const Entry& entityEntry = requiredEntry(raw.Values, "entity", where);
  const Entry& materialEntry = requiredEntry(raw.Values, "material", where);
  const Entry& roleEntry = requiredEntry(raw.Values, "role", where);

  spec.Entity = parseU16(entityEntry.Value, "entity",
                         location(path, entityEntry.Line));
  spec.Material = parseU16(materialEntry.Value, "material",
                           location(path, materialEntry.Line));

  if (roleEntry.Value == "domain") {
    spec.Role = flbgeo::VoxelRole::Domain;
  } else if (roleEntry.Value == "solid") {
    spec.Role = flbgeo::VoxelRole::Solid;
  } else {
    fail(location(path, roleEntry.Line) + ": role must be domain or solid, got '" +
         roleEntry.Value + "'");
  }

  switch (spec.Kind) {
    case SourceSpec::Shape::Box: {
      const Entry& boxEntry = requiredEntry(raw.Values, "box", where);
      const std::vector<T> box = parseRealList(
        boxEntry.Value, 6, location(path, boxEntry.Line) + ": box");
      spec.Box = AABB<T, 3>(Vector<T, 3>{box[0], box[1], box[2]},
                            Vector<T, 3>{box[3], box[4], box[5]});
      const Vector<T, 3>& extent = spec.Box.getExtension();
      for (int i = 0; i < 3; ++i) {
        if (!(extent[i] > T(0))) {
          fail(location(path, boxEntry.Line) + ": box min must be below max on "
               "each axis, got " + boxEntry.Value);
        }
      }
      break;
    }
    case SourceSpec::Shape::Sphere: {
      const Entry& centerEntry = requiredEntry(raw.Values, "center", where);
      const Entry& radiusEntry = requiredEntry(raw.Values, "radius", where);
      const std::vector<T> center = parseRealList(
        centerEntry.Value, 3, location(path, centerEntry.Line) + ": center");
      spec.Center = Vector<T, 3>{center[0], center[1], center[2]};
      spec.Radius =
        parseReal(radiusEntry.Value, location(path, radiusEntry.Line) + ": radius");
      if (!(spec.Radius > T(0))) {
        fail(location(path, radiusEntry.Line) + ": radius must be positive, got " +
             radiusEntry.Value);
      }
      break;
    }
    case SourceSpec::Shape::Stl: {
      const Entry& pathEntry = requiredEntry(raw.Values, "path", where);
      spec.Path = pathEntry.Value;
      if (const Entry* closedEntry = findEntry(raw.Values, "closed")) {
        if (closedEntry->Value == "0") {
          spec.Closed = false;
        } else if (closedEntry->Value == "1") {
          spec.Closed = true;
        } else {
          fail(location(path, closedEntry->Line) +
               ": closed must be 0 or 1, got '" + closedEntry->Value + "'");
        }
      }
      break;
    }
  }

  return spec;
}

struct PatchSpec {
  std::uint16_t Entity;
  flbgeo::PatchKind Kind;
  AABB<T, 3> Box;
};

PatchSpec parsePatch(const RawSource& raw, const std::string& path) {
  const std::string section = "[patch." + raw.Name + "]";
  const std::string where = location(path, raw.Line) + ": " + section;
  checkAllowedKeys(raw.Values, {"entity", "kind", "box"}, {}, {}, path, section);

  PatchSpec spec{};
  const Entry& entityEntry = requiredEntry(raw.Values, "entity", where);
  spec.Entity = parseU16(entityEntry.Value, "entity",
                         location(path, entityEntry.Line));

  const Entry& kindEntry = requiredEntry(raw.Values, "kind", where);
  if (kindEntry.Value == "wall") {
    spec.Kind = flbgeo::PatchKind::Wall;
  } else if (kindEntry.Value == "inlet") {
    spec.Kind = flbgeo::PatchKind::Inlet;
  } else if (kindEntry.Value == "outlet") {
    spec.Kind = flbgeo::PatchKind::Outlet;
  } else if (kindEntry.Value == "symmetry") {
    spec.Kind = flbgeo::PatchKind::Symmetry;
  } else {
    fail(location(path, kindEntry.Line) + ": " + section +
         ": kind must be wall, inlet, outlet or symmetry, got '" +
         kindEntry.Value + "'");
  }

  const Entry& boxEntry = requiredEntry(raw.Values, "box", where);
  const std::vector<T> box = parseRealList(
    boxEntry.Value, 6, location(path, boxEntry.Line) + ": box");
  spec.Box = AABB<T, 3>(Vector<T, 3>{box[0], box[1], box[2]},
                        Vector<T, 3>{box[3], box[4], box[5]});
  const Vector<T, 3>& extent = spec.Box.getExtension();
  for (int i = 0; i < 3; ++i) {
    if (!(extent[i] > T(0))) {
      fail(location(path, boxEntry.Line) +
           ": patch box min must be below max on each axis, got " +
           boxEntry.Value);
    }
  }
  return spec;
}

std::filesystem::path resolveStlPath(const std::string& recipePath,
                                     const std::string& sourcePath) {
  const std::filesystem::path base =
    std::filesystem::absolute(std::filesystem::path(recipePath)).parent_path();
  const std::filesystem::path path(sourcePath);
  return path.is_absolute() ? path : base / path;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 4 || std::string(argv[2]) != "-o" || std::string(argv[1]).empty() ||
      std::string(argv[3]).empty()) {
    std::cerr << "usage: geometry_build <recipe.geo> -o <out.flbgeo>" << std::endl;
    if (argc > 1 && std::string(argv[1]).rfind("-", 0) == 0) {
      std::cerr << "geometry_build: unknown argument '" << argv[1] << "'"
                << std::endl;
    } else if (argc > 2 && std::string(argv[2]) != "-o") {
      std::cerr << "geometry_build: unknown argument '" << argv[2] << "'"
                << std::endl;
    }
    return 1;
  }

  const std::string recipePath = argv[1];
  const std::string outputPath = argv[3];

  const RawRecipe raw = parseRecipe(recipePath);
  const GeometrySpec geometrySpec = parseGeometry(raw, recipePath);

  std::vector<SourceSpec> sources;
  sources.reserve(raw.Sources.size());
  for (const RawSource& source : raw.Sources) {
    sources.push_back(parseSource(source, recipePath));
  }

  const AABB<T, 3> domain(geometrySpec.DomainMin, geometrySpec.DomainMax);
  flbgeo::VoxelGeometryBuilder<T> builder(geometrySpec.Nx, geometrySpec.Ny,
                                         geometrySpec.Nz, domain,
                                         geometrySpec.CellLen,
                                         geometrySpec.BandWidth);

  flbgeo::GeometryRecipe<T> recipe;
  recipe.BandWidth = geometrySpec.BandWidth;
  recipe.Sources.reserve(sources.size());

  for (const SourceSpec& source : sources) {
    flbgeo::GeometrySource<T> recipeSource;
    recipeSource.Entity = source.Entity;
    recipeSource.Material = source.Material;
    recipeSource.Role = source.Role;

    switch (source.Kind) {
      case SourceSpec::Shape::Box:
        recipeSource.Kind = flbgeo::GeometrySourceKind::Box;
        recipeSource.Box = source.Box;
        builder.addBox(source.Box, source.Entity, source.Material, source.Role);
        break;
      case SourceSpec::Shape::Sphere:
        recipeSource.Kind = flbgeo::GeometrySourceKind::Sphere;
        recipeSource.Center = source.Center;
        recipeSource.Radius = source.Radius;
        builder.addSphere(source.Center, source.Radius, source.Entity,
                          source.Material, source.Role);
        break;
      case SourceSpec::Shape::Stl: {
        recipeSource.Kind = flbgeo::GeometrySourceKind::Stl;
        recipeSource.Path = source.Path;
        recipeSource.RequireClosed = source.Closed;

        const std::filesystem::path resolved =
          resolveStlPath(recipePath, source.Path);
        std::error_code error;
        if (!std::filesystem::is_regular_file(resolved, error)) {
          fail("source '" + source.Name + "': STL file '" + resolved.string() +
               "' does not exist or is not a regular file");
        }
        StlReader<T> reader(resolved.string(), geometrySpec.CellLen, T(1), 2);
        builder.addStl(reader, source.Entity, source.Material, source.Role,
                       source.Closed);
        break;
      }
    }

    recipe.Sources.push_back(std::move(recipeSource));
  }

  for (const RawSource& rawPatch : raw.Patches) {
    const PatchSpec patch = parsePatch(rawPatch, recipePath);
    builder.addPatch(patch.Entity, patch.Kind, patch.Box);
  }

  const flbgeo::VoxelGeometry<T> geometry = builder.build();
  flbgeo::writeGeometryFile(outputPath, geometry, recipe);

  std::cout << "[geometry_build] grid " << geometrySpec.Nx << "x"
            << geometrySpec.Ny << "x" << geometrySpec.Nz << ", sources "
            << recipe.Sources.size() << ", domain "
            << geometry.countRole(flbgeo::VoxelRole::Domain) << ", solid "
            << geometry.countRole(flbgeo::VoxelRole::Solid) << ", outside "
            << geometry.countRole(flbgeo::VoxelRole::Outside) << ", output "
            << outputPath << std::endl;

  return 0;
}
