#include <algorithm>
#include <iostream>

#include "Garfield/GeometrySimple.hh"

namespace Garfield {

GeometrySimple::GeometrySimple() : Geometry("GeometrySimple") {}

void GeometrySimple::AddSolid(Solid* solid, Medium* medium) {
  // Make sure the solid and the medium are defined.
  if (!solid || !medium) {
    std::cerr << m_className << "::AddSolid: Null pointer.\n";
    return;
  }

  // Update the bounding box ranges
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  if (!solid->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cerr << m_className << "::AddSolid: Solid has no bounding box.\n";
    return;
  }

  if (m_hasBoundingBox) {
    m_bbMin[0] = std::min(m_bbMin[0], xmin);
    m_bbMin[1] = std::min(m_bbMin[1], ymin);
    m_bbMin[2] = std::min(m_bbMin[2], zmin);
    m_bbMax[0] = std::max(m_bbMax[0], xmax);
    m_bbMax[1] = std::max(m_bbMax[1], ymax);
    m_bbMax[2] = std::max(m_bbMax[2], zmax);
  } else {
    m_bbMin[0] = xmin;
    m_bbMin[1] = ymin;
    m_bbMin[2] = zmin;
    m_bbMax[0] = xmax;
    m_bbMax[1] = ymax;
    m_bbMax[2] = zmax;
    m_hasBoundingBox = true;
  }

  // Add the new solid to the list.
  m_solids.emplace_back(std::make_pair(solid, medium));
}

Solid* GeometrySimple::GetSolid(const double x, const double y,
                                const double z, const bool tesselated) const {
  for (const auto& solid : m_solids) {
    if (solid.first->IsInside(x, y, z, tesselated)) return solid.first;
  }
  return nullptr;
}

Medium* GeometrySimple::GetMedium(
    const double x, const double y, const double z, 
    const bool tesselated) const {
  for (const auto& solid : m_solids) {
    if (solid.first->IsInside(x, y, z, tesselated)) {
      return solid.second;
    }
  }
  return m_medium;
}

Solid* GeometrySimple::GetSolid(const size_t i) const {
  if (i >= m_solids.size()) {
    std::cerr << m_className << "::GetSolid: Index out of range.\n";
    return nullptr;
  }
  return m_solids[i].first;
}

Solid* GeometrySimple::GetSolid(const size_t i, Medium*& medium) const {
  if (i >= m_solids.size()) {
    std::cerr << m_className << "::GetSolid: Index out of range.\n";
    return nullptr;
  }
  medium = m_solids[i].second;
  return m_solids[i].first;
}

void GeometrySimple::Clear() {
  m_solids.clear();
  m_medium = nullptr;
}

void GeometrySimple::PrintSolids() {
  std::cout << m_className << "::PrintSolids:\n";
  const auto nSolids = m_solids.size();
  if (nSolids == 1) {
    std::cout << "    1 solid\n";
  } else {
    std::cout << "    " << nSolids << " solids\n";
  }
  if (m_solids.empty()) return;
  std::cout << "      Index      Type    Medium\n";
  for (size_t i = 0; i < nSolids; ++i) {
    std::cout << "        " << i << "         ";
    if (m_solids[i].first->IsBox()) {
      std::cout << "box       ";
    } else if (m_solids[i].first->IsTube()) {
      std::cout << "tube      ";
    } else if (m_solids[i].first->IsSphere()) {
      std::cout << "sphere    ";
    } else if (m_solids[i].first->IsHole()) {
      std::cout << "hole      ";
    } else if (m_solids[i].first->IsRidge()) {
      std::cout << "ridge     ";
    } else if (m_solids[i].first->IsExtrusion()) {
      std::cout << "extrusion ";
    } else if (m_solids[i].first->IsWire()) {
      std::cout << "wire      ";
    } else {
      std::cout << "unknown  ";
    }
    if (m_solids[i].second) {
      std::cout << m_solids[i].second->GetName() << "\n";
    } else {
      std::cout << " ---\n";
    }
  }
}

bool GeometrySimple::IsInside(const double x, const double y,
                              const double z, const bool tesselated) const {
  if (!IsInBoundingBox(x, y, z)) return false;

  for (const auto& solid : m_solids) {
    if (solid.first->IsInside(x, y, z, tesselated)) return true;
  }
  return false;
}

bool GeometrySimple::IsInBoundingBox(const double x, const double y,
                                     const double z) const {
  if (!m_hasBoundingBox) {
    if (m_debug) {
      std::cerr << m_className << "::IsInBoundingBox:\n"
                << "    Bounding box is not defined.\n";
    }
    return true;
  }

  if (x >= m_bbMin[0] && x <= m_bbMax[0] &&
      y >= m_bbMin[1] && y <= m_bbMax[1] &&
      z >= m_bbMin[2] && z <= m_bbMax[2])
    return true;
  return false;
}
}
