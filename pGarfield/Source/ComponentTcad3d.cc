#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include "Garfield/ComponentTcad3d.hh"

namespace Garfield {

ComponentTcad3d::ComponentTcad3d() : ComponentTcadBase("Tcad3d") {}

void ComponentTcad3d::ElectricField(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& p, Medium*& m,
                                    int& status) {
  // Assume this will work.
  status = 0;
  ex = ey = ez = p = 0.;
  m = nullptr;
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::ElectricField:\n"
              << "    Field map is not available for interpolation.\n";
    status = -10;
    return;
  }
  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x)) {
    status = -6;
    return;
  }

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], x[2], w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    status = -6;
    return;
  }
  const Element& element = m_elements[i];
  const unsigned int nVertices = ElementVertices(element);
  for (unsigned int j = 0; j < nVertices; ++j) {
    const auto index = element.vertex[j];
    ex += w[j] * m_efield[index][0];
    ey += w[j] * m_efield[index][1];
    ez += w[j] * m_efield[index][2];
    p += w[j] * m_epot[index];
  }
  if (mirr[0]) ex = -ex;
  if (mirr[1]) ey = -ey;
  if (mirr[2]) ez = -ez;
  m = m_regions[element.region].medium;
  if (!m_regions[element.region].drift || !m) status = -5;
}

void ComponentTcad3d::ElectricField(const double x, const double y,
                                    const double z, double& ex, double& ey,
                                    double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

bool ComponentTcad3d::Interpolate(
    const double xin, const double yin, const double zin,
    const std::vector<std::array<double, 3> >& field,
    double& fx, double& fy, double& fz) {

  if (field.empty()) return false;
  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Make sure the point is inside the bounding box.
  if (!InBoundingBox(x)) return false;

  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], x[2], w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = ElementVertices(element);
  for (size_t j = 0; j < nVertices; ++j) {
    const auto index = element.vertex[j];
    fx += w[j] * field[index][0];
    fy += w[j] * field[index][1];
    fz += w[j] * field[index][2];
  }
  if (mirr[0]) fx = -fx;
  if (mirr[1]) fy = -fy;
  if (mirr[2]) fz = -fz;
  return true;
} 

bool ComponentTcad3d::Interpolate(
    const double xin, const double yin, const double zin,
    const std::vector<double>& field, double& f) {

  f = 0.;
  if (field.empty()) return false;
  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  // In case of periodicity, reduce to the cell volume.
  MapCoordinates(x, mirr);
  // Make sure the point is inside the bounding box.
  if (!InBoundingBox(x)) return false;

  std::array<double, nMaxVertices> w;
  const auto i = FindElement(x[0], x[1], x[2], w);
  // Stop if the point is outside the mesh.
  if (i >= m_elements.size()) return false;

  const Element& element = m_elements[i];
  const size_t nVertices = ElementVertices(element);
  for (size_t j = 0; j < nVertices; ++j) {
    f += w[j] * field[element.vertex[j]];
  }
  return true;
}

Medium* ComponentTcad3d::GetMedium(const double xin, const double yin,
                                   const double zin) {
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::GetMedium:\n"
              << "    Field map not available for interpolation.\n";
    return nullptr;
  }

  std::array<double, 3> x = {xin, yin, zin};
  std::array<bool, 3> mirr = {false, false, false};
  MapCoordinates(x, mirr);
  // Check if the point is inside the bounding box.
  if (!InBoundingBox(x)) return nullptr;

  // Determine the shape functions.
  std::array<double, nMaxVertices> w;
  const size_t i = FindElement(x[0], x[1], x[2], w);
  if (i >= m_elements.size()) {
    // Point is outside the mesh.
    return nullptr;
  }
  return m_regions[m_elements[i].region].medium;
}

void ComponentTcad3d::FillTree() {

  // Set up the octree.
  const float hx = 0.5 * (m_bbMax[0] - m_bbMin[0]);
  const float hy = 0.5 * (m_bbMax[1] - m_bbMin[1]);
  const float hz = 0.5 * (m_bbMax[2] - m_bbMin[2]);
  m_tree.reset(new TetrahedralTree(Vec3(m_bbMin[0] + hx, m_bbMin[1] + hy, m_bbMin[2] + hz),
                                   Vec3(hx, hy, hz)));

  // Insert the mesh nodes in the tree.
  const size_t nVertices = m_vertices.size();
  for (size_t i = 0; i < nVertices; ++i) {
    const auto& vtx = m_vertices[i];
    m_tree->InsertMeshNode(Vec3(vtx[0], vtx[1], vtx[2]), i);
  }

  // Insert the mesh elements in the tree.
  const size_t nElements = m_elements.size();
  for (size_t i = 0; i < nElements; ++i) {
    const Element& element = m_elements[i];
    const double bb[6] = {element.bbMin[0], element.bbMin[1], element.bbMin[2],
                          element.bbMax[0], element.bbMax[1], element.bbMax[2]};
    m_tree->InsertMeshElement(bb, i);
  }
}

bool ComponentTcad3d::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                     double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  xmin = m_bbMin[0];
  ymin = m_bbMin[1];
  zmin = m_bbMin[2];
  xmax = m_bbMax[0];
  ymax = m_bbMax[1];
  zmax = m_bbMax[2];
  if (m_periodic[0] || m_mirrorPeriodic[0]) {
    xmin = -std::numeric_limits<double>::infinity();
    xmax = +std::numeric_limits<double>::infinity();
  }
  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    ymin = -std::numeric_limits<double>::infinity();
    ymax = +std::numeric_limits<double>::infinity();
  }
  if (m_periodic[2] || m_mirrorPeriodic[2]) {
    zmin = -std::numeric_limits<double>::infinity();
    zmax = +std::numeric_limits<double>::infinity();
  }
  return true;
}

bool ComponentTcad3d::GetElementaryCell(
    double& xmin, double& ymin, double& zmin,
    double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  xmin = m_bbMin[0];
  ymin = m_bbMin[1];
  zmin = m_bbMin[2];
  xmax = m_bbMax[0];
  ymax = m_bbMax[1];
  zmax = m_bbMax[2];
  return true;
}

size_t ComponentTcad3d::FindElement(
    const double x, const double y, const double z,
    std::array<double, nMaxVertices>& w) const {

  w.fill(0.);

  std::vector<int> elementsToSearch;
  if (m_tree) {
    elementsToSearch = m_tree->GetElementsInBlock(Vec3(x, y, z));
  }
  const size_t nElementsToSearch = m_tree ? elementsToSearch.size() : m_elements.size(); 
  // Loop over the elements.
  for (size_t i = 0; i < nElementsToSearch; ++i) {
    const size_t idx = m_tree ? elementsToSearch[i] : i;
    if (InElement(x, y, z, m_elements[idx], w)) return idx;
  }

  if (m_debug) {
    std::cerr << m_className << "::FindElement:\n"
              << "    Point (" << x << ", " << y << ", " << z
              << ") is outside the mesh.\n";
  }
  return m_elements.size();
}

bool ComponentTcad3d::GetElement(const size_t i, double& vol,
                                 double& dmin, double& dmax, int& type,
                                 std::vector<size_t>& nodes, int& reg) const {
  nodes.clear();
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }

  const Element& element = m_elements[i];
  if (element.type == 2) {
    // Triangle
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const auto& v2 = m_vertices[element.vertex[2]];
    const double vx = (v1[1] - v0[1]) * (v2[2] - v0[2]) - 
                      (v1[2] - v0[2]) * (v2[1] - v0[1]);
    const double vy = (v1[2] - v0[2]) * (v2[0] - v0[0]) - 
                      (v1[0] - v0[0]) * (v2[2] - v0[2]);
    const double vz = (v1[0] - v0[0]) * (v2[1] - v0[1]) - 
                      (v1[1] - v0[1]) * (v2[0] - v0[0]);
    vol = sqrt(vx * vx + vy * vy + vz * vz);
    const double a = sqrt(pow(v1[0] - v0[0], 2) + pow(v1[1] - v0[1], 2) + 
                          pow(v1[2] - v0[2], 2));
    const double b = sqrt(pow(v2[0] - v0[0], 2) + pow(v2[1] - v0[1], 2) + 
                          pow(v2[2] - v0[2], 2));
    const double c = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) + 
                          pow(v1[2] - v2[2], 2));
    dmin = std::min({a, b, c});
    dmax = std::max({a, b, c});
  } else if (element.type == 5) {
    // Tetrahedron
    const auto& v0 = m_vertices[element.vertex[0]];
    const auto& v1 = m_vertices[element.vertex[1]];
    const auto& v2 = m_vertices[element.vertex[2]];
    const auto& v3 = m_vertices[element.vertex[3]];
    vol = fabs((v3[0] - v0[0]) * ((v1[1] - v0[1]) * (v2[2] - v0[2]) -
                                  (v2[1] - v0[1]) * (v1[2] - v0[2])) +
               (v3[1] - v0[1]) * ((v1[2] - v0[2]) * (v2[0] - v0[0]) -
                                  (v2[2] - v0[2]) * (v1[0] - v0[0])) +
               (v3[2] - v0[2]) * ((v1[0] - v0[0]) * (v2[1] - v0[1]) -
                                  (v3[0] - v0[0]) * (v1[1] - v0[1]))) /
          6.;
    // Loop over all pairs of m_vertices.
    for (size_t j = 0; j < nMaxVertices - 1; ++j) {
      const auto& vj = m_vertices[element.vertex[j]];
      for (size_t k = j + 1; k < nMaxVertices; ++k) {
        const auto& vk = m_vertices[element.vertex[k]];
        // Compute distance.
        const double dist = sqrt(pow(vj[0] - vk[0], 2) + pow(vj[1] - vk[1], 2) +
                                 pow(vj[2] - vk[2], 2));
        if (k == 1) {
          dmin = dmax = dist;
        } else {
          if (dist < dmin) dmin = dist;
          if (dist > dmax) dmax = dist;
        }
      }
    }
  } else {
    std::cerr << m_className << "::GetElement:\n"
              << "    Unexpected element type (" << type << ").\n";
    return false;
  }
  const size_t nVertices = ElementVertices(element);
  for (size_t j = 0; j < nVertices; ++j) {
    nodes.push_back(element.vertex[j]);
  }
  reg = element.region;
  return true;
}

bool ComponentTcad3d::GetNode(const size_t i, double& x, double& y,
                              double& z, double& v, double& ex, double& ey,
                              double& ez) const {
  if (i >= m_vertices.size()) {
    std::cerr << m_className << "::GetNode: Index out of range.\n";
    return false;
  }

  x = m_vertices[i][0];
  y = m_vertices[i][1];
  z = m_vertices[i][2];
  if (!m_epot.empty()) v = m_epot[i];
  if (!m_efield.empty()) {
    ex = m_efield[i][0];
    ey = m_efield[i][1];
    ez = m_efield[i][2];
  }
  return true;
}

bool ComponentTcad3d::InTetrahedron(
    const double x, const double y, const double z, const Element& element,
    std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  const auto& v2 = m_vertices[element.vertex[2]];
  const auto& v3 = m_vertices[element.vertex[3]];
  const double x10 = v1[0] - v0[0];
  const double y10 = v1[1] - v0[1];
  const double z10 = v1[2] - v0[2];

  const double x20 = v2[0] - v0[0];
  const double y20 = v2[1] - v0[1];
  const double z20 = v2[2] - v0[2];

  const double x30 = v3[0] - v0[0];
  const double y30 = v3[1] - v0[1];
  const double z30 = v3[2] - v0[2];

  const double x21 = v2[0] - v1[0];
  const double y21 = v2[1] - v1[1];
  const double z21 = v2[2] - v1[2];

  const double x31 = v3[0] - v1[0];
  const double y31 = v3[1] - v1[1];
  const double z31 = v3[2] - v1[2];

  const double x32 = v3[0] - v2[0];
  const double y32 = v3[1] - v2[1];
  const double z32 = v3[2] - v2[2];

  w[0] = (x - v1[0]) * (y21 * z31 - y31 * z21) +
         (y - v1[1]) * (z21 * x31 - z31 * x21) +
         (z - v1[2]) * (x21 * y31 - x31 * y21);

  w[0] /= x10 * (y31 * z21 - y21 * z31) + y10 * (z31 * x21 - z21 * x31) +
          z10 * (x31 * y21 - x21 * y31);
  if (w[0] < 0.) return false;

  w[1] = (x - v2[0]) * (-y20 * z32 + y32 * z20) +
         (y - v2[1]) * (-z20 * x32 + z32 * x20) +
         (z - v2[2]) * (-x20 * y32 + x32 * y20);

  w[1] /= x21 * (y20 * z32 - y32 * z20) + y21 * (z20 * x32 - z32 * x20) +
          z21 * (x20 * y32 - x32 * y20);
  if (w[1] < 0.) return false;

  w[2] = (x - v3[0]) * (y30 * z31 - y31 * z30) +
         (y - v3[1]) * (z30 * x31 - z31 * x30) +
         (z - v3[2]) * (x30 * y31 - x31 * y30);

  w[2] /= x32 * (y31 * z30 - y30 * z31) + y32 * (z31 * x30 - z30 * x31) +
          z32 * (x31 * y30 - x30 * y31);
  if (w[2] < 0.) return false;

  w[3] = (x - v0[0]) * (y20 * z10 - y10 * z20) +
         (y - v0[1]) * (z20 * x10 - z10 * x20) +
         (z - v0[2]) * (x20 * y10 - x10 * y20);

  w[3] /= x30 * (y20 * z10 - y10 * z20) + y30 * (z20 * x10 - z10 * x20) +
          z30 * (x20 * y10 - x10 * y20);
  if (w[3] < 0.) return false;

  if (m_debug) {
    // Reconstruct the point from the local coordinates.
    const double xr = w[0] * v0[0] + w[1] * v1[0] + w[2] * v2[0] + w[3] * v3[0];
    const double yr = w[0] * v0[1] + w[1] * v1[1] + w[2] * v2[1] + w[3] * v3[1];
    const double zr = w[0] * v0[2] + w[1] * v1[2] + w[2] * v2[2] + w[3] * v3[2];
    std::cout << m_className << "::InTetrahedron:\n"
              << "    Original coordinates:      (" << x << ", " << y << ", "
              << z << ")\n"
              << "    Local coordinates:         (" << w[0] << ", " << w[1]
              << ", " << w[2] << ", " << w[3] << ")\n"
              << "    Reconstructed coordinates: (" << xr << ", " << yr << ", "
              << zr << ")\n"
              << "    Checksum: " << w[0] + w[1] + w[2] + w[3] - 1. << "\n";
  }

  return true;
}

bool ComponentTcad3d::InTriangle(
    const double x, const double y, const double z, const Element& element,
    std::array<double, nMaxVertices>& w) const {
  const auto& v0 = m_vertices[element.vertex[0]];
  const auto& v1 = m_vertices[element.vertex[1]];
  const auto& v2 = m_vertices[element.vertex[2]];

  const double v1x = v1[0] - v0[0];
  const double v2x = v2[0] - v0[0];
  const double v1y = v1[1] - v0[1];
  const double v2y = v2[1] - v0[1];
  const double v1z = v1[2] - v0[2];
  const double v2z = v2[2] - v0[2];

  // Check whether the point lies in the plane of the triangle.
  // Compute the coefficients of the plane equation.
  const double a = v1y * v2z - v2y * v1z;
  const double b = v1z * v2x - v2z * v1x;
  const double c = v1x * v2y - v2x * v1y;
  const double d = a * v0[0] + b * v0[1] + c * v0[2];
  // Check if the point satisfies the plane equation.
  if (a * x + b * y + c * z != d) return false;

  // Map (x, y) onto local variables (b, c) such that
  // P = A + b * (B - A) + c * (C - A)
  // A point P is inside the triangle ABC if b, c > 0 and b + c < 1;
  // b, c are also weighting factors for points B, C
  w[1] = ((x - v0[0]) * v2y - (y - v0[1]) * v2x) / (v1x * v2y - v1y * v2x);
  if (w[1] < 0. || w[1] > 1.) return false;
  w[2] = ((v0[0] - x) * v1y - (v0[1] - y) * v1x) / (v1x * v2y - v1y * v2x);
  if (w[2] < 0. || w[1] + w[2] > 1.) return false;

  // Weighting factor for point A
  w[0] = 1. - w[1] - w[2];

  return true;
}

}
