#include <cmath>
#include <algorithm>
#include <iostream>

#include "Garfield/Polygon.hh"
#include "Garfield/SolidExtrusion.hh"

namespace Garfield {

SolidExtrusion::SolidExtrusion(const double lz,
                               const std::vector<double>& xp,
                               const std::vector<double>& yp)
    : Solid(0, 0, 0, "SolidExtrusion"), m_lZ(lz) {
  SetProfile(xp, yp);
}

SolidExtrusion::SolidExtrusion(const double lz,
                               const std::vector<double>& xp,
                               const std::vector<double>& yp,
                               const double cx, const double cy,  
                               const double cz,
                               const double dx, const double dy, 
                               const double dz)
    : SolidExtrusion(lz, xp, yp) {
  m_cX = cx;
  m_cY = cy;
  m_cZ = cz;
  SetDirection(dx, dy, dz);
}

bool SolidExtrusion::IsInside(const double x, const double y, const double z,
                              const bool /*tesselated*/) const {

  if (m_xp.empty()) return false;
  // Transform the point to local coordinates.
  double u = x, v = y, w = z;
  ToLocal(x, y, z, u, v, w);
  if (fabs(w) > m_lZ) return false;
  bool inside = false, edge = false;
  Polygon::Inside(m_xp, m_yp, u, v, inside, edge);
  return inside;
}

bool SolidExtrusion::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                              double& xmax, double& ymax, double& zmax) const {

  if (m_xp.empty()) return false;
  const double x0 = *std::min_element(m_xp.begin(), m_xp.end());
  const double x1 = *std::max_element(m_xp.begin(), m_xp.end());
  const double y0 = *std::min_element(m_yp.begin(), m_yp.end());
  const double y1 = *std::max_element(m_yp.begin(), m_yp.end());
  // Take the margins wide.
  const double r = std::max({fabs(x0), fabs(y0), fabs(x1), fabs(y1)});
  const double d = sqrt(r * r + m_lZ * m_lZ); 
  xmin = m_cX - d;
  xmax = m_cX + d;
  ymin = m_cY - d;
  ymax = m_cY + d;
  zmin = m_cZ - d;
  zmax = m_cZ + d;
  return true;
}

void SolidExtrusion::SetHalfLengthZ(const double lz) {
  if (lz > 0.) {
    m_lZ = lz;
  } else {
    std::cerr << "SolidExtrusion::SetHalfLengthZ: Half-length must be > 0.\n";
  }
}

void SolidExtrusion::SetProfile(const std::vector<double>& xp,
                                const std::vector<double>& yp) {

  if (xp.size() != yp.size()) {
    std::cerr << "SolidExtrusion::SetProfile:\n"
              << "    Mismatch between number of x and y coordinates.\n";
    return;
  }
  const auto np = xp.size();
  if (np < 3) {
    std::cerr << "SolidExtrusion::SetProfile: Too few points; rejected.\n";
    return;
  }
  if (!Polygon::NonTrivial(xp, yp)) {
    std::cerr << "SolidExtrusion::SetProfile: Not a valid polygon.\n";
    return;
  }
  const auto it = std::max_element(xp.begin(), xp.end());
  const unsigned int i0 = std::distance(xp.begin(), it);
  const unsigned int i1 = i0 < np - 1 ? i0 + 1 : 0;
  const unsigned int i2 = i1 < np - 1 ? i1 + 1 : 0; 
  const double det = (xp[i1] - xp[i0]) * (yp[i2] - yp[i0]) -
                     (xp[i2] - xp[i0]) * (yp[i1] - yp[i0]);
  if (det < 0.) {
    m_clockwise = true;
  } else if (det > 0.) {
    m_clockwise = false;
  } else {
    std::cerr << "SolidExtrusion::SetProfile:\n"
              << "    Unable to determine profile orientation;" 
              << "    assuming it is clockwise.\n";
    m_clockwise = true;
  }
  m_xp = xp;
  m_yp = yp;
}

bool SolidExtrusion::SolidPanels(std::vector<Panel>& panels) {
  // -----------------------------------------------------------------------
  //   PLAEXP - Generates a table of polygons for an extrusion.
  // -----------------------------------------------------------------------
  const auto id = GetId();
  const auto nPanels = panels.size();
  if (m_xp.empty()) {
    std::cerr << "SolidExtrusion::SolidPanels: Profile is not defined.\n";
    return false;
  }
  // Direction vector.
  const double fnorm = sqrt(m_dX * m_dX + m_dY * m_dY + m_dZ * m_dZ);
  if (fnorm <= 0) {
    std::cerr << "SolidExtrusion::SolidPanels:\n"
              << "    Zero norm direction vector; no panels generated.\n";
    return false;
  }
  const double a = m_dX / fnorm;
  const double b = m_dY / fnorm;
  const double c = m_dZ / fnorm;
  // Number of points
  const unsigned int np = m_xp.size();
  // Create the top lid.
  if (m_toplid) {
    Panel panel;
    panel.a = a;
    panel.b = b;
    panel.c = c;
    for (unsigned int i = 0; i < np; ++i) {
      // Rotate into place.
      double x, y, z;
      ToGlobal(m_xp[i], m_yp[i], m_lZ, x, y, z);
      panel.xv.push_back(x);
      panel.yv.push_back(y);
      panel.zv.push_back(z);
    }
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // Create the bottom lid.
  if (m_botlid) {
    Panel panel;
    panel.a = -a;
    panel.b = -b;
    panel.c = -c;
    for (unsigned int i = 0; i < np; ++i) {
      // Rotate into place.
      double x, y, z;
      ToGlobal(m_xp[i], m_yp[i], -m_lZ, x, y, z);
      panel.xv.push_back(x);
      panel.yv.push_back(y);
      panel.zv.push_back(z);
    }
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // Create the side panels.
  if (m_lZ > 0) {
    double x0, y0, z0;
    ToGlobal(m_xp.back(), m_yp.back(), -m_lZ, x0, y0, z0);
    double x1, y1, z1;
    ToGlobal(m_xp.back(), m_yp.back(), +m_lZ, x1, y1, z1);
    // Go around the extrusion.
    for (unsigned int i = 0; i < np; ++i) {
      // Bottom and top of the line along the axis of the extrusion.
      double x2, y2, z2;
      ToGlobal(m_xp[i], m_yp[i], +m_lZ, x2, y2, z2);
      double x3, y3, z3;
      ToGlobal(m_xp[i], m_yp[i], -m_lZ, x3, y3, z3);
      // Compute the normal vector.
      const unsigned int k = i == 0 ? np - 1 : i - 1;
      double xn = m_yp[k] - m_yp[i];
      double yn = m_xp[i] - m_xp[k];
      const double fn = sqrt(xn * xn + yn * yn);
      if (fn <= 0) {
        std::cerr << "SolidExtrusion::SolidPanels: Zero norm edge (warning).\n";
        continue;
      }
      if (m_clockwise) {
        xn = xn / fn;
        yn = yn / fn;
      } else {
        xn = -xn / fn;
        yn = -yn / fn;
      }
      Panel panel;
      panel.a = m_cPhi * m_cTheta * xn - m_sPhi * yn;
      panel.b = m_sPhi * m_cTheta * xn + m_cPhi * yn;
      panel.c = -m_sTheta * xn;
      panel.xv = {x0, x1, x2, x3};
      panel.yv = {y0, y1, y2, y3};
      panel.zv = {z0, z1, z2, z3};
      panel.colour = m_colour;
      panel.volume = id;
      panels.push_back(std::move(panel));
      // Shift the points.
      x0 = x3;
      y0 = y3;
      z0 = z3;
      x1 = x2;
      y1 = y2;
      z1 = z2;
    }
  }
  // Done, check panel count.
  std::cout << "SolidExtrusion::SolidPanels: " << panels.size() - nPanels
            << " panels.\n";
  return true;
}

double SolidExtrusion::GetDiscretisationLevel(const Panel& panel) {

  // Transform the normal vector to local coordinates.
  double u = 0., v = 0., w = 0.;
  VectorToLocal(panel.a, panel.b, panel.c, u, v, w);
  // Identify the vector.
  if (w > std::max(fabs(u), fabs(v))) {
    return m_dis[0];
  } else if (w < -std::max(fabs(u), fabs(v))) {
    return m_dis[1];
  }
  return m_dis[2];
}

void SolidExtrusion::Cut(const double x0, const double y0, const double z0,
                         const double xn, const double yn, const double zn,
                         std::vector<Panel>& panels) {

  //-----------------------------------------------------------------------
  //   PLAEXC - Cuts extrusion with a plane.
  //-----------------------------------------------------------------------

  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  const unsigned int np = m_xp.size(); 
  // Go through the lines of the top lid, first point.
  double x1, y1, z1;
  ToGlobal(m_xp.back(), m_yp.back(), m_lZ, x1, y1, z1);
  // Loop over the points.
  for (unsigned int i = 0; i < np; ++i) { 
    double x2, y2, z2;
    ToGlobal(m_xp[i], m_yp[i], m_lZ, x2, y2, z2);
    // Cut with the plane.
    double xc, yc, zc;
    if (Intersect(x1, y1, z1, x2, y2, z2, 
                  x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
      xv.push_back(xc);
      yv.push_back(yc);
      zv.push_back(zc);
    }
    // Shift the coordinates.
    x1 = x2;
    y1 = y2;
    z1 = z2;
  }

  if (m_lZ > 0.) {
    // Go through the lines of the bottom lid, first point.
    ToGlobal(m_xp.back(), m_yp.back(), -m_lZ, x1, y1, z1);
    // Loop over the points.
    for (unsigned int i = 0; i < np; ++i) {
      double x2, y2, z2;
      ToGlobal(m_xp[i], m_yp[i], -m_lZ, x2, y2, z2);
      double xc, yc, zc;
      if (Intersect(x1, y1, z1, x2, y2, z2,
                    x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
        xv.push_back(xc);
        yv.push_back(yc);
        zv.push_back(zc);
      }
      // Shift the coordinates.
      x1 = x2;
      y1 = y2;
      z1 = z2;
    }
    // Go through the ribs.
    for (unsigned int i = 0; i < np; ++i) {
      // Bottom and top of the line along the axis of the extrusion.
      ToGlobal(m_xp[i], m_yp[i], +m_lZ, x1, y1, z1);
      double x2, y2, z2;
      ToGlobal(m_xp[i], m_yp[i], -m_lZ, x2, y2, z2);
      double xc, yc, zc;
      if (Intersect(x1, y1, z1, x2, y2, z2, 
                    x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
        xv.push_back(xc);
        yv.push_back(yc);
        zv.push_back(zc);
      }
    }
  }
  // Get rid of butterflies.
  Polygon::EliminateButterflies(xv, yv, zv);
  if (xv.size() >= 3) {
    Panel panel;
    panel.a = xn;
    panel.b = yn;
    panel.c = zn;
    panel.xv = xv;
    panel.yv = yv;
    panel.zv = zv;
    panel.colour = m_colour;
    panel.volume = GetId();
    panels.push_back(std::move(panel));
  }
}
 
}
