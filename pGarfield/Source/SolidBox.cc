#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Polygon.hh"
#include "Garfield/SolidBox.hh"

namespace Garfield {

SolidBox::SolidBox(const double cx, const double cy, const double cz,
                   const double lx, const double ly, const double lz)
    : Solid(cx, cy, cz, "SolidBox"), m_lX(lx), m_lY(ly), m_lZ(lz) {}

SolidBox::SolidBox(const double cx, const double cy, const double cz,
                   const double lx, const double ly, const double lz,
                   const double dx, const double dy, const double dz)
    : SolidBox(cx, cy, cz, lx, ly, lz) {
  SetDirection(dx, dy, dz);
}

bool SolidBox::IsInside(const double x, const double y, const double z,
                        const bool /*tesselated*/) const {
  // Transform the point to local coordinates.
  double u = x, v = y, w = z;
  ToLocal(x, y, z, u, v, w);

  // See whether the point is inside.
  if (fabs(u) > m_lX || fabs(v) > m_lY || fabs(w) > m_lZ) {
    return false;
  }
  return true;
}

bool SolidBox::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                              double& xmax, double& ymax, double& zmax) const {
  if (m_cTheta == 1. && m_cPhi == 1.) {
    xmin = m_cX - m_lX;
    xmax = m_cX + m_lX;
    ymin = m_cY - m_lY;
    ymax = m_cY + m_lY;
    zmin = m_cZ - m_lZ;
    zmax = m_cZ + m_lZ;
    return true;
  }

  const double dd = sqrt(m_lX * m_lX + m_lY * m_lY + m_lZ * m_lZ);
  xmin = m_cX - dd;
  xmax = m_cX + dd;
  ymin = m_cY - dd;
  ymax = m_cY + dd;
  zmin = m_cZ - dd;
  zmax = m_cZ + dd;
  return true;
}

void SolidBox::SetHalfLengthX(const double lx) {
  if (lx > 0.) {
    m_lX = lx;
  } else {
    std::cerr << "SolidBox::SetHalfLengthX: Half-length must be > 0.\n";
  }
}

void SolidBox::SetHalfLengthY(const double ly) {
  if (ly > 0.) {
    m_lY = ly;
  } else {
    std::cerr << "SolidBox::SetHalfLengthY: Half-length must be > 0.\n";
  }
}

void SolidBox::SetHalfLengthZ(const double lz) {
  if (lz > 0.) {
    m_lZ = lz;
  } else {
    std::cerr << "SolidBox::SetHalfLengthZ: Half-length must be > 0.\n";
  }
}

bool SolidBox::SolidPanels(std::vector<Panel>& panels) {
  const auto id = GetId();
  const auto nPanels = panels.size();
  double xv0, yv0, zv0;
  double xv1, yv1, zv1;
  double xv2, yv2, zv2;
  double xv3, yv3, zv3;
  // Draw the 6 sides of the box, start with the x = xmin face.
  if (m_lY > 0 && m_lZ > 0) {
    ToGlobal(-m_lX, -m_lY, -m_lZ, xv0, yv0, zv0);
    ToGlobal(-m_lX, +m_lY, -m_lZ, xv1, yv1, zv1);
    ToGlobal(-m_lX, +m_lY, +m_lZ, xv2, yv2, zv2);
    ToGlobal(-m_lX, -m_lY, +m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = -m_cPhi * m_cTheta;
    panel.b = -m_sPhi * m_cTheta;
    panel.c = +m_sTheta;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // The x = xmax face.
  if (m_lX > 0 && m_lY > 0 && m_lZ > 0) {
    ToGlobal(+m_lX, -m_lY, -m_lZ, xv0, yv0, zv0);
    ToGlobal(+m_lX, +m_lY, -m_lZ, xv1, yv1, zv1);
    ToGlobal(+m_lX, +m_lY, +m_lZ, xv2, yv2, zv2);
    ToGlobal(+m_lX, -m_lY, +m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = m_cPhi * m_cTheta;
    panel.b = m_sPhi * m_cTheta;
    panel.c = -m_sTheta;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // The y = ymin face.
  if (m_lX > 0 && m_lZ > 0) {
    ToGlobal(-m_lX, -m_lY, -m_lZ, xv0, yv0, zv0);
    ToGlobal(+m_lX, -m_lY, -m_lZ, xv1, yv1, zv1);
    ToGlobal(+m_lX, -m_lY, +m_lZ, xv2, yv2, zv2);
    ToGlobal(-m_lX, -m_lY, +m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = m_sPhi;
    panel.b = -m_cPhi;
    panel.c = 0;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // The y = ymax face.
  if (m_lX > 0 && m_lY > 0 && m_lZ > 0) {
    ToGlobal(-m_lX, +m_lY, -m_lZ, xv0, yv0, zv0);
    ToGlobal(+m_lX, +m_lY, -m_lZ, xv1, yv1, zv1);
    ToGlobal(+m_lX, +m_lY, +m_lZ, xv2, yv2, zv2);
    ToGlobal(-m_lX, +m_lY, +m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = -m_sPhi;
    panel.b = +m_cPhi;
    panel.c = 0;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // The z = zmin face.
  if (m_lX > 0 && m_lY > 0) {
    ToGlobal(-m_lX, -m_lY, -m_lZ, xv0, yv0, zv0);
    ToGlobal(-m_lX, +m_lY, -m_lZ, xv1, yv1, zv1);
    ToGlobal(+m_lX, +m_lY, -m_lZ, xv2, yv2, zv2);
    ToGlobal(+m_lX, -m_lY, -m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = -m_cPhi * m_sTheta;
    panel.b = -m_sPhi * m_sTheta;
    panel.c = -m_cTheta;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // The z = zmax face.
  if (m_lX > 0 && m_lY > 0 && m_lZ > 0) {
    ToGlobal(-m_lX, -m_lY, +m_lZ, xv0, yv0, zv0);
    ToGlobal(-m_lX, +m_lY, +m_lZ, xv1, yv1, zv1);
    ToGlobal(+m_lX, +m_lY, +m_lZ, xv2, yv2, zv2);
    ToGlobal(+m_lX, -m_lY, +m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = +m_cPhi * m_sTheta;
    panel.b = +m_sPhi * m_sTheta;
    panel.c = +m_cTheta;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  // Done, check panel count.
  std::cout << "SolidBox::SolidPanels: " << panels.size() - nPanels
            << " panels.\n";
  return true;
}

double SolidBox::GetDiscretisationLevel(const Panel& panel) {

  // Transform the normal vector to local coordinates.
  double u = 0., v = 0., w = 0.;
  VectorToLocal(panel.a, panel.b, panel.c, u, v, w);
  // Identify the vector.
  if (u > std::max(std::abs(v), std::abs(w))) {
    return m_dis[0];
  } else if (u < -std::max(std::abs(v), std::abs(w))) {
    return m_dis[1];
  } else if (v > std::max(std::abs(u), std::abs(w))) {
    return m_dis[2];
  } else if (v < -std::max(std::abs(u), std::abs(w))) {
    return m_dis[3];
  } else if (w > std::max(std::abs(u), std::abs(v))) {
    return m_dis[4];
  } else if (w < -std::max(std::abs(u), std::abs(v))) {
    return m_dis[5];
  }
  if (m_debug) {
    std::cout << m_className << "::GetDiscretisationLevel:\n"
              << "    Found no match for the panel; return first value.\n";
  }
  return m_dis[0];
}

void SolidBox::Cut(const double x0, const double y0, const double z0,
                   const double xn, const double yn, const double zn,
                   std::vector<Panel>& panels) { 

  //-----------------------------------------------------------------------
  //   PLABXC - Cuts box with a plane.
  //-----------------------------------------------------------------------

  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  // Draw all 12 lines and cut.
  // The line (xmin,ymin,zmin) to (xmax,ymin,zmin).
  double x1, y1, z1;
  ToGlobal(-m_lX, -m_lY, -m_lZ, x1, y1, z1);
  double x2, y2, z2;
  ToGlobal(+m_lX, -m_lY, -m_lZ, x2, y2, z2);
  double xc, yc, zc;
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmin,ymax,zmin).
  ToGlobal(-m_lX, +m_lY, -m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmin,ymin,zmax).
  ToGlobal(-m_lX, -m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // The line (xmax,ymax,zmin) to (xmin,ymax,zmin).
  ToGlobal(+m_lX, +m_lY, -m_lZ, x1, y1, z1);
  ToGlobal(-m_lX, +m_lY, -m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmax,ymin,zmin).
  ToGlobal(+m_lX, -m_lY, -m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmax,ymax,zmax).
  ToGlobal(+m_lX, +m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // The line (xmin,ymax,zmax) to (xmax,ymax,zmax).
  ToGlobal(-m_lX, +m_lY, +m_lZ, x1, y1, z1);
  ToGlobal(+m_lX, +m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmin,ymin,zmax).
  ToGlobal(-m_lX, -m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmin,ymax,zmin).
  ToGlobal(-m_lX, +m_lY, -m_lZ, x1, y1, z1);
  ToGlobal(-m_lX, +m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // The line (xmax,ymin,zmax) to (xmin,ymin,zmax).
  ToGlobal(+m_lX, -m_lY, +m_lZ, x1, y1, z1);
  ToGlobal(-m_lX, -m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmax,ymax,zmax).
  ToGlobal(+m_lX, +m_lY, +m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmax,ymin,zmin).
  ToGlobal(+m_lX, -m_lY, -m_lZ, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
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
