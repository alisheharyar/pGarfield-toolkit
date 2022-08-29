#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Polygon.hh"
#include "Garfield/SolidRidge.hh"

namespace Garfield {

SolidRidge::SolidRidge(const double cx, const double cy, const double cz,
                       const double lx, const double ly, const double hz, 
                       const double hx)
    : Solid(cx, cy, cz, "SolidRidge"),
      m_lX(lx), 
      m_lY(ly),
      m_hz(hz),
      m_hx(hx) {}

SolidRidge::SolidRidge(const double cx, const double cy, const double cz,
                       const double lx, const double ly, const double hz,
                       const double hx,
                       const double dx, const double dy, const double dz)
    : SolidRidge(cx, cy, cz, lx, ly, hz, hx) {
  SetDirection(dx, dy, dz);
}

bool SolidRidge::IsInside(const double x, const double y, const double z,
                          const bool /*tesselated*/) const {
  // Transform the point to local coordinates.
  double u = x, v = y, w = z;
  ToLocal(x, y, z, u, v, w);

  bool inside = true;
  if (fabs(u) > m_lX || fabs(v) > m_lY || w < 0. || w > m_hz) {
    inside = false;
  } else if (u >= m_hx &&  m_hz * u + (m_lX - m_hx) * v > m_hz * m_lX) {
    inside = false;
  } else if (u <= m_hx && -m_hz * u + (m_lX + m_hx) * v > m_hz * m_lX) {
    inside = false; 
  }
  return inside;
}

bool SolidRidge::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                double& xmax, double& ymax, double& zmax) const {
  if (m_cTheta == 1. && m_cPhi == 1.) {
    xmin = m_cX - m_lX;
    xmax = m_cX + m_lX;
    ymin = m_cY - m_lY;
    ymax = m_cY + m_lY;
    zmin = m_cZ;
    zmax = m_cZ + m_hz;
    return true;
  }

  const double dd = sqrt(m_lX * m_lX + m_lY * m_lY + m_hz * m_hz);
  xmin = m_cX - dd;
  xmax = m_cX + dd;
  ymin = m_cY - dd;
  ymax = m_cY + dd;
  zmin = m_cZ - dd;
  zmax = m_cZ + dd;
  return true;
}

void SolidRidge::SetHalfLengthX(const double lx) {
  if (lx <= 0.) {
    std::cerr << "SolidRidge::SetHalfLengthX: Half-length must be > 0.\n";
    return;
  }
  m_lX = lx;
}

void SolidRidge::SetHalfLengthY(const double ly) {
  if (ly <= 0.) {
    std::cerr << "SolidRidge::SetHalfLengthY: Half-length must be > 0.\n";
    return;
  }
  m_lY = ly;
}

void SolidRidge::SetRidgeHeight(const double hz) {
  if (hz <= 0.) {
    std::cerr << "SolidRidge::SetRidgeHeight: Height must be > 0.\n";
    return;
  }
  m_hz = hz;
}

bool SolidRidge::SolidPanels(std::vector<Panel>& panels) {
  const auto id = GetId();
  const unsigned int nPanels = panels.size();
  // Direction vector.
  const double fnorm = sqrt(m_dX * m_dX + m_dY * m_dY + m_dZ * m_dZ);
  if (fnorm <= 0) {
    std::cerr << "SolidRidge::SolidPanels:\n"
              << "    Zero norm direction vector; no panels generated.\n";
    return false;
  }
  double xv0, yv0, zv0;
  double xv1, yv1, zv1;
  double xv2, yv2, zv2;
  double xv3, yv3, zv3;

  // Draw the 5 sides of the ridge, start with the floor.
  ToGlobal(-m_lX, -m_lY, 0, xv0, yv0, zv0);
  ToGlobal(-m_lX, +m_lY, 0, xv1, yv1, zv1);
  ToGlobal(+m_lX, +m_lY, 0, xv2, yv2, zv2);
  ToGlobal(+m_lX, -m_lY, 0, xv3, yv3, zv3);

  Panel base;
  base.a = -m_cPhi * m_sTheta;
  base.b = -m_sPhi * m_sTheta;
  base.c = -m_cTheta;
  base.xv = {xv0, xv1, xv2, xv3};
  base.yv = {yv0, yv1, yv2, yv3};
  base.zv = {zv0, zv1, zv2, zv3};
  base.colour = m_colour;
  base.volume = id;
  panels.push_back(std::move(base));

  // Side triangles at y=ymin and y=ymax.
  for (unsigned int i = 0; i < 2; ++i) {
    const double y = i == 0 ? -m_lY : +m_lY;
    ToGlobal(-m_lX, y, 0, xv0, yv0, zv0);
    ToGlobal(+m_lX, y, 0, xv1, yv1, zv1);
    ToGlobal(m_hx, y, m_hz, xv2, yv2, zv2);

    const double a = i == 0 ? +m_sPhi : -m_sPhi;
    const double b = i == 0 ? -m_cPhi : +m_cPhi;
    Panel side;
    side.a = a;
    side.b = b;
    side.c = 0.;
    side.xv = {xv0, xv1, xv2};
    side.yv = {yv0, yv1, yv2};
    side.zv = {zv0, zv1, zv2};
    side.colour = m_colour;
    side.volume = id;
    panels.push_back(std::move(side));
  }

  // The roof, parts at +x and -x.
  for (unsigned int i = 0; i < 2; ++i) {
    const double x = i == 0 ? +m_lX : -m_lX;
    ToGlobal(x, -m_lY, 0, xv0, yv0, zv0);
    ToGlobal(x, +m_lY, 0, xv1, yv1, zv1);
    ToGlobal(m_hx, +m_lY, m_hz, xv2, yv2, zv2);
    ToGlobal(m_hx, -m_lY, m_hz, xv3, yv3, zv3);
    const double dx = i == 0 ? m_lX - m_hx : m_lX + m_hx; 
    const double s = sqrt(m_hz * m_hz + dx * dx);
    const double xroof = i == 0 ? m_hz / s : -m_hz / s;
    const double zroof = dx / s;

    Panel roof;
    roof.a = m_cPhi * m_cTheta * xroof + m_cPhi * m_sTheta * zroof;
    roof.b = m_sPhi * m_cTheta * xroof + m_sPhi * m_sTheta * zroof;
    roof.c = -m_sTheta * xroof + m_cTheta * zroof;
    roof.xv = {xv0, xv1, xv2, xv3};
    roof.yv = {yv0, yv1, yv2, yv3};
    roof.zv = {zv0, zv1, zv2, zv3};
    roof.colour = m_colour;
    roof.volume = id;
    panels.push_back(std::move(roof));
  }
  std::cout << "SolidRidge::SolidPanels: " << panels.size() - nPanels
            << " panels.\n";
  return true;
}

double SolidRidge::GetDiscretisationLevel(const Panel& panel) {

  // Transform the normal vector to local coordinates.
  double u = 0., v = 0., w = 0.;
  VectorToLocal(panel.a, panel.b, panel.c, u, v, w);
  // Identify the vector.
  if (v > std::max(std::abs(u), std::abs(w))) {
    return m_dis[2];
  } else if (v < -std::max(std::abs(u), std::abs(w))) {
    return m_dis[3];
  } else if (w < -std::max(std::abs(u), std::abs(v))) {
    return m_dis[4];
  } else if (u > 0) {
    return m_dis[0];
  } else if (u < 0) {
    return m_dis[1];
  }
  if (m_debug) {
    std::cout << m_className << "::GetDiscretisationLevel:\n"
              << "    Found no match for the panel; return first value.\n";
  }
  return m_dis[0];
}

void SolidRidge::Cut(const double x0, const double y0, const double z0,
                     const double xn, const double yn, const double zn,
                     std::vector<Panel>& panels) {
  
  //-----------------------------------------------------------------------
  //   PLATBC - Cuts ridge with a plane.
  //-----------------------------------------------------------------------
 
  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  // Draw all 9 lines and cut.
  // The line (xmin,ymin,0) to (xmax,ymin,0).
  double x1, y1, z1;
  ToGlobal(-m_lX, -m_lY, 0., x1, y1, z1);
  double x2, y2, z2;
  ToGlobal(+m_lX, -m_lY, 0., x2, y2, z2); 
  double xc, yc, zc;
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xmin,ymax,0).
  ToGlobal(-m_lX, +m_lY, 0., x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  // ... to (xh,ymin,zh).
  ToGlobal(m_hx, -m_lY, m_hz, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }
  
  // The line (xmax,ymax,0) to (xmin,ymax,0).
  ToGlobal(+m_lX, +m_lY, 0., x1, y1, z1);
  ToGlobal(-m_lX, +m_lY, 0., x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }

  // ... to (xmax,ymin,0).
  ToGlobal(+m_lX, -m_lY, 0., x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }

  // ... to (xh,ymax,zh).
  ToGlobal(m_hx, +m_lY, m_hz, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }

  // The line (xmin,ymax,0) to (xh,ymax,zh).
  ToGlobal(-m_lX, +m_lY, 0., x1, y1, z1);
  ToGlobal(m_hx, +m_lY, m_hz, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }

  // The line (xh,ymax,zh) to (xh,ymin,zh)
  ToGlobal(m_hx, +m_lY, m_hz, x1, y1, z1);
  ToGlobal(m_hx, -m_lY, m_hz, x2, y2, z2);
  if (Intersect(x1, y1, z1, x2, y2, z2, x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
    xv.push_back(xc);
    yv.push_back(yc);
    zv.push_back(zc);
  }

  // The line (xh,ymin,zh) to (xmax,ymin,0)
  ToGlobal(m_hx, -m_lY, m_hz, x1, y1, z1);
  ToGlobal(+m_lX, -m_lY, 0., x2, y2, z2);
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
