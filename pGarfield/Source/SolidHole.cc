#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Polygon.hh"
#include "Garfield/SolidHole.hh"

namespace {

void CutBox(const std::array<double, 8>& xbox,
            const std::array<double, 8>& ybox,
            const std::array<double, 8>& zbox,
            std::vector<double>& xcut,
            std::vector<double>& ycut,
            std::vector<double>& zcut,
            const double x0, const double y0, const double z0,
            const double a, const double b, const double c) {

  //-----------------------------------------------------------------------
  //   PLABOX - Crossings between a box and a plane.
  //-----------------------------------------------------------------------
  xcut.clear();
  ycut.clear();
  zcut.clear();
  // Compute the, at most 6, crossings between plane and box.
  for (unsigned int i = 0; i < 8; ++i) {
    double xc, yc, zc;
    unsigned int j = i + 1;
    if (i == 3) {
      j = 0;
    } else if (i == 7) {
      j = 4;
    } 
    if (Garfield::Solid::Intersect(xbox[i], ybox[i], zbox[i], 
                                   xbox[j], ybox[j], zbox[j], 
                                   x0, y0, z0, a, b, c, xc, yc, zc)) {
      xcut.push_back(xc);
      ycut.push_back(yc);
      zcut.push_back(zc);
    }
  }
  for (unsigned int i = 0; i < 4; ++i) {
    double xc, yc, zc;
    const unsigned int j = i + 4;
    if (Garfield::Solid::Intersect(xbox[i], ybox[i], zbox[i],
                                   xbox[j], ybox[j], zbox[j],
                                   x0, y0, z0, a, b, c, xc, yc, zc)) {
      xcut.push_back(xc);
      ycut.push_back(yc);
      zcut.push_back(zc);
    }
  }

  // Eliminate the butterflies.
  Garfield::Polygon::EliminateButterflies(xcut, ycut, zcut);
}

}

namespace Garfield {

SolidHole::SolidHole(const double cx, const double cy, const double cz,
                     const double rup, const double rlow, 
                     const double lx, const double ly, const double lz)
    : Solid(cx, cy, cz, "SolidHole"),
      m_rUp(rup),
      m_rLow(rlow),
      m_lX(lx), 
      m_lY(ly),
      m_lZ(lz) {
  Update();
}

SolidHole::SolidHole(const double cx, const double cy, const double cz,
                     const double rup, const double rlow, 
                     const double lx, const double ly, const double lz,
                     const double dx, const double dy, const double dz)
    : SolidHole(cx, cy, cz, rup, rlow, lx, ly, lz) {
  SetDirection(dx, dy, dz);
}

void SolidHole::Update() {
  std::lock_guard<std::mutex> guard(m_mutex);
  const double alpha = Pi / (4. * (m_n - 1.));
  if (m_average) {
    m_fp = 2. / (1. + asinh(tan(alpha)) * cos(alpha) / tan(alpha));
  } else {
    m_fp = 1.;
  }
  m_fi = cos(alpha); 
}

bool SolidHole::IsInside(const double x, const double y, const double z,
                         const bool tesselated) const {
  // Transform the point to local coordinates.
  double u = x, v = y, w = z;
  ToLocal(x, y, z, u, v, w);

  if (fabs(u) > m_lX || fabs(v) > m_lY || fabs(w) > m_lZ) {
    return false;
  }
 
  double r = m_rLow + (w + m_lZ) * (m_rUp - m_rLow) / (2 * m_lZ);
  if (!tesselated) return (u * u + v * v >= r * r);
  const double rho = sqrt(u * u + v * v);
  if (m_average) r *= m_fp;
  if (rho > r) return true;
  if (rho < r * m_fi) return false;

  std::vector<double> xp;
  std::vector<double> yp;
  const double phi0 = -0.5 * HalfPi;
  const double dphi = HalfPi / double(m_n - 1);
  const unsigned int nP = 4 * (m_n - 1);
  for (unsigned int i = 0; i < nP; ++i) {
    // Bottom and top of the line along the axis of the cylinder.
    const double phi = phi0 + dphi * i;
    xp.push_back(r * cos(phi));
    yp.push_back(r * sin(phi));
  }
  bool inside = false;
  bool edge = false;
  Polygon::Inside(xp, yp, u, v, inside, edge);  
  return !inside;
}

bool SolidHole::GetBoundingBox(double& xmin, double& ymin, double& zmin,
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

void SolidHole::SetUpperRadius(const double r) {
  if (r <= 0.) {
    std::cerr << "SolidHole::SetUpperRadius: Radius must be > 0.\n";
    return;
  }
  m_rUp = r;
}

void SolidHole::SetLowerRadius(const double r) {
  if (r <= 0.) {
    std::cerr << "SolidHole::SetLowerRadius: Radius must be > 0.\n";
    return;
  }
  m_rLow = r;
}

void SolidHole::SetHalfLengthX(const double lx) {
  if (lx <= 0.) {
    std::cerr << "SolidHole::SetHalfLengthX: Half-length must be > 0.\n";
    return;
  }
  m_lX = lx;
}

void SolidHole::SetHalfLengthY(const double ly) {
  if (ly <= 0.) {
    std::cerr << "SolidHole::SetHalfLengthY: Half-length must be > 0.\n";
    return;
  }
  m_lY = ly;
}

void SolidHole::SetHalfLengthZ(const double lz) {
  if (lz <= 0.) {
    std::cerr << "SolidHole::SetHalfLengthZ: Half-length must be > 0.\n";
    return;
  }
  m_lZ = lz;
}

void SolidHole::SetSectors(const unsigned int n) {
  if (n < 1) {
    std::cerr << "SolidHole::SetSectors: Number must be > 0.\n";
    return;
  }
  m_n = n;
  Update();
}

bool SolidHole::SolidPanels(std::vector<Panel>& panels) {
  const auto id = GetId();
  const auto nPanels = panels.size();
  // Direction vector.
  const double fnorm = sqrt(m_dX * m_dX + m_dY * m_dY + m_dZ * m_dZ);
  if (fnorm <= 0) {
    std::cerr << "SolidHole::SolidPanels:\n"
              << "    Zero norm direction vector; no panels generated.\n";
    return false;
  }

  // Set the mean or the outer radius.
  double r1 = m_rLow;
  double r2 = m_rUp;
  if (m_average) {
    r1 *= m_fp;
    r2 *= m_fp;
  }

  double xv0, yv0, zv0;
  double xv1, yv1, zv1;
  double xv2, yv2, zv2;
  double xv3, yv3, zv3;
  // Draw the 6 sides of the box, start with the x=xmin face.
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
  // The x=xmax face.
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
  // The y=ymin face.
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
  // The y=ymax face.
  if (m_lX > 0 && m_lY > 0 && m_lZ > 0) {
    ToGlobal(-m_lX, +m_lY, -m_lZ, xv0, yv0, zv0);
    ToGlobal(+m_lX, +m_lY, -m_lZ, xv1, yv1, zv1);
    ToGlobal(+m_lX, +m_lY, +m_lZ, xv2, yv2, zv2);
    ToGlobal(-m_lX, +m_lY, +m_lZ, xv3, yv3, zv3);
    Panel panel;
    panel.a = -m_sPhi;
    panel.b = m_cPhi;
    panel.c = 0;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
  }
  const double phi0 = -0.5 * HalfPi;
  const double dphi = HalfPi / double(m_n - 1);
  // The faces at constant z have a hole, and are drawn in parts.
  for (int iside = -1; iside <= 1; iside += 2) {
    const double r = iside == -1 ? r1 : r2;
    const double w = m_lZ * iside;
    Panel panel;
    panel.a = iside * m_cPhi * m_sTheta;
    panel.b = iside * m_sPhi * m_sTheta;
    panel.c = iside * m_cTheta;
    panel.colour = m_colour;
    panel.volume = id;
    // Loop over the panels.
    for (unsigned int i = 0; i < m_n - 1; ++i) {
      // The panels for x=xmax.
      const double phi1 = phi0 + dphi * i;
      const double phi2 = phi1 + dphi;
      const double t1 = tan(phi1);
      const double t2 = tan(phi2);
      ToGlobal(r * cos(phi1), r * sin(phi1), w, xv0, yv0, zv0);
      ToGlobal(r * cos(phi2), r * sin(phi2), w, xv3, yv3, zv3);
      ToGlobal(m_lX, m_lY * t1, w, xv1, yv1, zv1);
      ToGlobal(m_lX, m_lY * t2, w, xv2, yv2, zv2);
      panel.xv = {xv0, xv1, xv2, xv3};
      panel.yv = {yv0, yv1, yv2, yv3};
      panel.zv = {zv0, zv1, zv2, zv3};
      panels.push_back(panel);
      // The panels for y=ymax.
      const double phi3 = phi1 + HalfPi;
      const double phi4 = phi3 + dphi;
      ToGlobal(r * cos(phi3), r * sin(phi3), w, xv0, yv0, zv0);
      ToGlobal(r * cos(phi4), r * sin(phi4), w, xv3, yv3, zv3);
      ToGlobal(-m_lX * t1, m_lY, w, xv1, yv1, zv1);
      ToGlobal(-m_lX * t2, m_lY, w, xv2, yv2, zv2);
      panel.xv = {xv0, xv1, xv2, xv3};
      panel.yv = {yv0, yv1, yv2, yv3};
      panel.zv = {zv0, zv1, zv2, zv3};
      panels.push_back(panel);
      // The panels for x=xmin.
      const double phi5 = phi3 + HalfPi;
      const double phi6 = phi5 + dphi;
      ToGlobal(r * cos(phi5), r * sin(phi5), w, xv0, yv0, zv0);
      ToGlobal(r * cos(phi6), r * sin(phi6), w, xv3, yv3, zv3);
      ToGlobal(-m_lX, -m_lY * t1, w, xv1, yv1, zv1);
      ToGlobal(-m_lX, -m_lY * t2, w, xv2, yv2, zv2);
      panel.xv = {xv0, xv1, xv2, xv3};
      panel.yv = {yv0, yv1, yv2, yv3};
      panel.zv = {zv0, zv1, zv2, zv3};
      panels.push_back(panel);
      // The panels for y=ymin.
      const double phi7 = phi5 + HalfPi;
      const double phi8 = phi7 + dphi;
      ToGlobal(r * cos(phi7), r * sin(phi7), w, xv0, yv0, zv0);
      ToGlobal(r * cos(phi8), r * sin(phi8), w, xv3, yv3, zv3);
      ToGlobal(m_lX * t1, -m_lY, w, xv1, yv1, zv1);
      ToGlobal(m_lX * t2, -m_lY, w, xv2, yv2, zv2);
      panel.xv = {xv0, xv1, xv2, xv3};
      panel.yv = {yv0, yv1, yv2, yv3};
      panel.zv = {zv0, zv1, zv2, zv3};
      panels.push_back(panel);
    }
  }
  // The panels of the central cylinder, compute the projection angles.
  const double alpha = atan2((r1 - r2) * cos(Pi / (4 * (m_n - 1))), 
                             2 * m_lZ);
  const double ci = cos(alpha);
  const double si = sin(alpha);
  // Initialise loop.
  ToGlobal(r1 * cos(phi0), r1 * sin(phi0), -m_lZ, xv0, yv0, zv0);
  ToGlobal(r2 * cos(phi0), r2 * sin(phi0), +m_lZ, xv1, yv1, zv1);
  // Go around the cylinder.
  const unsigned int nPoints = 4 * m_n - 3;
  for (unsigned int i = 1; i < nPoints; ++i) {
    // Bottom and top of the line along the axis of the cylinder.
    const double phi = phi0 + dphi * i;
    ToGlobal(r2 * cos(phi), r2 * sin(phi), +m_lZ, xv2, yv2, zv2);
    ToGlobal(r1 * cos(phi), r1 * sin(phi), -m_lZ, xv3, yv3, zv3);
    // Store the plane.
    Panel panel;
    const double phim = phi0 + dphi * (i - 0.5);
    const double cm = cos(phim);
    const double sm = sin(phim);
    panel.a = -m_cPhi * m_cTheta * cm * ci +
               m_sPhi *            sm * ci -
               m_cPhi * m_sTheta      * si;
    panel.b = -m_sPhi * m_cTheta * cm * ci -
               m_cPhi *            sm * ci -
               m_sPhi * m_sTheta      * si;
    panel.c = m_sTheta * cm * ci - m_cTheta * si;
    panel.xv = {xv0, xv1, xv2, xv3};
    panel.yv = {yv0, yv1, yv2, yv3};
    panel.zv = {zv0, zv1, zv2, zv3};
    panel.colour = m_colour;
    panel.volume = id;
    panels.push_back(std::move(panel));
    // Shift the points.
    xv0 = xv3;
    yv0 = yv3;
    zv0 = zv3;
    xv1 = xv2;
    yv1 = yv2;
    zv1 = zv2;
  }
  std::cout << "SolidHole::SolidPanels: " << panels.size() - nPanels
            << " panels.\n";
  return true;
}

double SolidHole::GetDiscretisationLevel(const Panel& panel) {
  
  // Transform the normal vector to local coordinates.
  double un = 0., vn = 0., wn = 0.;
  VectorToLocal(panel.a, panel.b, panel.c, un, vn, wn);
  // Transform one of the points (first).
  double u1 = 0., v1 = 0., w1 = 0.;
  ToLocal(panel.xv[0], panel.yv[0], panel.zv[0], u1, v1, w1);
  // Identify the vector.
  if (wn > std::max(std::abs(un), std::abs(vn))) {
    return m_dis[0];
  } else if (wn < -std::max(std::abs(un), std::abs(vn))) {
    return m_dis[1];
  } else if (un * u1 + vn * v1 + wn * w1 < 0) {
    return m_dis[2];
  } else if (un > std::max(std::abs(vn), std::abs(wn))) {
    return m_dis[3];
  } else if (un < -std::max(std::abs(vn), std::abs(wn))) {
    return m_dis[4];
  } else if (vn > std::max(std::abs(un), std::abs(wn))) {
    return m_dis[5];
  } else if (vn < -std::max(std::abs(un), std::abs(wn))) {
    return m_dis[6];
  }
  if (m_debug) {
    std::cout << m_className << "::GetDiscretisationLevel:\n"
              << "    Found no match for the panel; returning first value.\n";
  }
  return m_dis[0];
}

void SolidHole::Cut(const double x0, const double y0, const double z0,
                    const double xn, const double yn, const double zn,
                    std::vector<Panel>& panels) {

  //-----------------------------------------------------------------------
  //   PLACHC - Cuts a cylindrical hole with a plane.
  //-----------------------------------------------------------------------

  // Adjust the mean radii if requested.
  double r1 = m_rLow;
  double r2 = m_rUp;
  if (m_average) {
    r1 *= m_fp;
    r2 *= m_fp;
  }
  std::array<double, 8> xbox;
  std::array<double, 8> ybox;
  std::array<double, 8> zbox;
  // Loop over the boxes that make up the hole.
  const double dphi = HalfPi / (m_n - 1.);
  for (unsigned int i = 1; i <= m_n - 1; ++i) {
    const double phi1 = dphi * (i - 1.);
    const double phi2 = dphi * i;
    const double t1 = tan(-0.5 * HalfPi + phi1);
    const double t2 = tan(-0.5 * HalfPi + phi2);
    // The boxes ending at x=xmax.
    for (int iside : {-1, 1}) {
      const double r = iside < 0 ? r1 : r2;
      const double w = iside * m_lZ;
      const unsigned int k = iside < 0 ? 0 : 4;
      const double phi0 = -0.5 * HalfPi;
      ToGlobal(r * cos(phi0 + phi1), r * sin(phi0 + phi1), w, 
               xbox[k], ybox[k], zbox[k]);
      ToGlobal(r * cos(phi0 + phi2), r * sin(phi0 + phi2), w, 
               xbox[k + 3], ybox[k + 3], zbox[k + 3]);
      ToGlobal(m_lX, m_lY * t1, w, xbox[k + 1], ybox[k + 1], zbox[k + 1]);
      ToGlobal(m_lX, m_lY * t2, w, xbox[k + 2], ybox[k + 2], zbox[k + 2]);
    }
    std::vector<double> xv;
    std::vector<double> yv;
    std::vector<double> zv;
    CutBox(xbox, ybox, zbox, xv, yv, zv, x0, y0, z0, xn, yn, zn);
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
    xv.clear();
    yv.clear();
    zv.clear();
    // The panels for y=ymax.
    for (int iside : {-1, 1}) {
      const double r = iside < 0 ? r1 : r2;
      const double w = iside * m_lZ;
      const unsigned int k = iside < 0 ? 0 : 4;
      const double phi0 = 0.5 * HalfPi;
      ToGlobal(r * cos(phi0 + phi1), r * sin(phi0 + phi1), w, 
               xbox[k], ybox[k], zbox[k]);
      ToGlobal(r * cos(phi0 + phi2), r * sin(phi0 + phi2), w, 
               xbox[k + 3], ybox[k + 3], zbox[k + 3]);
      ToGlobal(-m_lX * t1, m_lY, w, xbox[k + 1], ybox[k + 1], zbox[k + 1]);
      ToGlobal(-m_lX * t2, m_lY, w, xbox[k + 2], ybox[k + 2], zbox[k + 2]);
    }
    CutBox(xbox, ybox, zbox, xv, yv, zv, x0, y0, z0, xn, yn, zn);
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
    xv.clear();
    yv.clear();
    zv.clear();
    // The panels for x=xmin.
    for (int iside : {-1, 1}) {
      const double r = iside < 0 ? r1 : r2;
      const double w = iside * m_lZ;
      const unsigned int k = iside < 0 ? 0 : 4;
      const double phi0 = 1.5 * HalfPi;
      ToGlobal(r * cos(phi0 + phi1), r * sin(phi0 + phi1), w, 
               xbox[k], ybox[k], zbox[k]);
      ToGlobal(r * cos(phi0 + phi2), r * sin(phi0 + phi2), w, 
               xbox[k + 3], ybox[k + 3], zbox[k + 3]);
      ToGlobal(-m_lX, -m_lY * t1, w, xbox[k + 1], ybox[k + 1], zbox[k + 1]);
      ToGlobal(-m_lX, -m_lY * t2, w, xbox[k + 2], ybox[k + 2], zbox[k + 2]);
    }
    CutBox(xbox, ybox, zbox, xv, yv, zv, x0, y0, z0, xn, yn, zn);
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
    xv.clear();
    yv.clear();
    zv.clear();
    // The panels for y=ymin.
    for (int iside : {-1, 1}) {
      const double r = iside < 0 ? r1 : r2;
      const double w = iside * m_lZ;
      const unsigned int k = iside < 0 ? 0 : 4;
      const double phi0 = -1.5 * HalfPi;
      ToGlobal(r * cos(phi0 + phi1), r * sin(phi0 + phi1), w, 
               xbox[k], ybox[k], zbox[k]);
      ToGlobal(r * cos(phi0 + phi2), r * sin(phi0 + phi2), w, 
               xbox[k + 3], ybox[k + 3], zbox[k + 3]);
      ToGlobal(m_lX * t1, -m_lY, w, xbox[k + 1], ybox[k + 1], zbox[k + 1]);
      ToGlobal(m_lX * t2, -m_lY, w, xbox[k + 2], ybox[k + 2], zbox[k + 2]);

    }
    CutBox(xbox, ybox, zbox, xv, yv, zv, x0, y0, z0, xn, yn, zn);
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

}
