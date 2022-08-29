#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Polygon.hh"
#include "Garfield/SolidSphere.hh"

namespace {

bool InPolyhedron(const std::vector<Garfield::Panel>& panels,
                  const double x, const double y, const double z,
                  const bool inv = false) {

  for (const auto& panel : panels) {
    double d = panel.a * (panel.xv[0] - x) + panel.b * (panel.yv[0] - y) + 
               panel.c * (panel.zv[0] - z);
    if (inv) d *= -1;
    if (d < 0.) return false;
  }
  return true;
}

}

namespace Garfield {

SolidSphere::SolidSphere(const double cx, const double cy, const double cz,
                         const double r)
    : Solid(cx, cy, cz, "SolidSphere") {
  SetRadius(r);
  UpdatePanels();
}

SolidSphere::SolidSphere(const double cx, const double cy, const double cz,
                         const double rmin, const double rmax)
    : Solid(cx, cy, cz, "SolidSphere") {
  SetRadii(rmin, rmax);
  UpdatePanels();
}

bool SolidSphere::IsInside(const double x, const double y, const double z,
                           const bool tesselated) const {
  // Transform the point to local coordinates.
  const double dx = x - m_cX;
  const double dy = y - m_cY;
  const double dz = z - m_cZ;

  if (fabs(dx) > m_rMax || fabs(dy) > m_rMax || fabs(dz) > m_rMax) {
    return false;
  }
  const double r = sqrt(dx * dx + dy * dy + dz * dz);
  if (!tesselated) return (r >= m_rMin && r <= m_rMax);
  if (r > m_rMax || !InPolyhedron(m_panelsO, dx, dy, dz)) return false;
  if (m_rMin > 0.) {
    return (r >= m_rMin || !InPolyhedron(m_panelsI, dx, dy, dz, true));
  }
  return true;
}

bool SolidSphere::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                 double& xmax, double& ymax,
                                 double& zmax) const {
  xmin = m_cX - m_rMax;
  xmax = m_cX + m_rMax;
  ymin = m_cY - m_rMax;
  ymax = m_cY + m_rMax;
  zmin = m_cZ - m_rMax;
  zmax = m_cZ + m_rMax;
  return true;
}

void SolidSphere::SetRadius(const double r) {
  if (r <= 0.) {
    std::cerr << "SolidSphere::SetRadius: Radius must be > 0.\n";
    return;
  }
  m_rMax = r;
  m_rMin = 0.;
  UpdatePanels();
}

void SolidSphere::SetRadii(const double rmin, const double rmax) {
  if (rmax <= 0.) {
    std::cerr << "SolidSphere::SetRadii: Outer radius must be > 0.\n";
    return;
  }
  if (rmin >= rmax) {
    std::cerr << "SolidSphere::SetRadii:\n"
              << "    Outer radius must be > inner radius.\n";
    return;
  }
  m_rMin = rmin;
  m_rMax = rmax;
  UpdatePanels();
}

void SolidSphere::SetMeridians(const unsigned int n) {
  if (n < 3) {
    std::cerr << "SolidSphere::SetMeridians: Number must be >= 3.\n";
    return;
  }
  m_n = n;
  UpdatePanels();
}

bool SolidSphere::SolidPanels(std::vector<Panel>& panels) {

  const auto nPanels = panels.size();
  panels.insert(panels.begin(), m_panelsO.begin(), m_panelsO.end());
  panels.insert(panels.begin(), m_panelsI.begin(), m_panelsI.end());
  std::cout << "SolidSphere::SolidPanels: " << panels.size() - nPanels
            << " panels.\n";
  return true;
}

void SolidSphere::UpdatePanels() {
  std::lock_guard<std::mutex> guard(m_mutex);
  m_panelsO.clear();
  m_panelsI.clear();
  const auto id = GetId();
  MakePanels(id, m_rMax, true, m_panelsO);
  if (m_rMin > 0.) {
    MakePanels(id, m_rMin, false, m_panelsI);
  } 
}

void SolidSphere::MakePanels(const int vol, const double r, const bool out,
                             std::vector<Panel>& panels) const {

  const double dphi = TwoPi / m_n;
  const double dtheta = Pi / m_n;
  // Loop over the sphere.
  for (unsigned int i = 1; i <= m_n; ++i) {
    const double phi0 = (i - 1.) * dphi;
    const double phi1 = phi0 + dphi;
    const double cphi0 = cos(phi0);
    const double sphi0 = sin(phi0);
    const double cphi1 = cos(phi1);
    const double sphi1 = sin(phi1);
    for (unsigned int j = 1; j <= m_n; ++j) {
      const double theta0 = -HalfPi + (j - 1.) * dtheta;
      const double theta1 = theta0 + dtheta;
      const double ctheta0 = cos(theta0);
      const double stheta0 = sin(theta0);
      const double ctheta1 = cos(theta1);
      const double stheta1 = sin(theta1);
      Panel panel;
      // Corners of this parcel.
      if (j == 1) {
        const double xv0 = m_cX + r * cphi0 * ctheta0;
        const double yv0 = m_cY + r * sphi0 * ctheta0;
        const double zv0 = m_cZ + r * stheta0;
        const double xv1 = m_cX + r * cphi1 * ctheta1;
        const double yv1 = m_cY + r * sphi1 * ctheta1;
        const double zv1 = m_cZ + r * stheta1;
        const double xv2 = m_cX + r * cphi0 * ctheta1;
        const double yv2 = m_cY + r * sphi0 * ctheta1;
        const double zv2 = m_cZ + r * stheta1;
        panel.xv = {xv0, xv1, xv2};
        panel.yv = {yv0, yv1, yv2};
        panel.zv = {zv0, zv1, zv2};
      } else if (j == m_n) {
        const double xv0 = m_cX + r * cphi0 * ctheta0;
        const double yv0 = m_cY + r * sphi0 * ctheta0;
        const double zv0 = m_cZ + r * stheta0;
        const double xv1 = m_cX + r * cphi1 * ctheta0;
        const double yv1 = m_cY + r * sphi1 * ctheta0;
        const double zv1 = m_cZ + r * stheta0;
        const double xv2 = m_cX + r * cphi1 * ctheta1;
        const double yv2 = m_cY + r * sphi1 * ctheta1;
        const double zv2 = m_cZ + r * stheta1;
        panel.xv = {xv0, xv1, xv2};
        panel.yv = {yv0, yv1, yv2};
        panel.zv = {zv0, zv1, zv2};
      } else {
        const double xv0 = m_cX + r * cphi0 * ctheta0;
        const double yv0 = m_cY + r * sphi0 * ctheta0;
        const double zv0 = m_cZ + r * stheta0;
        const double xv1 = m_cX + r * cphi1 * ctheta0;
        const double yv1 = m_cY + r * sphi1 * ctheta0;
        const double zv1 = m_cZ + r * stheta0;
        const double xv2 = m_cX + r * cphi1 * ctheta1;
        const double yv2 = m_cY + r * sphi1 * ctheta1;
        const double zv2 = m_cZ + r * stheta1;
        const double xv3 = m_cX + r * cphi0 * ctheta1;
        const double yv3 = m_cY + r * sphi0 * ctheta1;
        const double zv3 = m_cZ + r * stheta1;
        panel.xv = {xv0, xv1, xv2, xv3};
        panel.yv = {yv0, yv1, yv2, yv3};
        panel.zv = {zv0, zv1, zv2, zv3};
      }
      // Inclination angle in theta.
      const double alpha =
          atan2((ctheta0 - ctheta1) * sqrt((1. + cos(phi1 - phi0)) / 2),
                stheta1 - stheta0);
      const double calpha = cos(alpha);
      const double salpha = sin(alpha);
      // Store the panel.
      if (out) {
        panel.a = cos(0.5 * (phi0 + phi1)) * calpha;
        panel.b = sin(0.5 * (phi0 + phi1)) * calpha;
        panel.c = salpha;
      } else {
        panel.a = -cos(0.5 * (phi0 + phi1)) * calpha;
        panel.b = -sin(0.5 * (phi0 + phi1)) * calpha;
        panel.c = -salpha;
      }
      panel.volume = vol; 
      panels.push_back(std::move(panel));
    }
  }
}

double SolidSphere::GetDiscretisationLevel(const Panel& /*panel*/) {
  return m_dis;
}

void SolidSphere::Cut(const double x0, const double y0, const double z0,
                      const double xn, const double yn, const double zn,
                      std::vector<Panel>& panels) {

  //-----------------------------------------------------------------------
  //    PLASPC - Cuts sphere IVOL with a plane.
  //-----------------------------------------------------------------------
  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  // Loop over the sphere.
  const double r = m_rMax;
  for (unsigned int i = 1; i <= m_n; ++i) {
    // phi-Coordinates.
    const double phi0 = TwoPi * (i - 1.) / m_n;
    const double phi1 = TwoPi * i / m_n;
    for (unsigned int j = 1; j <= m_n; ++j) {
      // theta-Coordinates.
      const double theta0 = -HalfPi + Pi * (j - 1.) / m_n;
      const double theta1 = -HalfPi + Pi * j / m_n;
      // Reference point of this square.
      const double x1 = x0 + r * cos(phi0) * cos(theta0);
      const double y1 = y0 + r * sin(phi0) * cos(theta0);
      const double z1 = z0 + r * sin(theta0);
      // The meridian segment, doesn't exist at the S pole.
      if (j > 0) {
        const double x2 = x0 + r * cos(phi1) * cos(theta0);
        const double y2 = y0 + r * sin(phi1) * cos(theta0);
        const double z2 = z0 + r * sin(theta0);
        // Cut with the plane.
        double xc, yc, zc;
        if (Intersect(x1, y1, z1, x2, y2, z2, 
                      x0, y0, z0, xn, yn, zn, xc, yc, zc)) {
          xv.push_back(xc);
          yv.push_back(yc);
          zv.push_back(zc);
        }
      }
      // The latitude.
      const double x2 = x0 + r * cos(phi0) * cos(theta1);
      const double y2 = y0 + r * sin(phi0) * cos(theta1);
      const double z2 = z0 + r * sin(theta1);
      // Cut with the plane.
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
