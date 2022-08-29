#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/SolidWire.hh"

namespace Garfield {

SolidWire::SolidWire(const double cx, const double cy, const double cz,
                     const double rw, const double lz)
    : Solid(cx, cy, cz, "SolidWire"),
      m_r(rw),
      m_lZ(lz) {}

SolidWire::SolidWire(const double cx, const double cy, const double cz,
                     const double rw, const double lz,
                     const double dx, const double dy, const double dz)
    : SolidWire(cx, cy, cz, rw, lz) {
  SetDirection(dx, dy, dz);
}

bool SolidWire::IsInside(const double x, const double y, const double z,
                         const bool /*tesselated*/) const {
  // Transform the point to local coordinates.
  double u = x, v = y, w = z;
  ToLocal(x, y, z, u, v, w);

  if (fabs(w) > m_lZ) return false;

  const double r = sqrt(u * u + v * v);
  return r < m_r;
}

bool SolidWire::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                               double& xmax, double& ymax, double& zmax) const {
  if (m_cTheta == 1. && m_cPhi == 1.) {
    xmin = m_cX - m_r;
    xmax = m_cX + m_r;
    ymin = m_cY - m_r;
    ymax = m_cY + m_r;
    zmin = m_cZ - m_lZ;
    zmax = m_cZ + m_lZ;
    return true;
  }

  const double dd = sqrt(m_r * m_r + m_lZ * m_lZ);
  xmin = m_cX - dd;
  xmax = m_cX + dd;
  ymin = m_cY - dd;
  ymax = m_cY + dd;
  zmin = m_cZ - dd;
  zmax = m_cZ + dd;
  return true;
}

void SolidWire::SetRadius(const double r) {
  if (r <= 0.) {
    std::cerr << "SolidWire::SetRadius: Radius must be > 0.\n";
    return;
  }
  m_r = r;
}

void SolidWire::SetHalfLength(const double lz) {
  if (lz <= 0.) {
    std::cerr << "SolidWire::SetHalfLength: Half-length must be > 0.\n";
    return;
  }
  m_lZ = lz;
}

bool SolidWire::SolidPanels(std::vector<Panel>& /*panels*/) {
  return true;
}

double SolidWire::GetDiscretisationLevel(const Panel& /*panel*/) {
  return m_dis;
} 

void SolidWire::Cut(
    const double /*x0*/, const double /*y0*/, const double /*z0*/,
    const double /*xn*/, const double /*yn*/, const double /*zn*/,
    std::vector<Panel>& /*panels*/) {

}

}
