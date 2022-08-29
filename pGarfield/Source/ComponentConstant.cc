#include <iostream>
#include <numeric>

#include "Garfield/ComponentConstant.hh"
#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

ComponentConstant::ComponentConstant() : Component("Constant") {}

Medium* ComponentConstant::GetMedium(const double x, const double y, 
                                     const double z) {

  if (!m_hasArea) return Component::GetMedium(x, y, z);
  return InArea(x, y, z) ? m_medium : nullptr; 
}

void ComponentConstant::ElectricField(const double x, const double y,
                                      const double z, double& ex, double& ey,
                                      double& ez, Medium*& m, int& status) {
  ex = m_efield[0];
  ey = m_efield[1];
  ez = m_efield[2];
  m = GetMedium(x, y, z);
  if (!m) {
    // No medium at this point.
    status = -6;
    return;
  }

  if (m->IsDriftable()) {
    status = 0;
  } else {
    status = -5;
  }
}

void ComponentConstant::ElectricField(const double x, const double y,
                                      const double z, double& ex, double& ey,
                                      double& ez, double& v, Medium*& m,
                                      int& status) {
  ex = m_efield[0];
  ey = m_efield[1];
  ez = m_efield[2];
  if (m_hasPotential) {
    // Compute the potential at this point.
    const std::array<double, 3> d = {x - m_x0, y - m_y0, z - m_z0};
    v = m_v0 - std::inner_product(d.begin(), d.end(), m_efield.begin(), 0.);
  } else {
    v = 0.;
    if (m_debug) {
      std::cerr << m_className << "::ElectricField: Potential not defined.\n";
    }
  }
  m = GetMedium(x, y, z);
  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at ("
                << x << ", " << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (m->IsDriftable()) {
    status = 0;
  } else {
    status = -5;
  }
}

bool ComponentConstant::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_hasPotential) return false;

  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  if (!GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cerr << m_className << "::GetVoltageRange:\n"
              << "    Could not determine the bounding box.\n";
    return false;
  }
  // Calculate potentials at each corner
  const double pxmin = m_v0 - (xmin - m_x0) * m_efield[0];
  const double pxmax = m_v0 - (xmax - m_x0) * m_efield[0];
  const double pymin = m_v0 - (ymin - m_y0) * m_efield[1];
  const double pymax = m_v0 - (ymax - m_y0) * m_efield[1];
  const double pzmin = m_v0 - (zmin - m_z0) * m_efield[2];
  const double pzmax = m_v0 - (zmax - m_z0) * m_efield[2];
  double p[8];
  p[0] = pxmin + pymin + pzmin;
  p[1] = pxmin + pymin + pzmax;
  p[2] = pxmin + pymax + pzmin;
  p[3] = pxmin + pymax + pzmax;
  p[4] = pxmax + pymin + pzmin;
  p[5] = pxmax + pymin + pzmax;
  p[6] = pxmax + pymax + pzmin;
  p[7] = pxmax + pymax + pzmax;
  vmin = vmax = p[7];
  for (int i = 7; i--;) {
    if (p[i] > vmax) vmax = p[i];
    if (p[i] < vmin) vmin = p[i];
  }

  return true;
}

bool ComponentConstant::GetBoundingBox(
    double& xmin, double& ymin, double& zmin,
    double& xmax, double& ymax, double& zmax) {

  if (!m_hasArea) {
    return Component::GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  xmin = m_xmin[0];
  ymin = m_xmin[1];
  zmin = m_xmin[2];
  xmax = m_xmax[0];
  ymax = m_xmax[1];
  zmax = m_xmax[2];
  return true;
}

void ComponentConstant::WeightingField(const double x, const double y,
                                       const double z, double& wx, double& wy,
                                       double& wz, const std::string& label) {
  if (!m_hasWeightingField || label != m_label) return;

  Medium* m = GetMedium(x, y, z);
  if (!m) {
    wx = wy = wz = 0.;
    if (m_debug) {
      std::cout << m_className << "::WeightingField: No medium at ("
                << x << ", " << y << ", " << z << ")\n";
    }
    return;
  }
  wx = m_wfield[0];
  wy = m_wfield[1];
  wz = m_wfield[2];
}

double ComponentConstant::WeightingPotential(const double x, const double y,
                                             const double z,
                                             const std::string& label) {
  if (!m_hasWeightingPotential || label != m_label) return 0.;
  // Make sure we are in the active area.
  if (!GetMedium(x, y, z)) return 0.;
  // Compute the potential.
  const std::array<double, 3> d = {x - m_wx0, y - m_wy0, z - m_wz0};
  return m_w0 - std::inner_product(d.begin(), d.end(), m_wfield.begin(), 0.);
}

void ComponentConstant::SetElectricField(const double ex, const double ey,
                                         const double ez) {
  m_efield = {ex, ey, ez};
  if (ex * ex + ey * ey + ez * ez < Small) {
    std::cerr << m_className << "::SetElectricField: Field set to zero.\n";
  }
  m_ready = true;
}

void ComponentConstant::SetPotential(const double x, const double y,
                                     const double z, const double v) {
  m_x0 = x;
  m_y0 = y;
  m_z0 = z;
  m_v0 = v;
  m_hasPotential = true;
}

void ComponentConstant::SetWeightingField(const double wx, const double wy,
                                          const double wz,
                                          const std::string label) {
  m_label = label;
  m_wfield = {wx, wy, wz};
  m_hasWeightingField = true;
}

void ComponentConstant::SetWeightingPotential(const double x, const double y,
                                              const double z, const double v) {
  if (!m_hasWeightingField) {
    std::cerr << m_className << "::SetWeightingPotential:\n"
              << "    Set the weighting field first!\n";
    return;
  }
  m_wx0 = x;
  m_wy0 = y;
  m_wz0 = z;
  m_w0 = v;
  m_hasWeightingPotential = true;
}

void ComponentConstant::SetArea(
    const double xmin, const double ymin, const double zmin,
    const double xmax, const double ymax, const double zmax) {

  m_xmin[0] = std::min(xmin, xmax);
  m_xmin[1] = std::min(ymin, ymax);
  m_xmin[2] = std::min(zmin, zmax);
  m_xmax[0] = std::max(xmin, xmax);
  m_xmax[1] = std::max(ymin, ymax);
  m_xmax[2] = std::max(zmin, zmax);
  m_hasArea = true; 
}

void ComponentConstant::UnsetArea() {
  m_xmin.fill(0.);
  m_xmax.fill(0.);
  m_hasArea = false;
}

void ComponentConstant::Reset() {
  m_efield.fill(0.);
  m_hasPotential = false;
  m_wfield.fill(0.);
  m_hasWeightingField = false;
  m_hasWeightingPotential = false;
  m_label = "";
  m_ready = false;
  UnsetArea();
  m_medium = nullptr;
}

void ComponentConstant::UpdatePeriodicity() {
  if (m_debug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}
}
