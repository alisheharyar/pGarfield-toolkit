#include <iostream>

#include "Garfield/ComponentUser.hh"

namespace Garfield {

ComponentUser::ComponentUser() : Component("User") {}

void ComponentUser::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, Medium*& m, int& status) {
  if (!m_efield) {
    ex = ey = ez = 0.;
    m = nullptr;
    status = -10;
    return;
  }

  m_efield(x, y, z, ex, ey, ez);
  m = GetMedium(x, y, z);
  if (!m) {
    if (m_debug) {
      std::cerr << m_className << "::ElectricField:\n    (" << x << ", " << y
                << ", " << z << ") is not inside a medium.\n";
    }
    status = -6;
    return;
  }

  status = m->IsDriftable() ? 0 : -5;
}

void ComponentUser::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, double& v, Medium*& m,
                                  int& status) {
  if (!m_efield) {
    ex = ey = ez = v = 0.;
    m = nullptr;
    status = -10;
    return;
  }
  m_efield(x, y, z, ex, ey, ez);

  if (m_potential) {
    m_potential(x, y, z, v);
  } else {
    v = 0.;
  }

  m = GetMedium(x, y, z);
  if (!m) {
    if (m_debug) {
      std::cerr << m_className << "::ElectricField:\n    (" << x << ", " << y
                << ", " << z << ") is not inside a medium.\n";
    }
    status = -6;
    return;
  }

  status = m->IsDriftable() ? 0 : -5;
}

bool ComponentUser::GetVoltageRange(double& vmin, double& vmax) {
  vmin = vmax = 0.;
  return false;
}

void ComponentUser::MagneticField(const double x, const double y,
                                  const double z, double& bx, double& by,
                                  double& bz, int& status) {
  if (!m_bfield) {
    Component::MagneticField(x, y, z, bx, by, bz, status);
    return;
  }
  m_bfield(x, y, z, bx, by, bz);
  status = 0;
}

void ComponentUser::WeightingField(const double x, const double y,
                                   const double z, double& wx, double& wy,
                                   double& wz, const std::string& label) {
  wx = wy = wz = 0.;
  if (!m_wfield) return;
  m_wfield(x, y, z, wx, wy, wz, label);
}

double ComponentUser::WeightingPotential(const double x, const double y,
                                         const double z,
                                         const std::string& label) {
  double v = 0.;
  if (m_wpot) m_wpot(x, y, z, v, label);
  return v;
}

void ComponentUser::DelayedWeightingField(const double x, const double y, 
                                          const double z, const double t,
                                          double& wx, double& wy, double& wz,
                                          const std::string& label) {
  wx = wy = wz = 0.;
  if (m_dwfield) m_dwfield(x, y, z, t, wx, wy, wz, label);
}

bool ComponentUser::GetBoundingBox(
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

bool ComponentUser::HasMagneticField() const {
  return m_bfield ? true : Component::HasMagneticField();
}

void ComponentUser::SetElectricField(
    std::function<void(const double, const double, const double, 
                       double&, double&, double&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetElectricField: Function is empty.\n";
    return;
  }
  m_efield = f;
  m_ready = true;
}

void ComponentUser::SetPotential(
    std::function<void(const double, const double, const double, 
                       double&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetPotential: Function is empty.\n";
    return;
  }
  m_potential = f;
}

void ComponentUser::SetWeightingField(
    std::function<void(const double, const double, const double, 
                       double&, double&, double&, const std::string&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetWeightingField: Function is empty.\n";
    return;
  }
  m_wfield = f;
}

void ComponentUser::SetWeightingPotential(
    std::function<void(const double, const double, const double, 
                       double&, const std::string&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetWeightingPotential: Function is empty.\n";
    return;
  }
  m_wpot = f;
}

void ComponentUser::SetDelayedWeightingField(
    std::function<void(const double, const double, const double, const double,
                       double&, double&, double&, const std::string&)> f) {

  if (!f) {
    std::cerr << m_className << "::SetDelayedWeightingField: Function is empty.\n";
    return;
  }
  m_dwfield = f;
}

void ComponentUser::SetMagneticField(
    std::function<void(const double, const double, const double, 
                       double&, double&, double&)> f) {
  if (!f) {
    std::cerr << m_className << "::SetMagneticField: Function is empty.\n";
    return;
  }
  m_bfield = f;
}

void ComponentUser::SetArea(
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

void ComponentUser::UnsetArea() {
  m_xmin.fill(0.);
  m_xmax.fill(0.);
  m_hasArea = false;
}

void ComponentUser::Reset() {
  m_efield = nullptr;
  m_potential = nullptr;
  m_wfield = nullptr;
  m_wpot = nullptr;
  m_dwfield = nullptr;
  m_bfield = nullptr;
  m_ready = false;
  UnsetArea();
  m_medium = nullptr;
}

void ComponentUser::UpdatePeriodicity() {
  if (m_debug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}
}
