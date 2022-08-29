#include <cmath>
#include <iostream>

#include "Garfield/Component.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"

namespace Garfield {

Component::Component(const std::string& name) {
  m_className = "Component" + name;
}

void Component::SetGeometry(Geometry* geo) {
  // Make sure the geometry is defined
  if (!geo) {
    std::cerr << m_className << "::SetGeometry: Null pointer.\n";
    return;
  }
  m_geometry = geo;
}

Medium* Component::GetMedium(const double x, const double y, const double z) {
  if (!m_geometry) return nullptr;
  return m_geometry->GetMedium(x, y, z);
}

void Component::Clear() {
  m_geometry = nullptr;
  m_ready = false;
  // Reset periodicities.
  m_periodic.fill(false);
  m_mirrorPeriodic.fill(false);
  m_axiallyPeriodic.fill(false);
  m_rotationSymmetric.fill(false);
  // Reset the magnetic field.
  m_b0.fill(0.);
  Reset();
}

void Component::WeightingField(const double /*x*/, const double /*y*/,
                               const double /*z*/, double& wx, double& wy,
                               double& wz, const std::string& /*label*/) {
  if (m_debug) {
    std::cerr << m_className << "::WeightingField: Function not implemented.\n";
  }
  wx = wy = wz = 0.;
}

void Component::DelayedWeightingField(const double /*x*/, const double /*y*/,
                                      const double /*z*/, const double /*t*/,
                                      double& wx, double& wy, double& wz,
                                      const std::string& /*label*/) {
  if (m_debug) {
    std::cerr << m_className << "::DelayedWeightingField: Not implemented.\n";
  }
  wx = wy = wz = 0.;
}

double Component::WeightingPotential(const double /*x*/, const double /*y*/,
                                     const double /*z*/,
                                     const std::string& /*label*/) {
  if (m_debug) {
    std::cerr << m_className << "::WeightingPotential: Not implemented.\n";
  }
  return 0.;
}

double Component::DelayedWeightingPotential(const double /*x*/,
                                            const double /*y*/,
                                            const double /*z*/,
                                            const double /*t*/,
                                            const std::string& /*label*/) {
  if (m_debug) {
    std::cerr << m_className 
              << "::DelayedWeightingPotential: Not implemented.\n";
  }
  return 0.;
}


void Component::MagneticField(const double x, const double y, const double z,
                              double& bx, double& by, double& bz, int& status) {
  bx = m_b0[0];
  by = m_b0[1];
  bz = m_b0[2];
  if (m_debug) {
    std::cout << m_className << "::MagneticField: Field at (" << x << ", " << y
              << ", " << z << ") is (" << bx << ", " << by << ", " << bz
              << ")\n";
  }
  status = 0;
}

void Component::SetMagneticField(const double bx, const double by,
                                 const double bz) {
  m_b0 = {bx, by, bz};
}

bool Component::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                               double& xmax, double& ymax, double& zmax) {
  if (!m_geometry) return false;
  return m_geometry->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
}

bool Component::GetElementaryCell(double& xmin, double& ymin, double& zmin,
                                  double& xmax, double& ymax, double& zmax) {
  return GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
}

bool Component::CrossedWire(
    const double /*x0*/, const double /*y0*/, const double /*z0*/, 
    const double /*x1*/, const double /*y1*/, const double /*z1*/,
    double& /*xc*/, double& /*yc*/, double& /*zc*/, 
    const bool /*centre*/, double& /*rc*/) {
  return false;
}

bool Component::InTrapRadius(const double /*q0*/, const double x0,
                             const double y0, const double /*z0*/, double& xw,
                             double& yw, double& rw) {
  xw = x0;
  yw = y0;
  rw = 0.;
  return false;
}

bool Component::CrossedPlane(
    const double /*x0*/, const double /*y0*/, const double /*z0*/, 
    const double /*x1*/, const double /*y1*/, const double /*z1*/,
    double& /*xc*/, double& /*yc*/, double& /*zc*/) {
  return false;
} 

bool Component::HasMagneticField() const {
  return fabs(m_b0[0]) > Small || fabs(m_b0[1]) > Small || 
         fabs(m_b0[2]) > Small;
}

double Component::IntegrateFluxCircle(const double xc, const double yc,
                                      const double r, const unsigned int nI) {
  // FLDIN2, FCHK3
  if (nI == 0) {
    std::cerr << m_className << "::IntegrateFluxCircle:\n"
              << "    Number of intervals must be > 0.\n";
    return 0.;
  }
  // Number of Gaussian quadrature points per interval.
  constexpr size_t nG = 6;
  auto tg = Numerics::GaussLegendreNodes6();
  auto wg = Numerics::GaussLegendreWeights6();

  // Width and half-width of intervals.
  const double d = TwoPi / nI;
  const double h = 0.5 * d;
  // Arguments of ElectricField.
  double ex = 0., ey = 0., ez = 0.;
  Medium* m = nullptr;
  int status = 0;
  // Perform the integration.
  double s = 0.;
  for (size_t i = 0; i < nG; ++i) {
    const double phi0 = h * (1. + tg[i]);
    for (unsigned int k = 0; k < nI; ++k) {
      const double phi = phi0 + k * d;
      const double cp = cos(phi);
      const double sp = sin(phi);
      ElectricField(xc + cp * r, yc + sp * r, 0., ex, ey, ez, m, status);
      s += wg[i] * r * (ex * cp + ey * sp);
    }
  }
  return h * s * VacuumPermittivity;
}

double Component::IntegrateFluxSphere(const double xc, const double yc,
                                      const double zc, const double r,
                                      const unsigned int nI) {
  // FLDIN3, FCHK2, FCHK1
  if (nI == 0) {
    std::cerr << m_className << "::IntegrateFluxSphere:\n"
              << "    Number of intervals must be > 0.\n";
    return 0.;
  }
  // Number of Gaussian quadrature points.
  constexpr size_t nG = 6;
  auto tg = Numerics::GaussLegendreNodes6();
  auto wg = Numerics::GaussLegendreWeights6();

  const double r2 = r * r;
  // Width and half-width of theta intervals.
  const double dt = Pi / nI;
  const double ht = 0.5 * dt;
  // Width and half-width of phi intervals.
  const double dp = TwoPi / nI;
  const double hp = 0.5 * dp;
  // Arguments of ElectricField.
  double ex = 0., ey = 0., ez = 0.;
  Medium* m = nullptr;
  int status = 0;
  // Perform the integration.
  double s2 = 0.;
  // Loop over theta.
  for (size_t i = 0; i < nG; ++i) {
    const double theta0 = ht * (1. + tg[i]) - HalfPi;
    for (unsigned int k = 0; k < nI; ++k) {
      const double theta = theta0 + k * dt;
      const double ct = cos(theta);
      const double st = sin(theta);
      const double z = zc + st * r;
      double s1 = 0.;
      // Loop over phi.
      for (size_t ii = 0; ii < nG; ++ii) {
        const double phi0 = hp * (1. + tg[ii]);
        for (unsigned int kk = 0; kk < nI; ++kk) {
          const double phi = phi0 + kk * dp;
          const double cp = cos(phi);
          const double sp = sin(phi);
          const double x = xc + cp * ct * r;
          const double y = yc + sp * ct * r;
          ElectricField(x, y, z, ex, ey, ez, m, status);
          s1 += wg[ii] * ((ex * cp + ey * sp) * ct + ez * st);
        }
      }
      s2 += wg[i] * r2 * ct * hp * s1;
    }
  }
  return ht * s2 * VacuumPermittivity;
}

double Component::IntegrateFluxParallelogram(
    const double x0, const double y0, const double z0, const double dx1,
    const double dy1, const double dz1, const double dx2, const double dy2,
    const double dz2, const unsigned int nU, const unsigned int nV) {

  return IntegrateFluxParallelogram(
      x0, y0, z0, dx1, dy1, dz1, dx2, dy2, dz2, nU, nV, false, "");
}

double Component::IntegrateWeightingFluxParallelogram(const std::string& id,
    const double x0, const double y0, const double z0, const double dx1,
    const double dy1, const double dz1, const double dx2, const double dy2,
    const double dz2, const unsigned int nU, const unsigned int nV) {

  return IntegrateFluxParallelogram(
      x0, y0, z0, dx1, dy1, dz1, dx2, dy2, dz2, nU, nV, true, id);
}

double Component::IntegrateFluxParallelogram(
    const double x0, const double y0, const double z0, const double dx1,
    const double dy1, const double dz1, const double dx2, const double dy2,
    const double dz2, const unsigned int nU, const unsigned int nV,
    const bool wfield, const std::string& label) {

  // FLDIN4, FCHK4, FCHK5
  if (nU <= 1 || nV <= 1) {
    std::cerr << m_className << "::IntegrateFluxParallelogram:\n"
              << "    Number of points to integrate over must be > 1.\n";
    return 0.;
  }
  // Number of Gaussian quadrature points.
  constexpr size_t nG = 6;
  auto tg = Numerics::GaussLegendreNodes6();
  auto wg = Numerics::GaussLegendreWeights6();

  // Compute the normal vector.
  const double xn = dy1 * dz2 - dz1 * dy2;
  const double yn = dz1 * dx2 - dx1 * dz2;
  const double zn = dx1 * dy2 - dy1 * dx2;
  if (m_debug) {
    std::cout << m_className << "::IntegrateFluxParallelogram:\n"
              << "    Normal vector = " << xn << ", " << yn << ", " << zn 
              << ".\n";
  }
  // If this vector has zero norm, return 0 flux.
  const double d1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
  const double d2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
  if (xn * xn + yn * yn + zn * zn < 1.e-10 * sqrt(d1 * d2) ||
      d1 < 1.e-10 * d2 || d2 < 1.e-10 * d1) {
    std::cerr << m_className << "::IntegrateFluxParallelogram:\n"
              << "    Parallelogram does not have non-zero area.\n";
    return 0.;
  }

  // (Half-)step sizes in the two directions.
  const double du = 1. / nU;
  const double hu = 0.5 * du;
  const double dv = 1. / nV;
  const double hv = 0.5 * dv;
  // Arguments of ElectricField.
  double fx = 0., fy = 0., fz = 0.;
  Medium* m = nullptr;
  int status = 0;
  // Perform the integration.
  double s2 = 0.;
  for (size_t i = 0; i < nG; ++i) {
    const double v0 = hv * (1. + tg[i]);
    for (unsigned int k = 0; k < nV; ++k) {
      const double v = v0 + k * dv;
      double s1 = 0.;
      for (size_t ii = 0; ii < nG; ++ii) {
        const double u0 = hu * (1. + tg[ii]);
        for (unsigned int kk = 0; kk < nU; ++kk) {
          const double u = u0 + kk * du;
          const double x = x0 + u * dx1 + v * dx2;
          const double y = y0 + u * dy1 + v * dy2;
          const double z = z0 + u * dz1 + v * dz2;
          if (wfield) {
            WeightingField(x, y, z, fx, fy, fz, label);
          } else {
            ElectricField(x, y, z, fx, fy, fz, m, status);
          }
          s1 += wg[ii] * (fx * xn + fy * yn + fz * zn);
        }
      }
      s2 += wg[i] * hu * s1;
    }
  }
  return hv * s2;
}

double Component::IntegrateFluxLine(const double x0, const double y0,
                                    const double z0, const double x1,
                                    const double y1, const double z1,
                                    const double xp, const double yp,
                                    const double zp, const unsigned int nI,
                                    const int isign) {
  // FLDIN5
  // Normalise the norm vector.
  const double pmag2 = xp * xp + yp * yp + zp * zp;
  if (pmag2 <= 0.) {
    std::cerr << m_className << "::IntegrateFluxLine:\n"
              << "    Normal vector has zero length; flux set to 0.\n";
    return 0.;
  }
  const double pmag = sqrt(pmag2);
  const double xn = xp / pmag;
  const double yn = yp / pmag;
  const double zn = zp / pmag;

  // Check integration points.
  if (nI <= 1) {
    std::cerr << m_className << "::IntegrateFluxLine:\n"
              << "    Number of points to integrate over must be > 1.\n";
    return 0.;
  }
  // Ensure the segment has non-zero length.
  const double vx = x1 - x0;
  const double vy = y1 - y0;
  const double vz = z1 - z0;
  const double vmag2 = vx * vx + vy * vy + vz * vz;
  if (vmag2 <= 0.) {
    std::cerr << m_className << "::IntegrateFluxLine:\n"
              << "    Segment has zero length; flux set to 0.\n";
    return 0.;
  }
  const double vmag = sqrt(vmag2);
  // Segment should be perpendicular to the norm vector.
  if (fabs(vx * xn + vy * yn + vz * zn) > 1.e-4 * vmag) {
    std::cerr << m_className << "::IntegrateFluxLine:\n"
              << "    Segment is not perpendicular to norm vector.\n";
    return 0.;
  }

  // Perform the integration.
  constexpr size_t nG = 6;
  auto tg = Numerics::GaussLegendreNodes6();
  auto wg = Numerics::GaussLegendreWeights6();

  const double d = 1. / nI;
  const double h = 0.5 * d;
  // Arguments of ElectricField.
  double ex = 0., ey = 0., ez = 0.;
  Medium* m = nullptr;
  int status = 0;
  double s = 0.;
  for (size_t i = 0; i < nG; ++i) {
    const double u0 = h * (1. + tg[i]);
    for (unsigned int k = 0; k < nI; ++k) {
      const double u = u0 + k * d;
      const double x = x0 + u * vx;
      const double y = y0 + u * vy;
      const double z = z0 + u * vz;
      ElectricField(x, y, z, ex, ey, ez, m, status);
      double fn = ex * xn + ey * yn + ez * zn;
      if (isign != 0) {
        // TODO: -1?
        fn = isign * fn > 0 ? fabs(fn) : -1.;
      }
      s += wg[i] * fn;
    }
  }
  return s * vmag;
}

}  // namespace Garfield
