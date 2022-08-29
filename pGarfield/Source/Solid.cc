#include <cmath>
#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Solid.hh"

namespace Garfield {

unsigned int Solid::s_id = 0;

void Solid::SetDirection(const double dx, const double dy, const double dz) {
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d < Small) {
    std::cerr << m_className << ": Direction vector has zero norm.\n";
    return;
  }
  m_dX = dx / d;
  m_dY = dy / d;
  m_dZ = dz / d;
  double phi, theta;
  const double dt = sqrt(m_dX * m_dX + m_dY * m_dY);
  if (dt < Small) {
    phi = 0.;
    if (m_dZ > 0.) {
      theta = 0.;
    } else {
      theta = Pi;
    }
  } else {
    phi = atan2(m_dY, m_dX);
    theta = atan2(dt, m_dZ);
  }
  m_cTheta = cos(theta);
  m_sTheta = sin(theta);
  m_cPhi = cos(phi);
  m_sPhi = sin(phi);
}

bool Solid::GetProfile(std::vector<double>& /*xv*/,
                       std::vector<double>& /*yv*/) const {

  std::cerr << m_className << "::GetProfile: function not implemented.\n";
  return false;
}
 
double Solid::NotImplemented(const std::string& fcn) const {
  std::cerr << m_className << "::" << fcn << ": function not implemented.\n";
  return 0.;
}

bool Solid::Intersect(const double x1, const double y1, const double z1,
                      const double x2, const double y2, const double z2,
                      const double x0, const double y0, const double z0,
                      const double a, const double b, const double c,
                      double& xc, double& yc, double& zc) {
  //-----------------------------------------------------------------------
  //   PLALIN - Cuts an arbitrary plane with a line.
  //   Variables : (X1,Y1,Z1) : starting point of the line
  //               (X2,Y2,Z2) : end point of the line
  //               (X0,Y0,Z0) : point on the plane
  //               (A,B,C)    : parameters of the plane
  //-----------------------------------------------------------------------

  xc = yc = zc = 0.;
  bool on = false;
  // Form the two products.
  const double prod1 = (x0 - x1) * a + (y0 - y1) * b + (z0 - z1) * c;
  const double prod2 = (x2 - x1) * a + (y2 - y1) * b + (z2 - z1) * c;
  // Set a tolerance for lambda.
  constexpr double eps = 1.e-5;
  // Check the products are non-zero.
  const double dx = x2 - x1;
  const double dy = y2 - y1;
  const double dz = z2 - z1;
  const double d2 = dx * dx + dy * dy + dz * dz;
  if (std::abs(prod2) > 1.e-6 * sqrt((a * a + b * b + c * c) * d2)) {
    double s = prod1 / prod2;
    if (s >= -eps && s <= 1. + eps) on = true;
    s = std::max(0., std::min(1., s));
    xc = x1 + s * dx;
    yc = y1 + s * dy;
    zc = z1 + s * dz;
  }
  return on;
}

}
