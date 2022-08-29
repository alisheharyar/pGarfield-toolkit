#include "Garfield/ComponentNeBem3d.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Polygon.hh"
#include "NR.h"
#include "neBEM.h"
#include "neBEMInterface.h"

namespace {

unsigned int NextPoint(const unsigned int i, const unsigned int n) {
  const unsigned int j = i + 1;
  return j < n ? j : 0;
}

unsigned int PrevPoint(const unsigned int i, const unsigned int n) {
  return i > 0 ? i - 1 : n - 1;
}

double Mag(const std::array<double, 3>& a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

std::array<double, 3> UnitVector(const std::array<double, 3>& a) {
  const double mag = Mag(a);
  if (mag < 1.e-12) return a;
  const std::array<double, 3> b = {a[0] / mag, a[1] / mag, a[2] / mag};
  return b;
}

std::array<double, 3> CrossProduct(const std::array<double, 3>& u,
                                   const std::array<double, 3>& v) {
  const std::array<double, 3> w = {u[1] * v[2] - u[2] * v[1],
                                   u[2] * v[0] - u[0] * v[2],
                                   u[0] * v[1] - u[1] * v[0]};
  return w;
}

std::array<double, 3> LocalToGlobal(
    const double x, const double y, const double z,
    const std::array<std::array<double, 3>, 3>& dcos,
    const std::array<double, 3>& t) {
  // Initial vector.
  std::array<double, 4> u = {x, y, z, 1.};

  std::array<std::array<double, 4>, 4> rot;
  rot[0] = {dcos[0][0], dcos[1][0], dcos[2][0], 0.};
  rot[1] = {dcos[0][1], dcos[1][1], dcos[2][1], 0.};
  rot[2] = {dcos[0][2], dcos[1][2], dcos[2][2], 0.};
  rot[3] = {0., 0., 0., 1.};

  std::array<double, 4> v = {0., 0., 0., 0.};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      v[i] += rot[i][j] * u[j];
    }
  }
  const std::array<double, 3> a = {v[0] + t[0], v[1] + t[1], v[2] + t[2]};
  return a;
}

/// Compute lambda for a point on a line (0 = start, 1 = end).
double Lambda(const double x1, const double x0, const double x2,
              const double y1, const double y0, const double y2) {
  // Segment of zero length.
  if ((x1 - x2) == 0. && (y1 - y2) == 0.) {
    std::cerr << "ComponentNeBem3d::Lambda: Zero length segment.\n";
    return 2.;
  }

  double xl = 0.;
  const double dx1 = x0 - x1;
  const double dy1 = y0 - y1;
  const double dx2 = x0 - x2;
  const double dy2 = y0 - y2;
  if (dx1 * dx1 + dy1 * dy1 < dx2 * dx2 + dy2 * dy2) {
    // Point nearer to (x1, y1).
    if (fabs(y1 - y2) > fabs(x1 - x2)) {
      xl = dy1 / (y2 - y1);
    } else {
      xl = dx1 / (x2 - x1);
    }
  } else {
    // Point nearer to (x2, y2).
    if (fabs(y1 - y2) > fabs(x1 - x2)) {
      xl = 1. - dy2 / (y1 - y2);
    } else {
      xl = 1. - dx2 / (x1 - x2);
    }
  }
  return xl;
}

/// Determine whether a point (u, v) lies on a straight line
/// (x1, y1) to (x2, y2).
bool OnLine(const double x1, const double y1, const double x2, const double y2,
            const double u, const double v) {
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v)});
  epsx = std::max(1.e-10, epsx);
  epsy = std::max(1.e-10, epsy);

  if ((fabs(x1 - u) <= epsx && fabs(y1 - v) <= epsy) ||
      (fabs(x2 - u) <= epsx && fabs(y2 - v) <= epsy)) {
    // Point to be examined coincides with start or end.
    return true;
  } else if (fabs(x1 - x2) <= epsx && fabs(y1 - y2) <= epsy) {
    // The line (x1, y1) to (x2, y2) is in fact a point.
    return false;
  }
  double xc = 0., yc = 0.;
  if (fabs(u - x1) + fabs(v - y1) < fabs(u - x2) + fabs(v - y2)) {
    // (u, v) is nearer to (x1, y1).
    const double dx = (x2 - x1);
    const double dy = (y2 - y1);
    const double xl = ((u - x1) * dx + (v - y1) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x1;
      yc = y1;
    } else if (xl > 1.) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + xl * dx;
      yc = y1 + xl * dy;
    }
  } else {
    // (u, v) is nearer to (x2, y2).
    const double dx = (x1 - x2);
    const double dy = (y1 - y2);
    const double xl = ((u - x2) * dx + (v - y2) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x2;
      yc = y2;
    } else if (xl > 1.) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + xl * dx;
      yc = y2 + xl * dy;
    }
  }
  // See whether the point is on the line.
  if (fabs(u - xc) < epsx && fabs(v - yc) < epsy) {
    return true;
  }
  return false;
}

/// Determine whether the 2 straight lines (x1, y1) to (x2, y2)
/// and (u1, v1) to (u2, v2) cross at an intermediate point for both lines.
bool Crossing(const double x1, const double y1, const double x2,
              const double y2, const double u1, const double v1,
              const double u2, const double v2, double& xc, double& yc) {
  /// Matrix to compute the crossing point.
  std::array<std::array<double, 2>, 2> a;
  a[0][0] = y2 - y1;
  a[0][1] = v2 - v1;
  a[1][0] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  // Initial values.
  xc = 0.;
  yc = 0.;
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u1), fabs(u2)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v1), fabs(v2)});
  epsx = std::max(epsx, 1.e-10);
  epsy = std::max(epsy, 1.e-10);
  // Check for a point of one line located on the other line.
  if (OnLine(x1, y1, x2, y2, u1, v1)) {
    xc = u1;
    yc = v1;
    return true;
  } else if (OnLine(x1, y1, x2, y2, u2, v2)) {
    xc = u2;
    yc = v2;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x1, y1)) {
    xc = x1;
    yc = y1;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x2, y2)) {
    xc = x2;
    yc = y2;
    return true;
  } else if (fabs(det) < epsx * epsy) {
    // Parallel, non-touching.
    return false;
  }
  // Crossing, non-trivial lines: solve crossing equations.
  const double aux = a[1][1];
  a[1][1] = a[0][0] / det;
  a[0][0] = aux / det;
  a[1][0] = -a[1][0] / det;
  a[0][1] = -a[0][1] / det;
  // Compute crossing point.
  xc = a[0][0] * (x1 * y2 - x2 * y1) + a[1][0] * (u1 * v2 - u2 * v1);
  yc = a[0][1] * (x1 * y2 - x2 * y1) + a[1][1] * (u1 * v2 - u2 * v1);
  // See whether the crossing point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc)) {
    // Intersecting lines.
    return true;
  }
  // Crossing point not on both lines.
  return false;
}

void AddPoints(const std::vector<double>& xp1, const std::vector<double>& yp1,
               const std::vector<double>& xp2, const std::vector<double>& yp2,
               std::vector<double>& xl, std::vector<double>& yl,
               std::vector<int>& flags, std::vector<double>& qs,
               const double epsx, const double epsy) {
  struct Point {
    double x;
    double y;
    int flag;
    double q;
  };

  std::vector<Point> points;

  const unsigned int np1 = xp1.size();
  const unsigned int np2 = xp2.size();
  for (unsigned int i = 0; i < np1; ++i) {
    const double xi0 = xp1[i];
    const double yi0 = yp1[i];
    const double xi1 = xp1[NextPoint(i, np1)];
    const double yi1 = yp1[NextPoint(i, np1)];
    // Add the vertex.
    Point p1;
    p1.x = xi0;
    p1.y = yi0;
    p1.flag = 1;
    p1.q = 0.;
    // If also on 2 or vertex of 2, flag it as crossing or foreign.
    for (unsigned int j = 0; j < np2; ++j) {
      const double xj0 = xp2[j];
      const double yj0 = yp2[j];
      if (fabs(xj0 - xi0) < epsx && fabs(yj0 - yi0) < epsy) {
        p1.flag = 2;
      }
      const double xj1 = xp2[NextPoint(j, np2)];
      const double yj1 = yp2[NextPoint(j, np2)];
      if (OnLine(xj0, yj0, xj1, yj1, xi0, yi0) &&
          (fabs(xj0 - xi0) > epsx || fabs(yj0 - yi0) > epsy) &&
          (fabs(xj1 - xi0) > epsx || fabs(yj1 - yi0) > epsy)) {
        p1.flag = 3;
      }
    }
    points.push_back(std::move(p1));
    // Go over the line segments of the other polygon.
    std::vector<Point> pointsOther;
    for (unsigned int j = 0; j < np2; ++j) {
      const double xj0 = xp2[j];
      const double yj0 = yp2[j];
      // Add vertices of 2 that are on this line.
      if (OnLine(xi0, yi0, xi1, yi1, xj0, yj0) &&
          (fabs(xi0 - xj0) > epsx || fabs(yi0 - yj0) > epsy) &&
          (fabs(xi1 - xj0) > epsx || fabs(yi1 - yj0) > epsy)) {
        Point p2;
        p2.x = xj0;
        p2.y = yj0;
        p2.flag = 2;
        pointsOther.push_back(std::move(p2));
      }
      const double xj1 = xp2[NextPoint(j, np2)];
      const double yj1 = yp2[NextPoint(j, np2)];
      // Add crossing points.
      double xc = 0., yc = 0.;
      bool add = Crossing(xi0, yi0, xi1, yi1, xj0, yj0, xj1, yj1, xc, yc);
      if (add) {
        if ((fabs(xi0 - xc) < epsx && fabs(yi0 - yc) < epsy) ||
            (fabs(xi1 - xc) < epsx && fabs(yi1 - yc) < epsy) ||
            (fabs(xj0 - xc) < epsx && fabs(yj0 - yc) < epsy) ||
            (fabs(xj1 - xc) < epsx && fabs(yj1 - yc) < epsy)) {
          add = false;
        }
        if ((fabs(xi0 - xj0) < epsx && fabs(yi0 - yj0) < epsy) ||
            (fabs(xi0 - xj1) < epsx && fabs(yi0 - yj1) < epsy) ||
            (fabs(xi1 - xj0) < epsx && fabs(yi1 - yj0) < epsy) ||
            (fabs(xi1 - xj1) < epsx && fabs(yi1 - yj1) < epsy)) {
          add = false;
        }
      }
      if (add) {
        Point p2;
        p2.x = xc;
        p2.y = yc;
        p2.flag = 3;
        pointsOther.push_back(std::move(p2));
      }
    }
    // Compute the lambdas for these points.
    for (auto& p : pointsOther) {
      p.q = Lambda(xi0, p.x, xi1, yi0, p.y, yi1);
    }
    // Sort the list by using the lambdas.
    std::sort(
        pointsOther.begin(), pointsOther.end(),
        [](const Point& lhs, const Point& rhs) { return (lhs.q < rhs.q); });
    points.insert(points.end(), pointsOther.begin(), pointsOther.end());
  }

  for (const Point& p : points) {
    xl.push_back(p.x);
    yl.push_back(p.y);
    flags.push_back(p.flag);
    qs.push_back(p.q);
  }
}

/// Determine whether 2 panels are equal.
bool Equal(const Garfield::Panel& panel1, const Garfield::Panel& panel2,
           const double epsx, const double epsy) {
  const auto& xp1 = panel1.xv;
  const auto& yp1 = panel1.yv;
  const auto& xp2 = panel2.xv;
  const auto& yp2 = panel2.yv;
  if (xp1.empty() || xp2.empty()) return false;
  const unsigned int np1 = xp1.size();
  const unsigned int np2 = xp2.size();

  // Compare all points of 1 with all points of 2.
  for (unsigned int i = 0; i < np1; ++i) {
    // Loop over 2 until a match is found.
    bool match = false;
    for (unsigned int j = 0; j < np2; ++j) {
      if (fabs(xp2[j] - xp1[i]) < epsx && fabs(yp2[j] - yp1[i]) < epsy) {
        match = true;
        break;
      }
      const unsigned int jj = NextPoint(j, np2);
      if (OnLine(xp2[j], yp2[j], xp2[jj], yp2[jj], xp1[i], yp1[i])) {
        match = true;
        break;
      }
    }
    if (!match) return false;
  }

  // Compare all points of 2 with all points of 1.
  for (unsigned int i = 0; i < np2; ++i) {
    // Loop over 1 until a match is found.
    bool match = false;
    for (unsigned int j = 0; j < np1; ++j) {
      if (fabs(xp2[i] - xp1[j]) < epsx && fabs(yp2[i] - yp1[j]) < epsy) {
        match = true;
        break;
      }
      const unsigned int jj = NextPoint(j, np1);
      if (OnLine(xp1[j], yp1[j], xp1[jj], yp1[jj], xp2[i], yp2[i])) {
        match = true;
        break;
      }
    }
    if (!match) return false;
  }

  // If we get this far, the curves are the same.
  return true;
}

}  // namespace

namespace Garfield {

ComponentNeBem3d* gComponentNeBem3d = nullptr;

ComponentNeBem3d::ComponentNeBem3d() : Component("NeBem3d") {
  InitValues();
}

Medium* ComponentNeBem3d::GetMedium(const double x, const double y,
                                    const double z) {
  if (!m_geometry) return nullptr;
  return m_geometry->GetMedium(x, y, z, true);
}

void ComponentNeBem3d::InitValues() {

  // Set initial values of weighting field fast volume related parameters
  // Since MAXWtFld is not a variable, we do not need neBEM::
  for(int id = 1; id < MAXWtFld; ++id) {
    neBEM::OptFixedWtField[id] = 0;
    neBEM::OptWtFldFastVol[id] = 0;
    m_optWtFldFastVol[id] = 0;
    m_optCreateWtFldFastPF[id] = 0;
    m_optReadWtFldFastPF[id] = 0;
  }
}

void ComponentNeBem3d::ElectricField(const double x, const double y,
                                     const double z, double& ex, double& ey,
                                     double& ez, double& v, Medium*& m,
                                     int& status) {
  ex = ey = ez = v = 0.;
  status = 0;
  // Check if the requested point is inside a medium
  m = GetMedium(x, y, z);
  if (!m) {
    status = -6;
  } else if (!m->IsDriftable()) {
    status = -5;
  }

  if (!m_ready) {
    if (!Initialise()) {
      std::cerr << m_className << "::ElectricField: Initialisation failed.\n";
      status = -11;
      return;
    }
    m_ready = true;
  }

  // Construct a point.
  neBEM::Point3D point;
  point.X = 0.01 * x;
  point.Y = 0.01 * y;
  point.Z = 0.01 * z;

  // Compute the field.
  neBEM::Vector3D field;
  if (neBEM::neBEMPF(&point, &v, &field) != 0) {
    status = -10;
    return;
  }
  ex = 0.01 * field.X;
  ey = 0.01 * field.Y;
  ez = 0.01 * field.Z;
}

void ComponentNeBem3d::ElectricField(const double x, const double y,
                                     const double z, double& ex, double& ey,
                                     double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

bool ComponentNeBem3d::GetVoltageRange(double& vmin, double& vmax) {
  // Voltage and other bc have to come from the solids
  vmin = vmax = 0;
  return true;
}

void ComponentNeBem3d::WeightingField(const double x, const double y,
                                      const double z, double& wx, double& wy,
                                      double& wz, const std::string& label) {
  wx = wy = wz = 0.;
  if (m_wfields.count(label) == 0) return;
  const int id = m_wfields[label];
  neBEM::Point3D point;
  point.X = 0.01 * x;
  point.Y = 0.01 * y;
  point.Z = 0.01 * z;
  neBEM::Vector3D field;
  if (neBEM::neBEMWeightingField(&point, &field, id) == DBL_MAX) {
    std::cerr << m_className << "::WeightingField: Evaluation failed.\n";
    return;
  }
  wx = 0.01 * field.X;
  wy = 0.01 * field.Y;
  wz = 0.01 * field.Z;
}

double ComponentNeBem3d::WeightingPotential(const double x, const double y,
                                            const double z,
                                            const std::string& label) {
  if (m_wfields.count(label) == 0) return 0.;
  const int id = m_wfields[label];
  neBEM::Point3D point;
  point.X = 0.01 * x;
  point.Y = 0.01 * y;
  point.Z = 0.01 * z;
  neBEM::Vector3D field;
  const double v = neBEM::neBEMWeightingField(&point, &field, id);
  return v == DBL_MAX ? 0. : v;
}

void ComponentNeBem3d::AddPlaneX(const double x, const double v) {
  if (m_ynplan[0] && m_ynplan[1]) {
    std::cerr << m_className << "::AddPlaneX:\n"
              << "    Cannot have more than two planes at constant x.\n";
    return;
  }

  if (m_ynplan[0]) {
    m_ynplan[1] = true;
    if (x < m_coplan[0]) {
      m_coplan[1] = m_coplan[0];
      m_vtplan[1] = m_vtplan[0];
      m_coplan[0] = x;
      m_vtplan[0] = v;
    } else {
      m_coplan[1] = x;
      m_vtplan[1] = v;
    }
  } else {
    m_ynplan[0] = true;
    m_coplan[0] = x;
    m_vtplan[0] = v;
  }
  m_ready = false;
}

void ComponentNeBem3d::AddPlaneY(const double y, const double v) {
  if (m_ynplan[2] && m_ynplan[3]) {
    std::cerr << m_className << "::AddPlaneY:\n"
              << "    Cannot have more than two planes at constant y.\n";
    return;
  }

  if (m_ynplan[2]) {
    m_ynplan[3] = true;
    if (y < m_coplan[2]) {
      m_coplan[3] = m_coplan[2];
      m_vtplan[3] = m_vtplan[2];
      m_coplan[2] = y;
      m_vtplan[2] = v;
    } else {
      m_coplan[3] = y;
      m_vtplan[3] = v;
    }
  } else {
    m_ynplan[2] = true;
    m_coplan[2] = y;
    m_vtplan[2] = v;
  }
  m_ready = false;
}

void ComponentNeBem3d::AddPlaneZ(const double z, const double v) {
  if (m_ynplan[4] && m_ynplan[5]) {
    std::cerr << m_className << "::AddPlaneZ:\n"
              << "    Cannot have more than two planes at constant z.\n";
    return;
  }

  if (m_ynplan[4]) {
    m_ynplan[5] = true;
    if (z < m_coplan[4]) {
      m_coplan[5] = m_coplan[4];
      m_vtplan[5] = m_vtplan[4];
      m_coplan[4] = z;
      m_vtplan[4] = v;
    } else {
      m_coplan[5] = z;
      m_vtplan[5] = v;
    }
  } else {
    m_ynplan[4] = true;
    m_coplan[4] = z;
    m_vtplan[4] = v;
  }
  m_ready = false;
}

unsigned int ComponentNeBem3d::GetNumberOfPlanesX() const {
  if (m_ynplan[0] && m_ynplan[1]) {
    return 2;
  } else if (m_ynplan[0] || m_ynplan[1]) {
    return 1;
  }
  return 0;
}

unsigned int ComponentNeBem3d::GetNumberOfPlanesY() const {
  if (m_ynplan[2] && m_ynplan[3]) {
    return 2;
  } else if (m_ynplan[2] || m_ynplan[3]) {
    return 1;
  }
  return 0;
}

unsigned int ComponentNeBem3d::GetNumberOfPlanesZ() const {
  if (m_ynplan[4] && m_ynplan[5]) {
    return 2;
  } else if (m_ynplan[4] || m_ynplan[5]) {
    return 1;
  }
  return 0;
}

bool ComponentNeBem3d::GetPlaneX(const unsigned int i, double& x,
                                 double& v) const {
  if (i >= 2 || (i == 1 && !m_ynplan[1])) {
    std::cerr << m_className << "::GetPlaneX: Index out of range.\n";
    return false;
  }
  x = m_coplan[i];
  v = m_vtplan[i];
  return true;
}

bool ComponentNeBem3d::GetPlaneY(const unsigned int i, double& y,
                                 double& v) const {
  if (i >= 2 || (i == 1 && !m_ynplan[3])) {
    std::cerr << m_className << "::GetPlaneY: Index out of range.\n";
    return false;
  }
  y = m_coplan[i + 2];
  v = m_vtplan[i + 2];
  return true;
}

bool ComponentNeBem3d::GetPlaneZ(const unsigned int i, double& z,
                                 double& v) const {
  if (i >= 2 || (i == 1 && !m_ynplan[5])) {
    std::cerr << m_className << "::GetPlaneZ: Index out of range.\n";
    return false;
  }
  z = m_coplan[i + 4];
  v = m_vtplan[i + 4];
  return true;
}

void ComponentNeBem3d::SetTargetElementSize(const double length) {
  if (length < MinDist) {
    std::cerr << m_className << "::SetTargetElementSize: Value must be > "
              << MinDist << ".\n";
    return;
  }
  m_targetElementSize = length;
}

void ComponentNeBem3d::SetMinMaxNumberOfElements(const unsigned int nmin,
                                                 const unsigned int nmax) {
  if (nmin == 0 || nmax == 0) {
    std::cerr << m_className << "::SetMinMaxNumberOfElements:\n"
              << "    Values must be non-zero.\n";
    return;
  }
  m_minNbElementsOnLength = std::min(nmin, nmax);
  m_maxNbElementsOnLength = std::max(nmin, nmax);
}

void ComponentNeBem3d::SetNewModel(const unsigned int NewModel) {
  m_newModel = NewModel;
}

void ComponentNeBem3d::SetNewMesh(const unsigned int NewMesh) {
  m_newMesh = NewMesh;
}

void ComponentNeBem3d::SetNewBC(const unsigned int NewBC) { m_newBC = NewBC; }

void ComponentNeBem3d::SetNewPP(const unsigned int NewPP) { m_newPP = NewPP; }

void ComponentNeBem3d::SetModelOptions(const unsigned int NewModel,
                                       const unsigned int NewMesh,
                                       const unsigned int NewBC,
                                       const unsigned int NewPP) {
  m_newModel = NewModel;
  m_newMesh = NewMesh;
  m_newBC = NewBC;
  m_newPP = NewPP;
}

void ComponentNeBem3d::SetStoreInflMatrix(
    const unsigned int OptStoreInflMatrix) {
  m_optStoreInflMatrix = OptStoreInflMatrix;
}

void ComponentNeBem3d::SetReadInflMatrix(const unsigned int OptReadInflMatrix) {
  m_optReadInflMatrix = OptReadInflMatrix;
}

void ComponentNeBem3d::SetStoreInvMatrix(const unsigned int OptStoreInvMatrix) {
  m_optStoreInvMatrix = OptStoreInvMatrix;
}

void ComponentNeBem3d::SetReadInvMatrix(const unsigned int OptReadInvMatrix) {
  m_optReadInvMatrix = OptReadInvMatrix;
}

void ComponentNeBem3d::SetStorePrimitives(
    const unsigned int OptStorePrimitives) {
  m_optStorePrimitives = OptStorePrimitives;
}

void ComponentNeBem3d::SetReadPrimitives(const unsigned int OptReadPrimitives) {
  m_optReadPrimitives = OptReadPrimitives;
}

void ComponentNeBem3d::SetStoreElements(const unsigned int OptStoreElements) {
  m_optStoreElements = OptStoreElements;
}

void ComponentNeBem3d::SetReadElements(const unsigned int OptReadElements) {
  m_optReadElements = OptReadElements;
}

void ComponentNeBem3d::SetFormattedFile(const unsigned int OptFormattedFile) {
  m_optStoreFormatted = OptFormattedFile;
}

void ComponentNeBem3d::SetUnformattedFile(
    const unsigned int OptUnformattedFile) {
  m_optStoreUnformatted = OptUnformattedFile;
}

void ComponentNeBem3d::SetStoreReadOptions(
    const unsigned int OptStoreInflMatrix, const unsigned int OptReadInflMatrix,
    const unsigned int OptStoreInvMatrix, const unsigned int OptReadInvMatrix,
    const unsigned int OptStorePrimitives, const unsigned int OptReadPrimitives,
    const unsigned int OptStoreElements, const unsigned int OptReadElements,
    const unsigned int OptFormattedFile,
    const unsigned int OptUnformattedFile) {
  m_optStoreInflMatrix = OptStoreInflMatrix;
  m_optReadInflMatrix = OptReadInflMatrix;
  m_optStoreInvMatrix = OptStoreInvMatrix;
  m_optReadInvMatrix = OptReadInvMatrix;
  m_optStorePrimitives = OptStorePrimitives;
  m_optReadPrimitives = OptReadPrimitives;
  m_optStoreElements = OptStoreElements;
  m_optReadElements = OptReadElements;
  m_optStoreFormatted = OptFormattedFile;
  m_optStoreUnformatted = OptUnformattedFile;
}

void ComponentNeBem3d::SetReuseModel() {
  m_newModel = 0;
  m_newMesh = 0;
  m_newBC = 1;
  m_optReadInvMatrix = 1;
}

void ComponentNeBem3d::SetSystemChargeZero(
    const unsigned int OptSystemChargeZero) {
  m_optSystemChargeZero = OptSystemChargeZero;
}

void ComponentNeBem3d::SetValidateSolution(
    const unsigned int OptValidateSolution) {
  m_optValidateSolution = OptValidateSolution;
}

void ComponentNeBem3d::SetForceValidation(
    const unsigned int OptForceValidation) {
  m_optForceValidation = OptForceValidation;
}

void ComponentNeBem3d::SetRepeatLHMatrix(const unsigned int OptRepeatLHMatrix) {
  m_optRepeatLHMatrix = OptRepeatLHMatrix;
}

void ComponentNeBem3d::SetComputeOptions(const unsigned int OptSystemChargeZero,
                                         const unsigned int OptValidateSolution,
                                         const unsigned int OptForceVaildation,
                                         const unsigned int OptRepeatLHMatrix) {
  m_optSystemChargeZero = OptSystemChargeZero;
  m_optValidateSolution = OptValidateSolution;
  m_optForceValidation = OptForceVaildation;
  m_optRepeatLHMatrix = OptRepeatLHMatrix;
}

// Fast volume options (Physical potential and field)
void ComponentNeBem3d::SetFastVolOptions(const unsigned int OptFastVol,
                                         const unsigned int OptCreateFastPF,
                                         const unsigned int OptReadFastPF) {
  m_optFastVol = OptFastVol;
  m_optCreateFastPF = OptCreateFastPF;
  m_optReadFastPF = OptReadFastPF;
}

// Fast volume version
void ComponentNeBem3d::SetFastVolVersion(const unsigned int VersionFV) {
  m_versionFV = VersionFV;
}

// Fast volume blocks
void ComponentNeBem3d::SetFastVolBlocks(const unsigned int NbBlocksFV) {
  m_nbBlocksFV = NbBlocksFV;
}

// Needs to include IdWtField information for each of these WtFld functions
// Weighting potential and field related Fast volume options
void ComponentNeBem3d::SetWtFldFastVolOptions(
    const unsigned int IdWtField, const unsigned int OptWtFldFastVol,
    const unsigned int OptCreateWtFldFastPF,
    const unsigned int OptReadWtFldFastPF) {
  m_idWtField = IdWtField;
  m_optWtFldFastVol[IdWtField] = OptWtFldFastVol;
  m_optCreateWtFldFastPF[IdWtField] = OptCreateWtFldFastPF;
  m_optReadWtFldFastPF[IdWtField] = OptReadWtFldFastPF;
}

// Weighting field Fast volume version
void ComponentNeBem3d::SetWtFldFastVolVersion(
    const unsigned int IdWtField, const unsigned int VersionWtFldFV) {
  m_idWtField = IdWtField;
  m_versionWtFldFV[IdWtField] = VersionWtFldFV;
}

// Weighting field Fast volume blocks
void ComponentNeBem3d::SetWtFldFastVolBlocks(
    const unsigned int IdWtField, const unsigned int NbBlocksWtFldFV) {
  m_idWtField = IdWtField;
  m_nbBlocksWtFldFV[IdWtField] = NbBlocksWtFldFV;
}

// Known charge options
void ComponentNeBem3d::SetKnownChargeOptions(
    const unsigned int OptKnownCharge) {
  m_optKnownCharge = OptKnownCharge;
}

// Charging up options
void ComponentNeBem3d::SetChargingUpOptions(const unsigned int OptChargingUp) {
  m_optChargingUp = OptChargingUp;
}

void ComponentNeBem3d::SetPeriodicCopies(const unsigned int nx,
                                         const unsigned int ny,
                                         const unsigned int nz) {
  m_nCopiesX = nx;
  m_nCopiesY = ny;
  m_nCopiesZ = nz;
}

void ComponentNeBem3d::SetPeriodicityX(const double s) {
  if (s < Small) {
    std::cerr << m_className << "::SetPeriodicityX:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }
  m_periodicLength[0] = s;
  m_periodic[0] = true;
  m_mirrorPeriodic[0] = false;
  UpdatePeriodicity();
}

void ComponentNeBem3d::SetPeriodicityY(const double s) {
  if (s < Small) {
    std::cerr << m_className << "::SetPeriodicityY:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }
  m_periodicLength[1] = s;
  m_periodic[1] = true;
  m_mirrorPeriodic[1] = false;
  UpdatePeriodicity();
}

void ComponentNeBem3d::SetPeriodicityZ(const double s) {
  if (s < Small) {
    std::cerr << m_className << "::SetPeriodicityZ:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }
  m_periodicLength[2] = s;
  m_periodic[2] = true;
  m_mirrorPeriodic[2] = false;
  UpdatePeriodicity();
}

void ComponentNeBem3d::SetMirrorPeriodicityX(const double s) {
  if (s < Small) {
    std::cerr << m_className << "::SetMirrorPeriodicityX:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }
  m_periodicLength[0] = s;
  m_periodic[0] = false;
  m_mirrorPeriodic[0] = true;
  UpdatePeriodicity();
}

void ComponentNeBem3d::SetMirrorPeriodicityY(const double s) {
  if (s < Small) {
    std::cerr << m_className << "::SetMirrorPeriodicityY:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }
  m_periodicLength[1] = s;
  m_periodic[1] = false;
  m_mirrorPeriodic[1] = true;
  UpdatePeriodicity();
}

void ComponentNeBem3d::SetMirrorPeriodicityZ(const double s) {
  if (s < Small) {
    std::cerr << m_className << "::SetMirrorPeriodicityZ:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }
  m_periodicLength[2] = s;
  m_periodic[2] = false;
  m_mirrorPeriodic[2] = true;
  UpdatePeriodicity();
}

bool ComponentNeBem3d::GetPeriodicityX(double& s) const {
  if (!m_periodic[0] && !m_mirrorPeriodic[0]) {
    s = 0.;
    return false;
  }
  s = m_periodicLength[0];
  return true;
}

bool ComponentNeBem3d::GetPeriodicityY(double& s) const {
  if (!m_periodic[1] && !m_mirrorPeriodic[1]) {
    s = 0.;
    return false;
  }
  s = m_periodicLength[1];
  return true;
}

bool ComponentNeBem3d::GetPeriodicityZ(double& s) const {
  if (!m_periodic[2] && !m_mirrorPeriodic[2]) {
    s = 0.;
    return false;
  }
  s = m_periodicLength[2];
  return true;
}

bool ComponentNeBem3d::Initialise() {
  // Reset the lists.
  m_primitives.clear();
  m_elements.clear();

  if (!m_geometry) {
    std::cerr << m_className << "::Initialise: Geometry not set.\n";
    return false;
  }
  // Be sure we won't have intersections with the bounding box.
  // TODO! Loop over the solids and call PLACYE, PLACHE, PLABXE, PLASPE, ...
  // Loop over the solids.
  std::map<int, Solid*> solids;
  std::map<int, Solid::BoundaryCondition> bc;
  std::map<int, double> volt;
  std::map<int, double> eps;
  std::map<int, double> charge;
  const unsigned int nSolids = m_geometry->GetNumberOfSolids();
  std::vector<Panel> panelsIn;
  for (unsigned int i = 0; i < nSolids; ++i) {
    Medium* medium = nullptr;
    const auto solid = m_geometry->GetSolid(i, medium);
    if (!solid) continue;
    // Get the panels.
    solid->SolidPanels(panelsIn);
    // Get the boundary condition.
    const auto id = solid->GetId();
    solids[id] = solid;
    bc[id] = solid->GetBoundaryConditionType();
    volt[id] = solid->GetBoundaryPotential();
    if (bc[id] == Solid::Unknown) {
      std::cout << m_className << "::Initialise:\n"
                << "    Boundary conditions for solid " << id << " not set.\n";
      if (medium && medium->IsConductor()) {
        std::cout << "    Assuming the panels to be grounded.\n";
        bc[id] = Solid::Voltage;
        volt[id] = 0.;
      } else {
        std::cout << "    Assuming dielectric-dielectric interfaces.\n";
        bc[id] = Solid::Dielectric;
      }
    }
    charge[id] = solid->GetBoundaryChargeDensity();
    if (!medium) {
      eps[id] = 1.;
    } else {
      eps[id] = medium->GetDielectricConstant();
    }
  }
  // Apply cuts.
  // CALL CELSCT('APPLY')
  // Reduce to basic periodic copy.
  ShiftPanels(panelsIn);

  // Find contact panels and split into primitives.

  // *---------------------------------------------------------------------
  // * PLABEM - Prepares panels for BEM applications: removes the contacts
  // *          and cuts polygons to rectangles and right-angle triangles.
  // *---------------------------------------------------------------------

  // Establish tolerances.
  const double epsang = 1.e-6;  // BEMEPA
  const double epsxyz = 1.e-6;  // BEMEPD
  // CALL EPSSET('SET',EPSXYZ,EPSXYZ,EPSXYZ)

  const unsigned int nIn = panelsIn.size();
  if (m_debug) {
    std::cout << m_className << "::Initialise: Retrieved " << nIn
              << " panels from the solids.\n";
  }
  // Keep track of which panels have been processed.
  std::vector<bool> mark(nIn, false);
  // Count the number of interface panels that have been discarded.
  unsigned int nTrivial = 0;
  unsigned int nConflicting = 0;
  unsigned int nNotImplemented = 0;
  // Pick up panels which coincide potentially.
  for (unsigned int i = 0; i < nIn; ++i) {
    // Skip panels already done.
    if (mark[i]) continue;
    // Fetch panel parameters.
    const double a1 = panelsIn[i].a;
    const double b1 = panelsIn[i].b;
    const double c1 = panelsIn[i].c;
    const auto& xp1 = panelsIn[i].xv;
    const auto& yp1 = panelsIn[i].yv;
    const auto& zp1 = panelsIn[i].zv;
    const unsigned int np1 = xp1.size();
    // Establish its norm and offset.
    const double d1 = a1 * xp1[0] + b1 * yp1[0] + c1 * zp1[0];
    if (m_debug) {
      std::cout << "  Panel " << i << "\n    Norm vector: " << a1 << ", " << b1
                << ", " << c1 << ", " << d1 << ".\n";
    }
    // Rotation matrix.
    std::array<std::array<double, 3>, 3> rot;
    if (fabs(c1) <= fabs(a1) && fabs(c1) <= fabs(b1)) {
      // Rotation: removing C
      rot[0][0] = b1 / sqrt(a1 * a1 + b1 * b1);
      rot[0][1] = -a1 / sqrt(a1 * a1 + b1 * b1);
      rot[0][2] = 0.;
    } else if (fabs(b1) <= fabs(a1) && fabs(b1) <= fabs(c1)) {
      // Rotation: removing B
      rot[0][0] = c1 / sqrt(a1 * a1 + c1 * c1);
      rot[0][1] = 0.;
      rot[0][2] = -a1 / sqrt(a1 * a1 + c1 * c1);
    } else {
      // Rotation: removing A
      rot[0][0] = 0.;
      rot[0][1] = c1 / sqrt(b1 * b1 + c1 * c1);
      rot[0][2] = -b1 / sqrt(b1 * b1 + c1 * c1);
    }
    rot[2][0] = a1;
    rot[2][1] = b1;
    rot[2][2] = c1;
    rot[1][0] = rot[2][1] * rot[0][2] - rot[2][2] * rot[0][1];
    rot[1][1] = rot[2][2] * rot[0][0] - rot[2][0] * rot[0][2];
    rot[1][2] = rot[2][0] * rot[0][1] - rot[2][1] * rot[0][0];
    // Rotate to the x, y plane.
    std::vector<double> xp(np1, 0.);
    std::vector<double> yp(np1, 0.);
    std::vector<double> zp(np1, 0.);
    for (unsigned int k = 0; k < np1; ++k) {
      xp[k] = rot[0][0] * xp1[k] + rot[0][1] * yp1[k] + rot[0][2] * zp1[k];
      yp[k] = rot[1][0] * xp1[k] + rot[1][1] * yp1[k] + rot[1][2] * zp1[k];
      zp[k] = rot[2][0] * xp1[k] + rot[2][1] * yp1[k] + rot[2][2] * zp1[k];
    }
    // Store it.
    std::vector<Panel> newPanels;
    std::vector<int> vol1;
    std::vector<int> vol2;
    Panel panel1 = panelsIn[i];
    panel1.xv = xp;
    panel1.yv = yp;
    panel1.zv = zp;
    vol1.push_back(panel1.volume);
    vol2.push_back(-1);
    newPanels.push_back(std::move(panel1));
    // Pick up all matching planes.
    for (unsigned int j = i + 1; j < nIn; ++j) {
      if (mark[j]) continue;
      const double a2 = panelsIn[j].a;
      const double b2 = panelsIn[j].b;
      const double c2 = panelsIn[j].c;
      const auto& xp2 = panelsIn[j].xv;
      const auto& yp2 = panelsIn[j].yv;
      const auto& zp2 = panelsIn[j].zv;
      const unsigned int np2 = xp2.size();
      // See whether this matches the first.
      const double d2 = a2 * xp2[0] + b2 * yp2[0] + c2 * zp2[0];
      // Inner product.
      const double dot = a1 * a2 + b1 * b2 + c1 * c2;
      // Offset between the two planes.
      const double offset = d1 - d2 * dot;
      if (fabs(fabs(dot) - 1.) > epsang || fabs(offset) > epsxyz) continue;
      // Found a match.
      mark[j] = true;
      if (m_debug) std::cout << "    Match with panel " << j << ".\n";
      // Rotate this plane too.
      xp.assign(np2, 0.);
      yp.assign(np2, 0.);
      zp.assign(np2, 0.);
      for (unsigned int k = 0; k < np2; ++k) {
        xp[k] = rot[0][0] * xp2[k] + rot[0][1] * yp2[k] + rot[0][2] * zp2[k];
        yp[k] = rot[1][0] * xp2[k] + rot[1][1] * yp2[k] + rot[1][2] * zp2[k];
        zp[k] = rot[2][0] * xp2[k] + rot[2][1] * yp2[k] + rot[2][2] * zp2[k];
      }
      // Store it.
      Panel panel2 = panelsIn[j];
      panel2.xv = xp;
      panel2.yv = yp;
      panel2.zv = zp;
      vol1.push_back(panel2.volume);
      vol2.push_back(-1);
      newPanels.push_back(std::move(panel2));
    }
    std::vector<bool> obsolete(newPanels.size(), false);
    // Cut them as long as needed till no contacts remain.
    unsigned int jmin = 0;
    bool change = true;
    while (change) {
      change = false;
      const unsigned int n = newPanels.size();
      for (unsigned int j = 0; j < n; ++j) {
        if (obsolete[j] || j < jmin) continue;
        if (vol1[j] >= 0 && vol2[j] >= 0) continue;
        const auto& panelj = newPanels[j];
        for (unsigned int k = j + 1; k < n; ++k) {
          if (obsolete[k]) continue;
          if (vol1[k] >= 0 && vol2[k] >= 0) continue;
          const auto& panelk = newPanels[k];
          if (m_debug) std::cout << "    Cutting " << j << ", " << k << ".\n";
          // Separate contact and non-contact areas.
          std::vector<Panel> panelsOut;
          std::vector<int> itypo;
          EliminateOverlaps(panelj, panelk, panelsOut, itypo);
          const unsigned int nOut = panelsOut.size();
          if (nOut == 2) {
            // TODO: retrieve epsx, epsy from overlap finding?
            const double epsx = epsxyz;
            const double epsy = epsxyz;
            // If there are just 2 panels, see whether there is a new one.
            const bool equal1 = Equal(panelj, panelsOut[0], epsx, epsy);
            const bool equal2 = Equal(panelj, panelsOut[1], epsx, epsy);
            const bool equal3 = Equal(panelk, panelsOut[0], epsx, epsy);
            const bool equal4 = Equal(panelk, panelsOut[1], epsx, epsy);
            if ((equal1 || equal3) && (equal2 || equal4)) {
              if (m_debug) {
                std::cout << "    Original and new panels are identical.\n";
              }
            } else {
              change = true;
            }
          } else {
            change = true;
          }
          if (m_debug) std::cout << "    Change flag: " << change << ".\n";
          // If there is no change, keep the two panels and proceed.
          if (!change) continue;
          // Flag the existing panels as inactive.
          obsolete[j] = true;
          obsolete[k] = true;

          // Add the new panels.
          for (unsigned int l = 0; l < nOut; ++l) {
            if (itypo[l] == 1) {
              vol1.push_back(std::max(vol1[j], vol2[j]));
              vol2.push_back(-1);
            } else if (itypo[l] == 2) {
              vol1.push_back(std::max(vol1[k], vol2[k]));
              vol2.push_back(-1);
            } else {
              vol1.push_back(std::max(vol1[j], vol2[j]));
              vol2.push_back(std::max(vol1[k], vol2[k]));
            }
            newPanels.push_back(std::move(panelsOut[l]));
            obsolete.push_back(false);
          }
          jmin = j + 1;
          // Restart the loops.
          break;
        }
        if (change) break;
      }
    }
    // And rotate the panels back in place.
    const unsigned int nNew = newPanels.size();
    for (unsigned int j = 0; j < nNew; ++j) {
      if (obsolete[j]) continue;
      // Examine the boundary conditions.
      int interfaceType = 0;
      double potential = 0.;
      double lambda = 0.;
      double chargeDensity = 0.;
      if (m_debug) {
        std::cout << "    Volume 1: " << vol1[j] << ". Volume 2: " << vol2[j]
                  << ".\n";
      }
      if (vol1[j] < 0 && vol2[j] < 0) {
        // Shouldn't happen.
        continue;
      } else if (vol1[j] < 0 || vol2[j] < 0) {
        // Interface between a solid and vacuum/background.
        const auto vol = vol1[j] < 0 ? vol2[j] : vol1[j];
        interfaceType = InterfaceType(bc[vol]);
        if (bc[vol] == Solid::Dielectric ||
            bc[vol] == Solid::DielectricCharge) {
          if (fabs(eps[vol] - 1.) < 1.e-6) {
            // Same epsilon on both sides. Skip.
            interfaceType = 0;
          } else {
            lambda = (eps[vol] - 1.) / (eps[vol] + 1.);
          }
        } else if (bc[vol] == Solid::Voltage) {
          potential = volt[vol];
        }
        if (bc[vol] == Solid::Charge || bc[vol] == Solid::DielectricCharge) {
          chargeDensity = charge[vol];
        }
      } else {
        const auto bc1 = bc[vol1[j]];
        const auto bc2 = bc[vol2[j]];
        if (bc1 == Solid::Voltage || bc1 == Solid::Charge ||
            bc1 == Solid::Float) {
          interfaceType = InterfaceType(bc1);
          // First volume is a conductor. Other volume must be a dielectric.
          if (bc2 == Solid::Dielectric || bc2 == Solid::DielectricCharge) {
            if (bc1 == Solid::Voltage) {
              potential = volt[vol1[j]];
            } else if (bc1 == Solid::Charge) {
              chargeDensity = charge[vol1[j]];
            }
          } else {
            interfaceType = -1;
          }
          if (bc1 == Solid::Voltage && bc2 == Solid::Voltage) {
            const double v1 = volt[vol1[j]];
            const double v2 = volt[vol2[j]];
            if (fabs(v1 - v2) < 1.e-6 * (1. + fabs(v1) + fabs(v2))) {
              interfaceType = 0;
            }
          }
        } else if (bc1 == Solid::Dielectric || bc1 == Solid::DielectricCharge) {
          interfaceType = InterfaceType(bc2);
          // First volume is a dielectric.
          if (bc2 == Solid::Voltage) {
            potential = volt[vol2[j]];
          } else if (bc2 == Solid::Charge) {
            chargeDensity = charge[vol2[j]];
          } else if (bc2 == Solid::Dielectric ||
                     bc2 == Solid::DielectricCharge) {
            const double eps1 = eps[vol1[j]];
            const double eps2 = eps[vol2[j]];
            if (fabs(eps1 - eps2) < 1.e-6 * (1. + fabs(eps1) + fabs(eps2))) {
              // Same epsilon. Skip.
              interfaceType = 0;
            } else {
              lambda = (eps1 - eps2) / (eps1 + eps2);
            }
          }
        }
      }
      if (m_debug) {
        if (interfaceType < 0) {
          std::cout << "    Conflicting boundary conditions. Skip.\n";
        } else if (interfaceType < 1) {
          std::cout << "    Trivial interface. Skip.\n";
        } else if (interfaceType > 5) {
          std::cout << "    Interface type " << interfaceType
                    << " is not implemented. Skip.\n";
        }
      }
      if (interfaceType < 0) {
        ++nConflicting;
      } else if (interfaceType < 1) {
        ++nTrivial;
      } else if (interfaceType > 5) {
        ++nNotImplemented;
      }
      if (interfaceType < 1 || interfaceType > 5) continue;

      std::vector<Panel> panelsOut;
      // Reduce to rectangles and right-angle triangles.
      if (m_debug) std::cout << "    Creating primitives.\n";
      MakePrimitives(newPanels[j], panelsOut);
      // Loop over the rectangles and triangles.
      for (auto& panel : panelsOut) {
        const auto& up = panel.xv;
        const auto& vp = panel.yv;
        const auto& wp = panel.zv;
        const unsigned int np = up.size();
        // Rotate.
        xp.assign(np, 0.);
        yp.assign(np, 0.);
        zp.assign(np, 0.);
        for (unsigned int k = 0; k < np; ++k) {
          xp[k] = rot[0][0] * up[k] + rot[1][0] * vp[k] + rot[2][0] * wp[k];
          yp[k] = rot[0][1] * up[k] + rot[1][1] * vp[k] + rot[2][1] * wp[k];
          zp[k] = rot[0][2] * up[k] + rot[1][2] * vp[k] + rot[2][2] * wp[k];
        }
        Primitive primitive;
        primitive.a = panel.a;
        primitive.b = panel.b;
        primitive.c = panel.c;
        primitive.xv = xp;
        primitive.yv = yp;
        primitive.zv = zp;
        primitive.v = potential;
        primitive.q = chargeDensity;
        primitive.lambda = lambda;
        primitive.interface = interfaceType;
        // Set the requested discretization level (target element size).
        primitive.elementSize = -1.;
        if (solids.find(vol1[j]) != solids.end()) {
          const auto solid = solids[vol1[j]];
          if (solid) {
            primitive.elementSize = solid->GetDiscretisationLevel(panel);
          }
        }
        primitive.vol1 = vol1[j];
        primitive.vol2 = vol2[j];
        m_primitives.push_back(std::move(primitive));
      }
    }
  }

  // Add the wires.
  for (unsigned int i = 0; i < nSolids; ++i) {
    const auto solid = m_geometry->GetSolid(i);
    if (!solid) continue;
    if (!solid->IsWire()) continue;
    double x0 = 0., y0 = 0., z0 = 0.;
    solid->GetCentre(x0, y0, z0);
    double dx = 0., dy = 0., dz = 1.;
    solid->GetDirection(dx, dy, dz);
    const double dnorm = sqrt(dx * dx + dy * dy + dz * dz);
    if (dnorm < Small) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Wire has zero norm direction vector; skipped.\n";
      continue;
    }
    dx /= dnorm;
    dy /= dnorm;
    dz /= dnorm;
    const double h = solid->GetHalfLengthZ();
    Primitive primitive;
    primitive.a = solid->GetRadius();
    primitive.b = 0.;
    primitive.c = 0.;
    primitive.xv = {x0 - h * dx, x0 + h * dx};
    primitive.yv = {y0 - h * dy, y0 + h * dy};
    primitive.zv = {z0 - h * dz, z0 + h * dz};
    primitive.v = solid->GetBoundaryPotential();
    primitive.q = solid->GetBoundaryChargeDensity();
    primitive.lambda = 0.;
    primitive.interface = InterfaceType(solid->GetBoundaryConditionType());
    // Set the requested discretization level (target element size).
    Panel panel;
    primitive.elementSize = solid->GetDiscretisationLevel(panel);
    primitive.vol1 = solid->GetId();
    primitive.vol2 = -1;
    m_primitives.push_back(std::move(primitive));
  }
  // Print a warning if we have discarded some panels during the process.
  if (nTrivial > 0 || nConflicting > 0 || nNotImplemented > 0) {
    std::cerr << m_className << "::Initialise:\n";
    if (nConflicting > 0) {
      std::cerr << "    Skipped " << nConflicting
                << " panels with conflicting boundary conditions.\n";
    }
    if (nNotImplemented > 0) {
      std::cerr << "    Skipped " << nNotImplemented
                << " panels with not yet available boundary conditions.\n";
    }
    if (nTrivial > 0) {
      std::cerr << "    Skipped " << nTrivial
                << " panels with trivial boundary conditions.\n";
    }
  }
  if (m_debug) {
    std::cout << m_className << "::Initialise:\n"
              << "    Created " << m_primitives.size() << " primitives.\n";
  }

  // Discretize the primitives.
  for (const auto& primitive : m_primitives) {
    const auto nVertices = primitive.xv.size();
    if (nVertices < 2 || nVertices > 4) continue;
    std::vector<Element> elements;
    // Get the target element size.
    double targetSize = primitive.elementSize;
    if (targetSize < MinDist) targetSize = m_targetElementSize;
    if (nVertices == 2) {
      DiscretizeWire(primitive, targetSize, elements);
    } else if (nVertices == 3) {
      DiscretizeTriangle(primitive, targetSize, elements);
    } else if (nVertices == 4) {
      DiscretizeRectangle(primitive, targetSize, elements);
    }
    for (auto& element : elements) {
      element.interface = primitive.interface;
      element.lambda = primitive.lambda;
      element.bc = primitive.v;
      element.assigned = primitive.q;
      element.solution = 0.;
    }
    m_elements.insert(m_elements.end(),
                      std::make_move_iterator(elements.begin()),
                      std::make_move_iterator(elements.end()));
  }

  // Set the user options.
  // Number of threads and use of primary average properties
  neBEM::NbThreads = m_nThreads;
  neBEM::PrimAfter = m_primAfter;
  neBEM::WtFldPrimAfter = m_wtFldPrimAfter;
  neBEM::OptRmPrim = m_optRmPrim;

  // Fast volume details (physical potential and field related)
  neBEM::OptFastVol = m_optFastVol;
  neBEM::OptCreateFastPF = m_optCreateFastPF;
  neBEM::OptReadFastPF = m_optReadFastPF;
  neBEM::VersionFV = m_versionFV;
  neBEM::NbBlocksFV = m_nbBlocksFV;
  // Weighting potential and field related Fast volume details
  // Since MAXWtFld is not a variable, we do not need neBEM::
  for (int id = 1; id < MAXWtFld; ++id) {
    neBEM::OptWtFldFastVol[id] = m_optWtFldFastVol[id];
    neBEM::OptCreateWtFldFastPF[id] = m_optCreateWtFldFastPF[id];
    neBEM::OptReadWtFldFastPF[id] = m_optReadWtFldFastPF[id];
    neBEM::VersionWtFldFV[id] = m_versionWtFldFV[id];
    neBEM::NbBlocksWtFldFV[id] = m_nbBlocksWtFldFV[id];
  }
  // Known charge options
  neBEM::OptKnCh = m_optKnownCharge;
  // Charging up options
  neBEM::OptChargingUp = m_optChargingUp;

  // Initialize neBEM
  if (neBEM::neBEMInitialize() != 0) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Initialising neBEM failed.\n";
    return false;
  }
  gComponentNeBem3d = this;

  // Discretization related
  neBEM::MinNbElementsOnLength = m_minNbElementsOnLength;
  neBEM::MaxNbElementsOnLength = m_maxNbElementsOnLength;
  neBEM::ElementLengthRqstd = m_targetElementSize * 0.01;

  // New model / reuse existing model flag.
  neBEM::NewModel = m_newModel;
  neBEM::NewMesh = m_newMesh;
  neBEM::NewBC = m_newBC;
  neBEM::NewPP = m_newPP;

  // Store and read options.
  neBEM::OptStoreInflMatrix = m_optStoreInflMatrix;
  neBEM::OptReadInflMatrix = m_optReadInflMatrix;
  neBEM::OptStoreInvMatrix = m_optStoreInvMatrix;
  neBEM::OptReadInvMatrix = m_optReadInvMatrix;
  neBEM::OptStorePrimitives = m_optStorePrimitives;
  neBEM::OptReadPrimitives = m_optReadPrimitives;
  neBEM::OptStoreElements = m_optStoreElements;
  neBEM::OptReadElements = m_optReadElements;
  neBEM::OptFormattedFile = m_optStoreFormatted;
  neBEM::OptUnformattedFile = m_optStoreUnformatted;
  neBEM::OptSystemChargeZero = m_optSystemChargeZero;
  neBEM::OptValidateSolution = m_optValidateSolution;
  neBEM::OptForceValidation = m_optForceValidation;
  neBEM::OptRepeatLHMatrix = m_optRepeatLHMatrix;

  // Compute options
  neBEM::OptSystemChargeZero = m_optSystemChargeZero;
  neBEM::OptValidateSolution = m_optValidateSolution;
  neBEM::OptForceValidation = m_optForceValidation;
  neBEM::OptRepeatLHMatrix = m_optRepeatLHMatrix;

  // Pass debug level.
  neBEM::DebugLevel = m_debug ? 101 : 0;

  // Matrix inversion method (LU or SVD).
  if (m_inversion == Inversion::LU) {
    neBEM::OptInvMatProc = 0;
  } else {
    neBEM::OptInvMatProc = 1;
  }
  // Delete existing weighting fields (if any).
  if (neBEM::WtFieldChDen != NULL) {
    neBEM::neBEMDeleteAllWeightingFields();
  }
  // Transfer the geometry.
  if (neBEM::neBEMReadGeometry() != 0) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Transferring the geometry to neBEM failed.\n";
    return false;
  }

  // Discretization.
  int** elementNbs = neBEM::imatrix(1, neBEM::NbPrimitives, 1, 2);
  for (int i = 1; i <= neBEM::NbPrimitives; ++i) {
    const int vol1 = neBEM::VolRef1[i];
    double size1 = -1.;
    if (solids.find(vol1) != solids.end()) {
      const auto solid = solids[vol1];
      if (solid) {
        Panel panel;
        panel.a = neBEM::XNorm[i];
        panel.b = neBEM::YNorm[i];
        panel.c = neBEM::ZNorm[i];
        const int nv = neBEM::NbVertices[i];
        panel.xv.resize(nv);
        panel.yv.resize(nv);
        panel.zv.resize(nv);
        for (int j = 0; j < nv; ++j) {
          panel.xv[j] = neBEM::XVertex[i][j];
          panel.yv[j] = neBEM::YVertex[i][j];
          panel.zv[j] = neBEM::ZVertex[i][j];
        }
        size1 = solid->GetDiscretisationLevel(panel);
      }
    }
    int vol2 = neBEM::VolRef2[i];
    double size2 = -1.;
    if (solids.find(vol2) != solids.end()) {
      const auto solid = solids[vol2];
      if (solid) {
        Panel panel;
        panel.a = -neBEM::XNorm[i];
        panel.b = -neBEM::YNorm[i];
        panel.c = -neBEM::ZNorm[i];
        const int nv = neBEM::NbVertices[i];
        panel.xv.resize(nv);
        panel.yv.resize(nv);
        panel.zv.resize(nv);
        for (int j = 0; j < nv; ++j) {
          panel.xv[j] = neBEM::XVertex[i][j];
          panel.yv[j] = neBEM::YVertex[i][j];
          panel.zv[j] = neBEM::ZVertex[i][j];
        }
        size2 = solid->GetDiscretisationLevel(panel);
      }
    }
    double size = m_targetElementSize;
    if (size1 > 0. && size2 > 0.) {
      size = std::min(size1, size2);
    } else if (size1 > 0.) {
      size = size1;
    } else if (size2 > 0.) {
      size = size2;
    }
    // Convert from cm to m.
    size *= 0.01;

    // Work out the element dimensions.
    const double dx1 = neBEM::XVertex[i][0] - neBEM::XVertex[i][1];
    const double dy1 = neBEM::YVertex[i][0] - neBEM::YVertex[i][1];
    const double dz1 = neBEM::ZVertex[i][0] - neBEM::ZVertex[i][1];
    const double l1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
    int nb1 = (int)(sqrt(l1) / size);
    // Truncate to desired range.
    if (nb1 < neBEM::MinNbElementsOnLength) {
      nb1 = neBEM::MinNbElementsOnLength;
    } else if (nb1 > neBEM::MaxNbElementsOnLength) {
      nb1 = neBEM::MaxNbElementsOnLength;
    }
    int nb2 = 0;
    if (neBEM::NbVertices[i] > 2) {
      const double dx2 = neBEM::XVertex[i][2] - neBEM::XVertex[i][1];
      const double dy2 = neBEM::YVertex[i][2] - neBEM::YVertex[i][1];
      const double dz2 = neBEM::ZVertex[i][2] - neBEM::ZVertex[i][1];
      const double l2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
      nb2 = (int)(sqrt(l2) / size);
      if (nb2 < neBEM::MinNbElementsOnLength) {
        nb2 = neBEM::MinNbElementsOnLength;
      } else if (nb2 > neBEM::MaxNbElementsOnLength) {
        nb2 = neBEM::MaxNbElementsOnLength;
      }
    }
    elementNbs[i][1] = nb1;
    elementNbs[i][2] = nb2;
  }

  if (neBEM::neBEMDiscretize(elementNbs) != 0) {
    std::cerr << m_className << "::Initialise: Discretization failed.\n";
    neBEM::free_imatrix(elementNbs, 1, neBEM::NbPrimitives, 1, 2);
    return false;
  }
  neBEM::free_imatrix(elementNbs, 1, neBEM::NbPrimitives, 1, 2);
  if (neBEM::neBEMBoundaryInitialConditions() != 0) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Setting the boundary and initial conditions failed.\n";
    return false;
  }
  if (neBEM::neBEMSolve() != 0) {
    std::cerr << m_className << "::Initialise: Solution failed.\n";
    return false;
  }
  // Now the weighting fields.
  std::set<std::string> labels;
  for (unsigned int i = 0; i < nSolids; ++i) {
    const auto solid = m_geometry->GetSolid(i);
    if (!solid) continue;
    const std::string label = solid->GetLabel();
    if (!label.empty()) labels.insert(label);
  }
  for (const auto& label : labels) {
    std::vector<int> primitives;
    for (unsigned int i = 0; i < nSolids; ++i) {
      const auto solid = m_geometry->GetSolid(i);
      if (!solid) continue;
      if (solid->GetLabel() != label) continue;
      const int id = solid->GetId();
      // Add the primitives associated to this solid to the list.
      for (int j = 1; j <= neBEM::NbPrimitives; ++j) {
        if (neBEM::VolRef1[j] == id || neBEM::VolRef2[j] == id) {
          primitives.push_back(j);
        }
      }
    }
    // Request the weighting field for this list of primitives.
    const int np = primitives.size();
    const int id = neBEM::neBEMPrepareWeightingField(np, primitives.data());
    if (id < 0) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Weighting field calculation for readout group \""
                << label << "\" failed.\n";
      continue;
    } else {
      std::cout << m_className << "::Initialise:\n"
                << "    Prepared weighting field for readout group \"" << label
                << "\".\n";
      m_wfields[label] = id;
    }
  }
  // TODO! Not sure if we should call this here.
  // neBEM::neBEMEnd();
  m_ready = true;
  return true;
}

void ComponentNeBem3d::ShiftPanels(std::vector<Panel>& panels) const {
  // *---------------------------------------------------------------------
  // *   BEMBAS - Reduces panels to the basic period.
  // *---------------------------------------------------------------------
  const bool perx = m_periodic[0] || m_mirrorPeriodic[0];
  const bool pery = m_periodic[1] || m_mirrorPeriodic[1];
  const bool perz = m_periodic[2] || m_mirrorPeriodic[2];
  // Nothing to do if there is no periodicity.
  if (!perx && !pery && !perz) return;

  // Loop over the panels.
  for (auto& panel : panels) {
    const auto nv = panel.xv.size();
    if (nv == 0) continue;
    // Determine the centre of gravity.
    double xc = std::accumulate(panel.xv.begin(), panel.xv.end(), 0.);
    double yc = std::accumulate(panel.yv.begin(), panel.yv.end(), 0.);
    double zc = std::accumulate(panel.zv.begin(), panel.zv.end(), 0.);
    xc /= nv;
    yc /= nv;
    zc /= nv;
    // Any change ?
    constexpr double eps = 1.e-6;
    double rx = 0.;
    int nx = 0;
    if (perx && m_periodicLength[0] > 0.) {
      rx = xc / m_periodicLength[0];
      nx = std::round(rx);
      if (std::abs(rx - nx - 0.5) < eps) ++nx;
    }
    double ry = 0.;
    int ny = 0;
    if (pery && m_periodicLength[1] > 0.) {
      ry = yc / m_periodicLength[1];
      ny = std::round(ry);
      if (std::abs(ry - ny - 0.5) < eps) ++ny;
    }
    double rz = 0.;
    int nz = 0;
    if (perz && m_periodicLength[2] > 0.) {
      rz = zc / m_periodicLength[2];
      nz = std::round(rz);
      if (std::abs(rz - nz - 0.5) < eps) ++nz;
    }
    // Skip if there is no shift.
    if (nx == 0 && ny == 0 && nz == 0) continue;
    // Shift for x-periodicity.
    if (nx != 0) {
      const double shift = nx * m_periodicLength[0];
      for (auto& x : panel.xv) x -= shift;
    }
    // Shift for y-periodicity.
    if (ny != 0) {
      const double shift = ny * m_periodicLength[1];
      for (auto& y : panel.yv) y -= shift;
    }
    // Shift for z-periodicity.
    if (nz != 0) {
      const double shift = nz * m_periodicLength[2];
      for (auto& z : panel.zv) z -= shift;
    }
  }
}

int ComponentNeBem3d::InterfaceType(const Solid::BoundaryCondition bc) const {
  // 0: To be skipped
  // 1: Conductor-dielectric
  // 2: Conductor with known charge
  // 3: Conductor at floating potential
  // 4: Dielectric-dielectric
  // 5: Dielectric with given surface charge
  //
  // Check dielectric-dielectric formulation in
  // Jaydeep P. Bardhan,
  // Numerical solution of boundary-integral equations for molecular
  // electrostatics, J. Chem. Phys. 130, 094102 (2009)
  // https://doi.org/10.1063/1.3080769

  switch (bc) {
    case Solid::Voltage:
      return 1;
      break;
    case Solid::Charge:
      return 2;
      break;
    case Solid::Float:
      return 3;
      break;
    case Solid::Dielectric:
      return 4;
      break;
    case Solid::DielectricCharge:
      return 5;
      break;
    case Solid::ParallelField:
      // Symmetry boundary, E parallel.
      return 6;
      break;
    case Solid::PerpendicularField:
      // Symmetry boundary, E perpendicular.
      return 7;
      break;
    default:
      break;
  }
  return 0;
}

bool ComponentNeBem3d::DiscretizeWire(const Primitive& primitive,
                                      const double targetSize,
                                      std::vector<Element>& elements) const {
  const double dx = primitive.xv[1] - primitive.xv[0];
  const double dy = primitive.yv[1] - primitive.yv[0];
  const double dz = primitive.zv[1] - primitive.zv[0];
  const double lw = sqrt(dx * dx + dy * dy + dz * dz);
  // Determine the number of segments along the wire.
  unsigned int nSegments = NbOfSegments(lw, targetSize);
  const double elementSize = lw / nSegments;

  // Determine the direction cosines.
  // The z axis of the local coordinate system is along the wire.
  const std::array<double, 3> zu = {dx / lw, dy / lw, dz / lw};
  // Make the x axis orthogonal in the two largest components.
  std::array<double, 3> xu;
  if (fabs(zu[0]) >= fabs(zu[2]) && fabs(zu[1]) >= fabs(zu[2])) {
    xu = {-zu[1], zu[0], 0.};
  } else if (fabs(zu[0]) >= fabs(zu[1]) && fabs(zu[2]) >= fabs(zu[1])) {
    xu = {-zu[2], 0., zu[0]};
  } else {
    xu = {0., zu[2], -zu[1]};
  }
  xu = UnitVector(xu);
  // The y axis is given by the vectorial product of z and x.
  const std::array<double, 3> yu = UnitVector(CrossProduct(zu, xu));
  const std::array<std::array<double, 3>, 3> dcos = {xu, yu, zu};

  const double xincr = dx / nSegments;
  const double yincr = dy / nSegments;
  const double zincr = dz / nSegments;

  // TODO!
  const double radius = 1.;
  const double dA = TwoPi * radius * elementSize;
  for (unsigned int i = 0; i < nSegments; ++i) {
    const double x0 = primitive.xv[0] + i * xincr;
    const double y0 = primitive.yv[0] + i * yincr;
    const double z0 = primitive.zv[0] + i * zincr;
    Element element;
    element.origin = {x0 + 0.5 * xincr, y0 + 0.5 * yincr, z0 + 0.5 * zincr};
    element.xv = {x0, x0 + xincr};
    element.yv = {y0, y0 + yincr};
    element.zv = {z0, z0 + zincr};
    element.lx = radius;
    element.lz = elementSize;
    element.dA = dA;
    element.dcos = dcos;
    // Modify to be on surface?
    element.collocationPoint = element.origin;
    elements.push_back(std::move(element));
  }
  return true;
}

bool ComponentNeBem3d::DiscretizeTriangle(
    const Primitive& primitive, const double targetSize,
    std::vector<Element>& elements) const {
  // Origin of the local coordinate system is at the right angle corner.
  std::array<double, 3> corner = {primitive.xv[1], primitive.yv[1],
                                  primitive.zv[1]};
  // Determine the direction cosines.
  const double dx1 = primitive.xv[0] - primitive.xv[1];
  const double dy1 = primitive.yv[0] - primitive.yv[1];
  const double dz1 = primitive.zv[0] - primitive.zv[1];
  const double dx2 = primitive.xv[2] - primitive.xv[1];
  const double dy2 = primitive.yv[2] - primitive.yv[1];
  const double dz2 = primitive.zv[2] - primitive.zv[1];
  // Normal vector of the primitive.
  std::array<double, 3> nu = {primitive.a, primitive.b, primitive.c};
  nu = UnitVector(nu);
  // int flagDC = 1;
  // We begin with trial 1: one of the possible orientations.
  double lx = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  double lz = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
  std::array<double, 3> xu = {dx1 / lx, dy1 / lx, dz1 / lx};
  std::array<double, 3> zu = {dx2 / lz, dy2 / lz, dz2 / lz};
  std::array<double, 3> yu = CrossProduct(zu, xu);
  constexpr double tol = 1.e-3;
  if ((fabs(yu[0] - nu[0]) > tol) || (fabs(yu[1] - nu[1]) > tol) ||
      (fabs(yu[2] - nu[2]) > tol)) {
    // flagDC = 2;
    // Try the other orientation.
    std::swap(lx, lz);
    xu = {dx2 / lx, dy2 / lx, dz2 / lx};
    zu = {dx1 / lz, dy1 / lz, dz1 / lz};
    yu = CrossProduct(zu, xu);
    if ((fabs(yu[0] - nu[0]) > tol) || (fabs(yu[1] - nu[1]) > tol) ||
        (fabs(yu[2] - nu[2]) > tol)) {
      // No other possibility, search failed.
      std::cerr << m_className << "::DiscretizeTriangle:\n"
                << "    Could not establish direction vectors.\n";
      return false;
    }
  }
  std::array<std::array<double, 3>, 3> dcos = {xu, yu, zu};

  // TODO!
  /*
  double eps1 = primitive.eps1;
  double eps2 = primitive.eps2;
  if (flagDC == 1) std::swap(eps1, eps2);
  */

  // Determine the number of elements.
  unsigned int nx = NbOfSegments(lx, targetSize);
  unsigned int nz = NbOfSegments(lz, targetSize);
  double elementSizeX = lx / nx;
  double elementSizeZ = lz / nz;

  // Analyse the aspect ratio.
  double ar = elementSizeX / elementSizeZ;
  if (ar > 10.) {
    // Reduce nz.
    elementSizeZ = 0.1 * elementSizeX;
    nz = std::max(static_cast<int>(lz / elementSizeZ), 1);
    elementSizeZ = lz / nz;
  }
  if (ar < 0.1) {
    // Reduce nx.
    elementSizeX = 0.1 * elementSizeZ;
    nx = std::max(static_cast<int>(lx / elementSizeX), 1);
    elementSizeX = lx / nx;
  }
  const double dxdz = lx / lz;
  for (unsigned int k = 0; k < nz; ++k) {
    // Consider the k-th row.
    const double zlo = k * elementSizeZ;
    const double zhi = (k + 1) * elementSizeZ;
    double xlo = (lz - zlo) * dxdz;
    double xhi = (lz - zhi) * dxdz;
    // Triangular element on the k-th row.
    Element triangle;
    triangle.origin = LocalToGlobal(xhi, 0., zlo, dcos, corner);
    // Assign element values.
    const double triangleSizeX = xlo - xhi;
    triangle.lx = triangleSizeX;
    triangle.lz = elementSizeZ;
    triangle.dA = 0.5 * triangleSizeX * elementSizeZ;
    triangle.dcos = dcos;
    // Calculate the element vertices.
    const double xv0 = triangle.origin[0];
    const double yv0 = triangle.origin[1];
    const double zv0 = triangle.origin[2];
    const double xv1 = xv0 + triangleSizeX * xu[0];
    const double yv1 = yv0 + triangleSizeX * xu[1];
    const double zv1 = zv0 + triangleSizeX * xu[2];
    const double xv2 = xv0 + elementSizeZ * zu[0];
    const double yv2 = yv0 + elementSizeZ * zu[1];
    const double zv2 = zv0 + elementSizeZ * zu[2];
    // Assign vertices of the element.
    triangle.xv = {xv0, xv1, xv2};
    triangle.yv = {yv0, yv1, yv2};
    triangle.zv = {zv0, zv1, zv2};

    // Boundary condition is applied at the barycentre, not at the origin
    // of the element coordinate system which is at the right-angled corner.
    const double xb = triangleSizeX / 3.;
    const double yb = 0.;
    const double zb = elementSizeZ / 3.;
    triangle.collocationPoint =
        LocalToGlobal(xb, yb, zb, dcos, triangle.origin);
    elements.push_back(std::move(triangle));
    // No rectangular element on the last row.
    if (k == nz - 1) continue;

    // Determine number and size in x of the rectangular elements.
    const int nRect = xhi <= elementSizeX ? 1 : (int)(xhi / elementSizeX);
    const double rectSizeX = xhi / nRect;
    const double zc = 0.5 * (zlo + zhi);
    for (int i = 0; i < nRect; ++i) {
      // Centroid of the rectangular element.
      const double xc = (i + 0.5) * rectSizeX;
      const double yc = 0.;
      std::array<double, 3> centre = LocalToGlobal(xc, yc, zc, dcos, corner);
      // Assign element values.
      Element rect;
      rect.origin = centre;
      rect.lx = rectSizeX;
      rect.lz = elementSizeZ;
      rect.dA = rectSizeX * elementSizeZ;
      rect.dcos = dcos;
      // Boundary condition is applied at the origin of the rectangular
      // element coordinate system.
      rect.collocationPoint = centre;
      // Calculate the element vertices.
      const double hx = 0.5 * rectSizeX;
      const double hz = 0.5 * elementSizeZ;
      std::array<double, 3> p0 = LocalToGlobal(-hx, 0., -hz, dcos, centre);
      std::array<double, 3> p1 = LocalToGlobal(hx, 0., -hz, dcos, centre);
      std::array<double, 3> p2 = LocalToGlobal(hx, 0., hz, dcos, centre);
      std::array<double, 3> p3 = LocalToGlobal(-hx, 0., hz, dcos, centre);
      // Assign the vertices of the element.
      rect.xv = {p0[0], p1[0], p2[0], p3[0]};
      rect.yv = {p0[1], p1[1], p2[1], p3[1]};
      rect.zv = {p0[2], p1[2], p2[2], p3[2]};
      elements.push_back(std::move(rect));
    }
  }
  return true;
}

bool ComponentNeBem3d::DiscretizeRectangle(
    const Primitive& primitive, const double targetSize,
    std::vector<Element>& elements) const {
  std::array<double, 3> origin = {
      0.25 * std::accumulate(primitive.xv.begin(), primitive.xv.end(), 0.),
      0.25 * std::accumulate(primitive.yv.begin(), primitive.yv.end(), 0.),
      0.25 * std::accumulate(primitive.zv.begin(), primitive.zv.end(), 0.)};

  // Determine the direction cosines.
  const double dx1 = primitive.xv[1] - primitive.xv[0];
  const double dy1 = primitive.yv[1] - primitive.yv[0];
  const double dz1 = primitive.zv[1] - primitive.zv[0];
  const double dx2 = primitive.xv[2] - primitive.xv[1];
  const double dy2 = primitive.yv[2] - primitive.yv[1];
  const double dz2 = primitive.zv[2] - primitive.zv[1];
  double lx = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  double lz = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
  const std::array<double, 3> xu = {dx1 / lx, dy1 / lx, dz1 / lx};
  const std::array<double, 3> yu = {primitive.a, primitive.b, primitive.c};
  const std::array<double, 3> zu = CrossProduct(xu, yu);
  const std::array<std::array<double, 3>, 3> dcos = {xu, yu, zu};

  // Determine the number of elements.
  unsigned int nx = NbOfSegments(lx, targetSize);
  unsigned int nz = NbOfSegments(lz, targetSize);

  double elementSizeX = lx / nx;
  double elementSizeZ = lz / nz;

  // Analyze the element aspect ratio.
  double ar = elementSizeX / elementSizeZ;
  if (ar > 10.) {
    // Reduce nz.
    elementSizeZ = 0.1 * elementSizeX;
    nz = std::max(static_cast<int>(lz / elementSizeZ), 1);
    elementSizeZ = lz / nz;
  }
  if (ar < 0.1) {
    // Reduce nx.
    elementSizeX = 0.1 * elementSizeZ;
    nx = std::max(static_cast<int>(lx / elementSizeX), 1);
    elementSizeX = lx / nx;
  }

  const double dA = elementSizeX * elementSizeZ;
  for (unsigned int i = 0; i < nx; ++i) {
    // Centroid of the element.
    const double xav = -0.5 * lx + (i + 0.5) * elementSizeX;
    for (unsigned int k = 0; k < nz; ++k) {
      const double zav = -0.5 * lz + (k + 0.5) * elementSizeZ;

      std::array<double, 3> centre = LocalToGlobal(xav, 0., zav, dcos, origin);
      Element element;
      element.origin = centre;
      element.lx = elementSizeX;
      element.lz = elementSizeZ;
      element.dA = dA;
      element.dcos = dcos;
      element.collocationPoint = centre;
      // Calculate the element vertices.
      const double hx = 0.5 * elementSizeX;
      const double hz = 0.5 * elementSizeZ;
      std::array<double, 3> p0 = LocalToGlobal(-hx, 0., -hz, dcos, centre);
      std::array<double, 3> p1 = LocalToGlobal(hx, 0., -hz, dcos, centre);
      std::array<double, 3> p2 = LocalToGlobal(hx, 0., hz, dcos, centre);
      std::array<double, 3> p3 = LocalToGlobal(-hx, 0., hz, dcos, centre);
      // Assign the vertices of the element.
      element.xv = {p0[0], p1[0], p2[0], p3[0]};
      element.yv = {p0[1], p1[1], p2[1], p3[1]};
      element.zv = {p0[2], p1[2], p2[2], p3[2]};
      elements.push_back(std::move(element));
    }
  }
  return true;
}

unsigned int ComponentNeBem3d::NbOfSegments(const double length,
                                            const double target) const {
  // Check whether the length of the primitive is long enough.
  if (length < MinDist) return 1;
  unsigned int n = static_cast<unsigned int>(length / target);
  if (n < m_minNbElementsOnLength) {
    // Need to have a minimum number of elements per primitive...
    n = m_minNbElementsOnLength;
    if (length < n * MinDist) {
      // ...which may not be possible if the length is small.
      n = static_cast<unsigned int>(length / MinDist);
      if (n < 1) {
        // However, it is necessary to have at least one element!
        n = 1;
      }
    }
  }
  return std::min(n, m_maxNbElementsOnLength);
}

bool ComponentNeBem3d::EliminateOverlaps(const Panel& panel1,
                                         const Panel& panel2,
                                         std::vector<Panel>& panelsOut,
                                         std::vector<int>& itypo) {
  // *-----------------------------------------------------------------------
  // *   PLAOVL - Isolates the parts of plane 1 that are not hidden by 2.
  // *-----------------------------------------------------------------------

  const auto& xp1 = panel1.xv;
  const auto& yp1 = panel1.yv;
  const auto& zp1 = panel1.zv;
  const auto& xp2 = panel2.xv;
  const auto& yp2 = panel2.yv;
  const auto& zp2 = panel2.zv;
  // If the size of either is less than 3, simply return.
  if (xp1.size() <= 2 || xp2.size() <= 2) {
    return true;
  }
  // Compute the various tolerances.
  const double xmin1 = *std::min_element(std::begin(xp1), std::end(xp1));
  const double ymin1 = *std::min_element(std::begin(yp1), std::end(yp1));
  const double xmax1 = *std::max_element(std::begin(xp1), std::end(xp1));
  const double ymax1 = *std::max_element(std::begin(yp1), std::end(yp1));

  const double xmin2 = *std::min_element(std::begin(xp2), std::end(xp2));
  const double ymin2 = *std::min_element(std::begin(yp2), std::end(yp2));
  const double xmax2 = *std::max_element(std::begin(xp2), std::end(xp2));
  const double ymax2 = *std::max_element(std::begin(yp2), std::end(yp2));

  const double xmin = std::min(xmin1, xmin2);
  const double ymin = std::min(ymin1, ymin2);
  const double xmax = std::max(xmax1, xmax2);
  const double ymax = std::max(ymax1, ymax2);

  const double epsx = 1.e-6 * std::max(std::abs(xmax), std::abs(xmin));
  const double epsy = 1.e-6 * std::max(std::abs(ymax), std::abs(ymin));

  const double zsum1 = std::accumulate(std::begin(zp1), std::end(zp1), 0.);
  const double zsum2 = std::accumulate(std::begin(zp2), std::end(zp2), 0.);
  const double zmean = (zsum1 + zsum2) / (zp1.size() + zp2.size());

  std::array<std::vector<double>, 2> xl;
  std::array<std::vector<double>, 2> yl;
  std::array<std::vector<int>, 2> flags;
  std::array<std::vector<double>, 2> qs;
  // Establish the list of special points around polygon 1.
  AddPoints(xp1, yp1, xp2, yp2, xl[0], yl[0], flags[0], qs[0], epsx, epsy);
  // Establish the list of special points around polygon 2.
  AddPoints(xp2, yp2, xp1, yp1, xl[1], yl[1], flags[1], qs[1], epsx, epsy);

  bool ok = true;
  // Look up the cross-links: from plane 1 (2) to plane 2 (1).
  std::array<std::vector<int>, 2> links;
  for (unsigned int ic = 0; ic < 2; ++ic) {
    const unsigned int n1 = xl[ic].size();
    links[ic].assign(n1, -1);
    const unsigned int jc = ic == 0 ? 1 : 0;
    const unsigned int n2 = xl[jc].size();
    for (unsigned int i = 0; i < n1; ++i) {
      unsigned int nFound = 0;
      for (unsigned int j = 0; j < n2; ++j) {
        if (fabs(xl[ic][i] - xl[jc][j]) < epsx &&
            fabs(yl[ic][i] - yl[jc][j]) < epsy) {
          ++nFound;
          links[ic][i] = j;
        }
      }
      if (nFound == 0 && (flags[ic][i] == 2 || flags[ic][i] == 3)) {
        std::cerr << m_className << "::EliminateOverlaps: "
                  << "Warning. Expected match not found (" << ic << "-" << jc
                  << ").\n";
        links[ic][i] = -1;
        ok = false;
      } else if (nFound > 1) {
        std::cerr << m_className << "::EliminateOverlaps: "
                  << "Warning. More than 1 match found (" << ic << "-" << jc
                  << ").\n";
        links[ic][i] = -1;
        ok = false;
      }
    }
  }

  // List the points for debugging.
  if (m_debug) {
    for (unsigned int j = 0; j < 2; ++j) {
      std::cout << "      Polygon " << j << "\n      "
                << " No Type            x            y        Q   links\n";
      const unsigned int n = xl[j].size();
      for (unsigned int i = 0; i < n; ++i) {
        printf("        %3d %5d %13.6f %13.6f %5.3f %3d\n", i, flags[j][i],
               xl[j][i], yl[j][i], qs[j][i], links[j][i]);
      }
    }
  }
  if (!ok) return false;

  for (unsigned int ic = 0; ic < 2; ++ic) {
    // See whether all of 1 (2) is inside 2 (1).
    bool allInside = true;
    const unsigned int np = xl[ic].size();
    for (unsigned int i = 0; i < np; ++i) {
      if (flags[ic][i] != 1) {
        allInside = false;
        break;
      }
      bool inside = false, edge = false;
      if (ic == 0) {
        Polygon::Inside(xp2, yp2, xl[ic][i], yl[ic][i], inside, edge);
      } else {
        Polygon::Inside(xp1, yp1, xl[ic][i], yl[ic][i], inside, edge);
      }
      if (!(inside || edge)) {
        allInside = false;
        break;
      }
    }
    if (allInside) {
      // Apparently 1 (2) really is fully inside 2 (1).
      if (ic == 0) {
        if (m_debug) std::cout << "      Curve 0 fully inside 1.\n";
        // Write out curve 1.
        panelsOut.push_back(panel1);
      } else {
        if (m_debug) std::cout << "      Curve 1 fully inside 0.\n";
        // Write out curve 2.
        panelsOut.push_back(panel2);
      }
      panelsOut.back().zv.assign(panelsOut.back().xv.size(), zmean);
      itypo.push_back(3);
      std::vector<Panel> newPanels;
      if (ic == 0) {
        if (!TraceEnclosed(xl[0], yl[0], xl[1], yl[1], panel2, newPanels)) {
          return false;
        }
      } else {
        if (!TraceEnclosed(xl[1], yl[1], xl[0], yl[0], panel1, newPanels)) {
          return false;
        }
      }
      for (auto& panel : newPanels) {
        panel.zv.assign(panel.xv.size(), zmean);
        if (ic == 0) {
          itypo.push_back(2);
        } else {
          itypo.push_back(1);
        }
      }
      panelsOut.insert(panelsOut.end(), newPanels.begin(), newPanels.end());
      return true;
    }
  }

  for (unsigned int ic = 0; ic < 2; ++ic) {
    std::vector<Panel> newPanels;
    const unsigned int n = xl[ic].size();
    // Identify the parts of 1 (2) that are not overlapped, first mark.
    std::vector<bool> mark(n, false);
    bool done = false;
    while (!done) {
      if (m_debug) {
        std::cout << "      Searching for starting point on " << ic << ".\n";
      }
      done = true;
      // Try and find a new starting point
      for (unsigned int i = 0; i < n; ++i) {
        const unsigned int ii = NextPoint(i, n);
        // Skip parts already processed.
        if (mark[i] || mark[ii]) continue;
        // Skip if mid point is inside other volume.
        bool inside = false, edge = false;
        const double xm = 0.5 * (xl[ic][i] + xl[ic][ii]);
        const double ym = 0.5 * (yl[ic][i] + yl[ic][ii]);
        if (ic == 0) {
          Polygon::Inside(xp2, yp2, xm, ym, inside, edge);
        } else {
          Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
        }
        if (inside || edge) continue;
        // Found one.
        done = false;

        if (ic == 0) {
          // Trace this part of 1 outside 2.
          TraceNonOverlap(xp1, yp1, xl[0], yl[0], xl[1], yl[1], flags[0],
                          flags[1], links[0], links[1], mark, i, panel1,
                          newPanels);
        } else {
          // Trace this part of 2 outside 1.
          TraceNonOverlap(xp2, yp2, xl[1], yl[1], xl[0], yl[0], flags[1],
                          flags[0], links[1], links[0], mark, i, panel2,
                          newPanels);
        }
        break;
      }
    }
    for (auto& panel : newPanels) {
      panel.zv.assign(panel.xv.size(), zmean);
      itypo.push_back(ic + 1);
    }
    panelsOut.insert(panelsOut.end(), newPanels.begin(), newPanels.end());
    if (m_debug) {
      std::cout << "      No further non-overlapped areas of " << ic << ".\n";
    }
  }

  // Look for the overlapped parts.
  std::vector<Panel> newPanels;
  const unsigned int n1 = xl[0].size();
  std::vector<bool> mark1(n1, false);
  bool done = false;
  while (!done) {
    done = true;
    if (m_debug) {
      std::cout << "      Searching for starting point on overlap.\n";
    }
    for (unsigned int i = 0; i < n1; ++i) {
      // Skip points already processed.
      if (mark1[i]) continue;
      // Skip if not an edge point on both 1 and 2 or internal in 2.
      int ip1 = i;
      int ip2 = links[0][ip1];
      if (ip2 < 0 || flags[0][ip1] == 1) {
        bool inside = false, edge = false;
        Polygon::Inside(xp2, yp2, xl[0][ip1], yl[0][ip1], inside, edge);
        if (!(inside || edge)) continue;
      } else if (flags[1][ip2] == 1) {
        continue;
      }
      // Found one.
      done = false;
      TraceOverlap(xp1, yp1, xp2, yp2, xl[0], yl[0], xl[1], yl[1], flags[0],
                   links[0], links[1], mark1, ip1, ip2, panel1, newPanels);
      break;
    }
  }
  for (auto& panel : newPanels) {
    panel.zv.assign(panel.xv.size(), zmean);
    itypo.push_back(3);
  }
  panelsOut.insert(panelsOut.end(), newPanels.begin(), newPanels.end());
  // Finished
  if (m_debug) std::cout << "      No further overlapped areas.\n";
  return true;
}

bool ComponentNeBem3d::TraceEnclosed(const std::vector<double>& xl1,
                                     const std::vector<double>& yl1,
                                     const std::vector<double>& xl2,
                                     const std::vector<double>& yl2,
                                     const Panel& panel2,
                                     std::vector<Panel>& panelsOut) const {
  const int n1 = xl1.size();
  const int n2 = xl2.size();
  // Find 2 non-crossing connections: JP1-JP2 and KP1-KP2.
  unsigned int nFound = 0;
  int jp1 = 0, jp2 = 0;
  int kp1 = 0, kp2 = 0;
  for (int ip1 = 0; ip1 < n1; ++ip1) {
    const double x1 = xl1[ip1];
    const double y1 = yl1[ip1];
    for (int ip2 = 0; ip2 < n2; ++ip2) {
      if (nFound > 0 && ip2 == jp2) continue;
      const double x2 = xl2[ip2];
      const double y2 = yl2[ip2];
      bool cross = false;
      for (int k = 0; k < n1; ++k) {
        const int kk = NextPoint(k, n1);
        if (k == ip1 || kk == ip1) continue;
        double xc = 0., yc = 0.;
        cross =
            Crossing(x1, y1, x2, y2, xl1[k], yl1[k], xl1[kk], yl1[kk], xc, yc);
        if (cross) break;
      }
      if (cross) continue;
      if (m_debug) std::cout << "      No crossing with 1.\n";
      for (int k = 0; k < n2; ++k) {
        const int kk = NextPoint(k, n2);
        if (k == ip2 || kk == ip2) continue;
        double xc = 0., yc = 0.;
        cross =
            Crossing(x1, y1, x2, y2, xl2[k], yl2[k], xl2[kk], yl2[kk], xc, yc);
        if (cross) break;
      }
      if (cross) continue;
      if (nFound == 0) {
        jp1 = ip1;
        jp2 = ip2;
        if (m_debug) {
          std::cout << "      1st junction: " << jp1 << ", " << jp2 << ".\n";
        }
        ++nFound;
        break;
      } else {
        kp1 = ip1;
        kp2 = ip2;
        double xc = 0., yc = 0.;
        cross = Crossing(x1, y1, x2, y2, xl1[jp1], yl1[jp1], xl2[jp2], yl2[jp2],
                         xc, yc);
        if (!cross) {
          if (m_debug) {
            std::cout << "      2nd junction: " << kp1 << ", " << kp2 << ".\n";
          }
          ++nFound;
          break;
        }
      }
    }
    if (nFound > 1) break;
  }
  if (nFound < 2) {
    std::cerr << m_className << "::TraceEnclosed: Found no cut-out.\n";
    return false;
  }

  // Create part 1 of area 2.
  std::vector<double> xpl;
  std::vector<double> ypl;
  if (m_debug) std::cout << "      Creating part 1 of area 2.\n";
  for (int ip1 = jp1; ip1 <= kp1; ++ip1) {
    if (m_debug) std::cout << "        Adding " << ip1 << " on 1.\n";
    xpl.push_back(xl1[ip1]);
    ypl.push_back(yl1[ip1]);
  }
  // Try one way.
  int imax = jp2 < kp2 ? jp2 + n2 : jp2;
  int dir = +1;
  for (int i = kp2; i <= imax; ++i) {
    int ip2 = i % n2;
    if (m_debug) std::cout << "        Adding " << ip2 << " on 2.\n";
    xpl.push_back(xl2[ip2]);
    ypl.push_back(yl2[ip2]);
  }
  // Check for undesirable crossings.
  bool ok = true;
  for (int ip1 = 0; ip1 < n1; ++ip1) {
    if (ip1 == jp1 || ip1 == kp1) continue;
    bool inside = false, edge = false;
    Polygon::Inside(xpl, ypl, xl1[ip1], yl1[ip1], inside, edge);
    if (inside) {
      ok = false;
      break;
    }
  }
  if (!ok) {
    // Use the other way if this failed
    if (m_debug) std::cout << "      Trying the other direction.\n";
    xpl.resize(kp1 - jp1 + 1);
    ypl.resize(kp1 - jp1 + 1);
    imax = jp2 < kp2 ? kp2 : kp2 + n2;
    dir = -1;
    for (int i = imax; i >= jp2; --i) {
      const int ip2 = i % n2;
      if (m_debug) std::cout << "        Adding " << ip2 << " on 2.\n";
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
    }
  }

  // Save this part.
  Panel newPanel1 = panel2;
  newPanel1.xv = xpl;
  newPanel1.yv = ypl;
  panelsOut.push_back(std::move(newPanel1));
  if (m_debug) std::cout << "      Part 1 has " << xpl.size() << " nodes.\n";

  // Create part 2 of area 2.
  xpl.clear();
  ypl.clear();
  if (m_debug) std::cout << "      Creating part 2 of area 2.\n";
  imax = jp1 + n1;
  for (int i = kp1; i <= imax; ++i) {
    const int ip1 = i % n1;
    if (m_debug) std::cout << "        Adding " << ip1 << " on 1.\n";
    xpl.push_back(xl1[ip1]);
    ypl.push_back(yl1[ip1]);
  }
  // Add the part over area 2.
  if (dir == -1) {
    imax = jp2 > kp2 ? jp2 : jp2 + n2;
    for (int i = imax; i >= kp2; --i) {
      const int ip2 = i % n2;
      if (m_debug) std::cout << "        Adding " << ip2 << " on 2.\n";
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
    }
  } else {
    imax = jp2 > kp2 ? kp2 + n2 : kp2;
    for (int i = jp2; i <= imax; ++i) {
      const int ip2 = i % n2;
      if (m_debug) std::cout << "        Adding " << ip2 << " on 2.\n";
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
    }
  }
  // Save this part.
  Panel newPanel2 = panel2;
  newPanel2.xv = xpl;
  newPanel2.yv = ypl;
  panelsOut.push_back(std::move(newPanel2));
  if (m_debug) std::cout << "      Part 1 has " << xpl.size() << " nodes.\n";
  return true;
}

void ComponentNeBem3d::TraceNonOverlap(
    const std::vector<double>& xp1, const std::vector<double>& yp1,
    const std::vector<double>& xl1, const std::vector<double>& yl1,
    const std::vector<double>& xl2, const std::vector<double>& yl2,
    const std::vector<int>& flags1, const std::vector<int>& flags2,
    const std::vector<int>& links1, const std::vector<int>& links2,
    std::vector<bool>& mark1, int ip1, const Panel& panel1,
    std::vector<Panel>& panelsOut) const {
  const unsigned int n1 = xl1.size();
  const unsigned int n2 = xl2.size();

  // Remember the starting point.
  const int is1 = ip1;
  std::vector<double> xpl;
  std::vector<double> ypl;
  // Add the starting point.
  xpl.push_back(xl1[ip1]);
  ypl.push_back(yl1[ip1]);
  mark1[ip1] = true;
  if (m_debug) {
    std::cout << "      Start from point " << ip1 << " on curve 1.\n";
  }
  // Next point.
  ip1 = NextPoint(ip1, n1);
  xpl.push_back(xl1[ip1]);
  ypl.push_back(yl1[ip1]);
  mark1[ip1] = true;
  if (m_debug) {
    std::cout << "      Next point is " << ip1 << " on curve 1.\n";
  }

  // Keep track of the curve we are currently following.
  unsigned int il = 1;
  // Direction flag (-1: backward, 0: not set, 1: forward).
  int dir = 0;
  // End-of-curve flag.
  bool eoc = false;
  int ip2 = 0;
  while (!eoc) {
    if (il == 1 && flags1[std::max(ip1, 0)] == 1) {
      // On curve 1 and not on the edge of curve 2?
      ip1 = NextPoint(ip1, n1);
      if (ip1 == is1) {
        eoc = true;
        continue;
      }
      mark1[ip1] = true;
      xpl.push_back(xl1[ip1]);
      ypl.push_back(yl1[ip1]);
      if (m_debug) {
        std::cout << "      Went to point " << ip1 << " on curve 1.\n";
      }
    } else if (il == 1) {
      // On curve 1 and on the edge of curve 2?
      ip2 = links1[ip1];
      bool added = false;
      if (dir == +1 || dir == 0) {
        const double xm = 0.5 * (xl2[ip2] + xl2[NextPoint(ip2, n2)]);
        const double ym = 0.5 * (yl2[ip2] + yl2[NextPoint(ip2, n2)]);
        bool inside = false, edge = false;
        Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
        if (inside) {
          ip2 = NextPoint(ip2, n2);
          il = 2;
          dir = +1;
          ip1 = links2[ip2];
          if (ip1 == is1) {
            eoc = true;
            continue;
          } else if (ip1 >= 0) {
            mark1[ip1] = true;
          }
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          added = true;
          if (m_debug) {
            std::cout << "      Added point " << ip2 << " along 2 +.\n";
          }
        }
      }
      if (dir == -1 || dir == 0) {
        const double xm = 0.5 * (xl2[ip2] + xl2[PrevPoint(ip2, n2)]);
        const double ym = 0.5 * (yl2[ip2] + yl2[PrevPoint(ip2, n2)]);
        bool inside = false, edge = false;
        Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
        if (inside) {
          ip2 = PrevPoint(ip2, n2);
          il = 2;
          dir = -1;
          ip1 = links2[ip2];
          if (ip1 == is1) {
            eoc = true;
            continue;
          } else if (ip1 >= 0) {
            mark1[ip1] = true;
          }
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          added = true;
          if (m_debug) {
            std::cout << "      Added point " << ip2 << " along 2 -\n";
          }
        }
      }
      if (!added) {
        ip1 = NextPoint(ip1, n1);
        if (ip1 == is1) {
          eoc = true;
          continue;
        } else if (ip1 >= 0) {
          mark1[ip1] = true;
        }
        xpl.push_back(xl1[ip1]);
        ypl.push_back(yl1[ip1]);
        if (m_debug) std::cout << "      Continued over 1.\n";
      }
    } else if (il == 2 && flags2[std::max(ip2, 0)] == 1) {
      // On curve 2 normal vertex (outside 1 hopefully).
      ip2 = dir > 0 ? NextPoint(ip2, n2) : PrevPoint(ip2, n2);
      ip1 = links2[ip2];
      if (ip1 == is1) {
        eoc = true;
        continue;
      } else if (ip1 >= 0) {
        mark1[ip1] = true;
      }
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
      if (m_debug) {
        std::cout << "      Went to point " << ip2 << " on 2.\n";
      }
    } else if (il == 2) {
      // On curve 2 and on edge of 1.
      ip1 = links2[ip2];
      ip1 = NextPoint(ip1, n1);
      il = 1;
      if (ip1 == is1) {
        eoc = true;
        continue;
      }
      xpl.push_back(xl1[ip1]);
      ypl.push_back(yl1[ip1]);
      if (m_debug) std::cout << "      Resumed 1 at point " << ip1 << ".\n";
    } else {
      // Other cases should not occur.
      std::cerr << m_className << "::TraceNonOverlap: Unexpected case.\n";
      return;
    }
  }

  Panel newPanel = panel1;
  newPanel.xv = xpl;
  newPanel.yv = ypl;
  panelsOut.push_back(std::move(newPanel));
  if (m_debug) {
    std::cout << "      End of curve reached, " << xpl.size() << " points.\n";
  }
}

void ComponentNeBem3d::TraceOverlap(
    const std::vector<double>& xp1, const std::vector<double>& yp1,
    const std::vector<double>& xp2, const std::vector<double>& yp2,
    const std::vector<double>& xl1, const std::vector<double>& yl1,
    const std::vector<double>& xl2, const std::vector<double>& yl2,
    const std::vector<int>& flags1, const std::vector<int>& links1,
    const std::vector<int>& links2, std::vector<bool>& mark1, int ip1, int ip2,
    const Panel& panel1, std::vector<Panel>& panelsOut) const {
  int ip1L = -1;
  int ip1LL = -1;

  const unsigned int n1 = xl1.size();
  const unsigned int n2 = xl2.size();

  // Remember the starting points.
  const int is1 = ip1;
  const int is2 = ip2;
  std::vector<double> xpl;
  std::vector<double> ypl;
  xpl.push_back(xl1[ip1]);
  ypl.push_back(yl1[ip1]);
  mark1[ip1] = true;

  // Keep track of the curve we are currently following.
  unsigned int il = 1;
  // Direction flag (-1: backward, 0: not set, 1: forward).
  int dir = 0;
  // End-of-curve flag.
  bool eoc = false;

  if (m_debug) {
    std::cout << "      Start from point " << ip1 << " on curve " << il << "\n";
  }
  while (!eoc) {
    ip1LL = ip1L;
    ip1L = ip1;
    if (il == 1) {
      // We are on curve 1. Move to the next point.
      const int ii = NextPoint(ip1, n1);
      // Maybe finished over line 1?
      if (ii == is1) {
        eoc = true;
        continue;
      }
      // See whether the next point of 1 is on the edge or inside of 2.
      bool inside = false, edge = false;
      if (links1[ii] >= 0) {
        edge = true;
      } else if (flags1[ii] == 1) {
        Polygon::Inside(xp2, yp2, xl1[ii], yl1[ii], inside, edge);
      }
      // If it is, check that it doesn't leave 2 at any stage.
      if (inside || edge) {
        const double xm = 0.5 * (xl1[ip1] + xl1[ii]);
        const double ym = 0.5 * (yl1[ip1] + yl1[ii]);
        Polygon::Inside(xp2, yp2, xm, ym, inside, edge);
      }
      // If it is, continue over 1.
      if (inside || edge) {
        ip1 = ii;
        if (m_debug) {
          std::cout << "      Continued to point " << ip1 << " on " << il
                    << "\n";
        }
        xpl.push_back(xl1[ip1]);
        ypl.push_back(yl1[ip1]);
        mark1[ip1] = true;
        continue;
      }
      // Else we have to continue over 2, ensure we really are on curve 2.
      ip2 = links1[ip1];
      if (ip2 < 0) {
        std::cerr << m_className << "::TraceOverlap: "
                  << "No point 2 reference found; abandoned.\n";
        return;
      }
      // Impose a direction on 2 to avoid returning.
      if (dir == 0) {
        if (links2[NextPoint(ip2, n2)] == ip1LL &&
            links2[PrevPoint(ip2, n2)] == ip1LL) {
          std::cerr << m_className << "::TraceOverlap: "
                    << "Both 2+ and 2- return on 1; not stored.\n";
          return;
        } else if (links2[NextPoint(ip2, n2)] == ip1LL) {
          if (m_debug) {
            std::cout << "      2+ is a return to previous point on 1.\n";
          }
          dir = -1;
        } else if (links2[PrevPoint(ip2, n2)] == ip1LL) {
          if (m_debug) {
            std::cout << "      2- is a return to previous point on 1.\n";
          }
          dir = +1;
        } else {
          if (m_debug) std::cout << "      Both ways are OK.\n";
        }
      }
      // If not, try to continue over 2 in the + direction.
      if (dir == +1 || dir == 0) {
        ip2 = NextPoint(ip2, n2);
        if (ip2 == is2) {
          if (m_debug) std::cout << "      Return to start over 2+.\n";
          eoc = true;
          continue;
        }
        Polygon::Inside(xp1, yp1, xl2[ip2], yl2[ip2], inside, edge);
        if (inside || edge) {
          if (m_debug) {
            std::cout << "      Going to 2+ (point " << ip2 << " of 2).\n";
          }
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          dir = +1;
          if (links2[ip2] >= 0) {
            ip1 = links2[ip2];
            mark1[ip1] = true;
            il = 1;
            if (m_debug) {
              std::cout << "      This point is also on curve 1: " << ip1
                        << "\n";
            }
          } else {
            il = 2;
          }
          continue;
        }
        // Continuing in the + direction didn't work so go back a step.
        ip2 = PrevPoint(ip2, n2);
      }
      // Or if this still fails, try 2 in the - direction.
      if (dir == -1 || dir == 0) {
        ip2 = PrevPoint(ip2, n2);
        if (ip2 == is2) {
          if (m_debug) std::cout << "      Return to start over 2-\n";
          eoc = true;
          continue;
        }
        Polygon::Inside(xp1, yp1, xl2[ip2], yl2[ip2], inside, edge);
        if (inside || edge) {
          if (m_debug) {
            std::cout << "      Going to 2- (point " << ip2 << " of 2).\n";
          }
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          dir = -1;
          if (links2[ip2] >= 0) {
            ip1 = links2[ip2];
            mark1[ip1] = true;
            il = 1;
            if (m_debug) {
              std::cout << "      This point is also on 1: " << ip1 << ".\n";
            }
          } else {
            il = 2;
          }
          continue;
        }
      }
      // Should not get here.
      if (m_debug) std::cout << "      Dead end.\n";
      return;
    } else if (il == 2) {
      // We are on curve 2. Ensure the direction is set.
      if (dir == 0) {
        std::cerr << m_className << "::TraceOverlap: "
                  << "Direction not set; abandoned.\n";
        return;
      }
      // Move to the next point.
      ip2 = dir > 0 ? NextPoint(ip2, n2) : PrevPoint(ip2, n2);
      // Maybe finished over line 2?
      if (ip2 == is2) {
        // Reached the end.
        eoc = true;
        continue;
      }
      // Next step over 2.
      if (m_debug) {
        std::cout << "      Stepped over 2 to point " << ip2 << " of 2.\n";
      }
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
      if (links2[ip2] >= 0) {
        ip1 = links2[ip2];
        mark1[ip1] = true;
        il = 1;
        if (m_debug) {
          std::cout << "      This point is also on curve 1: " << ip1 << ".\n";
        }
      } else {
        il = 2;
      }
    }
  }

  if (xpl.size() <= 2) {
    if (m_debug) std::cout << "      Too few points.\n";
  } else {
    Panel newPanel = panel1;
    newPanel.xv = xpl;
    newPanel.yv = ypl;
    panelsOut.push_back(std::move(newPanel));
  }

  if (m_debug) {
    std::cout << "      End of curve reached, " << xpl.size() << " points.\n";
  }
}

bool ComponentNeBem3d::MakePrimitives(const Panel& panelIn,
                                      std::vector<Panel>& panelsOut) const {
  // *-----------------------------------------------------------------------
  // *   PLATRC - Cuts a polygon into right-angled triangles.
  // *-----------------------------------------------------------------------

  // Establish tolerances.
  // TODO! Class member?
  const double epsang = 1.e-6;  // BEMEPA

  if (panelIn.xv.empty() || panelIn.yv.empty() || panelIn.zv.empty()) {
    return false;
  }
  // Determine the mean z value.
  const double zsum =
      std::accumulate(std::begin(panelIn.zv), std::end(panelIn.zv), 0.);
  const double zmean = zsum / panelIn.zv.size();

  std::vector<Panel> stack;
  stack.push_back(panelIn);
  stack.back().zv.clear();
  for (unsigned int k = 0; k < stack.size(); ++k) {
    // Next polygon.
    const auto& xp1 = stack[k].xv;
    const auto& yp1 = stack[k].yv;
    const unsigned int np = xp1.size();
    if (m_debug) {
      std::cout << "    Polygon " << k << " with " << np << " nodes.\n";
    }
    if (np <= 2) {
      // Too few nodes.
      if (m_debug) std::cout << "      Too few points.\n";
      continue;
    }

    // See whether this is a right-angled triangle.
    if (np == 3) {
      const double x12 = xp1[0] - xp1[1];
      const double y12 = yp1[0] - yp1[1];
      const double x23 = xp1[1] - xp1[2];
      const double y23 = yp1[1] - yp1[2];
      const double x31 = xp1[2] - xp1[0];
      const double y31 = yp1[2] - yp1[0];
      const double x32 = xp1[2] - xp1[1];
      const double y32 = yp1[2] - yp1[1];
      const double x13 = xp1[0] - xp1[2];
      const double y13 = yp1[0] - yp1[2];
      const double x21 = xp1[1] - xp1[0];
      const double y21 = yp1[1] - yp1[0];
      const double s12 = x12 * x12 + y12 * y12;
      const double s32 = x32 * x32 + y32 * y32;
      const double s13 = x13 * x13 + y13 * y13;
      const double s23 = x23 * x23 + y23 * y23;
      const double s31 = x31 * x31 + y31 * y31;
      const double s21 = x21 * x21 + y21 * y21;
      if (fabs(x12 * x32 + y12 * y32) < epsang * sqrt(s12 * s32)) {
        if (m_debug) {
          std::cout << "      Right-angled triangle node 2 - done.\n";
        }
        panelsOut.push_back(stack[k]);
        continue;
      } else if (fabs(x13 * x23 + y13 * y23) < epsang * sqrt(s13 * s23)) {
        if (m_debug) {
          std::cout << "      Right-angled triangle node 3 - rearrange.\n";
        }
        Panel panel = stack[k];
        panel.xv = {xp1[1], xp1[2], xp1[0]};
        panel.yv = {yp1[1], yp1[2], yp1[0]};
        panelsOut.push_back(std::move(panel));
        continue;
      } else if (fabs(x31 * x21 + y31 * y21) < epsang * sqrt(s31 * s21)) {
        if (m_debug) {
          std::cout << "      Right-angled triangle node 1 - rearrange.\n";
        }
        Panel panel = stack[k];
        panel.xv = {xp1[2], xp1[0], xp1[1]};
        panel.yv = {yp1[2], yp1[0], yp1[1]};
        panelsOut.push_back(std::move(panel));
        continue;
      }
    }
    // See whether this is a rectangle.
    if (np == 4) {
      const double x12 = xp1[0] - xp1[1];
      const double y12 = yp1[0] - yp1[1];
      const double x23 = xp1[1] - xp1[2];
      const double y23 = yp1[1] - yp1[2];
      const double x34 = xp1[2] - xp1[3];
      const double y34 = yp1[2] - yp1[3];
      const double x43 = xp1[3] - xp1[2];
      const double y43 = yp1[3] - yp1[2];
      const double x14 = xp1[0] - xp1[3];
      const double y14 = yp1[0] - yp1[3];
      const double x32 = xp1[2] - xp1[1];
      const double y32 = yp1[2] - yp1[1];

      const double s12 = x12 * x12 + y12 * y12;
      const double s23 = x23 * x23 + y23 * y23;
      const double s32 = x32 * x32 + y32 * y32;
      const double s34 = x34 * x34 + y34 * y34;
      const double s43 = x43 * x43 + y43 * y43;
      const double s14 = x14 * x14 + y14 * y14;
      if (fabs(x12 * x32 + y12 * y32) < epsang * sqrt(s12 * s32) &&
          fabs(x23 * x43 + y23 * y43) < epsang * sqrt(s23 * s43) &&
          fabs(x14 * x34 + y14 * y34) < epsang * sqrt(s14 * s34)) {
        if (m_debug) std::cout << "      Rectangle.\n";
        panelsOut.push_back(stack[k]);
        continue;
      }
    }
    // See whether there are parallel sides, e.g. a trapezium (UK English).
    if (np >= 4 && SplitTrapezium(stack[k], stack, panelsOut, epsang)) continue;

    // Find a right-angled corner we can cut off.
    if (m_debug) std::cout << "      Trying to find a right-angle\n";
    bool corner = false;
    for (unsigned int ip = 0; ip < np; ++ip) {
      // Take only right angles.
      const unsigned int inext = NextPoint(ip, np);
      const unsigned int iprev = PrevPoint(ip, np);

      const double dxprev = xp1[iprev] - xp1[ip];
      const double dyprev = yp1[iprev] - yp1[ip];
      const double dxnext = xp1[inext] - xp1[ip];
      const double dynext = yp1[inext] - yp1[ip];
      if (fabs(dxprev * dxnext + dyprev * dynext) >
          epsang * sqrt((dxprev * dxprev + dyprev * dyprev) *
                        (dxnext * dxnext + dynext * dynext))) {
        continue;
      }
      // Ensure the midpoint is internal.
      if (np > 3) {
        const double xm = 0.5 * (xp1[iprev] + xp1[inext]);
        const double ym = 0.5 * (yp1[iprev] + yp1[inext]);
        bool inside = false, edge = false;
        Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
        if (!inside) continue;
      }
      // Check all vertex crossings.
      bool cross = false;
      for (unsigned int jp = 0; jp < np; ++jp) {
        // Accept immediate contact.
        const unsigned int jnext = NextPoint(jp, np);
        if (jp == iprev || jp == ip || jp == inext || jnext == iprev ||
            jnext == ip || jnext == inext)
          continue;
        // Check crossing.
        double xc = 0., yc = 0.;
        if (Crossing(xp1[iprev], yp1[iprev], xp1[inext], yp1[inext], xp1[jp],
                     yp1[jp], xp1[jnext], yp1[jnext], xc, yc)) {
          cross = true;
          break;
        }
      }
      if (cross) continue;
      // Found a triangle, introduce shorthand node references.
      if (m_debug) {
        std::cout << "      Cutting at right-angled corner " << ip << "\n";
      }
      corner = true;
      Panel panel = stack[k];
      panel.xv = {xp1[iprev], xp1[ip], xp1[inext]};
      panel.yv = {yp1[iprev], yp1[ip], yp1[inext]};
      stack.push_back(std::move(panel));
      // Eliminate this node from the polygon.
      stack[k].xv.erase(stack[k].xv.begin() + ip);
      stack[k].yv.erase(stack[k].yv.begin() + ip);
      stack.push_back(std::move(stack[k]));
      if (m_debug) {
        std::cout << "      Going for another pass, NP = " << np << "\n";
      }
      break;
    }
    if (corner) continue;

    // Find any corner we can cut off.
    if (m_debug) std::cout << "      Trying to find a corner\n";
    corner = false;
    for (unsigned int ip = 0; ip < np; ++ip) {  // 20
      const unsigned int iprev = PrevPoint(ip, np);
      const unsigned int inext = NextPoint(ip, np);
      // Ensure the midpoint is internal.
      if (np > 3) {
        const double xm = 0.5 * (xp1[iprev] + xp1[inext]);
        const double ym = 0.5 * (yp1[iprev] + yp1[inext]);
        bool inside = false, edge = false;
        Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
        if (!inside) continue;
      }
      // Check all vertex crossings.
      bool cross = false;
      for (unsigned int jp = 0; jp < np; ++jp) {
        const unsigned int jj = NextPoint(jp, np);
        // Accept immediate contact.
        if (jp == iprev || jp == ip || jp == inext || jj == iprev || jj == ip ||
            jj == inext)
          continue;
        // Check crossing.
        double xc = 0., yc = 0.;
        if (Crossing(xp1[iprev], yp1[iprev], xp1[inext], yp1[inext], xp1[jp],
                     yp1[jp], xp1[jj], yp1[jj], xc, yc)) {
          cross = true;
          break;
        }
      }
      if (cross) continue;
      // Found a triangle, introduce shorthand node references.
      if (m_debug) std::cout << "      Cutting at corner " << ip << "\n";
      corner = true;
      const double x1 = xp1[iprev];
      const double x2 = xp1[ip];
      const double x3 = xp1[inext];
      const double y1 = yp1[iprev];
      const double y2 = yp1[ip];
      const double y3 = yp1[inext];
      const double s21 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
      const double s31 = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1);
      const double s32 = (x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2);
      const double s12 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
      const double s13 = (x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3);
      const double s23 = (x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3);
      // Find the biggest opening angle.
      const double a1 =
          ((x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1)) / sqrt(s21 * s31);
      const double a2 =
          ((x3 - x2) * (x1 - x2) + (y3 - y2) * (y1 - y2)) / sqrt(s32 * s12);
      const double a3 =
          ((x1 - x3) * (x2 - x3) + (y1 - y3) * (y2 - y3)) / sqrt(s13 * s23);
      if (m_debug) {
        const double phi1 = acos(a1);
        const double phi2 = acos(a2);
        const double phi3 = acos(a3);
        std::cout << "      Angles: " << RadToDegree * phi1 << ", "
                  << RadToDegree * phi2 << ", " << RadToDegree * phi3 << "\n";
        const double sumphi = phi1 + phi2 + phi3;
        std::cout << "      Sum = " << RadToDegree * sumphi << "\n";
      }
      // See whether one angle is more or less right-angled.
      if (fabs(a1) < epsang || fabs(a2) < epsang || fabs(a3) < epsang) {
        if (m_debug) std::cout << "      Right-angled corner cut off.\n";
        Panel panel = stack[k];
        if (fabs(a1) < epsang) {
          panel.xv = {x3, x1, x2};
          panel.yv = {y3, y1, y2};
        } else if (fabs(a2) < epsang) {
          panel.xv = {x1, x2, x3};
          panel.yv = {y1, y2, y3};
        } else {
          panel.xv = {x2, x3, x1};
          panel.yv = {y2, y3, y1};
        }
        stack.push_back(std::move(panel));
      } else if (a1 <= a2 && a1 <= a3) {
        if (m_debug) std::cout << "      A1 < A2, A3 - adding 2 triangles.\n";
        const double xc = x2 + a2 * (x3 - x2) * sqrt(s12 / s32);
        const double yc = y2 + a2 * (y3 - y2) * sqrt(s12 / s32);
        Panel panel1 = stack[k];
        panel1.xv = {x3, xc, x1};
        panel1.yv = {y3, yc, y1};
        stack.push_back(std::move(panel1));
        Panel panel2 = stack[k];
        panel2.xv = {x2, xc, x1};
        panel2.yv = {y2, yc, y1};
        stack.push_back(std::move(panel2));
      } else if (a2 <= a1 && a2 <= a3) {
        if (m_debug) std::cout << "      A2 < A1, A3 - adding 2 triangles.\n";
        const double xc = x3 + a3 * (x1 - x3) * sqrt(s23 / s13);
        const double yc = y3 + a3 * (y1 - y3) * sqrt(s23 / s13);
        Panel panel1 = stack[k];
        panel1.xv = {x1, xc, x2};
        panel1.yv = {y1, yc, y2};
        stack.push_back(std::move(panel1));
        Panel panel2 = stack[k];
        panel2.xv = {x3, xc, x2};
        panel2.yv = {y3, yc, y2};
        stack.push_back(std::move(panel2));
      } else {
        if (m_debug) std::cout << "      A3 < A1, A2 - adding 2 triangles.\n";
        const double xc = x1 + a1 * (x2 - x1) * sqrt(s31 / s21);
        const double yc = y1 + a1 * (y2 - y1) * sqrt(s31 / s21);
        Panel panel1 = stack[k];
        panel1.xv = {x1, xc, x3};
        panel1.yv = {y1, yc, y3};
        stack.push_back(std::move(panel1));
        Panel panel2 = stack[k];
        panel2.xv = {x2, xc, x3};
        panel2.yv = {y2, yc, y3};
        stack.push_back(std::move(panel2));
      }
      // Eliminate this node from the polygon.
      stack[k].xv.erase(stack[k].xv.begin() + ip);
      stack[k].yv.erase(stack[k].yv.begin() + ip);
      stack.push_back(std::move(stack[k]));
      if (m_debug) {
        std::cout << "      Going for another pass, NP = " << np << "\n";
      }
      break;
    }
    if (corner) continue;
    std::cerr << m_className << "::MakePrimitives:\n    "
              << "Unable to identify a corner to cut, probably a degenerate "
                 "polygon.\n";
    // Next stack element.
  }
  for (auto& panel : panelsOut) {
    panel.zv.assign(panel.xv.size(), zmean);
  }
  return true;
}

bool ComponentNeBem3d::SplitTrapezium(const Panel panelIn,
                                      std::vector<Panel>& stack,
                                      std::vector<Panel>& panelsOut,
                                      const double epsang) const {
  const auto xp1 = panelIn.xv;
  const auto yp1 = panelIn.yv;
  const unsigned int np = xp1.size();
  for (unsigned int ip = 0; ip < np; ++ip) {
    const unsigned int inext = NextPoint(ip, np);
    const double xi0 = xp1[ip];
    const double yi0 = yp1[ip];
    const double xi1 = xp1[inext];
    const double yi1 = yp1[inext];
    const double dxi = xi0 - xi1;
    const double dyi = yi0 - yi1;
    const double si2 = dxi * dxi + dyi * dyi;
    for (unsigned int jp = ip + 2; jp < np; ++jp) {
      const unsigned int jnext = NextPoint(jp, np);
      // Skip adjacent segments.
      if (ip == jp || ip == jnext || inext == jp || inext == jnext) {
        continue;
      }
      const double xj0 = xp1[jp];
      const double yj0 = yp1[jp];
      const double xj1 = xp1[jnext];
      const double yj1 = yp1[jnext];
      const double dxj = xj0 - xj1;
      const double dyj = yj0 - yj1;
      const double sj2 = dxj * dxj + dyj * dyj;
      // Require parallelism.
      const double sij = sqrt(si2 * sj2);
      if (fabs(dxi * dxj + dyi * dyj + sij) > epsang * sij) continue;
      if (m_debug) {
        std::cout << "      Found parallel sections: " << ip << ", " << jp
                  << "\n";
      }
      // Avoid division by zero.
      if (sj2 <= 0 || si2 <= 0) {
        std::cerr << m_className << "::SplitTrapezium:\n    "
                  << "Zero norm segment found; skipped.\n";
        continue;
      }
      // Establish the cutting lines.
      const double xl1 =
          ((xi0 - xj0) * (xj1 - xj0) + (yi0 - yj0) * (yj1 - yj0)) / sj2;
      const double xl2 =
          ((xi1 - xj0) * (xj1 - xj0) + (yi1 - yj0) * (yj1 - yj0)) / sj2;
      const double xl3 =
          ((xj0 - xi0) * (xi1 - xi0) + (yj0 - yi0) * (yi1 - yi0)) / si2;
      const double xl4 =
          ((xj1 - xi0) * (xi1 - xi0) + (yj1 - yi0) * (yi1 - yi0)) / si2;
      if (m_debug) {
        std::cout << "      xl1 = " << xl1 << ", xl2 = " << xl2 << ", "
                  << "xl3 = " << xl3 << ", xl4 = " << xl4 << "\n";
      }
      // Check that there is at all a rectangle.
      const double r1 = (xl1 + epsang) * (1. + epsang - xl1);
      const double r2 = (xl2 + epsang) * (1. + epsang - xl2);
      const double r3 = (xl3 + epsang) * (1. + epsang - xl3);
      const double r4 = (xl4 + epsang) * (1. + epsang - xl4);
      if ((r1 < 0 && r4 < 0) || (r2 < 0 && r3 < 0)) {
        if (m_debug) std::cout << "      No rectangle.\n";
        continue;
      }
      // Determine the rectangular part.
      std::vector<double> xpl(4, 0.);
      std::vector<double> ypl(4, 0.);
      if (r1 >= 0) {
        xpl[0] = xi0;
        ypl[0] = yi0;
        xpl[1] = xj0 + xl1 * (xj1 - xj0);
        ypl[1] = yj0 + xl1 * (yj1 - yj0);
      } else if (r4 >= 0) {
        xpl[0] = xi0 + xl4 * (xi1 - xi0);
        ypl[0] = yi0 + xl4 * (yi1 - yi0);
        xpl[1] = xj1;
        ypl[1] = yj1;
      }
      if (r2 >= 0) {
        xpl[2] = xj0 + xl2 * (xj1 - xj0);
        ypl[2] = yj0 + xl2 * (yj1 - yj0);
        xpl[3] = xi1;
        ypl[3] = yi1;
      } else if (r3 >= 0) {
        xpl[2] = xj0;
        ypl[2] = yj0;
        xpl[3] = xi0 + xl3 * (xi1 - xi0);
        ypl[3] = yi0 + xl3 * (yi1 - yi0);
      }
      // Verify that the midpoints of these lines are internal.
      double xm = 0.5 * (xpl[0] + xpl[1]);
      double ym = 0.5 * (ypl[0] + ypl[1]);
      bool inside = false, edge = false;
      Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
      if (!(inside || edge)) {
        if (m_debug) std::cout << "      Midpoint 1 not internal.\n";
        continue;
      }
      xm = 0.5 * (xpl[2] + xpl[3]);
      ym = 0.5 * (ypl[2] + ypl[3]);
      Polygon::Inside(xp1, yp1, xm, ym, inside, edge);
      if (!(inside || edge)) {
        if (m_debug) std::cout << "      Midpoint 2 not internal.\n";
        continue;
      }

      const unsigned int iprev = PrevPoint(ip, np);
      const unsigned int jprev = PrevPoint(jp, np);
      // Ensure there are no crossings, accepting contact.
      bool cross = false;
      for (unsigned int i = 0; i < np; ++i) {
        if ((i == iprev && r1 >= 0) || i == ip || (i == inext && r2 >= 0) ||
            (i == jprev && r3 >= 0) || i == jp || (i == jnext && r4 >= 0))
          continue;
        const unsigned int ii = NextPoint(i, np);
        double xc = 0., yc = 0.;
        if (Crossing(xp1[i], yp1[i], xp1[ii], yp1[ii], xpl[0], ypl[0], xpl[1],
                     ypl[1], xc, yc) ||
            Crossing(xp1[i], yp1[i], xp1[ii], yp1[ii], xpl[2], ypl[2], xpl[3],
                     ypl[3], xc, yc)) {
          if (m_debug) {
            std::cout << "      Crossing (edge " << i << ", " << ii
                      << "), ip = " << ip << ", jp = " << jp << "\n";
          }
          cross = true;
          break;
        }
      }
      if (cross) continue;
      // Add the rectangular part.
      if ((fabs(xl1) < epsang && fabs(xl3) < epsang) ||
          (fabs(1. - xl2) < epsang && fabs(1. - xl4) < epsang)) {
        if (m_debug) std::cout << "      Not stored, degenerate.\n";
      } else {
        if (m_debug) std::cout << "      Adding rectangle.\n";
        Panel panel = panelIn;
        panel.xv = xpl;
        panel.yv = ypl;
        panelsOut.push_back(std::move(panel));
      }
      // First non-rectangular section.
      xpl.clear();
      ypl.clear();
      if (m_debug) std::cout << "      First non-rectangular section.\n";
      for (unsigned int i = jp + 1; i <= ip + np; ++i) {
        const unsigned int ii = i % np;
        xpl.push_back(xp1[ii]);
        ypl.push_back(yp1[ii]);
      }
      if (r1 >= 0 && r4 >= 0) {
        if (m_debug) std::cout << "      1-4 degenerate\n";
      } else if (r1 >= 0) {
        if (m_debug) std::cout << "      Using 1\n";
        xpl.push_back(xj0 + xl1 * (xj1 - xj0));
        ypl.push_back(yj0 + xl1 * (yj1 - yj0));
      } else if (r4 >= 0) {
        if (m_debug) std::cout << "      Using 4\n";
        xpl.push_back(xi0 + xl4 * (xi1 - xi0));
        ypl.push_back(yi0 + xl4 * (yi1 - yi0));
      } else {
        if (m_debug) std::cout << "      Neither 1 nor 4, should not happen\n";
      }
      if (xpl.size() < 3) {
        if (m_debug) {
          std::cout << "      Not stored, only " << xpl.size()
                    << " vertices.\n";
        }
      } else {
        if (m_debug) {
          std::cout << "      Adding non-rectangular part with " << xpl.size()
                    << " vertices.\n";
        }
        Panel panel = panelIn;
        panel.xv = xpl;
        panel.yv = ypl;
        stack.push_back(std::move(panel));
      }
      // Second non-rectangular section.
      xpl.clear();
      ypl.clear();
      if (m_debug) std::cout << "      Second non-rectangular section.\n";
      for (unsigned int i = ip + 1; i <= jp; ++i) {
        const unsigned int ii = i % np;
        xpl.push_back(xp1[ii]);
        ypl.push_back(yp1[ii]);
      }
      if (r2 >= 0 && r3 >= 0) {
        if (m_debug) std::cout << "      2-3 degenerate\n";
      } else if (r2 >= 0) {
        if (m_debug) std::cout << "      Using 2\n";
        xpl.push_back(xj0 + xl2 * (xj1 - xj0));
        ypl.push_back(yj0 + xl2 * (yj1 - yj0));
      } else if (r3 >= 0) {
        if (m_debug) std::cout << "      Using 3\n";
        xpl.push_back(xi0 + xl3 * (xi1 - xi0));
        ypl.push_back(yi0 + xl3 * (yi1 - yi0));
      } else {
        if (m_debug) std::cout << "      Neither 2 nor 3, should not happen\n";
      }
      if (xpl.size() < 3) {
        if (m_debug) {
          std::cout << "      Not stored, only " << xpl.size()
                    << " vertices.\n";
        }
      } else {
        if (m_debug) {
          std::cout << "      Adding non-rectangular part with " << xpl.size()
                    << " vertices.\n";
        }
        Panel panel = panelIn;
        panel.xv = xpl;
        panel.yv = ypl;
        stack.push_back(std::move(panel));
      }
      return true;
    }
  }
  return false;
}

bool ComponentNeBem3d::GetPrimitive(const unsigned int i, double& a, double& b,
                                    double& c, std::vector<double>& xv,
                                    std::vector<double>& yv,
                                    std::vector<double>& zv, int& interface,
                                    double& v, double& q,
                                    double& lambda) const {
  if (i >= m_primitives.size()) {
    std::cerr << m_className << "::GetPrimitive: Index out of range.\n";
    return false;
  }
  const auto& primitive = m_primitives[i];
  a = primitive.a;
  b = primitive.b;
  c = primitive.c;
  xv = primitive.xv;
  yv = primitive.yv;
  zv = primitive.zv;
  interface = primitive.interface;
  v = primitive.v;
  q = primitive.q;
  lambda = primitive.lambda;
  return true;
}

bool ComponentNeBem3d::GetPrimitive(const unsigned int i, double& a, double& b,
                                    double& c, std::vector<double>& xv,
                                    std::vector<double>& yv,
                                    std::vector<double>& zv, int& vol1,
                                    int& vol2) const {
  if (i >= m_primitives.size()) {
    std::cerr << m_className << "::GetPrimitive: Index out of range.\n";
    return false;
  }
  const auto& primitive = m_primitives[i];
  a = primitive.a;
  b = primitive.b;
  c = primitive.c;
  xv = primitive.xv;
  yv = primitive.yv;
  zv = primitive.zv;
  vol1 = primitive.vol1;
  vol2 = primitive.vol2;
  return true;
}

bool ComponentNeBem3d::GetVolume(const unsigned int vol, int& shape,
                                 int& material, double& epsilon,
                                 double& potential, double& charge, int& bc) {
  if (!m_geometry) return false;
  const unsigned int nSolids = m_geometry->GetNumberOfSolids();
  for (unsigned int i = 0; i < nSolids; ++i) {
    Medium* medium = nullptr;
    const auto solid = m_geometry->GetSolid(i, medium);
    if (!solid) continue;
    if (solid->GetId() != vol) continue;
    if (solid->IsTube() || solid->IsWire()) {
      shape = 1;
    } else if (solid->IsHole()) {
      shape = 2;
    } else if (solid->IsBox()) {
      shape = 3;
    } else if (solid->IsSphere()) {
      shape = 4;
    } else if (solid->IsRidge()) {
      shape = 5;
    } else if (solid->IsExtrusion()) {
      shape = 6;
    } else {
      std::cerr << m_className << "::GetVolume: Unknown solid shape.\n";
      return false;
    }
    material = medium ? medium->GetId() : 11;
    epsilon = medium ? medium->GetDielectricConstant() : 1.;
    potential = solid->GetBoundaryPotential();
    charge = solid->GetBoundaryChargeDensity();
    bc = solid->GetBoundaryConditionType();
    return true;
  }
  return false;
}

int ComponentNeBem3d::GetVolume(const double x, const double y,
                                const double z) {
  if (!m_geometry) return -1;
  const unsigned int nSolids = m_geometry->GetNumberOfSolids();
  for (unsigned int i = 0; i < nSolids; ++i) {
    Medium* medium = nullptr;
    const auto solid = m_geometry->GetSolid(i, medium);
    if (!solid) continue;
    if (solid->IsInside(x, y, z)) return solid->GetId();
  }
  return -1;
}

bool ComponentNeBem3d::GetElement(const unsigned int i, std::vector<double>& xv,
                                  std::vector<double>& yv,
                                  std::vector<double>& zv, int& interface,
                                  double& bc, double& lambda) const {
  if (i >= m_elements.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }
  const auto& element = m_elements[i];
  xv = element.xv;
  yv = element.yv;
  zv = element.zv;
  interface = element.interface;
  bc = element.bc;
  lambda = element.lambda;
  return true;
}

void ComponentNeBem3d::Reset() {
  m_primitives.clear();
  m_elements.clear();
  m_ynplan.fill(false);
  m_coplan.fill(0.);
  m_vtplan.fill(0.);
  m_ready = false;
}

void ComponentNeBem3d::UpdatePeriodicity() {
  for (unsigned int i = 0; i < 3; ++i) {
    // Cannot have regular and mirror periodicity at the same time.
    if (m_periodic[i] && m_mirrorPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Both simple and mirror periodicity requested. Reset.\n";
      m_periodic[i] = false;
      m_mirrorPeriodic[i] = false;
      continue;
    }
    if ((m_periodic[i] || m_mirrorPeriodic[i]) && m_periodicLength[i] < Small) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Periodic length is not set. Reset.\n";
      m_periodic[i] = m_mirrorPeriodic[i] = false;
    }
  }

  if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Axial periodicity is not available.\n";
    m_axiallyPeriodic.fill(false);
  }

  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Rotation symmetry is not available.\n";
    m_rotationSymmetric.fill(false);
  }
}
}  // namespace Garfield
