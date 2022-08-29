#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <complex>

#include "Garfield/ComponentNeBem2d.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Polygon.hh"
#include "Garfield/Random.hh"

namespace {

double Mag2(const std::array<double, 2>& x) {
  return x[0] * x[0] + x[1] * x[1];
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
  // Initial values.
  xc = 0.;
  yc = 0.;
  /// Matrix to compute the crossing point.
  std::array<std::array<double, 2>, 2> a;
  a[0][0] = y2 - y1;
  a[0][1] = v2 - v1;
  a[1][0] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  // Set tolerances.
  const double epsx = std::max(1.e-10, 1.e-10 * std::max({fabs(x1), fabs(x2), 
                                                          fabs(u1), fabs(u2)}));
  const double epsy = std::max(1.e-10, 1.e-10 * std::max({fabs(y1), fabs(y2),
                                                          fabs(v1), fabs(v2)}));
  // Check if the lines are parallel.
  if (fabs(det) < epsx * epsy) return false;
 
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

/// Determine whether two polygons are overlapping.
bool Intersecting(const std::vector<double>& xp1,
                  const std::vector<double>& yp1,
                  const std::vector<double>& xp2,
                  const std::vector<double>& yp2) {
  
  const double xmin1 = *std::min_element(std::begin(xp1), std::end(xp1));
  const double ymin1 = *std::min_element(std::begin(yp1), std::end(yp1));
  const double xmax1 = *std::max_element(std::begin(xp1), std::end(xp1));
  const double ymax1 = *std::max_element(std::begin(yp1), std::end(yp1));

  const double xmin2 = *std::min_element(std::begin(xp2), std::end(xp2));
  const double ymin2 = *std::min_element(std::begin(yp2), std::end(yp2));
  const double xmax2 = *std::max_element(std::begin(xp2), std::end(xp2));
  const double ymax2 = *std::max_element(std::begin(yp2), std::end(yp2));

  const double epsx = 1.e-6 * std::max({std::abs(xmax1), std::abs(xmin1),
                                        std::abs(xmax2), std::abs(xmin2)});
  const double epsy = 1.e-6 * std::max({std::abs(ymax1), std::abs(ymin1),
                                        std::abs(ymax2), std::abs(ymin2)});
  // Check if the bounding boxes overlap.
  if (xmax1 + epsx < xmin2 || xmax2 + epsx < xmin1) return false;
  if (ymax1 + epsy < ymin2 || ymax2 + epsy < ymin1) return false;

  const unsigned int n1 = xp1.size();
  const unsigned int n2 = xp2.size();
  for (unsigned int i = 0; i < n1; ++i) {
    const double x0 = xp1[i];
    const double y0 = yp1[i];
    const unsigned int ii = i < n1 - 1 ? i + 1 : 0;
    const double x1 = xp1[ii];
    const double y1 = yp1[ii];
    for (unsigned int j = 0; j < n2; ++j) {
      const unsigned int jj = j < n2 - 1 ? j + 1 : 0;
      const double u0 = xp2[j];
      const double v0 = yp2[j];
      const double u1 = xp2[jj];
      const double v1 = yp2[jj];
      double xc = 0., yc = 0.;
      if (!Crossing(x0, y0, x1, y1, u0, v0, u1, v1, xc, yc)) continue;
      if ((OnLine(x0, y0, x1, y1, u0, v0) || OnLine(x0, y0, x1, y1, u1, v1)) &&
          (OnLine(u0, v0, u1, v1, x0, y0) || OnLine(u0, v0, u1, v1, x1, y1))) {
        continue;
      }
      return true;
    }
  }
  return false;
}

/// Determine whether polygon 1 is fully enclosed by polygon 2.
bool Enclosed(const std::vector<double>& xp1,
              const std::vector<double>& yp1,
              const std::vector<double>& xp2,
              const std::vector<double>& yp2) {
  
  const double xmin1 = *std::min_element(std::begin(xp1), std::end(xp1));
  const double ymin1 = *std::min_element(std::begin(yp1), std::end(yp1));
  const double xmax1 = *std::max_element(std::begin(xp1), std::end(xp1));
  const double ymax1 = *std::max_element(std::begin(yp1), std::end(yp1));

  const double xmin2 = *std::min_element(std::begin(xp2), std::end(xp2));
  const double ymin2 = *std::min_element(std::begin(yp2), std::end(yp2));
  const double xmax2 = *std::max_element(std::begin(xp2), std::end(xp2));
  const double ymax2 = *std::max_element(std::begin(yp2), std::end(yp2));

  const double epsx = 1.e-6 * std::max({std::abs(xmax1), std::abs(xmin1),
                                        std::abs(xmax2), std::abs(xmin2)});
  const double epsy = 1.e-6 * std::max({std::abs(ymax1), std::abs(ymin1),
                                        std::abs(ymax2), std::abs(ymin2)});
  // Check the bounding boxes.
  if (xmax1 + epsx < xmin2 || xmax2 + epsx < xmin1) return false;
  if (ymax1 + epsy < ymin2 || ymax2 + epsy < ymin1) return false;
  if (xmin1 + epsx < xmin2 || xmax1 > xmax2 + epsx) return false; 
  if (ymin1 + epsy < ymin2 || ymax1 > ymax2 + epsy) return false; 
 
  const unsigned int n1 = xp1.size();
  for (unsigned int i = 0; i < n1; ++i) {
    bool inside = false, edge = false;
    Garfield::Polygon::Inside(xp2, yp2, xp1[i], yp1[i], inside, edge);
    if (!inside) return false;
  }
  return true;
}

}

namespace Garfield {

const double ComponentNeBem2d::InvEpsilon0 = 1. / VacuumPermittivity;
const double ComponentNeBem2d::InvTwoPiEpsilon0 = 1. / TwoPiEpsilon0;

ComponentNeBem2d::ComponentNeBem2d() : Component("NeBem2d") {}

void ComponentNeBem2d::ElectricField(const double x, const double y,
                                     const double z, double& ex, double& ey,
                                     double& ez, double& v, Medium*& m,
                                     int& status) {
  status = Field(x, y, z, ex, ey, ez, v, m, true);
}

void ComponentNeBem2d::ElectricField(const double x, const double y,
                                     const double z, double& ex, double& ey,
                                     double& ez, Medium*& m, int& status) {
  double v = 0.;
  status = Field(x, y, z, ex, ey, ez, v, m, false);
}

int ComponentNeBem2d::Field(const double x, const double y, const double z,                            double& ex, double& ey, double& ez, double& v,
                            Medium*& m, const bool opt) {

  ex = ey = ez = 0.;
  // Check if the requested point is inside the z-range.
  if (m_useRangeZ && (z < m_zmin || z > m_zmax)) return -6;

  // Check if the requested point is inside a medium.
  m = GetMedium(x, y, z);
  if (!m) return -6;

  // Inside a conductor?
  if (m->IsConductor()) {
    if (!opt) return -5;
    // Find the potential.
    for (const auto& region : m_regions) {
      bool inside = false, edge = false;
      Garfield::Polygon::Inside(region.xv, region.yv, x, y, inside, edge);
      if (inside || edge) {
        v = region.bc.second;
        break;
      }
    }
    return -5;
  }

  if (!m_ready) {
    if (!Initialise()) {
      std::cerr << m_className << "::ElectricField: Initialisation failed.\n";
      return -11;
    }
  }
  
  // See whether we are inside a wire.
  const unsigned int nWires = m_wires.size();
  for (unsigned int i = 0; i < nWires; ++i) {
    const double dx = x - m_wires[i].x;
    const double dy = y - m_wires[i].y;
    if (dx * dx + dy * dy < m_wires[i].r * m_wires[i].r) {
      v = m_wires[i].v;
      return i + 1;
    }
  }

  // Sum up the contributions from all straight-line elements.
  for (const auto& element : m_elements) {
    const double cphi = element.cphi;
    const double sphi = element.sphi;
    // Transform to local coordinates.
    double xL = 0., yL = 0.;
    ToLocal(x - element.x, y - element.y, cphi, sphi, xL, yL);
    if (opt) {
      // Compute the potential.
      v += LinePotential(element.a, xL, yL) * element.q;
    }
    // Compute the field in local coordinates.
    double fx = 0., fy = 0.;
    LineField(element.a, xL, yL, fx, fy);
    // Rotate to the global frame.
    ToGlobal(fx, fy, cphi, sphi, fx, fy);
    ex += fx * element.q;
    ey += fy * element.q;
  }

  // Add the contributions from the wires.
  for (const auto& wire : m_wires) {
    if (opt) {
      // Compute the potential.
      v += WirePotential(wire.r, x - wire.x, y - wire.y) * wire.q;
    }
    // Compute the field.
    double fx = 0., fy = 0.;
    WireField(wire.x, x - wire.x, y - wire.y, fx, fy);
    ex += fx * wire.q;
    ey += fy * wire.q;
  }

  for (const auto& box : m_spaceCharge) {
    if (opt) {
      v += BoxPotential(box.a, box.b, x - box.x, y - box.y, box.v0) * box.q;
    }
    double fx = 0., fy = 0.;
    BoxField(box.a, box.b, x - box.x, y - box.y, fx, fy);
    ex += fx * box.q;
    ey += fy * box.q;
  }
  return 0;
}

bool ComponentNeBem2d::GetVoltageRange(double& vmin, double& vmax) {
  bool gotValue = false;
  for (const auto& region : m_regions) {
    if (region.bc.first != Voltage) continue;
    if (!gotValue) {
      vmin = vmax = region.bc.second;
      gotValue = true;
    } else {
      vmin = std::min(vmin, region.bc.second);
      vmax = std::max(vmax, region.bc.second);
    }
  }

  for (const auto& segment : m_segments) {
    if (segment.bc.first != Voltage) continue;
    if (!gotValue) {
      vmin = vmax = segment.bc.second;
      gotValue = true;
    } else {
      vmin = std::min(vmin, segment.bc.second);
      vmax = std::max(vmax, segment.bc.second);
    }
  }

  for (const auto& wire : m_wires) {
    if (!gotValue) {
      vmin = vmax = wire.v;
      gotValue = true;
    } else {
      vmin = std::min(vmin, wire.v);
      vmax = std::max(vmax, wire.v);
    }
  }
  return gotValue;
}

Medium* ComponentNeBem2d::GetMedium(const double x, const double y, 
                                    const double z) {

  if (m_geometry) return m_geometry->GetMedium(x, y, z);
  for (const auto& region : m_regions) {
    bool inside = false, edge = false;
    Garfield::Polygon::Inside(region.xv, region.yv, x, y, inside, edge);
    if (inside || edge) return region.medium;
  }
  return m_medium;
}

bool ComponentNeBem2d::GetBoundingBox(
    double& xmin, double& ymin, double& zmin,
    double& xmax, double& ymax, double& zmax) {

  if (m_geometry) {
    return m_geometry->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  return GetElementaryCell(xmin, ymin, zmin, xmax, ymax, zmax);
}

bool ComponentNeBem2d::GetElementaryCell(
    double& xmin, double& ymin, double& zmin,
    double& xmax, double& ymax, double& zmax) {

  zmin = m_zmin;
  zmax = m_zmax;
  bool gotValue = false;
  for (const auto& region : m_regions) {
    const auto& xv = region.xv;
    const auto& yv = region.yv;
    if (!gotValue) {
      xmin = *std::min_element(std::begin(xv), std::end(xv));
      ymin = *std::min_element(std::begin(yv), std::end(yv));
      xmax = *std::max_element(std::begin(xv), std::end(xv));
      ymax = *std::max_element(std::begin(yv), std::end(yv));
      gotValue = true;
    } else {
      xmin = std::min(*std::min_element(std::begin(xv), std::end(xv)), xmin);
      ymin = std::min(*std::min_element(std::begin(yv), std::end(yv)), ymin);
      xmax = std::max(*std::max_element(std::begin(xv), std::end(xv)), xmax);
      ymax = std::max(*std::max_element(std::begin(yv), std::end(yv)), ymax);
    }
  }
  for (const auto& seg : m_segments) {
    if (!gotValue) {
      xmin = std::min(seg.x0[0], seg.x1[0]);
      xmax = std::max(seg.x0[0], seg.x1[0]);
      ymin = std::min(seg.x0[1], seg.x1[1]);
      ymax = std::max(seg.x0[1], seg.x1[1]);
      gotValue = true;
    } else {
      xmin = std::min({xmin, seg.x0[0], seg.x1[0]});
      xmax = std::max({xmax, seg.x0[0], seg.x1[0]});
      ymin = std::min({ymin, seg.x0[1], seg.x1[1]});
      ymax = std::max({ymax, seg.x0[1], seg.x1[1]});
    }
  }
  for (const auto& wire : m_wires) {
    if (!gotValue) {
      xmin = xmax = wire.x;
      ymin = ymax = wire.y;
    } else {
      xmin = std::min(xmin, wire.x);
      xmax = std::max(xmax, wire.x);
      ymin = std::min(ymin, wire.y);
      ymax = std::max(ymax, wire.y);
    }
  }
  return gotValue; 
}

bool ComponentNeBem2d::CrossedWire(
    const double x0, const double y0, const double z0, 
    const double x1, const double y1, const double z1,
    double& xc, double& yc, double& zc, const bool centre, double& rc) {
  xc = x0;
  yc = y0;
  zc = z0;

  if (m_wires.empty()) return false;

  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double d2 = dx * dx + dy * dy;
  // Make sure the step length is non-zero.
  if (d2 < Small) return false;
  const double invd2 = 1. / d2;

  // Both coordinates are assumed to be located inside
  // the drift area and inside a drift medium.
  // This should have been checked before this call.

  double dMin2 = 0.;
  for (const auto& wire : m_wires) {
    const double xw = wire.x;
    const double yw = wire.y;
    // Calculate the smallest distance between track and wire.
    const double xIn0 = dx * (xw - x0) + dy * (yw - y0);
    // Check if the minimum is located before (x0, y0).
    if (xIn0 < 0.) continue;
    const double xIn1 = -(dx * (xw - x1) + dy * (yw - y1));
    // Check if the minimum is located behind (x1, y1).
    if (xIn1 < 0.) continue;
    // Minimum is located between (x0, y0) and (x1, y1).
    const double xw0 = xw - x0;
    const double xw1 = xw - x1;
    const double yw0 = yw - y0;
    const double yw1 = yw - y1;
    const double dw02 = xw0 * xw0 + yw0 * yw0;
    const double dw12 = xw1 * xw1 + yw1 * yw1;
    if (xIn1 * xIn1 * dw02 > xIn0 * xIn0 * dw12) {
      dMin2 = dw02 - xIn0 * xIn0 * invd2;
    } else {
      dMin2 = dw12 - xIn1 * xIn1 * invd2;
    }
    const double r2 = wire.r * wire.r;
    if (dMin2 < r2) {
      // Wire has been crossed.
      if (centre) {
        xc = xw;
        yc = yw;
      } else {
        // Find the point of intersection.
        const double p = -xIn0 * invd2;
        const double q = (dw02 - r2) * invd2;
        const double s = sqrt(p * p - q);
        const double t = std::min(-p + s, -p - s);
        xc = x0 + t * dx;
        yc = y0 + t * dy;
        zc = z0 + t * (z1 - z0);
      }
      rc = wire.r;
      return true;
    }
  }
  return false;
}

bool ComponentNeBem2d::InTrapRadius(const double q, const double x,
                                    const double y, const double /*z*/,
                                    double& xw, double& yw, double& rw) {

  for (const auto& wire : m_wires) {
    // Skip wires with the wrong charge.
    if (q * wire.q > 0.) continue;
    const double dx = wire.x - x;
    const double dy = wire.y - y;
    const double rt = wire.r * wire.ntrap;
    if (dx * dx + dy * dy < rt * rt) {
      xw = wire.x;
      yw = wire.y;
      rw = wire.r;
      return true;
    }
  }
  return false;
}

void ComponentNeBem2d::SetRangeZ(const double zmin, const double zmax) {

   if (fabs(zmax - zmin) <= 0.) {
     std::cerr << m_className << "::SetRangeZ: Zero range is not permitted.\n";
     return;
   }
   m_zmin = std::min(zmin, zmax);
   m_zmax = std::max(zmin, zmax);
   m_useRangeZ = true;
}

bool ComponentNeBem2d::AddSegment(const double x0, const double y0,
                                  const double x1, const double y1,
                                  const double v, const int ndiv) {
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  if (dx * dx + dy * dy < Small) {
    std::cerr << m_className << "::AddSegment: Length must be > 0.\n";
    return false;
  }

  Segment segment;
  segment.x0 = {x0, y0};
  segment.x1 = {x1, y1};
  segment.bc = std::make_pair(Voltage, v);
  segment.region1 = -1;
  segment.region2 = -1;
  segment.ndiv = ndiv;
  m_segments.push_back(std::move(segment));

  if (m_debug) {
    std::cout << m_className << "::AddSegment:\n    (" 
              << x0 << ", " << y0 << ") - (" << x1 << ", " << y1 << ")\n"
              << "    Potential: " << v << " V\n";
  }

  m_ready = false;
  return true;
}

bool ComponentNeBem2d::AddWire(const double x, const double y, const double d,
                               const double v, const int ntrap) {
  if (d < Small) {
    std::cerr << m_className << "::AddWire: Diameter must be > 0.\n";
    return false;
  }
  if (ntrap <= 0) {
    std::cerr << m_className << "::AddWire: Nbr. of trap radii must be > 0.\n";
    return false;
  }

  Wire wire;
  wire.x = x;
  wire.y = y;
  wire.r = 0.5 * d;
  wire.v = v;
  wire.q = 0.;
  wire.ntrap = ntrap;
  m_wires.push_back(std::move(wire));

  if (m_debug) {
    std::cout << m_className << "::AddWire:\n"
              << "    Centre: (" << x << ", " << y << ")\n"
              << "    Diameter: " << d << " cm\n"
              << "    Potential: " << v << " V\n";
  }

  m_ready = false;
  return true;
}

bool ComponentNeBem2d::AddRegion(const std::vector<double>& xp,
                                 const std::vector<double>& yp, 
                                 Medium* medium, const unsigned int bctype,
                                 const double v, const int ndiv) {

  if (xp.size() != yp.size()) {
    std::cerr << m_className << "::AddRegion:\n"
              << "    Mismatch between number of x- and y-coordinates.\n";
    return false;
  }
  if (xp.size() < 3) {
    std::cerr << m_className << "::AddRegion: Too few points.\n";
    return false;
  }
  if (bctype != 1 && bctype != 4) {
    std::cerr << m_className << "::AddRegion: Invalid boundary condition.\n";
    return false;
  }

  // Check if this is a valid polygon (no self-crossing).
  const unsigned int np = xp.size();
  if (np > 3) {
    for (unsigned int i0 = 0; i0 < np; ++i0) {
      const unsigned int i1 = i0 < np - 1 ? i0 + 1 : 0;
      for (unsigned int j = 0; j < np - 3; ++j) {
        const unsigned int j0 = i1 < np - 1 ? i1 + 1 : 0;
        const unsigned int j1 = j0 < np - 1 ? j0 + 1 : 0; 
        double xc = 0., yc = 0.;
        if (Crossing(xp[i0], yp[i0], xp[i1], yp[i1], 
                     xp[j0], yp[j0], xp[j1], yp[j1], xc, yc)) {
          std::cerr << m_className << "::AddRegion: Edges cross each other.\n";
          return false;
        }
      }
    }
  } 
  std::vector<double> xv = xp;
  std::vector<double> yv = yp;
  const double xmin = *std::min_element(std::begin(xv), std::end(xv));
  const double ymin = *std::min_element(std::begin(yv), std::end(yv));
  const double xmax = *std::max_element(std::begin(xv), std::end(xv));
  const double ymax = *std::max_element(std::begin(yv), std::end(yv));

  const double epsx = 1.e-6 * std::max(std::abs(xmax), std::abs(xmin));
  const double epsy = 1.e-6 * std::max(std::abs(ymax), std::abs(ymin));

  const double f = Polygon::Area(xp, yp);
  if (std::abs(f) < std::max(1.e-10, epsx * epsy)) {
    std::cerr << m_className << "::AddRegion: Degenerate polygon.\n";
    return false;
  } else if (f > 0.) {
    // Make sure all polygons have the same "handedness".
    if (m_debug) {
      std::cout << m_className << "::AddRegion: Reversing orientation.\n";
    }
    std::reverse(xv.begin(), xv.end());
    std::reverse(yv.begin(), yv.end());
  }
  for (const auto& region : m_regions) {
    if (Intersecting(xv, yv, region.xv, region.yv)) {
      std::cerr << m_className << "::AddRegion:\n"
                << "    Polygon intersects an existing region.\n";
      return false;
    }
  }
  Region region;
  region.xv = xv;
  region.yv = yv;
  region.medium = medium;
  if (bctype == 1) {
    region.bc = std::make_pair(Voltage, v);
  } else if (bctype == 4) {
    region.bc = std::make_pair(Dielectric, v);
  }
  region.depth = 0;
  region.ndiv = ndiv;
  m_regions.push_back(std::move(region));
  return true;
}

void ComponentNeBem2d::AddChargeDistribution(const double x, const double y,
                                             const double a, const double b,
                                             const double rho) {
  if (a < Small || b < Small) {
    std::cerr << m_className << "::AddChargeDistribution:\n"
              << "    Lengths must be > 0.\n";
    return;
  }
  const double a2 = a * a;
  const double b2 = b * b;
  const double v0 = -2 * (Pi * b2 + 2 * atan(b / a) * (a2 - b2));
  SpaceCharge box;
  box.x = x;
  box.y = y;
  box.a = a;
  box.b = b;
  box.q = rho;
  box.v0 = v0;
  m_spaceCharge.push_back(std::move(box));
}

void ComponentNeBem2d::SetNumberOfDivisions(const unsigned int ndiv) {
  if (ndiv == 0) {
    std::cerr << m_className << "::SetNumberOfDivisions:\n"
              << "    Number of divisions must be greater than zero.\n";
    return;
  }

  m_nDivisions = ndiv;
  m_ready = false;
}

void ComponentNeBem2d::SetNumberOfCollocationPoints(const unsigned int ncoll) {
  if (ncoll == 0) {
    std::cerr << m_className << "::SetNumberOfCollocationPoints:\n"
              << "    Number of points must be greater than zero.\n";
    return;
  }

  m_nCollocationPoints = ncoll;
  m_ready = false;
}

void ComponentNeBem2d::SetMaxNumberOfIterations(const unsigned int niter) {
  if (niter == 0) {
    std::cerr << m_className << "::SetMaxNumberOfIterations:\n"
              << "    Number of iterations must be greater than zero.\n";
    return;
  }
  m_nMaxIterations = niter;
}

bool ComponentNeBem2d::GetRegion(const unsigned int i,
                                 std::vector<double>& xv, 
                                 std::vector<double>& yv,
                                 Medium*& medium, unsigned int& bctype, 
                                 double& v) {
  if (i >= m_regions.size()) return false;
  if (!m_ready) {
    if (!Initialise()) return false;
  }
  const auto& region = m_regions[i];
  xv = region.xv;
  yv = region.yv;
  medium = region.medium;
  bctype = region.bc.first == Voltage ? 1 : 4;
  v = region.bc.second;
  return true; 
}

bool ComponentNeBem2d::GetSegment(const unsigned int i, 
    double& x0, double& y0, double& x1, double& y1, double& v) const {
 
  if (i >= m_segments.size()) return false;
  const auto& seg = m_segments[i];
  x0 = seg.x0[0];
  y0 = seg.x0[1];
  x1 = seg.x1[0];
  y1 = seg.x1[1];
  v = seg.bc.second;
  return true;
}

bool ComponentNeBem2d::GetWire(const unsigned int i,
  double& x, double& y, double& d, double& v, double& q) const {
 
  if (i >= m_wires.size()) return false;
  const auto& wire = m_wires[i];
  x = wire.x;
  y = wire.y;
  d = 2 * wire.r;
  v = wire.v;
  q = wire.q;
  return true;
}

bool ComponentNeBem2d::GetElement(const unsigned int i,
  double& x0, double& y0, double& x1, double& y1, double& q) const {

  if (i >= m_elements.size()) return false;
  const auto& element = m_elements[i];
  ToGlobal(-element.a, 0., element.cphi, element.sphi, x0, y0);
  ToGlobal( element.a, 0., element.cphi, element.sphi, x1, y1);
  x0 += element.x;
  y0 += element.y;
  x1 += element.x;
  y1 += element.y;
  q = element.q;
  return true;
}

bool ComponentNeBem2d::Initialise() {

  m_ready = false;
  m_elements.clear();

  double vmin = 0., vmax = 0.; 
  if (!GetVoltageRange(vmin, vmax)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Could not determine the voltage range.\n";
    return false;
  }
  if (fabs(vmin - vmax) < 1.e-6 * (vmin + vmax)) {
    std::cerr << m_className << "::Initialise:\n"
              << "    All potentials are the same.\n";
    return false;
  }

  if (m_debug) std::cout << m_className << "::Initialise:\n";
  // Loop over the regions.
  const unsigned int nRegions = m_regions.size();
  if (m_debug) std::cout << "    " << nRegions << " regions.\n";
  std::vector<std::vector<unsigned int> > motherRegions(nRegions);
  for (unsigned int i = 0; i < nRegions; ++i) {
    auto& region = m_regions[i];
    // Check if the region is fully enclosed by other ones.
    for (unsigned int j = 0; j < nRegions; ++j) {
      if (i == j) continue;
      const auto& other = m_regions[j];
      if (!Enclosed(region.xv, region.yv, other.xv, other.yv)) continue;
      motherRegions[i].push_back(j);
      if (m_debug) {
        std::cout << "    Region " << i << " is enclosed by region "
                  << j << ".\n";
      }
    }
    region.depth = motherRegions[i].size();
  }

  std::vector<Segment> segments;
  for (unsigned int i = 0; i < nRegions; ++i) {
    const auto& region = m_regions[i];
    int outerRegion = -1;
    for (const unsigned int k : motherRegions[i]) {
      if (outerRegion < 0) {
        outerRegion = k;
      } else if (m_regions[outerRegion].depth < m_regions[k].depth) {
        outerRegion = k;
      }
    }
    // Add the segments bounding this region. 
    const unsigned int n = region.xv.size();
    for (unsigned int j = 0; j < n; ++j) {
      const unsigned int k = j < n - 1 ? j + 1 : 0;
      Segment seg;
      seg.x0 = {region.xv[j], region.yv[j]};
      seg.x1 = {region.xv[k], region.yv[k]};
      seg.region1 = i;
      seg.region2 = outerRegion;
      seg.bc = region.bc;
      seg.ndiv = region.ndiv;
      segments.push_back(std::move(seg));
    }
  }
  // Add the segments specified by the user.
  segments.insert(segments.end(), m_segments.begin(), m_segments.end());
  const unsigned int nSegments = segments.size();
  if (m_debug) std::cout << "    " << nSegments << " segments.\n";
  std::vector<bool> done(nSegments, false);
  // Look for overlaps.
  for (unsigned int i = 0; i < nSegments; ++i) {
    if (done[i]) continue;
    if (m_debug) {
      std::cout << "    Segment " << i << ". (" 
                << segments[i].x0[0] << ", " << segments[i].x0[1] << ") - (" 
                << segments[i].x1[0] << ", " << segments[i].x1[1] << ")\n";
    }
    const double x0 = segments[i].x0[0];
    const double x1 = segments[i].x1[0];
    const double y0 = segments[i].x0[1];
    const double y1 = segments[i].x1[1];
    // Pick up all collinear segments.
    std::vector<Segment> newSegments;
    for (unsigned int j = i + 1; j < nSegments; ++j) {
      const double u0 = segments[j].x0[0];
      const double u1 = segments[j].x1[0];
      const double v0 = segments[j].x0[1];
      const double v1 = segments[j].x1[1];
      const double epsx = std::max(1.e-10 * std::max({fabs(x0), fabs(x1), 
                                                      fabs(u0), fabs(u1)}), 
                                   1.e-10);
      const double epsy = std::max(1.e-10 * std::max({fabs(y0), fabs(y1), 
                                                      fabs(v0), fabs(v1)}), 
                                   1.e-10);
      const double a00 = y1 - y0;
      const double a01 = v1 - v0;
      const double a10 = x0 - x1;
      const double a11 = u0 - u1;
      const double det = a00 * a11 - a10 * a01;
      const double tol = epsx * epsy;
      // Skip non-parallel segments.
      if (std::abs(det) > tol) continue;

      if (std::abs(x0 * (y1 - v0) + x1 * (v0 - y0) + u0 * (y0 - y1)) > tol) {
        continue;
      }
      newSegments.push_back(segments[j]);
      done[j] = true;
    }
    newSegments.push_back(segments[i]);
    if (newSegments.size() > 1) {
      if (m_debug) {
        std::cout << "      Determining overlaps of " << newSegments.size() 
                  << " collinear segments.\n";
      }
      EliminateOverlaps(newSegments);
      if (m_debug) {
        std::cout << "      " << newSegments.size() 
                  << " segments after splitting/merging.\n";
      }
    } 
    for (const auto& segment : newSegments) {
      double lambda = 0.;
      if (segment.bc.first == Dielectric) {
        // Dielectric-dielectric interface.
        const int reg1 = segment.region1;
        const int reg2 = segment.region2;
        double eps1 = 1.;
        if (reg1 < 0) {
          if (m_medium) eps1 = m_medium->GetDielectricConstant();
        } else if (m_regions[reg1].medium) {
          eps1 = m_regions[reg1].medium->GetDielectricConstant();
        }
        double eps2 = 1.;
        if (reg2 < 0) {
          if (m_medium) eps2 = m_medium->GetDielectricConstant();
        } else if (m_regions[reg2].medium) {
          eps2 = m_regions[reg2].medium->GetDielectricConstant();
        }
        if (fabs(eps1 - eps2) < 1.e-6 * (1. + fabs(eps1) + fabs(eps2))) {
          if (m_debug) std::cout << "      Same epsilon. Skip.\n";
          continue;
        }
        // Compute lambda.
        lambda = (eps2 - eps1) / (eps1 + eps2);
        if (m_debug) std::cout << "      Lambda = " << lambda << "\n";
      }
      const int ndiv = segment.ndiv <= 0 ? m_nDivisions : segment.ndiv;
      Discretise(segment, m_elements, lambda, ndiv);
    }
  }  
  std::vector<std::vector<double> > influenceMatrix;
  std::vector<std::vector<double> > inverseMatrix;

  bool converged = false;
  unsigned int nIter = 0;
  while (!converged) {
    ++nIter;
    if (m_autoSize) {
      std::cout << m_className << "::Initialise: Iteration " << nIter << "\n";
    }
    const unsigned int nElements = m_elements.size();
    const unsigned int nEntries = nElements + m_wires.size() + 1;
    if (m_debug) {
      std::cout << "    " << nElements << " elements.\n"
                << "    Matrix has " << nEntries << " rows/columns.\n";
    }
    // Compute the influence matrix.
    influenceMatrix.assign(nEntries, std::vector<double>(nEntries, 0.));
    if (!ComputeInfluenceMatrix(influenceMatrix)) {
      std::cerr << m_className << "::Initialise:\n"
                << "     Error computing the influence matrix.\n";
      return false;
    }

    // Invert the influence matrix.
    inverseMatrix.assign(nEntries, std::vector<double>(nEntries, 0.));
    if (!InvertMatrix(influenceMatrix, inverseMatrix)) {
      std::cerr << m_className << "::Initialise: Matrix inversion failed.\n";
      return false;
    }
    if (m_debug) std::cout << "    Matrix inversion ok.\n";

    // Compute the right hand side vector (boundary conditions).
    std::vector<double> boundaryConditions(nEntries, 0.);
    for (unsigned int i = 0; i < nElements; ++i) {
      if (m_elements[i].bc.first == Voltage) {
        boundaryConditions[i] = m_elements[i].bc.second;
        for (const auto& box : m_spaceCharge) {
          const double x = m_elements[i].x - box.x;
          const double y = m_elements[i].y - box.y;
          const double vs = BoxPotential(box.a, box.b, x, y, box.v0) * box.q;
          boundaryConditions[i] -= vs;
        }
      } else {
        for (const auto& box : m_spaceCharge) {
          const double x = m_elements[i].x - box.x;
          const double y = m_elements[i].y - box.y;
          double fx = 0., fy = 0.;
          BoxField(box.a, box.b, x, y, fx, fy);
          // Rotate to the local frame of the target element.
          ToLocal(fx, fy, m_elements[i].cphi, m_elements[i].sphi, fx, fy);
          boundaryConditions[i] -= box.q * fy;
        }
      }
    }
    const unsigned int nWires = m_wires.size();
    for (unsigned int i = 0; i < nWires; ++i) {
      boundaryConditions[nElements + i] = m_wires[i].v;
      for (const auto& box : m_spaceCharge) {
        const double x = m_wires[i].x - box.x;
        const double y = m_wires[i].y - box.y;
        const double vs = BoxPotential(box.a, box.b, x, y, box.v0) * box.q;
        boundaryConditions[nElements + i] -= vs;
      }
    }

    double qsum = 0.;
    for (const auto& box : m_spaceCharge) {
      qsum += 4 * box.q * box.a * box.b;
    }
    boundaryConditions.back() = -qsum;

    // Solve for the charge distribution.
    if (!Solve(inverseMatrix, boundaryConditions)) {
      std::cerr << m_className << "::Initialise: Solution failed.\n";
      return false;
    }
    if (m_debug) std::cout << "    Solution ok.\n";
    const double tol = 1.e-6 * fabs(vmax - vmin);
    std::vector<bool> ok(nElements, true);
    converged = CheckConvergence(tol, ok);
    if (!m_autoSize) break;
    if (nIter >= m_nMaxIterations) break;
    for (unsigned int j = 0; j < nElements; ++j) {
      if (!ok[j]) {
        SplitElement(m_elements[j], m_elements);
        if (m_debug) std::cout << "    Splitting element " << j << ".\n";
      }
    }
  }
  // Sort the regions by depth (innermost first).
  std::sort(m_regions.begin(), m_regions.end(),  
            [](const Region& lhs, const Region& rhs) { 
              return (lhs.depth > rhs.depth);
            });
  m_ready = true;
  return true;
}

void ComponentNeBem2d::EliminateOverlaps(std::vector<Segment>& segments) {

  if (segments.empty()) return;
  const unsigned int nIn = segments.size();
  // Find the first/last point along the line.
  std::array<double, 2> x0 = segments[0].x0;
  std::array<double, 2> x1 = segments[0].x1;
  // Use x or y coordinate depending on the orientation of the line.
  const unsigned int ic = fabs(x1[1] - x0[1]) > fabs(x1[0] - x0[0]) ? 1 : 0;
  std::vector<bool> swapped(nIn, false);
  for (unsigned int i = 0; i < nIn; ++i) {
    const auto& seg = segments[i];
    std::array<double, 2> u0 = seg.x0;
    std::array<double, 2> u1 = seg.x1;
    if (u0[ic] > u1[ic]) {
      // Swap points.
      std::swap(u0, u1);
      swapped[i] = true;
    }
    if (u0[ic] < x0[ic]) x0 = u0;
    if (u1[ic] > x1[ic]) x1 = u1;
  }
  const std::array<double, 2> d = {x1[0] - x0[0], x1[1] - x0[1]};

  // Make a list of all points and their linear coordinate.
  std::vector<std::pair<double, std::vector<unsigned int> > > points;
  for (unsigned int i = 0; i < nIn; ++i) {
    for (const auto& xl : {segments[i].x0, segments[i].x1}) {
      const std::array<double, 2> d0 = {xl[0] - x0[0], xl[1] - x0[1]};
      const std::array<double, 2> d1 = {x1[0] - xl[0], x1[1] - xl[1]};
      double lambda = 0.;
      if (Mag2(d0) < Mag2(d1)) {
        // Point nearer to x0.
        lambda = d0[ic] / d[ic];
      } else {
        // Point nearer to p1.
        lambda = 1. - d1[ic] / d[ic];
      }
      // Add the point to the list.
      bool found = false;
      for (auto& point : points) {
        if (fabs(point.first - lambda) < 1.e-6) {
          found = true;
          point.second.push_back(i);
          break;
        }
      }
      if (found) continue;
      points.push_back(std::make_pair(lambda, 
                                      std::vector<unsigned int>({i})));
    }
  }
  // Sort the points by linear coordinate.
  std::sort(std::begin(points), std::end(points));
  
  std::vector<Segment> newSegments;
  const unsigned int nPoints = points.size();
  std::array<double, 2> xl = {x0[0] + points[0].first * d[0], 
                              x0[1] + points[0].first * d[1]};
  std::vector<unsigned int> left = points[0].second;
  for (unsigned int i = 1; i < nPoints; ++i) {
    Segment seg = segments[left.front()];
    seg.x0 = xl;
    xl = {x0[0] + points[i].first * d[0], x0[1] + points[i].first * d[1]};
    seg.x1 = xl;
    if (swapped[left.front()]) std::swap(seg.x0, seg.x1);
    // Sort out the boundary conditions.
    if (left.size() > 1) {
      for (unsigned int j = 1; j < left.size(); ++j) {
        const auto& other = segments[left[j]];
        if (seg.bc.first == Dielectric) {
          if (other.bc.first == Dielectric) {
            // Dielectric-dielectric interface.
            if ((seg.x1[ic] - seg.x0[ic]) * (other.x1[ic] - other.x0[ic]) > 0) {
              // Same orientation.
              continue;
            }
            seg.region2 = other.region1;
          } else {
            seg.bc = other.bc;
          }
        } else if (seg.bc.first != other.bc.first) {
          std::cerr << m_className << "::EliminateOverlaps:\n"
                    << "    Warning: conflicting boundary conditions.\n";
        }
      }
    }
    newSegments.push_back(std::move(seg));
    for (unsigned int k : points[i].second) {
      const auto it = std::find(left.begin(), left.end(), k);
      if (it == left.end()) {
        left.push_back(k);
      } else {
        left.erase(it);
      }
    }
  }
  segments.swap(newSegments); 
}

bool ComponentNeBem2d::Discretise(const Segment& seg,
                                  std::vector<Element>& elements,
                                  const double lambda, 
                                  const unsigned int ndiv) {

  if (ndiv < 1) {
    std::cerr << m_className << "::Discretise: Number of elements < 1.\n";
    return false;
  }
  const double phi = atan2(seg.x1[1] - seg.x0[1], seg.x1[0] - seg.x0[0]);
  const double cphi = cos(phi);
  const double sphi = sin(phi);
  const double dx = (seg.x1[0] - seg.x0[0]) / ndiv;
  const double dy = (seg.x1[1] - seg.x0[1]) / ndiv;
  const double a = 0.5 * sqrt(dx * dx + dy * dy);
  double x = seg.x0[0] - 0.5 * dx;
  double y = seg.x0[1] - 0.5 * dy;
  for (unsigned int i = 0; i < ndiv; ++i) {
    x += dx;
    y += dy;
    Element element;
    element.cphi = cphi;
    element.sphi = sphi;
    element.x = x;
    element.y = y;
    element.a = a;
    element.bc = seg.bc;
    element.lambda = lambda;
    elements.push_back(std::move(element));
  }
  return true;
}
 
bool ComponentNeBem2d::ComputeInfluenceMatrix(
    std::vector<std::vector<double> >& infmat) const {

  const unsigned int nL = m_elements.size();
  const unsigned int nE = nL + m_wires.size();
  // Loop over the target elements (F).
  for (unsigned int iF = 0; iF < nE; ++iF) {
    const auto bcF = iF < nL ? m_elements[iF].bc.first : Voltage;
    const double cphiF = iF < nL ? m_elements[iF].cphi : 1.;
    const double sphiF = iF < nL ? m_elements[iF].sphi : 0.;
    // Collocation point.
    const double xF = iF < nL ? m_elements[iF].x : m_wires[iF - nL].x;
    const double yF = iF < nL ? m_elements[iF].y : m_wires[iF - nL].y;
    
    // Loop over the source elements (S).
    for (unsigned int jS = 0; jS < nE; ++jS) {
      // Calculate the influence coefficient.
      double infCoeff = 0.;
      if (jS < nL) { 
        // Straight line element.
        const auto& src = m_elements[jS];
        double xL = 0., yL = 0.;
        ToLocal(xF - src.x, yF - src.y, src.cphi, src.sphi, xL, yL);
        if (bcF == Voltage) {
          infCoeff = LinePotential(src.a, xL, yL);
        } else if (bcF == Dielectric) {
          // Dielectric-dielectric interface.
          // Normal component of the displacement vector is continuous.
          if (iF == jS) {
            // Self-influence.
            infCoeff = 1. / (2. * src.lambda * VacuumPermittivity);
          } else {
            // Compute flux at the collocation point.
            double fx = 0., fy = 0.;
            LineField(src.a, xL, yL, fx, fy);
            // Rotate to the global frame.
            ToGlobal(fx, fy, src.cphi, src.sphi, fx, fy);
            // Rotate to the local frame of the target element.
            ToLocal(fx, fy, cphiF, sphiF, fx, fy);
            infCoeff = fy;
          }
        }
      } else {
        // Wire.
        const auto& src = m_wires[jS - nL];
        if (bcF == Voltage) {
          infCoeff = WirePotential(src.r, xF - src.x, yF - src.y);
        } else if (bcF == Dielectric) {
          double fx = 0., fy = 0.;
          WireField(src.r, xF - src.x, yF - src.y, fx, fy);
          ToLocal(fx, fy, cphiF, sphiF, fx, fy);
          infCoeff = fy;
        }
      }
      infmat[iF][jS] = infCoeff;
    }
  }

  // Add charge neutrality condition.
  for (unsigned int i = 0; i < nE; ++i) {
    if (i < nL) {
      infmat[nE][i] = m_elements[i].a;
    } else {
      infmat[nE][i] = m_wires[i - nL].r;
    }
    infmat[i][nE] = 0.;
  }
  infmat[nE][nE] = 0.;

  return true;
}

void ComponentNeBem2d::SplitElement(Element& oldElement, 
  std::vector<Element>& elements) {
  oldElement.a *= 0.5;

  Element newElement = oldElement;
  double dx = 0., dy = 0.;
  ToGlobal(newElement.a, 0., newElement.cphi, newElement.sphi, dx, dy);
  oldElement.x += dx;
  oldElement.y += dy;
  newElement.x -= dx;
  newElement.y -= dx;

  elements.push_back(std::move(newElement));
}

bool ComponentNeBem2d::InvertMatrix(
    std::vector<std::vector<double> >& influenceMatrix,
    std::vector<std::vector<double> >& inverseMatrix) const {

  const unsigned int nEntries = influenceMatrix.size();

  // Temporary arrays for LU decomposition/substitution
  std::vector<double> col(nEntries, 0.);
  std::vector<int> index(nEntries, 0);

  // Decompose the influence matrix
  if (!LUDecomposition(influenceMatrix, index)) {
    std::cerr << m_className << "::InvertMatrix: LU decomposition failed.\n";
    return false;
  }

  // Initialise the inverse influence matrix
  inverseMatrix.assign(nEntries, std::vector<double>(nEntries, 0.));
  // Invert the matrix.
  for (unsigned int j = 0; j < nEntries; ++j) {
    col.assign(nEntries, 0.);
    col[j] = 1.;
    LUSubstitution(influenceMatrix, index, col);
    for (unsigned int i = 0; i < nEntries; ++i) inverseMatrix[i][j] = col[i];
  }

  // Clear the influence matrix.
  influenceMatrix.clear();

  return true;
}

bool ComponentNeBem2d::LUDecomposition(std::vector<std::vector<double> >& mat,
                                       std::vector<int>& index) const {
  // The influence matrix is replaced by the LU decomposition of a rowwise
  // permutation of itself. The implementation is based on:
  // W. H. Press,
  // Numerical recipes in C++: the Art of Scientific Computing (version 2.11)

  const unsigned int n = m_elements.size() + m_wires.size();
  // v stores the implicit scaling of each row
  std::vector<double> v(n, 0.);

  // Loop over rows to get the implicit scaling information.
  for (unsigned int i = 0; i < n; ++i) {
    double big = 0.;
    for (unsigned int j = 0; j < n; ++j) {
      big = std::max(big, fabs(mat[i][j]));
    }
    if (big == 0.) return false;
    // Save the scaling
    v[i] = 1. / big;
  }

  // Loop over columns
  unsigned int imax = 0;
  for (unsigned int j = 0; j < n; ++j) {
    for (unsigned int i = 0; i < j; ++i) {
      double sum = mat[i][j];
      for (unsigned int k = 0; k < i; ++k) {
        sum -= mat[i][k] * mat[k][j];
      }
      mat[i][j] = sum;
    }
    // Initialise for the search for the largest pivot element
    double big = 0.;
    for (unsigned int i = j; i < n; ++i) {
      double sum = mat[i][j];
      for (unsigned int k = 0; k < j; ++k) {
        sum -= mat[i][k] * mat[k][j];
      }
      mat[i][j] = sum;
      // Is the figure of merit for the pivot better than the best so far?
      const double dum = v[i] * fabs(sum);
      if (dum >= big) {
        big = dum;
        imax = i;
      }
    }
    // Do we need to interchange rows?
    if (j != imax) {
      for (unsigned k = 0; k < n; ++k) {
        const double dum = mat[imax][k];
        mat[imax][k] = mat[j][k];
        mat[j][k] = dum;
      }
      // Interchange the scale factor
      v[imax] = v[j];
    }
    index[j] = imax;
    if (mat[j][j] == 0.) mat[j][j] = Small;
    if (j != n - 1) {
      // Divide by the pivot element
      const double dum = 1. / mat[j][j];
      for (unsigned int i = j + 1; i < n; ++i) {
        mat[i][j] *= dum;
      }
    }
  }

  return true;
}

void ComponentNeBem2d::LUSubstitution(
    const std::vector<std::vector<double> >& mat, const std::vector<int>& index,
    std::vector<double>& col) const {

  const unsigned int n = m_elements.size() + m_wires.size();
  unsigned int ii = 0;
  // Forward substitution
  for (unsigned i = 0; i < n; ++i) {
    const unsigned int ip = index[i];
    double sum = col[ip];
    col[ip] = col[i];
    if (ii != 0) {
      for (unsigned j = ii - 1; j < i; ++j) {
        sum -= mat[i][j] * col[j];
      }
    } else if (sum != 0.) {
      ii = i + 1;
    }
    col[i] = sum;
  }

  // Backsubstitution
  for (int i = n - 1; i >= 0; i--) {
    double sum = col[i];
    for (unsigned j = i + 1; j < n; ++j) {
      sum -= mat[i][j] * col[j];
    }
    col[i] = sum / mat[i][i];
  }
}

bool ComponentNeBem2d::Solve(const std::vector<std::vector<double> >& invmat,
                             const std::vector<double>& bc) {
  const unsigned int nEntries = bc.size();
  const unsigned int nElements = m_elements.size();
  for (unsigned int i = 0; i < nElements; ++i) {
    double solution = 0.;
    for (unsigned int j = 0; j < nEntries; ++j) {
      solution += invmat[i][j] * bc[j];
    }
    m_elements[i].q = solution;
  }
  const unsigned int nWires = m_wires.size();
  for (unsigned int i = 0; i < nWires; ++i) {
    double solution = 0.;
    for (unsigned int j = 0; j < nEntries; ++j) {
      solution += invmat[nElements + i][j] * bc[j];
    }
    m_wires[i].q = solution;
  }

  if (m_debug) {
    std::cout << m_className << "::Solve:\n  Element  Solution\n";
    for (unsigned int i = 0; i < nElements; ++i) {
      std::printf(" %8u   %15.5f\n", i, m_elements[i].q);
    }
    if (!m_wires.empty()) {
      std::cout << "   Wire    Solution\n";
      for (unsigned int i = 0; i < nWires; ++i) {
        std::printf("  %8u   %15.5f\n", i, m_wires[i].q);
      }
    }
  }
  return true;
}

bool ComponentNeBem2d::CheckConvergence(const double tol,
                                        std::vector<bool>& ok) {

  // Potential and normal component of the electric field
  // evaluated at the collocation points.
  std::vector<double> v(m_nCollocationPoints, 0.);
  std::vector<double> n(m_nCollocationPoints, 0.);

  if (m_debug) {
    std::cout << m_className << "::CheckConvergence:\n"
              << "  element #  type          LHS              RHS\n";
  }
  const double scale = 1. / m_nCollocationPoints;
  unsigned int i = 0;
  for (const auto& tgt : m_elements) {
    v.assign(m_nCollocationPoints, 0.);
    n.assign(m_nCollocationPoints, 0.);
    double dx = 0., dy = 0.;
    ToGlobal(2 * tgt.a, 0., tgt.cphi, tgt.sphi, dx, dy);
    const double x0 = tgt.x - 0.5 * dx;
    const double y0 = tgt.y - 0.5 * dy;
    // Loop over the collocation points.
    for (unsigned int k = 0; k < m_nCollocationPoints; ++k) {
      double xG = x0;
      double yG = y0;
      if (m_randomCollocation) {
        const double r = RndmUniformPos();
        xG += r * dx;
        yG += r * dy;
      } else {
        const double s = (k + 1.) / (m_nCollocationPoints + 1.);
        xG += s * dx;
        yG += s * dy;
      }
      // Sum up the contributions from all boundary elements.
      for (const auto& src : m_elements) {
        double xL = 0., yL = 0.;
        // Transform to local coordinate system.
        ToLocal(xG - src.x, yG - src.y, src.cphi, src.sphi, xL, yL);
        // Compute the potential.
        v[k] += LinePotential(src.a, xL, yL) * src.q;
        // Compute the field.
        double fx = 0., fy = 0.;
        LineField(src.a, xL, yL, fx, fy);
        // Rotate to the global frame.
        ToGlobal(fx, fy, src.cphi, src.sphi, fx, fy);
        // Rotate to the local frame of the test element.
        ToLocal(fx, fy, tgt.cphi, tgt.sphi, fx, fy);
        n[k] += fy * src.q;
      }

      for (const auto& src : m_wires) {
        // Compute the potential.
        v[k] += WirePotential(src.r, xG - src.x, yG - src.y) * src.q;
        // Compute the field.
        double fx = 0., fy = 0.;
        WireField(src.r, xG - src.x, yG - src.y, fx, fy);
        // Rotate to the local frame of the test element.
        ToLocal(fx, fy, tgt.cphi, tgt.sphi, fx, fy);
        n[k] += fy * src.q;
      }
      
      for (const auto& box : m_spaceCharge) {
        const double xL = xG - box.x;
        const double yL = yG - box.y;
        v[k] += BoxPotential(box.a, box.b, xL, yL, box.v0) * box.q;
        double fx = 0., fy = 0.;
        BoxField(box.a, box.b, xL, yL, fx, fy);
        ToLocal(fx, fy, tgt.cphi, tgt.sphi, fx, fy);
        n[k] += fy * box.q; 
      }
    }
    const double v0 = scale * std::accumulate(v.begin(), v.end(), 0.);
    const double n0 = scale * std::accumulate(n.begin(), n.end(), 0.);
    double n1 = 0.;
    if (tgt.bc.first == Voltage) {
      const double dv = v0 - tgt.bc.second;
      if (fabs(dv) > tol) ok[i] = false;
      if (m_debug) {
        std::printf(" %8u  cond.  %15.5f  %15.5f %15.5f\n", 
                    i, v0, tgt.bc.second, dv);
      }
    } else if (tgt.bc.first == Dielectric) {
      // Dielectric-dielectric interface
      // TODO.
      n1 = n0 + 0.5 * InvEpsilon0 * tgt.q / tgt.lambda;
      if (m_debug) std::printf(" %8u  diel.  %15.5f  %15.5f\n", i, n0, n1);
    }
    ++i;
  }

  for (const auto& tgt : m_wires) {
    v.assign(m_nCollocationPoints, 0.);
    double x0 = tgt.x;
    double y0 = tgt.y;
    // Loop over the collocation points.
    for (unsigned int k = 0; k < m_nCollocationPoints; ++k) {
      const double phi = TwoPi * RndmUniform();
      const double xG = x0 + tgt.r * cos(phi);
      const double yG = y0 + tgt.r * sin(phi);
      // Sum up the contributions from all boundary elements.
      for (const auto& src : m_elements) {
        // Transform to local coordinate system.
        double xL = 0., yL = 0.;
        ToLocal(xG - src.x, yG - src.y, src.cphi, src.sphi, xL, yL);
        // Compute the potential.
        v[k] += LinePotential(src.a, xL, yL) * src.q;
      }
      for (const auto& src : m_wires) {
        v[k] += WirePotential(src.r, xG - src.x, yG - src.y) * src.q;
      }
      for (const auto& box : m_spaceCharge) {
        const double xL = xG - box.x;
        const double yL = yG - box.y;
        v[k] += BoxPotential(box.a, box.b, xL, yL, box.v0) * box.q;
      }
    }
    const double v0 = scale * std::accumulate(v.begin(), v.end(), 0.);
    if (m_debug) {
      std::printf(" %8u  wire   %15.5f  %15.5f\n", i, v0, tgt.v);
    }
    ++i;
  }

  return true;
}

double ComponentNeBem2d::LinePotential(const double a, const double x,
                                       const double y) const {
  double p = 0.;
  const double amx = a - x;
  const double apx = a + x;
  if (fabs(y) > Small) {
    const double y2 = y * y;
    p = 2. * a - y * (atan(amx / y) + atan(apx / y)) -
        0.5 * amx * log(amx * amx + y2) - 0.5 * apx * log(apx * apx + y2);
  } else if (fabs(x) != a) {
    p = 2. * a - 0.5 * amx * log(amx * amx) - 0.5 * apx * log(apx * apx);
  } else {
    p = 2. * a * (1. - log(2. * a));
  }

  return InvTwoPiEpsilon0 * p;
}

double ComponentNeBem2d::WirePotential(const double r0, const double x,
                                       const double y) const {
  const double r = sqrt(x * x + y * y);
  if (r >= r0) {
    return -log(r) * r0 * InvEpsilon0;
  }

  // Inside the wire the potential is constant
  return -log(r0) * r0 * InvEpsilon0;
}

void ComponentNeBem2d::LineField(const double a, 
                                 const double x, const double y,
                                 double& ex, double& ey) const {
  const double amx = a - x;
  const double apx = a + x;
  if (fabs(y) > 0.) {
    const double y2 = y * y;
    ex = 0.5 * log((apx * apx + y2) / (amx * amx + y2));
    ey = atan(amx / y) + atan(apx / y);
  } else if (fabs(x) != a) {
    ex = 0.5 * log(apx * apx / (amx * amx));
    ey = 0.;
  } else {
    // Singularity at the end points of the line
    constexpr double eps2 = 1.e-24;
    ex = 0.25 * log(pow(apx * apx - eps2, 2) / pow(amx * amx - eps2, 2));
    ey = 0.;
  }
  ex *= InvTwoPiEpsilon0;
  ey *= InvTwoPiEpsilon0;
}

void ComponentNeBem2d::WireField(const double r0, 
                                 const double x, const double y,
                                 double& ex, double& ey) const {
  const double r02 = r0 * r0;
  const double r2 = x * x + y * y;
  if (r2 > r02) {
    ex = x * r0 / r2;
    ey = y * r0 / r2;
  } else if (r2 == r02) {
    ex = 0.5 * x / r0;
    ey = 0.5 * y / r0;
  } else {
    // Inside the wire the field is zero.
    ex = ey = 0.;
    return;
  }
  ex *= InvEpsilon0;
  ey *= InvEpsilon0;
}

double ComponentNeBem2d::BoxPotential(const double a, const double b, 
                                      const double x, const double y,
                                      const double v0) const {

  const double invc2 = 1. / (a * a + b * b);
  double v1 = 0., v2 = 0.;
  if (fabs(x) > a || fabs(y) > b) {
    // Outside the rectangle.
    const std::array<double, 2> dx = {x - a, x + a};
    const std::array<double, 2> dy = {y - b, y + b};
    const std::array<double, 2> dx2 = {dx[0] * dx[0], dx[1] * dx[1]};
    const std::array<double, 2> dy2 = {dy[0] * dy[0], dy[1] * dy[1]};
    v1 = dx[0] * dy[0] * log((dx2[0] + dy2[0]) * invc2) 
       - dx[1] * dy[0] * log((dx2[1] + dy2[0]) * invc2)
       + dx[1] * dy[1] * log((dx2[1] + dy2[1]) * invc2) 
       - dx[0] * dy[1] * log((dx2[0] + dy2[1]) * invc2);
    std::array<double, 4> alpha = {atan2(dy[0], dx[0]), atan2(dy[0], dx[1]),
                                   atan2(dy[1], dx[1]), atan2(dy[1], dx[0])};
    if (x < 0.) {
      for (size_t i = 0; i < 4; ++i) if (alpha[i] < 0.) alpha[i] += TwoPi;
    }
    v2 = dx2[0] * (alpha[0] - alpha[3]) + dx2[1] * (alpha[2] - alpha[1]) +
         dy2[0] * (alpha[1] - alpha[0]) + dy2[1] * (alpha[3] - alpha[2]);
  } else {
    // Inside the rectangle.
    const std::array<double, 2> dx = {a - x, a + x};
    const std::array<double, 2> dy = {b - y, b + y};
    const std::array<double, 2> dx2 = {dx[0] * dx[0], dx[1] * dx[1]};
    const std::array<double, 2> dy2 = {dy[0] * dy[0], dy[1] * dy[1]};
    v1 = dx[0] * dy[0] * log((dx2[0] + dy2[0]) * invc2) +
         dy[0] * dx[1] * log((dx2[1] + dy2[0]) * invc2) +
         dx[1] * dy[1] * log((dx2[1] + dy2[1]) * invc2) + 
         dy[1] * dx[0] * log((dx2[0] + dy2[1]) * invc2);
    const double beta1 = atan2(dy[0], dx[0]);
    const double beta2 = atan2(dx[1], dy[0]);
    const double beta3 = atan2(dy[1], dx[1]);
    const double beta4 = atan2(dx[0], dy[1]); 
    v2 = dx2[0] * (beta1 - beta4) + dy2[0] * (beta2 - beta1) +
         dx2[1] * (beta3 - beta2) + dy2[1] * (beta4 - beta3);
    v2 += HalfPi * (dx2[0] + dy2[0] + dx2[1] + dy2[1]); 
  }
  return -InvTwoPiEpsilon0 * 0.5 * (v1 + v2 - v0);
}  

void ComponentNeBem2d::BoxField(const double a, const double b,
                                const double x, const double y,
                                double& ex, double& ey) const {
  const std::array<double, 2> dx = {x - a, x + a};
  const std::array<double, 2> dy = {y - b, y + b};
  const std::array<double, 2> dx2 = {dx[0] * dx[0], dx[1] * dx[1]};
  const std::array<double, 2> dy2 = {dy[0] * dy[0], dy[1] * dy[1]};
  const double r1 = dx2[0] + dy2[0];
  const double r2 = dx2[1] + dy2[0];
  const double r3 = dx2[1] + dy2[1];
  const double r4 = dx2[0] + dy2[1];
  ex = 0.5 * (dy[0] * log(r1 / r2) + dy[1] * log(r3 / r4));
  ey = 0.5 * (dx[0] * log(r1 / r4) + dx[1] * log(r3 / r2));
  if (fabs(x) > a || fabs(y) > b) {
    std::array<double, 4> alpha = {atan2(dy[0], dx[0]), atan2(dy[0], dx[1]),
                                   atan2(dy[1], dx[1]), atan2(dy[1], dx[0])};
    if (x < 0.) {
      for (size_t i = 0; i < 4; ++i) if (alpha[i] < 0.) alpha[i] += TwoPi;
    }
    ex += dx[0] * (alpha[0] - alpha[3]) + dx[1] * (alpha[2] - alpha[1]);
    ey -= dy[0] * (alpha[0] - alpha[1]) + dy[1] * (alpha[2] - alpha[3]);
  } else {
    const double beta1 = atan2(-dy[0], -dx[0]);
    const double beta2 = atan2( dx[1], -dy[0]);
    const double beta3 = atan2( dy[1],  dx[1]);
    const double beta4 = atan2(-dx[0],  dy[1]); 
    ex += dx[0] * (HalfPi + beta1 - beta4) + dx[1] * (HalfPi + beta3 - beta2);
    ey -= dy[0] * (beta1 - beta2 - HalfPi) + dy[1] * (beta3 - beta4 - HalfPi);
  }
  ex *= InvTwoPiEpsilon0;
  ey *= InvTwoPiEpsilon0;
}

void ComponentNeBem2d::Reset() {
  m_regions.clear();
  m_segments.clear();
  m_wires.clear();
  m_elements.clear();
  m_spaceCharge.clear();
  m_ready = false;
}

void ComponentNeBem2d::UpdatePeriodicity() {
  std::cerr << m_className << "::UpdatePeriodicity:\n"
            << "    Periodicities are not supported.\n";
}

void ComponentNeBem2d::ToLocal(const double xIn, const double yIn, 
                               const double cphi, const double sphi,
                               double& xOut, double& yOut) const {
  xOut = +cphi * xIn + sphi * yIn;
  yOut = -sphi * xIn + cphi * yIn;
}

void ComponentNeBem2d::ToGlobal(const double xIn, const double yIn,
                                const double cphi, const double sphi,
                                double& xOut, double& yOut) const {
  xOut = cphi * xIn - sphi * yIn;
  yOut = sphi * xIn + cphi * yIn;
}

}
