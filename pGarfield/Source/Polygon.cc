#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "Garfield/Polygon.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"

namespace {

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

/// Determines whether the 2 straight lines (x1, y1) to (x2, y2)
/// and (u1, u2) to (v1, v2) cross at an intermediate point for both lines.
bool Crossing(const double x1, const double y1, 
              const double x2, const double y2,
              const double u1, const double v1,
              const double u2, const double v2) {

  std::array<std::array<double, 2>, 2> a;
  // Matrix to compute the crossing point.
  a[0][0] = y2 - y1;
  a[1][0] = v2 - v1;
  a[0][1] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  // Set tolerances.
  const double epsx = std::max(1.e-10 * std::max({std::abs(x1), std::abs(x2),
                                                  std::abs(u1), std::abs(u2)}),
                               1.e-10);
  const double epsy = std::max(1.e-10 * std::max({std::abs(y1), std::abs(y2),
                                                  std::abs(v1), std::abs(v2)}),
                               1.e-10);
  // Check for a point of one line located on the other line.
  if (OnLine(x1, y1, x2, y2, u1, v1) || OnLine(x1, y1, x2, y2, u2, v2) || 
      OnLine(u1, v1, u2, v2, x1, y1) || OnLine(u1, v1, u2, v2, x2, y2)) {
    // Point on other line.
    return true;
  }
  if (std::abs(det) < epsx * epsy) {
    // Parallel, non-touching
    return false;
  }
  // Crossing, non-trivial lines: solve crossing equations.
  std::swap(a[0][0], a[1][1]);
  const double invdet = 1. / det;
  a[0][0] *= invdet;
  a[1][1] *= invdet;
  a[0][1] = -a[0][1] * invdet;
  a[1][0] = -a[1][0] * invdet;
  // Compute crossing point.
  const double xc = a[0][0] * (x1 * y2 - x2 * y1) + 
                    a[0][1] * (u1 * v2 - u2 * v1);
  const double yc = a[1][0] * (x1 * y2 - x2 * y1) +
                    a[1][1] * (u1 * v2 - u2 * v1);
  // See whether the crossing point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc)) {
    return true;
  }
  // Crossing point not on both lines.
  return false;
}

}

namespace Garfield {

namespace Polygon {

void Inside(const std::vector<double>& xpl, const std::vector<double>& ypl,
            const double x, const double y, bool& inside, bool& edge) {
  // Initial settings.
  inside = false;
  edge = false;
  const unsigned int npl = xpl.size();
  if (ypl.size() != npl) return;
  // Special treatment for few points.
  if (npl < 2) {
    return;
  } else if (npl == 2) {
    edge = OnLine(xpl[0], ypl[0], xpl[1], ypl[1], x, y);
    return;
  }
  // Determine the range of the data.
  const double xmin = *std::min_element(std::begin(xpl), std::end(xpl));
  const double xmax = *std::max_element(std::begin(xpl), std::end(xpl));
  const double ymin = *std::min_element(std::begin(ypl), std::end(ypl));
  const double ymax = *std::max_element(std::begin(ypl), std::end(ypl));

  // Set tolerances.
  double epsx = 1.e-8 * std::max(fabs(xmin), fabs(xmax));
  double epsy = 1.e-8 * std::max(fabs(ymin), fabs(ymax));
  epsx = std::max(epsx, 1.e-8);
  epsy = std::max(epsy, 1.e-8);

  // Ensure that we have a range.
  if (fabs(xmax - xmin) <= epsx) {
    if (y >= ymin - epsy && y <= ymax + epsy &&
        fabs(xmax + xmin - 2 * x) <= epsx) {
      edge = true;
    } else {
      edge = false;
    }
  } else if (fabs(ymax - ymin) <= epsy) {
    if (x >= xmin - epsx && x <= xmax + epsx &&
        fabs(ymax + ymin - 2 * y) <= epsy) {
      edge = true;
    } else {
      edge = false;
    }
  }
  // Choose a point at "infinity".
  double xinf = xmin - fabs(xmax - xmin);
  double yinf = ymin - fabs(ymax - ymin);

  unsigned int nIter = 0;
  bool ok = false;
  while (!ok && nIter < 100) {
    ok = true;
    // Loop over the edges counting intersections.
    unsigned int nCross = 0;
    for (unsigned int j = 0; j < npl; ++j) {
      const unsigned int jj = j < npl - 1 ? j + 1 : 0;
      // Flag points located on one of the edges.
      if (OnLine(xpl[j], ypl[j], xpl[jj], ypl[jj], x, y)) {
        edge = true;
        return;
      }
      // Count mid-line intersects.
      double xc = 0., yc = 0.;
      if (Crossing(x, y, xinf, yinf, xpl[j], ypl[j], xpl[jj], ypl[jj], xc,
                   yc)) {
        ++nCross;
      }
      // Ensure that the testing line doesn't cross a corner.
      if (OnLine(x, y, xinf, yinf, xpl[j], ypl[j])) {
        xinf = xmin - ::Garfield::RndmUniform() * fabs(xmax - xinf);
        yinf = ymin - ::Garfield::RndmUniform() * fabs(ymax - yinf);
        ok = false;
        break;
      }
    }
    if (ok) {
      // Set the INSIDE flag.
      if (nCross != 2 * (nCross / 2)) inside = true;
      return;
    }
    ++nIter;
  }

  std::cerr << "Polygon::Inside:\n    Warning. Unable to verify "
            << "whether a point is internal; setting to edge.\n";
  inside = false;
  edge = true;
}

double Area(const std::vector<double>& xp, const std::vector<double>& yp) {

  double f = 0.;
  const unsigned int n = xp.size();
  for (unsigned int i = 0; i < n; ++i) {
    const unsigned int ii = i < n - 1 ? i + 1 : 0;
    f += xp[i] * yp[ii] - xp[ii] * yp[i];  
  }
  return 0.5 * f; 
}

bool NonTrivial(const std::vector<double>& xp, 
                const std::vector<double>& yp) {

  // -----------------------------------------------------------------------
  // PLACHK - Checks whether a set of points builds a non-trivial 
  //           polygon in the (x,y) plane.
  // -----------------------------------------------------------------------

  // First check number of points.
  if (xp.size() != yp.size()) return false;
  const unsigned int np = xp.size();
  if (xp.size() < 3) return false;
  // Find a second point at maximum distance of the first.
  double d1 = 0.;
  double x1 = 0.;
  double y1 = 0.;
  unsigned int i1 = 0;
  double xmin = xp[0];
  double ymin = yp[0];
  double xmax = xp[0];
  double ymax = yp[0];
  for (unsigned int i = 1; i < np; ++i) {
    xmin = std::min(xmin, xp[i]);
    ymin = std::min(ymin, yp[i]);
    xmax = std::max(xmax, xp[i]);
    ymax = std::max(ymax, yp[i]);
    const double dx = xp[i] - xp[0];
    const double dy = yp[i] - yp[0];
    const double d = dx * dx + dy * dy;
    if (d > d1) {
      x1 = dx;
      y1 = dy;
      i1 = i;
      d1 = d;
    }
  }
  // Set tolerances.
  double epsx = 1.e-6 * (std::abs(xmax) + std::abs(xmin));
  double epsy = 1.e-6 * (std::abs(ymax) + std::abs(ymin));
  if (epsx <= 0) epsx = 1.e-6;
  if (epsy <= 0) epsy = 1.e-6;
  // See whether there is a range at all.
  if (std::abs(xmax - xmin) <= epsx && std::abs(ymax - ymin) <= epsy) {
    // Single point.
    return false;
  }
  // See whether there is a second point.
  if (d1 <= epsx * epsx + epsy * epsy || i1 == 0) return false;

  // Find a third point maximising the external product.
  double d2 = 0.; 
  unsigned int i2 = 0;
  for (unsigned int i = 1; i < np; ++i) {
    if (i == i1) continue;
    const double dx = xp[i] - xp[0];
    const double dy = yp[i] - yp[0];
    const double d = std::abs(x1 * dy - y1 * dx);
    if (d > d2) {
      d2 = d;
      i2 = i;
    }
  }
  if (d2 <= epsx * epsy || i2 <= 0) {
    // No third point.
    return false;
  }
  // Seems to be OK.
  return true;
}

void EliminateButterflies(std::vector<double>& xp, std::vector<double>& yp,
                          std::vector<double>& zp) {
  //----------------------------------------------------------------------
  //   BUTFLD - Tries to eliminate "butterflies", i.e. the crossing of 2
  //            adjacent segments of a polygon, by point exchanges.
  //----------------------------------------------------------------------
  unsigned int np = xp.size();
  if (np <= 3) return;
  // Compute range.
  const double xmin = *std::min_element(std::begin(xp), std::end(xp));
  const double xmax = *std::max_element(std::begin(xp), std::end(xp));
  const double ymin = *std::min_element(std::begin(yp), std::end(yp));
  const double ymax = *std::max_element(std::begin(yp), std::end(yp));
  const double zmin = *std::min_element(std::begin(zp), std::end(zp));
  const double zmax = *std::max_element(std::begin(zp), std::end(zp));
  // Set tolerances.
  const double epsx = std::max(1.e-10 * std::abs(xmax - xmin), 1.e-6);
  const double epsy = std::max(1.e-10 * std::abs(ymax - ymin), 1.e-6);
  const double epsz = std::max(1.e-10 * std::abs(zmax - zmin), 1.e-6);
  double xsurf = 0.;
  double ysurf = 0.;
  double zsurf = 0.;
  for (unsigned int i = 2; i < np; ++i) {
    const double x1 = xp[i - 1] - xp[0];
    const double y1 = yp[i - 1] - yp[0];
    const double z1 = zp[i - 1] - zp[0];
    xsurf += std::abs((yp[i] - yp[0]) * z1 - y1 * (zp[i] - zp[0]));
    ysurf += std::abs((xp[i] - xp[0]) * z1 - x1 * (zp[i] - zp[0]));
    zsurf += std::abs((xp[i] - xp[0]) * y1 - x1 * (yp[i] - yp[0]));
  }
  // Eliminate points appearing twice, initialise marks.
  std::vector<bool> mark(np, false);
  // Scan the list.
  for (unsigned int i = 0; i < np; ++i) {
    if (mark[i]) continue;
    for (unsigned int j = i + 1; j < np; ++j) {
      if (std::abs(xp[i] - xp[j]) <= epsx && std::abs(yp[i] - yp[j]) <= epsy &&
          std::abs(zp[i] - zp[j]) <= epsz) {
        mark[j] = true;
      }
    }
  }
  // And remove the duplicate points.
  unsigned int nNew = 0;
  for (unsigned int i = 0; i < np; ++i) {
    if (mark[i]) continue;
    xp[nNew] = xp[i];
    yp[nNew] = yp[i];
    zp[nNew] = zp[i];
    ++nNew;
  }
  // std::cout << "ElminateButterflies: old/new number of points: " << np 
  //           << "/" << nNew << "\n";
  // Update the number of points.
  np = nNew;
  xp.resize(np);
  yp.resize(np);
  zp.resize(np);
  // No risk of having a butterfly with less than 4 points.
  if (np <= 3) return;
  // Select the axis with the largest norm.
  unsigned int iaxis = 0;
  if (xsurf > ysurf && xsurf > zsurf) {
    iaxis = 1;
  } else if (ysurf > zsurf) {
    iaxis = 2; 
  } else {
    iaxis = 3;
  }
  // Set number of passes to avoid endless loop.
  unsigned int nPass = 0;
  bool repass = true;
  while (repass) {
    // Make a pass.
    ++nPass;
    repass = false;
    for (unsigned int i = 0; i < np; ++i) {
      const unsigned int ii = (i + 1) % np;
      for (unsigned int j = i + 2; j < np; ++j) {
        const unsigned int jj = (j + 1) % np;
        if (j + 1 >= np && jj >= i) continue;
        // Check for a crossing.
        if (iaxis == 1 && !Crossing(yp[i], zp[i], yp[ii], zp[ii],
                                    yp[j], zp[j], yp[jj], zp[jj])) {
          continue;
        } else if (iaxis == 2 && !Crossing(xp[i], zp[i], xp[ii], zp[ii],
                                           xp[j], zp[j], xp[jj], zp[jj])) {
          continue;
        } else if (iaxis == 3 && !Crossing(xp[i], yp[i], xp[ii], yp[ii],
                                           xp[j], yp[j], xp[jj], yp[jj])) {
          continue;
        }
        // If there is a crossing, exchange the portion in between.
        for (unsigned int k = 0; k < (j - i) / 2; ++k) {
          const unsigned int k1 = (i + k + 1) % np;
          const unsigned int k2 = (j - k) % np;
          std::swap(xp[k1], xp[k2]);
          std::swap(yp[k1], yp[k2]);
          std::swap(zp[k1], zp[k2]);
        }
        // And remember we have to do another pass after this.
        repass = true;
      }
    }
    // See whether we have to do another pass.
    if (repass && nPass > np) {
      std::cerr << "Butterfly: unable to eliminate the internal crossings;\n"
                << "    plot probably incorrect.\n";
      break;
    }
  }
}

}

}
