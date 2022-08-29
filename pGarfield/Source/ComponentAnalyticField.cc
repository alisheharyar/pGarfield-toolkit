#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <TCanvas.h>
#include <TGraph.h>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/ViewBase.hh"

namespace {

constexpr double Internal2Newton = 1.e-15 * Garfield::TwoPiEpsilon0 * 100.;

std::pair<std::complex<double>, std::complex<double> > Th1(
    const std::complex<double>& zeta, const double p1, const double p2) {
  const std::complex<double> zsin = sin(zeta);
  const std::complex<double> zcof = 4. * zsin * zsin - 2.;
  std::complex<double> zu = -p1 - zcof * p2;
  std::complex<double> zunew = 1. - zcof * zu - p2;
  const std::complex<double> zterm1 = (zunew + zu) * zsin;
  zu = -3. * p1 - zcof * 5. * p2;
  zunew = 1. - zcof * zu - 5. * p2;
  const std::complex<double> zterm2 = (zunew - zu) * cos(zeta);
  return std::make_pair(std::move(zterm1), std::move(zterm2));
}

// Transformation from Cartesian and polar coordinates.
void Cartesian2Polar(const double x, const double y, double& r,
                     double& theta) {
  if (x == 0. && y == 0.) {
    r = theta = 0.;
    return;
  }
  r = sqrt(x * x + y * y);
  theta = atan2(y, x) * Garfield::RadToDegree;
}

// Transformation from polar to Cartesian coordinates.
void Polar2Cartesian(const double r, const double theta, double& x,
                     double& y) {
  const double thetap = theta * Garfield::DegreeToRad;
  x = r * cos(thetap);
  y = r * sin(thetap);
}

// Transformation (rho, phi) to (r, theta) via the map
// (r, theta) = (exp(rho), 180 * phi / Pi).
void Internal2Polar(const double rho, const double phi, double& r,
                    double& theta) {
  // CFMRTP
  r = exp(rho);
  theta = Garfield::RadToDegree * phi;
}

// Transformation (r, theta) to (rho, phi) via the map
// (rho, phi) = (log(r), Pi * theta / 180).
void Polar2Internal(const double r, const double theta, double& rho,
                    double& phi) {
  // CFMPTR
  rho = r > 0. ? log(r) : -25.;
  phi = Garfield::DegreeToRad * theta;
}

// Transformation (x, y) to (rho, phi) via the conformal map
// (x, y) = exp(rho, phi)
void Cartesian2Internal(const double x, const double y, double& rho,
                        double& phi) {
  // CFMCTR
  if (x == 0 && y == 0) {
    rho = -25.;
    phi = 0.;
    return;
  }
  // const std::complex<double> z = log(std::complex<double>(x, y));
  // rho = real(z);
  // phi = imag(z);
  rho = 0.5 * log(x * x + y * y);
  phi = atan2(y, x); 
}

// Transformation (rho, phi) to (x, y) via the conformal map
// (x, y) = exp(rho, phi)
void Internal2Cartesian(const double rho, const double phi, double& x, 
                        double& y) {
  // CFMRTC
  // const std::complex<double> z = exp(std::complex<double>(rho, phi));
  // x = real(z);
  // y = imag(z);
  const double r = exp(rho);
  x = r * cos(phi);
  y = r * sin(phi);
}

bool FitDipoleMoment(const std::vector<double>& angle, 
                     const std::vector<double>& volt,
                     double& ampdip, double& phidip, const bool dbg) {
  //-----------------------------------------------------------------------
  //   DIPFIT - Determines the dipole moment of a wire.
  //-----------------------------------------------------------------------

  // Initial values.
  phidip = 0.;
  ampdip = 0.;
  const unsigned int n = angle.size();
  // Initial search for a maximum.
  double phiMax = 0.;
  double sumMax = 0.;
  constexpr unsigned int nTry = 100;
  std::array<double, nTry> phiTry;
  std::array<double, nTry> sumTry;
  for (unsigned int i = 0; i < nTry; ++i) {
    // Make the internal product with a shifted cosine.
    phiTry[i] = i * Garfield::TwoPi / nTry;
    sumTry[i] = 0.;
    for (unsigned int j = 0; j < n; ++j) {
      sumTry[i] += volt[j] * cos(phiTry[i] - angle[j]);
    }
    sumTry[i] *= 2. / n;
    // See whether this one beats earlier.
    if (sumTry[i] > sumMax) {
      phiMax = phiTry[i];
      sumMax = sumTry[i];
    }
  }
  if (dbg) {
    std::printf("    Maximum of scan at phi = %12.5f, product = %12.5f\n",
                phiMax, sumMax);
    // CALL GRGRPH(phiTry, sumTry, nTry, 'Angle [radians]',
    //             'Cos projection', 'Search of maximum')
  }
  phidip = phiMax;
  ampdip = sumMax;
  // Scan in the neighbourbood
  constexpr double eps = 0.1;
  double x1 = phiMax - eps;
  double x2 = phiMax;
  double x3 = phiMax + eps;
  double f1 = 0.;
  double f2 = sumMax;
  double f3 = 0.;
  for (unsigned int j = 0; j < n; ++j) {
    f1 += volt[j] * cos(x1 - angle[j]);
    f3 += volt[j] * cos(x3 - angle[j]);
  }
  f1 *= 2. / n;
  f3 *= 2. / n;
  // Refine the estimate by parabolic extremum search.
  const double epsf = 1.e-3 * sumMax;
  constexpr double epsx = 1.e-3 * Garfield::TwoPi;
  constexpr unsigned int nMaxIter = 10;
  for (unsigned int i = 0; i < nMaxIter; ++i) {
    if (dbg) std::cout << "    Start of iteration " << i << ".\n";
    // Estimate parabolic extremum.
    const double det = (f1 - f2) * x3 + (f3 - f1) * x2 + (f2 - f3) * x1;
    if (std::abs(det) <= 0.) {
      std::cerr << "    Warning: Determinant = 0; parabolic search stopped.\n";
      phidip = x2;
      ampdip = f2;
      return false;
    }
    const double xp = ((f1 - f2) * x3 * x3 + (f3 - f1) * x2 * x2 + 
                       (f2 - f3) * x1 * x1) / (2 * det);
    double fp = 0.;
    for (unsigned int j = 0; j < n; ++j) {
      fp += volt[j] * cos(xp - angle[j]);
    }
    fp *= 2. / n;
    // Debugging output.
    if (dbg) {
      std::printf("    Point 1:  x = %15.8f f = %15.8f\n", x1, f1);
      std::printf("    Point 2:  x = %15.8f f = %15.8f\n", x2, f2);
      std::printf("    Point 3:  x = %15.8f f = %15.8f\n", x3, f3);
      std::printf("    Parabola: x = %15.8f f = %15.8f\n", xp, fp);
    }
    // Check that the new estimate doesn't coincide with an old point.
    const double tol = epsx * (epsx + std::abs(xp));
    if (fabs(xp - x1) < tol || fabs(xp - x2) < tol || fabs(xp - x3) < tol) {
      if (dbg) {
        std::cout << "    Location convergence criterion satisfied.\n";
      }
      phidip = xp;
      ampdip = fp;
      return true;
    }
    // Check convergence.
    if (std::abs(fp - f1) < epsf * (std::abs(fp) + std::abs(f1) + epsf)) {
      if (dbg) {
        std::cout << "    Function value convergence criterion satisfied.\n";
      }
      phidip = xp;
      ampdip = fp;
      return true;
    }
    // Store the value in the table.
    if (fp > f1) {
      f3 = f2;
      x3 = x2;
      f2 = f1;
      x2 = x1;
      f1 = fp;
      x1 = xp;
    } else if (fp > f2) {
      f3 = f2;
      x3 = x2;
      f2 = fp;
      x2 = xp;
    } else if (fp > f3) {
      f3 = fp;
      x3 = xp;
    } else {
      std::cerr << "    Warning: Parabolic extremum is worse "
                << "than current optimum; search stopped.\n";
      std::printf("    Point 1:  x = %15.8f f = %15.8f\n", x1, f1);
      std::printf("    Point 2:  x = %15.8f f = %15.8f\n", x2, f2);
      std::printf("    Point 3:  x = %15.8f f = %15.8f\n", x3, f3);
      std::printf("    Parabola: x = %15.8f f = %15.8f\n", xp, fp);
      phidip = x2;
      ampdip = f2;
      return false;
    }
  }
  // No convergence.
  std::cerr << "    Warning: No convergence after maximum number of steps.\n"
            << "    Current extremum f = " << f2 << "\n"
            << "    Found for        x = " << x2 << "\n";
  phidip = x2;
  ampdip = f2;
  return false;
}

double MirrorCoordinate(const double x, const double xp, const double xw,
                        const double sx) { 

  // Find the plane nearest to the wire.
  double cx = xp - sx * round((xp - xw) / sx);
  return 2 * cx - x - xw;
}

bool SameSide(const double x0, const double x1, const double xp) {

  return ((x0 <= xp && x1 <= xp) || (x0 >= xp && x1 >= xp));
}

}  // namespace

namespace Garfield {

using CMatrix = std::vector<std::vector<std::complex<double> > >;
using DMatrix = std::vector<std::vector<double> >;

ComponentAnalyticField::ComponentAnalyticField() : Component("AnalyticField") {
  CellInit();
}

Medium* ComponentAnalyticField::GetMedium(const double xin, const double yin,
                                          const double zin) {
  if (m_geometry) return m_geometry->GetMedium(xin, yin, zin);

  // Make sure the cell is prepared.
  if (!m_cellset && !Prepare()) return nullptr;

  double xpos = xin, ypos = yin;
  if (m_polar) Cartesian2Internal(xin, yin, xpos, ypos);
  // In case of periodicity, move the point into the basic cell.
  if (m_perx) xpos -= m_sx * round(xin / m_sx);

  double arot = 0.;
  if (m_pery && m_tube) {
    Cartesian2Polar(xin, yin, xpos, ypos);
    arot = RadToDegree * m_sy * round(DegreeToRad * ypos / m_sy);
    ypos -= arot;
    Polar2Cartesian(xpos, ypos, xpos, ypos);
  } else if (m_pery) {
    ypos -= m_sy * round(ypos / m_sy);
  }

  // Move the point to the correct side of the plane.
  if (m_perx && m_ynplan[0] && xpos <= m_coplan[0]) xpos += m_sx;
  if (m_perx && m_ynplan[1] && xpos >= m_coplan[1]) xpos -= m_sx;
  if (m_pery && m_ynplan[2] && ypos <= m_coplan[2]) ypos += m_sy;
  if (m_pery && m_ynplan[3] && ypos >= m_coplan[3]) ypos -= m_sy;

  // In case (xpos, ypos) is located behind a plane there is no field.
  if (m_tube) {
    if (!InTube(xpos, ypos, m_cotube, m_ntube)) return nullptr;
  } else {
    if ((m_ynplan[0] && xpos < m_coplan[0]) ||
        (m_ynplan[1] && xpos > m_coplan[1]) ||
        (m_ynplan[2] && ypos < m_coplan[2]) ||
        (m_ynplan[3] && ypos > m_coplan[3])) {
      return nullptr;
    }
  }

  // If (xpos, ypos) is within a wire, there is no field either.
  for (const auto& wire : m_w) {
    double dx = xpos - wire.x;
    double dy = ypos - wire.y;
    // Correct for periodicities.
    if (m_perx) dx -= m_sx * round(dx / m_sx);
    if (m_pery) dy -= m_sy * round(dy / m_sy);
    // Check the actual position.
    if (dx * dx + dy * dy < wire.r * wire.r) return nullptr;
  }
  return m_medium;
}

bool ComponentAnalyticField::GetVoltageRange(double& pmin, double& pmax) {
  // Make sure the cell is prepared.
  if (!m_cellset && !Prepare()) {
    std::cerr << m_className << "::GetVoltageRange: Cell not set up.\n";
    return false;
  }

  pmin = m_vmin;
  pmax = m_vmax;
  return true;
}

bool ComponentAnalyticField::GetBoundingBox(
    double& x0, double& y0, double& z0,
    double& x1, double& y1, double& z1) {
  // If a geometry is present, try to get the bounding box from there.
  if (m_geometry) {
    if (m_geometry->GetBoundingBox(x0, y0, z0, x1, y1, z1)) return true;
  }
  // Otherwise, return the cell dimensions.
  return GetElementaryCell(x0, y0, z0, x1, y1, z1);
}

bool ComponentAnalyticField::GetElementaryCell(
    double& x0, double& y0, double& z0,
    double& x1, double& y1, double& z1) {
  if (!m_cellset && !Prepare()) return false;
  if (m_polar) {
    double rmax, thetamax;
    Internal2Polar(m_xmax, m_ymax, rmax, thetamax);
    x0 = -rmax;
    y0 = -rmax;
    x1 = +rmax;
    y1 = +rmax;
  } else {
    x0 = m_xmin;
    y0 = m_ymin;
    x1 = m_xmax;
    y1 = m_ymax;
  }
  z0 = m_zmin;
  z1 = m_zmax;
  return true;
}

void ComponentAnalyticField::PrintCell() {
  //-----------------------------------------------------------------------
  //   CELPRT - Subroutine printing all available information on the cell.
  //-----------------------------------------------------------------------

  // Make sure the cell is prepared.
  if (!m_cellset && !Prepare()) {
    std::cerr << m_className << "::PrintCell: Cell not set up.\n";
    return;
  }
  std::cout << m_className
            << "::PrintCell: Cell identification: " << GetCellType() << "\n";
  // Print positions of wires, applied voltages and resulting charges.
  if (!m_w.empty()) {
    std::cout << "  Table of the wires\n";
    if (m_polar) {
      std::cout << "  Nr    Diameter     r       phi     Voltage";
    } else {
      std::cout << "  Nr    Diameter     x        y      Voltage";
    }
    std::cout << "      Charge    Tension   Length   Density  Label\n";
    if (m_polar) {
      std::cout << "       [micron]     [cm]     [deg]    [Volt]";
    } else {
      std::cout << "       [micron]     [cm]     [cm]     [Volt]";
    }
    std::cout << "     [pC/cm]       [g]      [cm]    [g/cm3]\n";
    for (unsigned int i = 0; i < m_nWires; ++i) {
      const auto& w = m_w[i];
      double xw = w.x;
      double yw = w.y;
      double dw = 2 * w.r;
      if (m_polar) {
        Internal2Polar(w.x, w.y, xw, yw);
        dw *= xw;
      }
      std::printf(
          "%4d %9.2f %9.4f %9.4f %9.3f %12.4f %9.2f %9.2f %9.2f \"%s\"\n", i,
          1.e4 * dw, xw, yw, w.v, w.e * TwoPiEpsilon0 * 1.e-3, w.tension,
          w.u, w.density, w.type.c_str());
    }
  }
  // Print information on the tube if present.
  if (m_tube) {
    std::string shape;
    if (m_ntube == 0) {
      shape = "Circular";
    } else if (m_ntube == 3) {
      shape = "Triangular";
    } else if (m_ntube == 4) {
      shape = "Square";
    } else if (m_ntube == 5) {
      shape = "Pentagonal";
    } else if (m_ntube == 6) {
      shape = "Hexagonal";
    } else if (m_ntube == 7) {
      shape = "Heptagonal";
    } else if (m_ntube == 8) {
      shape = "Octagonal";
    } else {
      shape = "Polygonal with " + std::to_string(m_ntube) + " corners";
    }
    std::cout << "  Enclosing tube\n"
              << "    Potential:  " << m_vttube << " V\n"
              << "    Radius:     " << m_cotube << " cm\n"
              << "    Shape:      " << shape << "\n"
              << "    Label:      " << m_planes[4].type << "\n";
  }
  // Print data on the equipotential planes.
  if (m_ynplan[0] || m_ynplan[1] || m_ynplan[2] || m_ynplan[3]) {
    std::cout << "  Equipotential planes\n";
    // First those at const x or r.
    const std::string xr = m_polar ? "r" : "x";
    if (m_ynplan[0] && m_ynplan[1]) {
      std::cout << "    There are two planes at constant " << xr << ":\n";
    } else if (m_ynplan[0] || m_ynplan[1]) {
      std::cout << "    There is one plane at constant " << xr << ":\n";
    }
    for (unsigned int i = 0; i < 2; ++i) {
      if (!m_ynplan[i]) continue;
      if (m_polar) {
        std::cout << "    r = " << exp(m_coplan[i]) << " cm, ";
      } else {
        std::cout << "    x = " << m_coplan[i] << " cm, ";
      }
      if (fabs(m_vtplan[i]) > 1.e-4) {
        std::cout << "potential = " << m_vtplan[i] << " V, ";
      } else {
        std::cout << "earthed, ";
      }
      const auto& plane = m_planes[i];
      if (plane.type.empty() && plane.type != "?") {
        std::cout << "label = " << plane.type << ", ";
      }
      const unsigned int nStrips = plane.strips1.size() + plane.strips2.size();
      const unsigned int nPixels = plane.pixels.size();
      if (nStrips == 0 && nPixels == 0) {
        std::cout << "no strips or pixels.\n";
      } else if (nPixels == 0) {
        std::cout << nStrips << " strips.\n";
      } else if (nStrips == 0) {
        std::cout << nPixels << " pixels.\n";
      } else {
        std::cout << nStrips << " strips, " << nPixels << " pixels.\n";
      }
      for (const auto& strip : plane.strips2) {
        std::cout << "      ";
        if (m_polar) {
          double gap = i == 0 ? expm1(strip.gap) : -expm1(-strip.gap);
          gap *= exp(m_coplan[i]);
          std::cout << RadToDegree * strip.smin << " < phi < "
                    << RadToDegree * strip.smax
                    << " degrees, gap = " << gap << " cm";
        } else {
          std::cout << strip.smin << " < y < " << strip.smax
                    << " cm, gap = " << strip.gap << " cm";
        }
        if (!strip.type.empty() && strip.type != "?") {
          std::cout << " (\"" << strip.type << "\")";
        }
        std::cout << "\n";
      }
      for (const auto& strip : plane.strips1) {
        std::cout << "      " << strip.smin << " < z < " << strip.smax;
        if (m_polar) {
          double gap = i == 0 ? expm1(strip.gap) : -expm1(-strip.gap);
          gap *= exp(m_coplan[i]);
          std::cout << " cm, gap = " << gap << " cm";
        } else {
          std::cout << " cm, gap = " << strip.gap << " cm";
        }
        if (!strip.type.empty() && strip.type != "?") {
          std::cout << " (\"" << strip.type << "\")";
        }
        std::cout << "\n";
      }
      for (const auto& pix : plane.pixels) {
        std::cout << "      ";
        if (m_polar) {
          std::cout << RadToDegree * pix.smin << " < phi < "
                    << RadToDegree * pix.smax << " degrees, "; 
        } else {
          std::cout << pix.smin << " < y < " << pix.smax << " cm, ";
        }
        std::cout << pix.zmin << " < z < " << pix.zmax << " cm, gap = ";
        if (m_polar) {
          double gap = i == 0 ? expm1(pix.gap) : -expm1(-pix.gap);
          gap *= exp(m_coplan[i]);
          std::cout << gap << " cm";
        } else {
          std::cout << pix.gap << " cm";
        }
        if (!pix.type.empty() && pix.type != "?") {
          std::cout << " (\"" << pix.type << "\")";
        }
        std::cout << "\n";
      }
    }
    // Next the planes at constant y or phi
    const std::string yphi = m_polar ? "phi" : "y";
    if (m_ynplan[2] && m_ynplan[3]) {
      std::cout << "    There are two planes at constant " << yphi << ":\n";
    } else if (m_ynplan[2] || m_ynplan[3]) {
      std::cout << "    There is one plane at constant " << yphi << ":\n";
    }
    for (unsigned int i = 2; i < 4; ++i) {
      if (!m_ynplan[i]) continue;
      if (m_polar) {
        std::cout << "    phi = " << RadToDegree * m_coplan[i] << " degrees, ";
      } else {
        std::cout << "    y = " << m_coplan[i] << " cm, ";
      }
      if (fabs(m_vtplan[i]) > 1.e-4) {
        std::cout << "potential = " << m_vtplan[i] << " V, ";
      } else {
        std::cout << "earthed, ";
      }
      const auto& plane = m_planes[i];
      if (plane.type.empty() && plane.type != "?") {
        std::cout << "label = " << plane.type << ", ";
      }
      const unsigned int nStrips = plane.strips1.size() + plane.strips2.size();
      const unsigned int nPixels = plane.pixels.size();
      if (nStrips == 0 && nPixels == 0) {
        std::cout << "no strips or pixels.\n";
      } else if (nPixels == 0) {
        std::cout << nStrips << " strips.\n";
      } else if (nStrips == 0) {
        std::cout << nPixels << " pixels.\n";
      } else {
        std::cout << nStrips << " strips, " << nPixels << " pixels.\n";
      }
      for (const auto& strip : plane.strips2) {
        std::cout << "      ";
        if (m_polar) {
          std::cout << exp(strip.smin) << " < r < " << exp(strip.smax) 
                    << " cm, gap = " << RadToDegree * strip.gap << " degrees";
        } else {
          std::cout << strip.smin << " < x < " << strip.smax
                    << " cm, gap = " << strip.gap << " cm";
        }
        if (!strip.type.empty() && strip.type != "?") {
          std::cout << " (\"" << strip.type << "\")";
        }
        std::cout << "\n";
      }
      for (const auto& strip : plane.strips1) {
        std::cout << "      " << strip.smin << " < z < " << strip.smax;
        if (m_polar) {
          std::cout << " cm, gap = " << RadToDegree * strip.gap << " degrees";
        } else {
          std::cout << " cm, gap = " << strip.gap << " cm";
        }
        if (!strip.type.empty() && strip.type != "?") {
          std::cout << " (\"" << strip.type << "\")";
        }
        std::cout << "\n";
      }
      for (const auto& pix : plane.pixels) {
        std::cout << "      ";
        if (m_polar) {
          std::cout << exp(pix.smin) << " < r < " << exp(pix.smax) << " cm, "; 
        } else {
          std::cout << pix.smin << " < x < " << pix.smax << " cm, ";
        }
        std::cout << pix.zmin << " < z < " << pix.zmax << " cm, gap = ";
        if (m_polar) {
          std::cout << RadToDegree * pix.gap << " degrees";
        } else {
          std::cout << pix.gap << " cm";
        }
        if (!pix.type.empty() && pix.type != "?") {
          std::cout << " (\"" << pix.type << "\")";
        }
        std::cout << "\n";
      }
    }
  }
  // Print the type of periodicity.
  std::cout << "  Periodicity\n";
  if (m_perx) {
    std::cout << "    The cell is repeated every ";
    if (m_polar) {
      std::cout << exp(m_sx) << " cm in r.\n";
    } else {
      std::cout << m_sx << " cm in x.\n";
    }
  } else {
    if (m_polar) {
      std::cout << "    The cell is not periodic in r.\n";
    } else {
      std::cout << "    The cell has no translation periodicity in x.\n";
    }
  }
  if (m_pery) {
    std::cout << "    The cell is repeated every ";
    if (m_polar) {
      std::cout << RadToDegree * m_sy << " degrees in phi.\n";
    } else {
      std::cout << m_sy << " cm in y.\n";
    }
  } else {
    if (m_polar) {
      std::cout << "    The cell is not periodic in phi.\n";
    } else {
      std::cout << "    The cell has no translation periodicity in y.\n";
    }
  }
  std::cout << "  Other data\n";
  std::cout << "    Gravity vector: (" << m_down[0] << ", " << m_down[1]
            << ", " << m_down[2] << ") [g].\n";
  std::cout << "  Cell dimensions:\n    ";
  if (!m_polar) {
    std::cout << m_xmin << " < x < " << m_xmax << " cm, " << m_ymin << " < y < "
              << m_ymax << " cm.\n";
  } else {
    double xminp, yminp;
    Internal2Polar(m_xmin, m_ymin, xminp, yminp);
    double xmaxp, ymaxp;
    Internal2Polar(m_xmax, m_ymax, xmaxp, ymaxp);
    std::cout << xminp << " < r < " << xmaxp << " cm, " << yminp << " < phi < "
              << ymaxp << " degrees.\n";
  }
  std::cout << "  Potential range:\n    " << m_vmin << " < V < " << m_vmax
            << " V.\n";
  // Print voltage shift in case no equipotential planes are present.
  if (!(m_ynplan[0] || m_ynplan[1] || m_ynplan[2] || m_ynplan[3] || m_tube)) {
    std::cout << "  All voltages have been shifted by " << m_v0
              << " V to avoid net wire charge.\n";
  } else {
    // Print the net charge on wires.
    double sum = 0.;
    for (const auto& w : m_w) sum += w.e;
    std::cout << "  The net charge on the wires is "
              << 1.e-3 * TwoPiEpsilon0 * sum << " pC/cm.\n";
  }
}

bool ComponentAnalyticField::CrossedWire(
    const double xx0, const double yy0, const double z0, 
    const double xx1, const double yy1, const double z1,
    double& xc, double& yc, double& zc, const bool centre, double& rc) {
  xc = xx0;
  yc = yy0;
  zc = z0;

  if (m_w.empty()) return false;

  double x0 = xx0;
  double y0 = yy0;
  double x1 = xx1;
  double y1 = yy1;
  if (m_polar) {
    Cartesian2Internal(xx0, yy0, x0, y0);
    Cartesian2Internal(xx1, yy1, x1, y1);
  }
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double d2 = dx * dx + dy * dy;
  // Check that the step length is non-zero.
  if (d2 < Small) return false;
  const double invd2 = 1. / d2;

  // Check if a whole period has been crossed.
  if ((m_perx && fabs(dx) >= m_sx) || (m_pery && fabs(dy) >= m_sy)) {
    std::cerr << m_className << "::CrossedWire:\n"
              << "    Particle crossed more than one period.\n";
    return false;
  }

  // Both coordinates are assumed to be located inside
  // the drift area and inside a drift medium.
  // This should have been checked before this call.

  const double xm = 0.5 * (x0 + x1);
  const double ym = 0.5 * (y0 + y1);
  double dMin2 = 0.;
  for (const auto& wire : m_w) {
    double xw = wire.x;
    if (m_perx) xw += m_sx * round((xm - xw) / m_sx);
    double yw = wire.y;
    if (m_pery) yw += m_sy * round((ym - yw) / m_sy);
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
    // Add in the times nTrap to account for the trap radius.
    const double r2 = wire.r * wire.r;
    if (dMin2 < r2) {
      // Wire has been crossed.
      if (centre) {
        if (m_polar) {
          Internal2Cartesian(xw, yw, xc, yc);
        } else {
          xc = xw;
          yc = yw;
        }
      } else {
        // Find the point of intersection.
        const double p = -xIn0 * invd2;
        const double q = (dw02 - r2) * invd2;
        const double s = sqrt(p * p - q);
        const double t = std::min(-p + s, -p - s);
        if (m_polar) {
          Internal2Cartesian(x0 + t * dx, y0 + t * dy, xc, yc);
        } else {
          xc = x0 + t * dx;
          yc = y0 + t * dy;
        }
        zc = z0 + t * (z1 - z0);
      }
      rc = wire.r;
      if (m_polar) rc *= exp(wire.x);
      return true;
    }
  }
  return false;
}

bool ComponentAnalyticField::InTrapRadius(const double qin, const double xin,
                                          const double yin, const double zin,
                                          double& xw, double& yw,
                                          double& rw) {

  double x0 = xin;
  double y0 = yin;
  if (m_polar) {
    Cartesian2Internal(xin, yin, x0, y0);
  }
  // In case of periodicity, move the point into the basic cell.
  int nX = 0, nY = 0, nPhi = 0;
  if (m_perx) {
    nX = int(round(x0 / m_sx));
    x0 -= m_sx * nX;
  }
  if (m_pery && m_tube) {
    Cartesian2Polar(xin, yin, x0, y0);
    nPhi = int(round(DegreeToRad * y0 / m_sy));
    y0 -= RadToDegree * m_sy * nPhi;
    Polar2Cartesian(x0, y0, x0, y0);
  } else if (m_pery) {
    nY = int(round(y0 / m_sy));
    y0 -= m_sy * nY;
  }
  // Move the point to the correct side of the plane.
  std::array<bool, 4> shift = {false, false, false, false};
  if (m_perx && m_ynplan[0] && x0 <= m_coplan[0]) {
    x0 += m_sx;
    shift[0] = true;
  }
  if (m_perx && m_ynplan[1] && x0 >= m_coplan[1]) {
    x0 -= m_sx;
    shift[1] = true;
  }
  if (m_pery && m_ynplan[2] && y0 <= m_coplan[2]) {
    y0 += m_sy;
    shift[2] = true;
  }
  if (m_pery && m_ynplan[3] && y0 >= m_coplan[3]) {
    y0 -= m_sy;
    shift[3] = true;
  }

  for (const auto& wire : m_w) {
    // Skip wires with the wrong charge.
    if (qin * wire.e > 0.) continue;
    const double dxw0 = wire.x - x0;
    const double dyw0 = wire.y - y0;
    const double r2 = dxw0 * dxw0 + dyw0 * dyw0;
    const double rTrap = wire.r * wire.nTrap;
    if (r2 < rTrap * rTrap) {
      xw = wire.x;
      yw = wire.y;
      rw = wire.r;
      if (shift[0]) xw -= m_sx;
      if (shift[1]) xw += m_sx;
      if (shift[2]) yw -= m_sy;
      if (shift[3]) yw += m_sy;
      if (m_pery && m_tube) {
        double rhow, phiw;
        Cartesian2Polar(xw, yw, rhow, phiw);
        phiw += RadToDegree * m_sy * nPhi;
        Polar2Cartesian(rhow, phiw, xw, yw);
      } else if (m_pery) {
        y0 += m_sy * nY;
      }
      if (m_perx) xw += m_sx * nX;
      if (m_polar) {
        Internal2Cartesian(xw, yw, xw, yw);
        rw *= exp(wire.x);
      }
      if (m_debug) {
        std::cout << m_className << "::InTrapRadius: (" << xin << ", "
                  << yin << ", " << zin << ")" << " within trap radius.\n";
      }
      return true;
    }
  }

  return false;
}

bool ComponentAnalyticField::CrossedPlane(
    const double xx0, const double yy0, const double z0, 
    const double xx1, const double yy1, const double z1,
    double& xc, double& yc, double& zc) {

  double x0 = xx0;
  double y0 = yy0;
  double x1 = xx1;
  double y1 = yy1;
  if (m_polar) {
    Cartesian2Internal(xx0, yy0, x0, y0);
    Cartesian2Internal(xx1, yy1, x1, y1);
  }

  double smin = -1.;
  for (size_t i = 0; i < 2; ++i) {
    if (!m_ynplan[i]) continue;
    if (SameSide(x0, x1, m_coplan[i])) continue;
    const double s = (m_coplan[i] - x0) / (x1 - x0);
    if (smin < 0. || s < smin) {
      smin = s;
    }
  }
  for (size_t i = 2; i < 4; ++i) {
    if (!m_ynplan[i]) continue;
    if (SameSide(y0, y1, m_coplan[i])) continue;
    const double s = (m_coplan[i] - y0) / (y1 - y0);
    if (smin < 0. || s < smin) {
      smin = s;
    }
  }
  if (smin < 0.) return false;
  xc = x0 + smin * (x1 - x0);
  yc = y0 + smin * (y1 - y0);
  zc = z0 + smin * (z1 - z0);
  if (m_polar) Internal2Cartesian(xc, yc, xc, yc);
  return true;
} 

void ComponentAnalyticField::AddWire(const double x, const double y,
                                     const double diameter,
                                     const double voltage,
                                     const std::string& label,
                                     const double length, const double tension,
                                     double rho, const int ntrap) {
  // Check if the provided parameters make sense.
  if (diameter <= 0.) {
    std::cerr << m_className << "::AddWire: Unphysical wire diameter.\n";
    return;
  }

  if (tension <= 0.) {
    std::cerr << m_className << "::AddWire: Unphysical wire tension.\n";
    return;
  }

  if (rho <= 0.) {
    std::cerr << m_className << "::AddWire: Unphysical wire density.\n";
    return;
  }

  if (length <= 0.) {
    std::cerr << m_className << "::AddWire: Unphysical wire length.\n";
    return;
  }

  if (ntrap <= 0) {
    std::cerr << m_className << "::AddWire: Nbr. of trap radii must be > 0.\n";
    return;
  }

  if (m_polar && x <= diameter) {
    std::cerr << m_className << "::AddWire: Wire is too close to the origin.\n";
    return;
  }

  // Create a new wire.
  Wire wire;
  if (m_polar) {
    double r = 0., phi = 0.;
    Polar2Internal(x, y, r, phi);
    wire.x = r;
    wire.y = phi;
    wire.r = 0.5 * diameter / x; 
  } else {
    wire.x = x;
    wire.y = y;
    wire.r = 0.5 * diameter;
  }
  wire.v = voltage;
  wire.u = length;
  wire.type = label;
  wire.e = 0.;
  wire.ind = -1;
  wire.nTrap = ntrap;
  wire.tension = tension;
  wire.density = rho;
  // Add the wire to the list
  m_w.push_back(std::move(wire));
  ++m_nWires;

  // Force recalculation of the capacitance and signal matrices.
  m_cellset = false;
  m_sigset = false;
}

void ComponentAnalyticField::AddTube(const double radius, const double voltage,
                                     const int nEdges,
                                     const std::string& label) {
  // Check if the provided parameters make sense.
  if (radius <= 0.0) {
    std::cerr << m_className << "::AddTube: Unphysical tube dimension.\n";
    return;
  }

  if (nEdges < 3 && nEdges != 0) {
    std::cerr << m_className << "::AddTube: Unphysical number of tube edges ("
              << nEdges << ")\n";
    return;
  }

  // If there is already a tube defined, print a warning message.
  if (m_tube) {
    std::cout << m_className << "::AddTube:\n"
              << "    Warning: Existing tube settings will be overwritten.\n";
  }

  // Set the coordinate system.
  m_tube = true;
  m_polar = false;

  // Set the tube parameters.
  m_cotube = radius;
  m_cotube2 = radius * radius;
  m_vttube = voltage;

  m_ntube = nEdges;

  m_planes[4].type = label;
  m_planes[4].ind = -1;

  // Force recalculation of the capacitance and signal matrices.
  m_cellset = false;
  m_sigset = false;
}

void ComponentAnalyticField::AddPlaneX(const double x, const double v,
                                       const std::string& label) {
  if (m_polar) {
    std::cerr << m_className << "::AddPlaneX:\n"
              << "    Not compatible with polar coordinates; ignored.\n";
    return;
  }
  if (m_ynplan[0] && m_ynplan[1]) {
    std::cerr << m_className << "::AddPlaneX:\n"
              << "    Cannot have more than two planes at constant x.\n";
    return;
  }

  if (m_ynplan[0]) {
    m_ynplan[1] = true;
    m_coplan[1] = x;
    m_vtplan[1] = v;
    m_planes[1].type = label;
    m_planes[1].ind = -1;
  } else {
    m_ynplan[0] = true;
    m_coplan[0] = x;
    m_vtplan[0] = v;
    m_planes[0].type = label;
    m_planes[0].ind = -1;
  }

  // Force recalculation of the capacitance and signal matrices.
  m_cellset = false;
  m_sigset = false;
}

void ComponentAnalyticField::AddPlaneY(const double y, const double v,
                                       const std::string& label) {
  if (m_polar) {
    std::cerr << m_className << "::AddPlaneY:\n"
              << "    Not compatible with polar coordinates; ignored.\n";
    return;
  }
  if (m_ynplan[2] && m_ynplan[3]) {
    std::cerr << m_className << "::AddPlaneY:\n"
              << "    Cannot have more than two planes at constant y.\n";
    return;
  }

  if (m_ynplan[2]) {
    m_ynplan[3] = true;
    m_coplan[3] = y;
    m_vtplan[3] = v;
    m_planes[3].type = label;
    m_planes[3].ind = -1;
  } else {
    m_ynplan[2] = true;
    m_coplan[2] = y;
    m_vtplan[2] = v;
    m_planes[2].type = label;
    m_planes[2].ind = -1;
  }

  // Force recalculation of the capacitance and signal matrices.
  m_cellset = false;
  m_sigset = false;
}

void ComponentAnalyticField::AddPlaneR(const double r, const double v,
                                       const std::string& label) {
  if (!m_polar) {
    std::cerr << m_className << "::AddPlaneR:\n"
              << "    Not compatible with Cartesian coordinates; ignored.\n";
    return;
  }
  if (r <= 0.) {
    std::cerr << m_className << "::AddPlaneR:\n"
              << "    Radius must be larger than zero; plane ignored.\n";
    return;
  }

  if (m_ynplan[0] && m_ynplan[1]) {
    std::cerr << m_className << "::AddPlaneR:\n"
              << "    Cannot have more than two circular planes.\n";
    return;
  }

  if (m_ynplan[0]) {
    m_ynplan[1] = true;
    m_coplan[1] = log(r);
    m_vtplan[1] = v;
    m_planes[1].type = label;
    m_planes[1].ind = -1;
  } else {
    m_ynplan[0] = true;
    m_coplan[0] = log(r);
    m_vtplan[0] = v;
    m_planes[0].type = label;
    m_planes[0].ind = -1;
  }

  // Force recalculation of the capacitance and signal matrices.
  m_cellset = false;
  m_sigset = false;
}

void ComponentAnalyticField::AddPlanePhi(const double phi, const double v,
                                         const std::string& label) {
  if (!m_polar) {
    std::cerr << m_className << "::AddPlanePhi:\n"
              << "    Not compatible with Cartesian coordinates; ignored.\n";
    return;
  }
  if (m_ynplan[2] && m_ynplan[3]) {
    std::cerr << m_className << "::AddPlanePhi:\n"
              << "    Cannot have more than two planes at constant phi.\n";
    return;
  }

  if (m_ynplan[2]) {
    m_ynplan[3] = true;
    m_coplan[3] = phi * DegreeToRad;
    m_vtplan[3] = v;
    m_planes[3].type = label;
    m_planes[3].ind = -1;
  } else {
    m_ynplan[2] = true;
    m_coplan[2] = phi * DegreeToRad;
    m_vtplan[2] = v;
    m_planes[2].type = label;
    m_planes[2].ind = -1;
    // Switch off default periodicity.
    if (m_pery && std::abs(m_sy - TwoPi) < 1.e-4) {
      m_pery = false;
    }
  }

  // Force recalculation of the capacitance and signal matrices.
  m_cellset = false;
  m_sigset = false;
}

void ComponentAnalyticField::AddStripOnPlaneX(const char dir, const double x,
                                              const double smin,
                                              const double smax,
                                              const std::string& label,
                                              const double gap) {
  if (m_polar || (!m_ynplan[0] && !m_ynplan[1])) {
    std::cerr << m_className << "::AddStripOnPlaneX:\n"
              << "    There are no planes at constant x.\n";
    return;
  }

  if (dir != 'y' && dir != 'Y' && dir != 'z' && dir != 'Z') {
    std::cerr << m_className << "::AddStripOnPlaneX:\n"
              << "    Invalid direction (" << dir << ").\n"
              << "    Only strips in y or z direction are possible.\n";
    return;
  }

  if (fabs(smax - smin) < Small) {
    std::cerr << m_className << "::AddStripOnPlaneX:\n"
              << "    Strip width must be greater than zero.\n";
    return;
  }

  Strip newStrip;
  newStrip.type = label;
  newStrip.ind = -1;
  newStrip.smin = std::min(smin, smax);
  newStrip.smax = std::max(smin, smax);
  newStrip.gap = gap > Small ? gap : -1.;

  int iplane = 0;
  if (m_ynplan[1]) {
    const double d0 = fabs(m_coplan[0] - x);
    const double d1 = fabs(m_coplan[1] - x);
    if (d1 < d0) iplane = 1;
  }

  if (dir == 'y' || dir == 'Y') {
    m_planes[iplane].strips1.push_back(std::move(newStrip));
  } else {
    m_planes[iplane].strips2.push_back(std::move(newStrip));
  }
}

void ComponentAnalyticField::AddStripOnPlaneY(const char dir, const double y,
                                              const double smin,
                                              const double smax,
                                              const std::string& label,
                                              const double gap) {
  if (m_polar || (!m_ynplan[2] && !m_ynplan[3])) {
    std::cerr << m_className << "::AddStripOnPlaneY:\n"
              << "    There are no planes at constant y.\n";
    return;
  }

  if (dir != 'x' && dir != 'X' && dir != 'z' && dir != 'Z') {
    std::cerr << m_className << "::AddStripOnPlaneY:\n"
              << "    Invalid direction (" << dir << ").\n"
              << "    Only strips in x or z direction are possible.\n";
    return;
  }

  if (fabs(smax - smin) < Small) {
    std::cerr << m_className << "::AddStripOnPlaneY:\n"
              << "    Strip width must be greater than zero.\n";
    return;
  }

  Strip newStrip;
  newStrip.type = label;
  newStrip.ind = -1;
  newStrip.smin = std::min(smin, smax);
  newStrip.smax = std::max(smin, smax);
  newStrip.gap = gap > Small ? gap : -1.;

  int iplane = 2;
  if (m_ynplan[3]) {
    const double d2 = fabs(m_coplan[2] - y);
    const double d3 = fabs(m_coplan[3] - y);
    if (d3 < d2) iplane = 3;
  }

  if (dir == 'x' || dir == 'X') {
    m_planes[iplane].strips1.push_back(std::move(newStrip));
  } else {
    m_planes[iplane].strips2.push_back(std::move(newStrip));
  }
}

void ComponentAnalyticField::AddStripOnPlaneR(const char dir, const double r,
                                              const double smin,
                                              const double smax,
                                              const std::string& label,
                                              const double gap) {
  if (!m_polar || (!m_ynplan[0] && !m_ynplan[1])) {
    std::cerr << m_className << "::AddStripOnPlaneR:\n"
              << "    There are no planes at constant r.\n";
    return;
  }

  if (dir != 'p' && dir != 'P' && dir != 'z' && dir != 'Z') {
    std::cerr << m_className << "::AddStripOnPlaneR:\n"
              << "    Invalid direction (" << dir << ").\n"
              << "    Only strips in p(hi) or z direction are possible.\n";
    return;
  }

  if (fabs(smax - smin) < Small) {
    std::cerr << m_className << "::AddStripOnPlaneR:\n"
              << "    Strip width must be greater than zero.\n";
    return;
  }

  Strip newStrip;
  newStrip.type = label;
  newStrip.ind = -1;
  if (dir == 'z' || dir == 'Z') {
    const double phimin = smin * DegreeToRad;
    const double phimax = smax * DegreeToRad;
    newStrip.smin = std::min(phimin, phimax);
    newStrip.smax = std::max(phimin, phimax);
  } else {
    newStrip.smin = std::min(smin, smax);
    newStrip.smax = std::max(smin, smax);
  }
  newStrip.gap = gap > Small ? gap : -1.;

  int iplane = 0;
  if (m_ynplan[1]) {
    const double rho = r > 0. ? log(r) : -25.;
    const double d0 = fabs(m_coplan[0] - rho);
    const double d1 = fabs(m_coplan[1] - rho);
    if (d1 < d0) iplane = 1;
  }

  if (dir == 'p' || dir == 'P') {
    m_planes[iplane].strips1.push_back(std::move(newStrip));
  } else {
    m_planes[iplane].strips2.push_back(std::move(newStrip));
  }
}

void ComponentAnalyticField::AddStripOnPlanePhi(const char dir, 
                                                const double phi,
                                                const double smin,
                                                const double smax,
                                                const std::string& label,
                                                const double gap) {
  if (!m_polar || (!m_ynplan[2] && !m_ynplan[3])) {
    std::cerr << m_className << "::AddStripOnPlanePhi:\n"
              << "    There are no planes at constant phi.\n";
    return;
  }

  if (dir != 'r' && dir != 'R' && dir != 'z' && dir != 'Z') {
    std::cerr << m_className << "::AddStripOnPlanePhi:\n"
              << "    Invalid direction (" << dir << ").\n"
              << "    Only strips in r or z direction are possible.\n";
    return;
  }

  if (fabs(smax - smin) < Small) {
    std::cerr << m_className << "::AddStripOnPlanePhi:\n"
              << "    Strip width must be greater than zero.\n";
    return;
  }

  Strip newStrip;
  newStrip.type = label;
  newStrip.ind = -1;
  if (dir== 'z' || dir == 'Z') {
    if (smin < Small || smax < Small) {
      std::cerr << m_className << "::AddStripOnPlanePhi:\n"
                << "    Radius must be greater than zero.\n";
      return;
    }
    const double rhomin = log(smin);
    const double rhomax = log(smax);
    newStrip.smin = std::min(rhomin, rhomax);
    newStrip.smax = std::max(rhomin, rhomax);
  } else {
    newStrip.smin = std::min(smin, smax);
    newStrip.smax = std::max(smin, smax);
  }
  newStrip.gap = gap > Small ? DegreeToRad * gap : -1.;

  int iplane = 2;
  if (m_ynplan[3]) {
    const double d2 = fabs(m_coplan[2] - phi * DegreeToRad);
    const double d3 = fabs(m_coplan[3] - phi * DegreeToRad);
    if (d3 < d2) iplane = 3;
  }

  if (dir == 'r' || dir == 'R') {
    m_planes[iplane].strips1.push_back(std::move(newStrip));
  } else {
    m_planes[iplane].strips2.push_back(std::move(newStrip));
  }
}


void ComponentAnalyticField::AddPixelOnPlaneX(
    const double x, const double ymin, const double ymax, const double zmin,
    const double zmax, const std::string& label, const double gap,
    const double rot) {
  if (m_polar || (!m_ynplan[0] && !m_ynplan[1])) {
    std::cerr << m_className << "::AddPixelOnPlaneX:\n"
              << "    There are no planes at constant x.\n";
    return;
  }

  if (fabs(ymax - ymin) < Small || fabs(zmax - zmin) < Small) {
    std::cerr << m_className << "::AddPixelOnPlaneX:\n"
              << "    Pixel width must be greater than zero.\n";
    return;
  }

  Pixel pixel;
  pixel.type = label;
  pixel.ind = -1;
  pixel.smin = std::min(ymin, ymax);
  pixel.smax = std::max(ymin, ymax);
  pixel.zmin = std::min(zmin, zmax);
  pixel.zmax = std::max(zmin, zmax);
  pixel.gap = gap > Small ? gap : -1.;
  if (fabs(rot) > 1.e-9) {
    pixel.cphi = cos(rot);
    pixel.sphi = sin(rot);
  }

  int iplane = 0;
  if (m_ynplan[1]) {
    const double d0 = fabs(m_coplan[0] - x);
    const double d1 = fabs(m_coplan[1] - x);
    if (d1 < d0) iplane = 1;
  }

  m_planes[iplane].pixels.push_back(std::move(pixel));
}

void ComponentAnalyticField::AddPixelOnPlaneY(
    const double y, const double xmin, const double xmax, const double zmin,
    const double zmax, const std::string& label, const double gap,
    const double rot) {
  if (m_polar || (!m_ynplan[2] && !m_ynplan[3])) {
    std::cerr << m_className << "::AddPixelOnPlaneY:\n"
              << "    There are no planes at constant y.\n";
    return;
  }

  if (fabs(xmax - xmin) < Small || fabs(zmax - zmin) < Small) {
    std::cerr << m_className << "::AddPixelOnPlaneY:\n"
              << "    Pixel width must be greater than zero.\n";
    return;
  }

  Pixel pixel;
  pixel.type = label;
  pixel.ind = -1;
  pixel.smin = std::min(xmin, xmax);
  pixel.smax = std::max(xmin, xmax);
  pixel.zmin = std::min(zmin, zmax);
  pixel.zmax = std::max(zmin, zmax);
  pixel.gap = gap > Small ? gap : -1.;
  if (fabs(rot) > 1.e-9) {
    pixel.cphi = cos(rot);
    pixel.sphi = sin(rot);
  }

  int iplane = 2;
  if (m_ynplan[3]) {
    const double d0 = fabs(m_coplan[2] - y);
    const double d1 = fabs(m_coplan[3] - y);
    if (d1 < d0) iplane = 3;
  }

  m_planes[iplane].pixels.push_back(std::move(pixel));
}

void ComponentAnalyticField::AddPixelOnPlaneR(
    const double r, const double phimin, const double phimax, 
    const double zmin, const double zmax, const std::string& label, 
    const double gap) {
  if (!m_polar || (!m_ynplan[0] && !m_ynplan[1])) {
    std::cerr << m_className << "::AddPixelOnPlaneR:\n"
              << "    There are no planes at constant r.\n";
    return;
  }

  if (fabs(phimax - phimin) < Small || fabs(zmax - zmin) < Small) {
    std::cerr << m_className << "::AddPixelOnPlaneR:\n"
              << "    Pixel width must be greater than zero.\n";
    return;
  }

  Pixel pixel;
  pixel.type = label;
  pixel.ind = -1;
  const double smin = phimin * DegreeToRad;
  const double smax = phimax * DegreeToRad;
  pixel.smin = std::min(smin, smax);
  pixel.smax = std::max(smin, smax);
  pixel.zmin = std::min(zmin, zmax);
  pixel.zmax = std::max(zmin, zmax);
  pixel.gap = gap > Small ? gap : -1.;

  int iplane = 0;
  if (m_ynplan[1]) {
    const double rho = r > 0. ? log(r) : -25.;
    const double d0 = fabs(m_coplan[0] - rho);
    const double d1 = fabs(m_coplan[1] - rho);
    if (d1 < d0) iplane = 1;
  }

  m_planes[iplane].pixels.push_back(std::move(pixel));
}

void ComponentAnalyticField::AddPixelOnPlanePhi(
    const double phi, const double rmin, const double rmax, 
    const double zmin, const double zmax, const std::string& label, 
    const double gap) {
  if (!m_polar || (!m_ynplan[2] && !m_ynplan[3])) {
    std::cerr << m_className << "::AddPixelOnPlanePhi:\n"
              << "    There are no planes at constant phi.\n";
    return;
  }

  if (fabs(rmax - rmin) < Small || fabs(zmax - zmin) < Small) {
    std::cerr << m_className << "::AddPixelOnPlanePhi:\n"
              << "    Pixel width must be greater than zero.\n";
    return;
  }
  if (rmin < Small || rmax < Small) {
    std::cerr << m_className << "::AddPixelOnPlanePhi:\n"
              << "    Radius must be greater than zero.\n";
    return;
  }
  Pixel pixel;
  pixel.type = label;
  pixel.ind = -1;
  const double smin = log(rmin);
  const double smax = log(rmax);
  pixel.smin = std::min(smin, smax);
  pixel.smax = std::max(smin, smax);
  pixel.zmin = std::min(zmin, zmax);
  pixel.zmax = std::max(zmin, zmax);
  pixel.gap = gap > Small ? DegreeToRad * gap : -1.;

  int iplane = 2;
  if (m_ynplan[3]) {
    const double d0 = fabs(m_coplan[2] - phi * DegreeToRad);
    const double d1 = fabs(m_coplan[3] - phi * DegreeToRad);
    if (d1 < d0) iplane = 3;
  }

  m_planes[iplane].pixels.push_back(std::move(pixel));
}

void ComponentAnalyticField::EnableDipoleTerms(const bool on) {

  m_cellset = false;
  m_sigset = false;
  m_dipole = on;
}

void ComponentAnalyticField::SetPeriodicityX(const double s) {
  if (m_polar) {
    std::cerr << m_className << "::SetPeriodicityX:\n"
              << "    Cannot use x-periodicity with polar coordinates.\n";
    return;
  }
  if (s < Small) {
    std::cerr << m_className << "::SetPeriodicityX:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }

  m_periodic[0] = true;
  m_sx = s;
  UpdatePeriodicity();
}

void ComponentAnalyticField::SetPeriodicityY(const double s) {
  if (m_polar) {
    std::cerr << m_className << "::SetPeriodicityY:\n"
              << "    Cannot use y-periodicity with polar coordinates.\n";
    return;
  }
  if (s < Small) {
    std::cerr << m_className << "::SetPeriodicityY:\n"
              << "    Periodic length must be greater than zero.\n";
    return;
  }

  m_periodic[1] = true;
  m_sy = s;
  UpdatePeriodicity();
}

void ComponentAnalyticField::SetPeriodicityPhi(const double s) {
  if (!m_polar && !m_tube) {
    std::cerr << m_className << "::SetPeriodicityPhi:\n"
              << "    Cannot use phi-periodicity with Cartesian coordinates.\n";
    return;
  }
  if (std::abs(360. - s * int(360. / s)) > 1.e-4) {
    std::cerr << m_className << "::SetPeriodicityPhi:\n"
              << "    Phi periods must divide 360; ignored.\n";
    return;
  }

  m_periodic[1] = true;
  m_sy = DegreeToRad * s;
  m_mtube = int(360. / s);
  UpdatePeriodicity();
}

bool ComponentAnalyticField::GetPeriodicityX(double& s) {
  if (!m_periodic[0] || m_polar) {
    s = 0.;
    return false;
  }

  s = m_sx;
  return true;
}

bool ComponentAnalyticField::GetPeriodicityY(double& s) {
  if (!m_periodic[1] || m_polar) {
    s = 0.;
    return false;
  }

  s = m_sy;
  return true;
}

bool ComponentAnalyticField::GetPeriodicityPhi(double& s) {
  if (!m_periodic[1] || (!m_polar && !m_tube)) {
    s = 0.;
    return false;
  }

  s = RadToDegree * m_sy;
  return true;
}

void ComponentAnalyticField::UpdatePeriodicity() {
  // Check if the settings have actually changed.
  if (m_perx && !m_periodic[0]) {
    m_perx = false;
    m_cellset = false;
    m_sigset = false;
  } else if (!m_perx && m_periodic[0]) {
    if (m_sx < Small) {
      std::cerr << m_className << "::UpdatePeriodicity:\n";
      std::cerr << "    Periodicity in x direction was enabled"
                << " but periodic length is not set.\n";
    } else {
      m_perx = true;
      m_cellset = false;
      m_sigset = false;
    }
  }

  if (m_pery && !m_periodic[1]) {
    m_pery = false;
    m_cellset = false;
    m_sigset = false;
  } else if (!m_pery && m_periodic[1]) {
    if (m_sy < Small) {
      std::cerr << m_className << "::UpdatePeriodicity:\n";
      std::cerr << "    Periodicity in y direction was enabled"
                << " but periodic length is not set.\n";
    } else {
      m_pery = true;
      m_cellset = false;
      m_sigset = false;
    }
  }

  // Check if symmetries other than x/y periodicity have been requested
  if (m_periodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicity in z is not possible.\n";
  }

  if (m_mirrorPeriodic[0] || m_mirrorPeriodic[1] || m_mirrorPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Mirror periodicity is not possible.\n";
  }

  if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Axial periodicity is not possible.\n";
  }

  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Rotation symmetry is not possible.\n";
  }
}

void ComponentAnalyticField::SetCartesianCoordinates() {
  if (m_polar) {
    std::cout << m_className << "::SetCartesianCoordinates:\n    "
              << "Switching to Cartesian coordinates; resetting the cell.\n";
    CellInit();
  }
  m_polar = false;
}

void ComponentAnalyticField::SetPolarCoordinates() {
  if (!m_polar) {
    std::cout << m_className << "::SetPolarCoordinates:\n    "
              << "Switching to polar coordinates; resetting the cell.\n";
    CellInit();
  }
  m_polar = true;
  // Set default phi period.
  m_pery = true;
  m_sy = TwoPi;
}

void ComponentAnalyticField::AddCharge(const double x, const double y,
                                       const double z, const double q) {
  // Convert from fC to internal units (division by 4 pi epsilon0).
  Charge3d charge;
  charge.x = x;
  charge.y = y;
  charge.z = z;
  charge.e = q / FourPiEpsilon0;
  m_ch3d.push_back(std::move(charge));
}

void ComponentAnalyticField::ClearCharges() {
  m_ch3d.clear();
  m_nTermBessel = 10;
  m_nTermPoly = 100;
}

void ComponentAnalyticField::PrintCharges() const {
  std::cout << m_className << "::PrintCharges:\n";
  if (m_ch3d.empty()) {
    std::cout << "    No charges present.\n";
    return;
  }
  std::cout << "      x [cm]      y [cm]      z [cm]      charge [fC]\n";
  for (const auto& charge : m_ch3d) {
    std::cout << "     " << std::setw(9) << charge.x << "   " << std::setw(9)
              << charge.y << "   " << std::setw(9) << charge.z << "   "
              << std::setw(11) << charge.e * FourPiEpsilon0 << "\n";
  }
}

unsigned int ComponentAnalyticField::GetNumberOfPlanesX() const {
  if (m_polar) {
    return 0;
  } else if (m_ynplan[0] && m_ynplan[1]) {
    return 2;
  } else if (m_ynplan[0] || m_ynplan[1]) {
    return 1;
  }
  return 0;
}

unsigned int ComponentAnalyticField::GetNumberOfPlanesY() const {
  if (m_polar) {
    return 0;
  } else if (m_ynplan[2] && m_ynplan[3]) {
    return 2;
  } else if (m_ynplan[2] || m_ynplan[3]) {
    return 1;
  }
  return 0;
}

unsigned int ComponentAnalyticField::GetNumberOfPlanesR() const {
  if (!m_polar) {
    return 0;
  } else if (m_ynplan[0] && m_ynplan[1]) {
    return 2;
  } else if (m_ynplan[0] || m_ynplan[1]) {
    return 1;
  }
  return 0;
}

unsigned int ComponentAnalyticField::GetNumberOfPlanesPhi() const {
  if (!m_polar) {
    return 0;
  } else if (m_ynplan[2] && m_ynplan[3]) {
    return 2;
  } else if (m_ynplan[2] || m_ynplan[3]) {
    return 1;
  }
  return 0;
}

bool ComponentAnalyticField::GetWire(const unsigned int i, double& x, double& y,
                                     double& diameter, double& voltage,
                                     std::string& label, double& length,
                                     double& charge, int& ntrap) const {
  if (i >= m_nWires) {
    std::cerr << m_className << "::GetWire: Index out of range.\n";
    return false;
  }

  if (m_polar) {
    double r = 0., theta = 0.;
    Internal2Polar(m_w[i].x, m_w[i].y, r, theta);
    x = r;
    y = theta;
    diameter = 2 * m_w[i].r * r;
  } else {
    x = m_w[i].x;
    y = m_w[i].y;
    diameter = 2 * m_w[i].r;
  }
  voltage = m_w[i].v;
  label = m_w[i].type;
  length = m_w[i].u;
  charge = m_w[i].e;
  ntrap = m_w[i].nTrap;
  return true;
}

bool ComponentAnalyticField::GetPlaneX(const unsigned int i, double& x,
                                       double& voltage,
                                       std::string& label) const {
  if (m_polar || i >= 2 || (i == 1 && !m_ynplan[1])) {
    std::cerr << m_className << "::GetPlaneX: Index out of range.\n";
    return false;
  }

  x = m_coplan[i];
  voltage = m_vtplan[i];
  label = m_planes[i].type;
  return true;
}

bool ComponentAnalyticField::GetPlaneY(const unsigned int i, double& y,
                                       double& voltage,
                                       std::string& label) const {
  if (m_polar || i >= 2 || (i == 1 && !m_ynplan[3])) {
    std::cerr << m_className << "::GetPlaneY: Index out of range.\n";
    return false;
  }

  y = m_coplan[i + 2];
  voltage = m_vtplan[i + 2];
  label = m_planes[i + 2].type;
  return true;
}

bool ComponentAnalyticField::GetPlaneR(const unsigned int i, double& r,
                                       double& voltage,
                                       std::string& label) const {
  if (!m_polar || i >= 2 || (i == 1 && !m_ynplan[1])) {
    std::cerr << m_className << "::GetPlaneR: Index out of range.\n";
    return false;
  }

  r = exp(m_coplan[i]);
  voltage = m_vtplan[i];
  label = m_planes[i].type;
  return true;
}

bool ComponentAnalyticField::GetPlanePhi(const unsigned int i, double& phi,
                                         double& voltage,
                                         std::string& label) const {
  if (!m_polar || i >= 2 || (i == 1 && !m_ynplan[3])) {
    std::cerr << m_className << "::GetPlanePhi: Index out of range.\n";
    return false;
  }

  phi = RadToDegree * m_coplan[i + 2];
  voltage = m_vtplan[i + 2];
  label = m_planes[i + 2].type;
  return true;
}

bool ComponentAnalyticField::GetTube(double& r, double& voltage, int& nEdges,
                                     std::string& label) const {
  if (!m_tube) return false;
  r = m_cotube;
  voltage = m_vttube;
  nEdges = m_ntube;
  label = m_planes[4].type;
  return true;
}

bool ComponentAnalyticField::ElectricFieldAtWire(const unsigned int iw,
                                                 double& ex, double& ey) {
  //-----------------------------------------------------------------------
  //   FFIELD - Subroutine calculating the electric field at a given wire
  //            position, as if the wire itself were not there but with
  //            the presence of its mirror images.
  //   VARIABLES : IW         : wire number
  //               EX, EY     : x- and y-component of the electric field.
  //   (Last changed on 27/ 1/96.)
  //-----------------------------------------------------------------------
  ex = ey = 0.;
  // Check the wire number.
  if (iw >= m_nWires) {
    std::cerr << m_className << "::ElectricFieldAtWire: Index out of range.\n";
    return false;
  }
  // Set the flags appropriately.
  std::vector<bool> cnalso(m_nWires, true);
  cnalso[iw] = false;
 
  const double xpos = m_w[iw].x;
  const double ypos = m_w[iw].y;
  // Call the appropriate function.
  switch (m_cellType) {
    case A00:
      FieldAtWireA00(xpos, ypos, ex, ey, cnalso);
      break;
    case B1X:
      FieldAtWireB1X(xpos, ypos, ex, ey, cnalso);
      break;
    case B1Y:
      FieldAtWireB1Y(xpos, ypos, ex, ey, cnalso);
      break;
    case B2X:
      FieldAtWireB2X(xpos, ypos, ex, ey, cnalso);
      break;
    case B2Y:
      FieldAtWireB2Y(xpos, ypos, ex, ey, cnalso);
      break;
    case C10:
      FieldAtWireC10(xpos, ypos, ex, ey, cnalso);
      break;
    case C2X:
      FieldAtWireC2X(xpos, ypos, ex, ey, cnalso);
      break;
    case C2Y:
      FieldAtWireC2Y(xpos, ypos, ex, ey, cnalso);
      break;
    case C30:
      FieldAtWireC30(xpos, ypos, ex, ey, cnalso);
      break;
    case D10:
      FieldAtWireD10(xpos, ypos, ex, ey, cnalso);
      break;
    case D20:
      FieldAtWireD20(xpos, ypos, ex, ey, cnalso);
      break;
    case D30:
      FieldAtWireD30(xpos, ypos, ex, ey, cnalso);
      break;
    default:
      std::cerr << m_className << "::ElectricFieldAtWire:\n"
                << "    Unknown cell type (id " << m_cellType << ")\n";
      return false;
  }
  // Correct for the equipotential planes.
  ex -= m_corvta;
  ey -= m_corvtb;
  if (m_polar) {
    const double r = exp(xpos);
    const double er = ex / r;
    const double ep = ey / r;
    const double ct = cos(ypos);
    const double st = sin(ypos);
    ex = +ct * er - st * ep;
    ey = +st * er + ct * ep;
  }
  return true;
}

void ComponentAnalyticField::SetScanningGrid(const unsigned int nX,
                                             const unsigned int nY) {
  if (nX < 2) {
    std::cerr << m_className << "::SetScanningGrid:\n"
              << "    Number of x-lines must be > 1.\n";
  } else {
    m_nScanX = nX;
  }
  if (nY < 2) {
    std::cerr << m_className << "::SetScanningGrid:\n"
              << "    Number of y-lines must be > 1.\n";
  } else {
    m_nScanY = nY;
  }
}

void ComponentAnalyticField::SetScanningArea(const double xmin,
                                             const double xmax,
                                             const double ymin,
                                             const double ymax) {
  if (std::abs(xmax - xmin) < Small || std::abs(ymax - ymin) < Small) {
    std::cerr << m_className << "::SetScanningArea:\n"
              << "    Zero range not permitted.\n";
    return;
  }
  m_scanRange = ScanningRange::User;
  m_xScanMin = std::min(xmin, xmax);
  m_xScanMax = std::max(xmin, xmax);
  m_yScanMin = std::min(ymin, ymax);
  m_yScanMax = std::max(ymin, ymax);
}

void ComponentAnalyticField::SetScanningAreaFirstOrder(const double scale) {
  m_scanRange = ScanningRange::FirstOrder;
  if (scale > 0.) {
    m_scaleRange = scale;
  } else {
    std::cerr << m_className << "::SetScanningAreaFirstOrder:\n"
              << "    Scaling factor must be > 0.\n";
  }
}

void ComponentAnalyticField::SetGravity(const double dx, const double dy, 
                                        const double dz) {

  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d > 0.) {
    m_down[0] = dx / d;
    m_down[1] = dy / d;
    m_down[2] = dz / d;
  } else {
    std::cerr << m_className << "::SetGravity:\n"
              << "    The gravity vector has zero norm ; ignored.\n";
  } 
} 

void ComponentAnalyticField::GetGravity(double& dx, double& dy, 
                                        double& dz) const {

  dx = m_down[0];
  dy = m_down[1];
  dz = m_down[2];
}

bool ComponentAnalyticField::ForcesOnWire(
    const unsigned int iw, std::vector<double>& xMap, std::vector<double>& yMap,
    std::vector<std::vector<double> >& fxMap,
    std::vector<std::vector<double> >& fyMap) {
  if (!m_cellset && !Prepare()) {
    std::cerr << m_className << "::ForcesOnWire: Cell not set up.\n";
    return false;
  }

  if (iw >= m_nWires) {
    std::cerr << m_className << "::ForcesOnWire: Wire index out of range.\n";
    return false;
  }
  if (m_polar) {
    std::cerr << m_className << "::ForcesOnWire: Cannot handle polar cells.\n";
    return false;
  }
  const auto& wire = m_w[iw];
  // Compute a 'safe box' around the wire.
  double bxmin = m_perx ? wire.x - 0.5 * m_sx : 2 * m_xmin - m_xmax;
  double bxmax = m_perx ? wire.x + 0.5 * m_sx : 2 * m_xmax - m_xmin;
  double bymin = m_pery ? wire.y - 0.5 * m_sy : 2 * m_ymin - m_ymax;
  double bymax = m_pery ? wire.y + 0.5 * m_sy : 2 * m_ymax - m_ymin;

  // If the initial area is almost zero in 1 direction, make it square.
  if (std::abs(bxmax - bxmin) < 0.1 * std::abs(bymax - bymin)) {
    bxmin = wire.x - 0.5 * std::abs(bymax - bymin);
    bxmax = wire.x + 0.5 * std::abs(bymax - bymin);
  } else if (std::abs(bymax - bymin) < 0.1 * std::abs(bxmax - bxmin)) {
    bymin = wire.y - 0.5 * std::abs(bxmax - bxmin);
    bymax = wire.y + 0.5 * std::abs(bxmax - bxmin);
  }
  const double dw = 2 * wire.r;
  // Scan the other wires.
  for (unsigned int j = 0; j < m_nWires; ++j) {
    if (j == iw) continue;
    const double xj = m_w[j].x;
    const double yj = m_w[j].y;
    const double dj = 2 * m_w[j].r;
    double xnear = m_perx ? xj - m_sx * round((xj - wire.x) / m_sx) : xj;
    double ynear = m_pery ? yj - m_sy * round((yj - wire.y) / m_sy) : yj;
    if (std::abs(xnear - wire.x) > std::abs(ynear - wire.y)) {
      if (xnear < wire.x) {
        bxmin = std::max(bxmin, xnear + dj + dw);
        if (m_perx) bxmax = std::min(bxmax, xnear + m_sx - dj - dw);
      } else {
        bxmax = std::min(bxmax, xnear - dj - dw);
        if (m_perx) bxmin = std::max(bxmin, xnear - m_sx + dj + dw);
      }
    } else {
      if (ynear < wire.y) {
        bymin = std::max({bymin, ynear - dj - dw, ynear + dj + dw});
        if (m_pery) bymax = std::min(bymax, ynear + m_sy - dj - dw);
      } else {
        bymax = std::min({bymax, ynear - dj - dw, ynear + dj + dw});
        if (m_pery) bymin = std::max(bymin, ynear - m_sy + dj + dw);
      }
    }
  }
  // Scan the planes.
  if (m_ynplan[0]) bxmin = std::max(bxmin, m_coplan[0] + dw);
  if (m_ynplan[1]) bxmax = std::min(bxmax, m_coplan[1] - dw);
  if (m_ynplan[2]) bymin = std::max(bymin, m_coplan[2] + dw);
  if (m_ynplan[3]) bymax = std::min(bymax, m_coplan[3] - dw);

  // If there is a tube, check all corners.
  if (m_tube) {
    const double d2 = m_cotube2 - dw * dw;
    if (d2 < Small) {
      std::cerr << m_className << "::ForcesOnWire:\n    Diameter of wire " << iw
                << " is too large compared to the tube.\n";
      return false;
    }

    double corr = sqrt((bxmin * bxmin + bymin * bymin) / d2);
    if (corr > 1.) {
      bxmin /= corr;
      bymin /= corr;
    }
    corr = sqrt((bxmin * bxmin + bymax * bymax) / d2);
    if (corr > 1.) {
      bxmin /= corr;
      bymax /= corr;
    }
    corr = sqrt((bxmax * bxmax + bymin * bymin) / d2);
    if (corr > 1.) {
      bxmax /= corr;
      bymin /= corr;
    }
    corr = sqrt((bxmax * bxmax + bymax * bymax) / d2);
    if (corr > 1) {
      bxmax /= corr;
      bymax /= corr;
    }
  }
  // Make sure we found a reasonable 'safe area'.
  if ((bxmin - wire.x) * (wire.x - bxmax) <= 0 ||
      (bymin - wire.y) * (wire.y - bymax) <= 0) {
    std::cerr << m_className << "::ForcesOnWire:\n    Unable to find an area "
              << "free of elements around wire " << iw << ".\n";
    return false;
  }
  // Now set a reasonable scanning range.
  double sxmin = bxmin;
  double sxmax = bxmax;
  double symin = bymin;
  double symax = bymax;
  if (m_scanRange == ScanningRange::User) {
    // User-specified range.
    sxmin = wire.x + m_xScanMin;
    symin = wire.y + m_yScanMin;
    sxmax = wire.x + m_xScanMax;
    symax = wire.y + m_yScanMax;
  } else if (m_scanRange == ScanningRange::FirstOrder) {
    // Get the field and force at the nominal position.
    double ex = 0., ey = 0.;
    ElectricFieldAtWire(iw, ex, ey);
    double fx = -ex * wire.e * Internal2Newton;
    double fy = -ey * wire.e * Internal2Newton;
    if (m_useGravitationalForce) {
      // Mass per unit length [kg / cm].
      const double m = 1.e-3 * wire.density * Pi * wire.r * wire.r;
      fx -= m_down[0] * GravitationalAcceleration * m;
      fy -= m_down[1] * GravitationalAcceleration * m;
    }
    const double u2 = wire.u * wire.u;
    const double shiftx =
        -125 * fx * u2 / (GravitationalAcceleration * wire.tension);
    const double shifty =
        -125 * fy * u2 / (GravitationalAcceleration * wire.tension);
    // If 0th order estimate of shift is not small.
    const double tol = 0.1 * wire.r;
    if (std::abs(shiftx) > tol || std::abs(shifty) > tol) {
      sxmin = std::max(bxmin, std::min(wire.x + m_scaleRange * shiftx,
                                       wire.x - shiftx / m_scaleRange));
      symin = std::max(bymin, std::min(wire.y + m_scaleRange * shifty,
                                       wire.y - shifty / m_scaleRange));
      sxmax = std::min(bxmax, std::max(wire.x + m_scaleRange * shiftx,
                                       wire.x - shiftx / m_scaleRange));
      symax = std::min(bymax, std::max(wire.y + m_scaleRange * shifty,
                                       wire.y - shifty / m_scaleRange));
      // If one is very small, make the area square within bounds.
      if (std::abs(sxmax - sxmin) < 0.1 * std::abs(symax - symin)) {
        sxmin = std::max(bxmin, wire.x - 0.5 * std::abs(symax - symin));
        sxmax = std::min(bxmax, wire.x + 0.5 * std::abs(symax - symin));
      } else if (std::abs(symax - symin) < 0.1 * std::abs(sxmax - sxmin)) {
        symin = std::max(bymin, wire.y - 0.5 * std::abs(sxmax - sxmin));
        symax = std::min(bymax, wire.y + 0.5 * std::abs(sxmax - sxmin));
      }
    }
  }
  if (m_debug) {
    std::cout << m_className << "::ForcesOnWire:\n";
    std::printf("    Free area %12.5e < x < %12.5e\n", bxmin, bxmax);
    std::printf("              %12.5e < y < %12.5e\n", bymin, bymax);
    std::printf("    Scan area %12.5e < x < %12.5e\n", sxmin, sxmax);
    std::printf("              %12.5e < y < %12.5e\n", symin, symax);
  }

  xMap.resize(m_nScanX);
  const double stepx = (sxmax - sxmin) / (m_nScanX - 1);
  for (unsigned int i = 0; i < m_nScanX; ++i) {
    xMap[i] = sxmin + i * stepx;
  }
  yMap.resize(m_nScanY);
  const double stepy = (symax - symin) / (m_nScanY - 1);
  for (unsigned int i = 0; i < m_nScanY; ++i) {
    yMap[i] = symin + i * stepy;
  }
  // Save the original coordinates of the wire.
  const double x0 = wire.x;
  const double y0 = wire.y;
  // Prepare interpolation tables.
  fxMap.assign(m_nScanX, std::vector<double>(m_nScanY, 0.));
  fyMap.assign(m_nScanX, std::vector<double>(m_nScanY, 0.));
  bool ok = true;
  for (unsigned int i = 0; i < m_nScanX; ++i) {
    for (unsigned int j = 0; j < m_nScanY; ++j) {
      // Get the wire position for this shift.
      m_w[iw].x = xMap[i];
      m_w[iw].y = yMap[j];
      // Verify the current situation.
      if (!WireCheck()) {
        std::cerr << m_className << "::ForcesOnWire: Wire " << iw << ".\n"
                  << "    Scan involves a disallowed wire position.\n";
        ok = false;
        continue;
      }
      // Recompute the charges for this configuration.
      if (!Setup()) {
        std::cerr << m_className << "::ForcesOnWire: Wire " << iw << ".\n"
                  << "    Failed to compute charges at a scan point.\n";
        ok = false;
        continue;
      }
      // Compute the forces.
      double ex = 0., ey = 0.;
      ElectricFieldAtWire(iw, ex, ey);
      fxMap[i][j] = -ex * wire.e * Internal2Newton;
      fyMap[i][j] = -ey * wire.e * Internal2Newton;
    }
  }
  // Place the wire back in its shifted position.
  m_w[iw].x = x0;
  m_w[iw].y = y0;
  Setup();
  return ok;
}

void ComponentAnalyticField::SetNumberOfSteps(const unsigned int n) {
  if (n == 0) {
    std::cerr << m_className << "::SetNumberOfSteps:\n"
              << "    Number of steps must be > 0.\n";
    return;
  }
  m_nSteps = n;
}

bool ComponentAnalyticField::WireDisplacement(
    const unsigned int iw, const bool detailed, std::vector<double>& csag,
    std::vector<double>& xsag, std::vector<double>& ysag, double& stretch,
    const bool print) {
  if (!m_cellset && !Prepare()) {
    std::cerr << m_className << "::WireDisplacement: Cell not set up.\n";
    return false;
  }
  if (iw >= m_nWires) {
    std::cerr << m_className
              << "::WireDisplacement: Wire index out of range.\n";
    return false;
  }
  if (m_polar) {
    std::cerr << m_className 
              << "::WireDisplacement: Cannot handle polar cells.\n";
    return false;
  }
  const auto& wire = m_w[iw];
  // Save the original coordinates.
  const double x0 = wire.x;
  const double y0 = wire.y;
  // First-order approximation.
  if (!Setup()) {
    std::cerr << m_className << "::WireDisplacement:\n"
              << "    Charge calculation failed at central position.\n";
    return false;
  }

  double fx = 0., fy = 0.;
  if (m_useElectrostaticForce) {
    double ex = 0., ey = 0.;
    ElectricFieldAtWire(iw, ex, ey);
    fx -= ex * wire.e * Internal2Newton;
    fy -= ey * wire.e * Internal2Newton;
  }
  if (m_useGravitationalForce) {
    // Mass per unit length [kg / cm].
    const double m = 1.e-3 * wire.density * Pi * wire.r * wire.r;
    fx -= m_down[0] * GravitationalAcceleration * m;
    fy -= m_down[1] * GravitationalAcceleration * m;
  }
  const double u2 = wire.u * wire.u;
  const double shiftx =
      -125 * fx * u2 / (GravitationalAcceleration * wire.tension);
  const double shifty =
      -125 * fy * u2 / (GravitationalAcceleration * wire.tension);
  // Get the elongation from this.
  const double s = 4 * sqrt(shiftx * shiftx + shifty * shifty) / wire.u;
  double length = wire.u;
  if (s > Small) {
    const double t = sqrt(1 + s * s);
    length *= 0.5 * (t + log(s + t) / s);
  }
  stretch = (length - wire.u) / wire.u;
  if (print) {
    std::cout << m_className
              << "::WireDisplacement:\n"
              << "    Forces and displacement in 0th order.\n"
              << "    Wire information: number   = " << iw << "\n"
              << "                      type     = " << wire.type << "\n"
              << "                      location = (" << wire.x << ", " << wire.y
              << ") cm\n"
              << "                      voltage  = " << wire.v << " V\n"
              << "                      length   = " << wire.u << " cm\n"
              << "                      tension  = " << wire.tension << " g\n"
              << "    In this position: Fx       = " << fx << " N/cm\n"
              << "                      Fy       = " << fy << " N/cm\n"
              << "                      x-shift  = " << shiftx << " cm\n"
              << "                      y-shift  = " << shifty << " cm\n"
              << "                      stretch  = " << 100. * stretch << "%\n";
  }
  if (!detailed) {
    csag = {0.};
    xsag = {shiftx};
    ysag = {shifty};
    return true;
  }
  std::vector<double> xMap(m_nScanX, 0.);
  std::vector<double> yMap(m_nScanY, 0.);
  std::vector<std::vector<double> > fxMap(m_nScanX,
                                          std::vector<double>(m_nScanY, 0.));
  std::vector<std::vector<double> > fyMap(m_nScanX,
                                          std::vector<double>(m_nScanY, 0.));
  if (!ForcesOnWire(iw, xMap, yMap, fxMap, fyMap)) {
    std::cerr << m_className << "::WireDisplacement:\n"
              << "    Could not compute the electrostatic force table.\n";
    return false;
  }
  // Compute the detailed wire shift.
  if (!SagDetailed(wire, xMap, yMap, fxMap, fyMap, csag, xsag, ysag)) {
    std::cerr << m_className << "::WireDisplacement: Wire " << iw << ".\n"
              << "    Computation of the wire sag failed.\n";
    return false;
  }
  // Verify that the wire is in range.
  const double sxmin = xMap.front();
  const double sxmax = xMap.back();
  const double symin = yMap.front();
  const double symax = yMap.back();
  const unsigned int nSag = xsag.size();
  bool outside = false;
  length = 0.;
  double xAvg = 0.;
  double yAvg = 0.;
  double xMax = 0.;
  double yMax = 0.;
  for (unsigned int i = 0; i < nSag; ++i) {
    if (x0 + xsag[i] < sxmin || x0 + xsag[i] > sxmax || y0 + ysag[i] < symin ||
        y0 + ysag[i] > symax) {
      outside = true;
    }
    xAvg += xsag[i];
    yAvg += ysag[i];
    xMax = std::max(xMax, std::abs(xsag[i]));
    yMax = std::max(yMax, std::abs(ysag[i]));
    if (i == 0) continue;
    const double dx = xsag[i] - xsag[i - 1];
    const double dy = ysag[i] - ysag[i - 1];
    const double dz = csag[i] - csag[i - 1];
    length += sqrt(dx * dx + dy * dy + dz * dz);
  }
  xAvg /= nSag;
  yAvg /= nSag;
  stretch = (length - wire.u) / wire.u;
  // Warn if a point outside the scanning area was found.
  if (outside) {
    std::cerr
        << m_className << "::WireDisplacement: Warning.\n    "
        << "The wire profile is located partially outside the scanning area.\n";
  }
  if (print) {
    std::cout << "    Sag profile for wire " << iw << ".\n"
              << " Point     z [cm]   x-sag [um]   y-sag [um]\n";
    for (unsigned int i = 0; i < nSag; ++i) {
      std::printf(" %3d   %10.4f  %10.4f  %10.4f\n", 
                  i, csag[i], xsag[i] * 1.e4, ysag[i] * 1.e4);
    }
    std::printf("    Average sag in x and y: %10.4f and %10.4f micron\n",
                1.e4 * xAvg, 1.e4 * yAvg);
    std::printf("    Maximum sag in x and y: %10.4f and %10.4f micron\n",
                1.e4 * xMax, 1.e4 * yMax);
    std::cout << "    Elongation:             " << 100. * stretch << "%\n";
  }
  return true;
}

int ComponentAnalyticField::Field(const double xin, const double yin,
                                  const double zin, double& ex, double& ey,
                                  double& ez, double& volt, const bool opt) {
  //-----------------------------------------------------------------------
  //   EFIELD - Subroutine calculating the electric field and the potential
  //            at a given place. It makes use of the routines POT...,
  //            depending on the type of the cell.
  //   VARIABLES : XPOS       : x-coordinate of the place where the field
  //                            is to be calculated.
  //               YPOS, ZPOS : y- and z-coordinates
  //               EX, EY, EZ : x-, y-, z-component of the electric field.
  //               VOLT       : potential at (XPOS,YPOS).
  //               IOPT       : 1 if both E and V are required, 0 if only E
  //                            is to be computed.
  //               ILOC       : Tells where the point is located (0: normal
  //                            I > 0: in wire I, -1: outside a plane,
  //                            -5: in a material, -6: outside the mesh,
  //                            -10: unknown potential).
  //   (Last changed on 28/ 9/07.)
  //-----------------------------------------------------------------------

  // Initialise the field for returns without actual calculations.
  ex = ey = ez = volt = 0.;

  // Make sure the charges have been calculated.
  if (!m_cellset && !Prepare()) return -11;

  double xpos = xin, ypos = yin;
  if (m_polar) Cartesian2Internal(xin, yin, xpos, ypos);
  // In case of periodicity, move the point into the basic cell.
  if (m_perx) xpos -= m_sx * round(xin / m_sx);
  double arot = 0.;
  if (m_pery && m_tube) {
    Cartesian2Polar(xin, yin, xpos, ypos);
    arot = RadToDegree * m_sy * round(DegreeToRad * ypos / m_sy);
    ypos -= arot;
    Polar2Cartesian(xpos, ypos, xpos, ypos);
  } else if (m_pery) {
    ypos -= m_sy * round(ypos / m_sy);
  }

  // Move the point to the correct side of the plane.
  if (m_perx && m_ynplan[0] && xpos <= m_coplan[0]) xpos += m_sx;
  if (m_perx && m_ynplan[1] && xpos >= m_coplan[1]) xpos -= m_sx;
  if (m_pery && m_ynplan[2] && ypos <= m_coplan[2]) ypos += m_sy;
  if (m_pery && m_ynplan[3] && ypos >= m_coplan[3]) ypos -= m_sy;

  // In case (XPOS,YPOS) is located behind a plane there is no field.
  if (m_tube) {
    if (!InTube(xpos, ypos, m_cotube, m_ntube)) {
      volt = m_vttube;
      return -4;
    }
  } else {
    if (m_ynplan[0] && xpos < m_coplan[0]) {
      volt = m_vtplan[0];
      return -4;
    }
    if (m_ynplan[1] && xpos > m_coplan[1]) {
      volt = m_vtplan[1];
      return -4;
    }
    if (m_ynplan[2] && ypos < m_coplan[2]) {
      volt = m_vtplan[2];
      return -4;
    }
    if (m_ynplan[3] && ypos > m_coplan[3]) {
      volt = m_vtplan[3];
      return -4;
    }
  }

  // If (xpos, ypos) is within a wire, there is no field either.
  for (int i = m_nWires; i--;) {
    double dx = xpos - m_w[i].x;
    double dy = ypos - m_w[i].y;
    // Correct for periodicities.
    if (m_perx) dx -= m_sx * round(dx / m_sx);
    if (m_pery) dy -= m_sy * round(dy / m_sy);
    // Check the actual position.
    if (dx * dx + dy * dy < m_w[i].r * m_w[i].r) {
      volt = m_w[i].v;
      return i + 1;
    }
  }

  // Call the appropriate potential calculation function.
  switch (m_cellType) {
    case A00:
      FieldA00(xpos, ypos, ex, ey, volt, opt);
      break;
    case B1X:
      FieldB1X(xpos, ypos, ex, ey, volt, opt);
      break;
    case B1Y:
      FieldB1Y(xpos, ypos, ex, ey, volt, opt);
      break;
    case B2X:
      FieldB2X(xpos, ypos, ex, ey, volt, opt);
      break;
    case B2Y:
      FieldB2Y(xpos, ypos, ex, ey, volt, opt);
      break;
    case C10:
      FieldC10(xpos, ypos, ex, ey, volt, opt);
      break;
    case C2X:
      FieldC2X(xpos, ypos, ex, ey, volt, opt);
      break;
    case C2Y:
      FieldC2Y(xpos, ypos, ex, ey, volt, opt);
      break;
    case C30:
      FieldC30(xpos, ypos, ex, ey, volt, opt);
      break;
    case D10:
      FieldD10(xpos, ypos, ex, ey, volt, opt);
      break;
    case D20:
      FieldD20(xpos, ypos, ex, ey, volt, opt);
      break;
    case D30:
      FieldD30(xpos, ypos, ex, ey, volt, opt);
      break;
    default:
      // Unknown cell type
      std::cerr << m_className << "::Field:\n";
      std::cerr << "    Unknown cell type (id " << m_cellType << ")\n";
      return -10;
  }

  // Add dipole terms if requested
  if (m_dipole) {
    double exd = 0., eyd = 0., voltd = 0.;
    switch (m_cellType) {
      case A00:
        DipoleFieldA00(xpos, ypos, exd, eyd, voltd, opt);
        break;
      case B1X:
        DipoleFieldB1X(xpos, ypos, exd, eyd, voltd, opt);
        break;
      case B1Y:
        DipoleFieldB1Y(xpos, ypos, exd, eyd, voltd, opt);
        break;
      case B2X:
        DipoleFieldB2X(xpos, ypos, exd, eyd, voltd, opt);
        break;
      case B2Y:
        DipoleFieldB2Y(xpos, ypos, exd, eyd, voltd, opt);
        break;
      default:
        break;
    }
    ex += exd;
    ey += eyd;
    volt += voltd;
  }

  // Rotate the field in some special cases.
  if (m_pery && m_tube) {
    double xaux, yaux;
    Cartesian2Polar(ex, ey, xaux, yaux);
    yaux += arot;
    Polar2Cartesian(xaux, yaux, ex, ey);
  }

  // Correct for the equipotential planes.
  ex -= m_corvta;
  ey -= m_corvtb;
  volt += m_corvta * xpos + m_corvtb * ypos + m_corvtc;

  // Add three dimensional point charges.
  if (!m_ch3d.empty()) {
    double ex3d = 0., ey3d = 0., ez3d = 0., volt3d = 0.;
    switch (m_cellType) {
      case A00:
      case B1X:
      case B1Y:
        Field3dA00(xin, yin, zin, ex3d, ey3d, ez3d, volt3d);
        break;
      case B2X:
        Field3dB2X(xin, yin, zin, ex3d, ey3d, ez3d, volt3d);
        break;
      case B2Y:
        Field3dB2Y(xin, yin, zin, ex3d, ey3d, ez3d, volt3d);
        break;
      case D10:
        Field3dD10(xin, yin, zin, ex3d, ey3d, ez3d, volt3d);
        break;
      default:
        Field3dA00(xin, yin, zin, ex3d, ey3d, ez3d, volt3d);
        break;
    }
    ex += ex3d;
    ey += ey3d;
    ez += ez3d;
    volt += volt3d;
  }

  if (m_polar) {
    const double r = exp(xpos);
    const double er = ex / r;
    const double ep = ey / r;
    const double theta = atan2(yin, xin);
    const double ct = cos(theta);
    const double st = sin(theta);
    ex = +ct * er - st * ep;
    ey = +st * er + ct * ep; 
  }
  return 0;
}

void ComponentAnalyticField::CellInit() {

  m_cellset = false;
  m_sigset = false;

  // Coordinate system
  m_polar = false;

  // Cell type
  m_cellType = A00;

  // Bounding box and voltage range.
  m_xmin = m_xmax = 0.;
  m_ymin = m_ymax = 0.;
  m_zmin = m_zmax = 0.;
  m_vmin = m_vmax = 0.;

  // Periodicities
  m_perx = m_pery = false;
  m_periodic[0] = m_periodic[1] = false;
  m_sx = m_sy = 1.;

  // Signals
  m_nFourier = 1;
  m_cellTypeFourier = A00;
  m_fperx = false;
  m_fpery = false;
  m_mxmin = 0;
  m_mxmax = 0;
  m_mymin = 0;
  m_mymax = 0;
  m_mfexp = 0;

  m_readout.clear();

  // Wires.
  m_nWires = 0;
  m_w.clear();

  // Dipole settings
  m_dipole = false;
  m_cosph2.clear();
  m_sinph2.clear();
  m_amp2.clear();

  // B2 type cells
  m_b2sin.clear();
  // C type cells
  m_mode = 0;
  m_zmult = std::complex<double>(0., 0.);
  m_p1 = m_p2 = m_c1 = 0.;
  // D3 type cells
  m_zw.clear();
  m_kappa = 0.;

  // Reference potential
  m_v0 = 0.;
  m_corvta = m_corvtb = m_corvtc = 0.;

  // Planes
  for (int i = 0; i < 5; ++i) {
    m_planes[i].type = '?';
    m_planes[i].ind = -1;
    m_planes[i].ewxcor = 0.;
    m_planes[i].ewycor = 0.;
    m_planes[i].strips1.clear();
    m_planes[i].strips2.clear();
    m_planes[i].pixels.clear();
  }
  for (int i = 0; i < 4; ++i) {
    m_ynplan[i] = false;
    m_coplan[i] = 0.;
    m_vtplan[i] = 0.;
  }
  // Plane shorthand
  m_ynplax = m_ynplay = false;
  m_coplax = m_coplay = 1.;

  // Tube properties
  m_tube = false;
  m_ntube = 0;
  m_mtube = 1;
  m_cotube = 1.;
  m_cotube2 = 1.;
  m_vttube = 0.;

  // Weighting charges
  m_qwire.clear();
  m_qplane.clear();

  // 3D charges
  m_ch3d.clear();

  // Gravity
  m_down = {0, 0, 1};
}

bool ComponentAnalyticField::Prepare() {
  std::lock_guard<std::mutex> guard(m_mutex);
  // Check that the cell makes sense.
  if (!CellCheck()) {
    std::cerr << m_className << "::Prepare:\n"
              << "    The cell does not meet the requirements.\n";
    return false;
  }
  if (m_debug) std::cout << m_className << "::Prepare: Cell check ok.\n";

  // Determine the cell type.
  if (!CellType()) {
    std::cerr << m_className << "::Prepare:\n"
              << "    Type identification of the cell failed.\n";
    return false;
  }
  if (m_debug) {
    std::cout << m_className << "::Prepare:\n"
              << "    Cell is of type " << CellType() << ".\n";
  }

  // Calculate the charges.
  if (!Setup()) {
    std::cerr << m_className << "::Prepare: Calculation of charges failed.\n";
    return false;
  }
  if (m_debug) {
    std::cout << m_className << "::Prepare:\n"
              << "    Calculation of charges was successful.\n";
  }

  // Assign default gaps for strips and pixels.
  if (!PrepareStrips()) {
    std::cerr << m_className << "::Prepare: Strip/pixel preparation failed.\n";
    return false;
  }

  m_cellset = true;

  // Add dipole terms if required
  if (m_dipole) {
    if (!SetupDipoleTerms()) {
      std::cerr << m_className << "::Prepare:\n"
                << "    Computing the dipole moments failed.\n";
      m_dipole = false;
    }
  }
  return true;
}

bool ComponentAnalyticField::CellCheck() {
  //-----------------------------------------------------------------------
  //   CELCHK - Subroutine checking the wire positions, The equipotential
  //            planes and the periodicity. Two planes having different
  //            voltages are not allowed to have a common line, wires are
  //            not allowed to be at the same position etc.
  //            This routine determines also the cell-dimensions.
  //   VARIABLE  : WRONG(I)   : .TRUE. if wire I will be removed
  //               IPLAN.     : Number of wires with coord > than plane .
  //   (Last changed on 16/ 2/05.)
  //-----------------------------------------------------------------------

  // Checks on the planes, first move the x planes to the basic cell.
  if (m_perx) {
    const std::string xr = m_polar ? "r" : "x";
    double conew1 = m_coplan[0] - m_sx * round(m_coplan[0] / m_sx);
    double conew2 = m_coplan[1] - m_sx * round(m_coplan[1] / m_sx);
    // Check that they are not one on top of the other.
    if (m_ynplan[0] && m_ynplan[1] && conew1 == conew2) {
      if (conew1 > 0.)
        conew1 -= m_sx;
      else
        conew2 += m_sx;
    }
    // Print some warnings if the planes have been moved.
    if ((conew1 != m_coplan[0] && m_ynplan[0]) ||
        (conew2 != m_coplan[1] && m_ynplan[1])) {
      std::cout << m_className << "::CellCheck:\n    The planes in "
                << xr << " are moved to the basic period.\n"
                << "    This should not affect the results.\n";
    }
    m_coplan[0] = conew1;
    m_coplan[1] = conew2;

    // Two planes should now be separated by SX, cancel PERX if not.
    if (m_ynplan[0] && m_ynplan[1] && fabs(m_coplan[1] - m_coplan[0]) != m_sx) {
      std::cerr << m_className << "::CellCheck:\n    The separation of the "
                << xr << " planes does not match the period.\n"
                << "    The periodicity is cancelled.\n";
      m_perx = false;
    }
    // If there are two planes left, they should have identical V's.
    if (m_ynplan[0] && m_ynplan[1] && m_vtplan[0] != m_vtplan[1]) {
      std::cerr << m_className << "::CellCheck:\n    The voltages of the two "
                << xr << " planes differ.\n"
                << "    The periodicity is cancelled.\n";
      m_perx = false;
    }
  }

  // Idem for the y or phi planes: move them to the basic period.
  if (m_pery) {
    const std::string yp = m_polar ? "phi" : "y";
    double conew3 = m_coplan[2] - m_sy * round(m_coplan[2] / m_sy);
    double conew4 = m_coplan[3] - m_sy * round(m_coplan[3] / m_sy);
    // Check that they are not one on top of the other.
    if (m_ynplan[2] && m_ynplan[3] && conew3 == conew4) {
      if (conew3 > 0.)
        conew3 -= m_sy;
      else
        conew4 += m_sy;
    }
    // Print some warnings if the planes have been moved.
    if ((conew3 != m_coplan[2] && m_ynplan[2]) ||
        (conew4 != m_coplan[3] && m_ynplan[3])) {
      std::cout << m_className << "::CellCheck:\n    The planes in "
                << yp << " are moved to the basic period.\n"
                << "    This should not affect the results.\n";
    }
    m_coplan[2] = conew3;
    m_coplan[3] = conew4;

    // Two planes should now be separated by SY, cancel PERY if not.
    if (m_ynplan[2] && m_ynplan[3] && fabs(m_coplan[3] - m_coplan[2]) != m_sy) {
      std::cerr << m_className << "::CellCheck:\n    The separation of the two "
                << yp << " planes does not match the period.\n"
                << "    The periodicity is cancelled.\n";
      m_pery = false;
    }
    // If there are two planes left, they should have identical V's.
    if (m_ynplan[2] && m_ynplan[3] && m_vtplan[2] != m_vtplan[3]) {
      std::cerr << m_className << "::CellCheck:\n    The voltages of the two "
                << yp << " planes differ.\n"
                << "    The periodicity is cancelled.\n";
      m_pery = false;
    }
  }

  // Check that there is no voltage conflict of crossing planes.
  for (int i = 0; i < 2; ++i) {
    for (int j = 2; j < 3; ++j) {
      if (m_ynplan[i] && m_ynplan[j] && m_vtplan[i] != m_vtplan[j]) {
        const std::string yp = m_polar ? "phi" : "y";
        std::cerr << m_className << "::CellCheck:\n"
                  << "    Conflicting potential of two crossing planes.\n"
                  << "    One " << yp << " plane is removed.\n";
        m_ynplan[j] = false;
      }
    }
  }

  // Make sure the coordinates of the planes are properly ordered.
  for (int i = 0; i < 3; i += 2) {
    if (m_ynplan[i] && m_ynplan[i + 1]) {
      if (m_coplan[i] == m_coplan[i + 1]) {
        std::cerr << m_className << "::CellCheck:\n"
                  << "    Two planes are on top of each other.\n"
                  << "    One of them is removed.\n";
        m_ynplan[i + 1] = false;
      }
      if (m_coplan[i] > m_coplan[i + 1]) {
        if (m_debug) {
          std::cout << m_className << "::CellCheck:\n    Planes "
                    << i << " and " << i + 1 << " are interchanged.\n";
        }
        // Interchange the two planes.
        const double cohlp = m_coplan[i];
        m_coplan[i] = m_coplan[i + 1];
        m_coplan[i + 1] = cohlp;

        const double vthlp = m_vtplan[i];
        m_vtplan[i] = m_vtplan[i + 1];
        m_vtplan[i + 1] = vthlp;

        Plane plahlp = m_planes[i];
        m_planes[i] = m_planes[i + 1];
        m_planes[i + 1] = plahlp;
      }
    }
  }

  // Checks on the wires, start moving them to the basic x period.
  if (m_perx) {
    for (auto& wire : m_w) {
      double xnew = wire.x - m_sx * round(wire.x / m_sx);
      if (m_ynplan[0] && xnew <= m_coplan[0]) xnew += m_sx;
      if (m_ynplan[1] && xnew >= m_coplan[1]) xnew -= m_sx;
      if (fabs(xnew - wire.x) > 1.e-8) {
        double xprt = wire.x;
        double yprt = wire.y;
        if (m_polar) Internal2Polar(wire.x, wire.y, xprt, yprt);
        const std::string xr = m_polar ? "r" : "x";
        std::cout << m_className << "::CellCheck:\n    The " << wire.type
                  << "-wire at (" << xprt << ", " << yprt
                  << ") is moved to the basic " << xr << " period.\n"
                  << "    This should not affect the results.\n";
      }
      wire.x = xnew;
    }
  }

  // In case of y-periodicity, all wires should be in the first y-period.
  if (m_tube && m_pery) {
    for (unsigned int i = 0; i < m_nWires; ++i) {
      double xnew = m_w[i].x;
      double ynew = m_w[i].y;
      Cartesian2Polar(xnew, ynew, xnew, ynew);
      if (int(round(DegreeToRad * ynew / m_sy)) != 0) {
        std::cout << m_className << "::CellCheck:\n";
        std::cout << "    The " << m_w[i].type << "-wire at (" << m_w[i].x
                  << ", " << m_w[i].y
                  << ") is moved to the basic phi period.\n";
        std::cout << "    This should not affect the results.\n";
        ynew -= RadToDegree * m_sy * round(DegreeToRad * ynew / m_sy);
        Polar2Cartesian(xnew, ynew, m_w[i].x, m_w[i].y);
      }
    }
  } else if (m_pery) {
    for (auto& wire : m_w) {
      double ynew = wire.y - m_sy * round(wire.y / m_sy);
      if (m_ynplan[2] && ynew <= m_coplan[2]) ynew += m_sy;
      if (m_ynplan[3] && ynew >= m_coplan[3]) ynew -= m_sy;
      if (fabs(ynew - wire.y) > 1.e-8) {
        double xprt = wire.x;
        double yprt = wire.y;
        if (m_polar) Internal2Polar(wire.x, wire.y, xprt, yprt);
        const std::string yp = m_polar ? "phi" : "y";
        std::cout << m_className << "::CellCheck:\n    The " << wire.type
                  << "-wire at (" << xprt << ", " << yprt
                  << ") is moved to the basic " << yp << " period.\n"
                  << "    This should not affect the results.\n";
      }
      wire.y = ynew;
    }
  }

  // Make sure the plane numbering is standard: P1 wires P2, P3 wires P4.
  int iplan1 = 0, iplan2 = 0, iplan3 = 0, iplan4 = 0;
  for (const auto& wire : m_w) {
    if (m_ynplan[0] && wire.x <= m_coplan[0]) ++iplan1;
    if (m_ynplan[1] && wire.x <= m_coplan[1]) ++iplan2;
    if (m_ynplan[2] && wire.y <= m_coplan[2]) ++iplan3;
    if (m_ynplan[3] && wire.y <= m_coplan[3]) ++iplan4;
  }

  // Find out whether smaller (-1) or larger (+1) coord. are to be kept.
  const int imid = int(m_nWires) / 2;
  if (m_ynplan[0] && m_ynplan[1]) {
    if (iplan1 > imid) {
      m_ynplan[1] = false;
      iplan1 = -1;
    } else {
      iplan1 = +1;
    }
    if (iplan2 < imid) {
      m_ynplan[0] = false;
      iplan2 = +1;
    } else {
      iplan2 = -1;
    }
  }
  if (m_ynplan[0] && !m_ynplan[1]) {
    if (iplan1 > imid) {
      iplan1 = -1;
    } else {
      iplan1 = +1;
    }
  }
  if (m_ynplan[1] && !m_ynplan[0]) {
    if (iplan2 < imid) {
      iplan2 = +1;
    } else {
      iplan2 = -1;
    }
  }

  if (m_ynplan[2] && m_ynplan[3]) {
    if (iplan3 > imid) {
      m_ynplan[3] = false;
      iplan3 = -1;
    } else {
      iplan3 = +1;
    }
    if (iplan4 < imid) {
      m_ynplan[2] = false;
      iplan4 = +1;
    } else {
      iplan4 = -1;
    }
  }
  if (m_ynplan[2] && !m_ynplan[3]) {
    if (iplan3 > imid) {
      iplan3 = -1;
    } else {
      iplan3 = +1;
    }
  }
  if (m_ynplan[3] && !m_ynplan[2]) {
    if (iplan4 < imid) {
      iplan4 = +1;
    } else {
      iplan4 = -1;
    }
  }

  // Adapt the numbering of the planes if necessary.
  if (iplan1 == -1) {
    m_ynplan[0] = false;
    m_ynplan[1] = true;
    m_coplan[1] = m_coplan[0];
    m_vtplan[1] = m_vtplan[0];
    m_planes[1] = m_planes[0];
  }

  if (iplan2 == +1) {
    m_ynplan[1] = false;
    m_ynplan[0] = true;
    m_coplan[0] = m_coplan[1];
    m_vtplan[0] = m_vtplan[1];
    m_planes[0] = m_planes[1];
  }

  if (iplan3 == -1) {
    m_ynplan[2] = false;
    m_ynplan[3] = true;
    m_coplan[3] = m_coplan[2];
    m_vtplan[3] = m_vtplan[2];
    m_planes[3] = m_planes[2];
  }

  if (iplan4 == +1) {
    m_ynplan[3] = false;
    m_ynplan[2] = true;
    m_coplan[2] = m_coplan[3];
    m_vtplan[2] = m_vtplan[3];
    m_planes[2] = m_planes[3];
  }

  std::vector<bool> wrong(m_nWires, false);
  // Second pass for the wires, check position relative to the planes.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const double rw = m_w[i].r;
    const double dw = 2. * rw;
    if (m_ynplan[0] && m_w[i].x - rw <= m_coplan[0]) wrong[i] = true;
    if (m_ynplan[1] && m_w[i].x + rw >= m_coplan[1]) wrong[i] = true;
    if (m_ynplan[2] && m_w[i].y - rw <= m_coplan[2]) wrong[i] = true;
    if (m_ynplan[3] && m_w[i].y + rw >= m_coplan[3]) wrong[i] = true;
    if (m_tube) {
      if (!InTube(m_w[i].x, m_w[i].y, m_cotube, m_ntube)) {
        std::cerr << m_className << "::CellCheck:\n";
        std::cerr << "    The " << m_w[i].type << "-wire at (" << m_w[i].x
                  << ", " << m_w[i].y << ") is located outside the tube.\n";
        std::cerr << "    This wire is removed.\n";
        wrong[i] = true;
      }
    } else if (wrong[i]) {
      double xprt = m_w[i].x;
      double yprt = m_w[i].y;
      if (m_polar) Internal2Polar(m_w[i].x, m_w[i].y, xprt, yprt);
      std::cerr << m_className << "::CellCheck:\n    The " << m_w[i].type
                << "-wire at (" << xprt << ", " << yprt << ") is located "
                << "outside the planes.\n    This wire is removed.\n";
    } else if ((m_perx && dw >= m_sx) || (m_pery && dw >= m_sy)) {
      double xprt = m_w[i].x;
      double yprt = m_w[i].y;
      if (m_polar) Internal2Polar(m_w[i].x, m_w[i].y, xprt, yprt);
      std::cerr << m_className << "::CellCheck:\n    The diameter of the "
                << m_w[i].type << "-wire at (" << xprt << ", " << yprt
                << ") exceeds 1 period.\n    This wire is removed.\n";
      wrong[i] = true;
    }
  }

  // Check the wire spacing.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    if (wrong[i]) continue;
    for (unsigned int j = i + 1; j < m_nWires; ++j) {
      if (wrong[j]) continue;
      double xsepar = 0.;
      double ysepar = 0.;
      if (m_tube) {
        if (m_pery) {
          double xaux1, xaux2, yaux1, yaux2;
          Cartesian2Polar(m_w[i].x, m_w[i].y, xaux1, yaux1);
          Cartesian2Polar(m_w[j].x, m_w[j].y, xaux2, yaux2);
          yaux1 -= m_sy * round(yaux1 / m_sy);
          yaux2 -= m_sy * round(yaux2 / m_sy);
          Polar2Cartesian(xaux1, yaux1, xaux1, yaux1);
          Polar2Cartesian(xaux2, yaux2, xaux2, yaux2);
          xsepar = xaux1 - xaux2;
          ysepar = yaux1 - yaux2;
        } else {
          xsepar = m_w[i].x - m_w[j].x;
          ysepar = m_w[i].y - m_w[j].y;
        }
      } else {
        xsepar = fabs(m_w[i].x - m_w[j].x);
        if (m_perx) xsepar -= m_sx * round(xsepar / m_sx);
        ysepar = fabs(m_w[i].y - m_w[j].y);
        if (m_pery) ysepar -= m_sy * round(ysepar / m_sy);
      }
      const double rij = m_w[i].r + m_w[j].r;
      if (xsepar * xsepar + ysepar * ysepar > rij * rij) continue;
      double xprti = m_w[i].x;
      double yprti = m_w[i].y;
      double xprtj = m_w[j].x;
      double yprtj = m_w[j].y;
      if (m_polar) {
        Internal2Polar(m_w[i].x, m_w[i].y, xprti, yprti);
        Internal2Polar(m_w[j].x, m_w[j].y, xprtj, yprtj);
      }
      std::cerr << m_className << "::CellCheck:\n    Wires " << m_w[i].type
                << " at (" << xprti << ", " << yprti << ") and " << m_w[j].type
                << " at (" << xprtj << ", " << yprtj
                << ") overlap at least partially.\n"
                << "    The latter wire is removed.\n";
      wrong[j] = true;
    }
  }

  // Remove the wires which are not acceptable for one reason or another.
  const int iWires = m_nWires;
  m_nWires = 0;
  for (int i = 0; i < iWires; ++i) {
    if (!wrong[i]) {
      m_w[m_nWires] = m_w[i];
      ++m_nWires;
    }
  }

  // Ensure that some elements are left.
  int nElements = m_nWires;
  if (m_ynplan[0]) ++nElements;
  if (m_ynplan[1]) ++nElements;
  if (m_ynplan[2]) ++nElements;
  if (m_ynplan[3]) ++nElements;
  if (m_tube) ++nElements;

  if (nElements < 2) {
    std::cerr << m_className << "::CellCheck:\n";
    std::cerr << "    At least 2 elements are necessary.\n";
    std::cerr << "    Cell rejected.\n";
    return false;
  }

  // Determine maximum and minimum coordinates and potentials.
  bool setx = false;
  bool sety = false;
  bool setz = false;
  bool setv = false;

  m_xmin = m_xmax = 0.;
  m_ymin = m_ymax = 0.;
  m_zmin = m_zmax = 0.;
  m_vmin = m_vmax = 0.;

  // Loop over the wires.
  for (const auto& wire : m_w) {
    const double rw = wire.r;
    if (setx) {
      m_xmin = std::min(m_xmin, wire.x - rw);
      m_xmax = std::max(m_xmax, wire.x + rw);
    } else {
      m_xmin = wire.x - rw;
      m_xmax = wire.x + rw;
      setx = true;
    }
    if (sety) {
      m_ymin = std::min(m_ymin, wire.y - rw);
      m_ymax = std::max(m_ymax, wire.y + rw);
    } else {
      m_ymin = wire.y - rw;
      m_ymax = wire.y + rw;
      sety = true;
    }
    if (setz) {
      m_zmin = std::min(m_zmin, -0.5 * wire.u);
      m_zmax = std::max(m_zmax, +0.5 * wire.u);
    } else {
      m_zmin = -0.5 * wire.u;
      m_zmax = +0.5 * wire.u;
      setz = true;
    }
    if (setv) {
      m_vmin = std::min(m_vmin, wire.v);
      m_vmax = std::max(m_vmax, wire.v);
    } else {
      m_vmin = m_vmax = wire.v;
      setv = true;
    }
  }
  // Consider the planes.
  for (int i = 0; i < 4; ++i) {
    if (!m_ynplan[i]) continue;
    if (i < 2) {
      if (setx) {
        m_xmin = std::min(m_xmin, m_coplan[i]);
        m_xmax = std::max(m_xmax, m_coplan[i]);
      } else {
        m_xmin = m_xmax = m_coplan[i];
        setx = true;
      }
    } else {
      if (sety) {
        m_ymin = std::min(m_ymin, m_coplan[i]);
        m_ymax = std::max(m_ymax, m_coplan[i]);
      } else {
        m_ymin = m_ymax = m_coplan[i];
        sety = true;
      }
    }
    if (setv) {
      m_vmin = std::min(m_vmin, m_vtplan[i]);
      m_vmax = std::max(m_vmax, m_vtplan[i]);
    } else {
      m_vmin = m_vmax = m_vtplan[i];
      setv = true;
    }
  }

  // Consider the tube.
  if (m_tube) {
    m_xmin = -1.1 * m_cotube;
    m_xmax = +1.1 * m_cotube;
    setx = true;
    m_ymin = -1.1 * m_cotube;
    m_ymax = +1.1 * m_cotube;
    sety = true;
    m_vmin = std::min(m_vmin, m_vttube);
    m_vmax = std::max(m_vmax, m_vttube);
    setv = true;
  }

  // In case of x-periodicity, XMAX-XMIN should be SX,
  if (m_perx && m_sx > (m_xmax - m_xmin)) {
    m_xmin = -0.5 * m_sx;
    m_xmax = +0.5 * m_sx;
    setx = true;
  }
  // in case of y-periodicity, YMAX-YMIN should be SY,
  if (m_pery && m_sy > (m_ymax - m_ymin)) {
    m_ymin = -0.5 * m_sy;
    m_ymax = +0.5 * m_sy;
    sety = true;
  }
  // in case the cell is polar, the y range should be < 2 pi.
  if (m_polar && (m_ymax - m_ymin) >= TwoPi) {
    m_ymin = -Pi;
    m_ymax = +Pi;
    sety = true;
  }

  // Fill in missing dimensions.
  if (setx && m_xmin != m_xmax && (m_ymin == m_ymax || !sety)) {
    m_ymin -= 0.5 * fabs(m_xmax - m_xmin);
    m_ymax += 0.5 * fabs(m_xmax - m_xmin);
    sety = true;
  }
  if (sety && m_ymin != m_ymax && (m_xmin == m_xmax || !setx)) {
    m_xmin -= 0.5 * fabs(m_ymax - m_ymin);
    m_xmax += 0.5 * fabs(m_ymax - m_ymin);
    setx = true;
  }

  if (!setz) {
    m_zmin = -0.25 * (fabs(m_xmax - m_xmin) + fabs(m_ymax - m_ymin));
    m_zmax = +0.25 * (fabs(m_xmax - m_xmin) + fabs(m_ymax - m_ymin));
    setz = true;
  }

  // Ensure that all dimensions are now set.
  if (!(setx && sety && setz)) {
    std::cerr << m_className << "::CellCheck:\n";
    std::cerr << "    Unable to establish"
              << " default dimensions in all directions.\n";
  }

  // Check that at least some different voltages are present.
  if (m_vmin == m_vmax || !setv) {
    std::cerr << m_className << "::CellCheck:\n";
    std::cerr << "    All potentials in the cell are the same.\n";
    std::cerr << "    There is no point in going on.\n";
    return false;
  }

  // Cell seems to be alright since it passed all critical tests.
  return true;
}

bool ComponentAnalyticField::WireCheck() const {
  //-----------------------------------------------------------------------
  //    CELWCH - Subroutine checking the wire positions only, contrary
  //             to CELCHK, this routine does not modify the cell.
  //-----------------------------------------------------------------------

  if (m_nWires == 0) return false;

  if (m_nWires == 1 &&
      !(m_ynplan[0] || m_ynplan[1] || m_ynplan[2] || m_ynplan[3]) && !m_tube) {
    return false;
  }
  // Check position relative to the planes.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    if (m_ynplan[0] && wire.x - wire.r <= m_coplan[0]) return false;
    if (m_ynplan[1] && wire.x + wire.r >= m_coplan[1]) return false;
    if (m_ynplan[2] && wire.y - wire.r <= m_coplan[2]) return false;
    if (m_ynplan[3] && wire.y + wire.r >= m_coplan[3]) return false;
    if (m_tube) {
      if (!InTube(wire.x, wire.y, m_cotube, m_ntube)) return false;
    } else if ((m_perx && 2 * wire.r >= m_sx) || 
               (m_pery && 2 * wire.r >= m_sy)) {
      return false;
    }
  }
  // Check the wire spacing.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const double xi = m_w[i].x;
    const double yi = m_w[i].y;
    for (unsigned int j = i + 1; j < m_nWires; ++j) {
      const double xj = m_w[j].x;
      const double yj = m_w[j].y;
      double xsepar = std::abs(xi - xj);
      double ysepar = std::abs(yi - yj);
      if (m_tube) {
        if (m_pery) {
          double xaux1 = 0., yaux1 = 0.;
          double xaux2 = 0., yaux2 = 0.;
          Cartesian2Polar(xi, yi, xaux1, yaux1);
          Cartesian2Polar(xj, yj, xaux2, yaux2);
          yaux1 -= m_sy * round(yaux1 / m_sy);
          yaux2 -= m_sy * round(yaux2 / m_sy);
          Polar2Cartesian(xaux1, yaux1, xaux1, yaux1);
          Polar2Cartesian(xaux2, yaux2, xaux2, yaux2);
          xsepar = xaux1 - xaux2;
          ysepar = yaux1 - yaux2;
        }
      } else {
        if (m_perx) xsepar -= m_sx * round(xsepar / m_sx);
        if (m_pery) ysepar -= m_sy * round(ysepar / m_sy);
      }
      const double rij = m_w[i].r + m_w[j].r;
      if (xsepar * xsepar + ysepar * ysepar < rij * rij) return false;
    }
  }
  return true;
}

bool ComponentAnalyticField::CellType() {
  // Tube geometries
  if (m_tube) {
    if (m_ntube == 0) {
      if (m_pery) {
        m_cellType = D20;
      } else {
        m_cellType = D10;
      }
    } else if (m_ntube >= 3 && m_ntube <= 8) {
      if (m_pery) {
        m_cellType = D40;
      } else {
        m_cellType = D30;
      }
    } else {
      std::cerr << m_className << "::CellType:\n"
                << "    Potentials for tube with " << m_ntube
                << " edges are not yet available.\n"
                << "    Using a round tube instead.\n";
      m_ntube = 0;
      m_cellType = D30;
    }
    return true;
  }

  // Find the 'A' type cell.
  if (!(m_perx || m_pery) && !(m_ynplan[0] && m_ynplan[1]) &&
      !(m_ynplan[2] && m_ynplan[3])) {
    m_cellType = A00;
    return true;
  }

  // Find the 'B1X' type cell.
  if (m_perx && !m_pery && !(m_ynplan[0] || m_ynplan[1]) &&
      !(m_ynplan[2] && m_ynplan[3])) {
    m_cellType = B1X;
    return true;
  }

  // Find the 'B1Y' type cell.
  if (m_pery && !m_perx && !(m_ynplan[0] && m_ynplan[1]) &&
      !(m_ynplan[2] || m_ynplan[3])) {
    m_cellType = B1Y;
    return true;
  }

  // Find the 'B2X' type cell.
  if (m_perx && !m_pery && !(m_ynplan[2] && m_ynplan[3])) {
    m_cellType = B2X;
    return true;
  }

  if (!(m_perx || m_pery) && !(m_ynplan[2] && m_ynplan[3]) &&
      (m_ynplan[0] && m_ynplan[1])) {
    m_sx = fabs(m_coplan[1] - m_coplan[0]);
    m_cellType = B2X;
    return true;
  }

  // Find the 'B2Y' type cell.
  if (m_pery && !m_perx && !(m_ynplan[0] && m_ynplan[1])) {
    m_cellType = B2Y;
    return true;
  }

  if (!(m_perx || m_pery) && !(m_ynplan[0] && m_ynplan[1]) &&
      (m_ynplan[2] && m_ynplan[3])) {
    m_sy = fabs(m_coplan[3] - m_coplan[2]);
    m_cellType = B2Y;
    return true;
  }

  // Find the 'C1 ' type cell.
  if (!(m_ynplan[0] || m_ynplan[1] || m_ynplan[2] || m_ynplan[3]) && m_perx &&
      m_pery) {
    m_cellType = C10;
    return true;
  }

  // Find the 'C2X' type cell.
  if (!((m_ynplan[2] && m_pery) || (m_ynplan[2] && m_ynplan[3]))) {
    if (m_ynplan[0] && m_ynplan[1]) {
      m_sx = fabs(m_coplan[1] - m_coplan[0]);
      m_cellType = C2X;
      return true;
    }
    if (m_perx && m_ynplan[0]) {
      m_cellType = C2X;
      return true;
    }
  }

  // Find the 'C2Y' type cell.
  if (!((m_ynplan[0] && m_perx) || (m_ynplan[0] && m_ynplan[1]))) {
    if (m_ynplan[2] && m_ynplan[3]) {
      m_sy = fabs(m_coplan[3] - m_coplan[2]);
      m_cellType = C2Y;
      return true;
    }
    if (m_pery && m_ynplan[2]) {
      m_cellType = C2Y;
      return true;
    }
  }

  // Find the 'C3 ' type cell.
  if (m_perx && m_pery) {
    m_cellType = C30;
    return true;
  }

  if (m_perx) {
    m_sy = fabs(m_coplan[3] - m_coplan[2]);
    m_cellType = C30;
    return true;
  }

  if (m_pery) {
    m_sx = fabs(m_coplan[1] - m_coplan[0]);
    m_cellType = C30;
    return true;
  }

  if (m_ynplan[0] && m_ynplan[1] && m_ynplan[2] && m_ynplan[3]) {
    m_sx = fabs(m_coplan[1] - m_coplan[0]);
    m_sy = fabs(m_coplan[3] - m_coplan[2]);
    m_cellType = C30;
    return true;
  }

  // Cell is not recognised.
  return false;
}

std::string ComponentAnalyticField::GetCellType(const Cell) const {
  switch (m_cellType) {
    case A00:
      return "A  ";
    case B1X:
      return "B1X";
    case B1Y:
      return "B1Y";
    case B2X:
      return "B2X";
    case B2Y:
      return "B2Y";
    case C10:
      return "C1 ";
    case C2X:
      return "C2X";
    case C2Y:
      return "C2Y";
    case C30:
      return "C3 ";
    case D10:
      return "D1 ";
    case D20:
      return "D2 ";
    case D30:
      return "D3 ";
    case D40:
      return "D4 ";
    default:
      break;
  }
  return "Unknown";
}

bool ComponentAnalyticField::PrepareStrips() {
  // -----------------------------------------------------------------------
  //    CELSTR - Assigns default anode-cathode gaps, if applicable.
  //    (Last changed on  7/12/00.)
  // -----------------------------------------------------------------------

  double gapDef[4] = {0., 0., 0., 0.};

  // Compute default gaps.
  if (m_ynplan[0]) {
    if (m_ynplan[1]) {
      gapDef[0] = m_coplan[1] - m_coplan[0];
    } else if (m_nWires <= 0) {
      gapDef[0] = -1.;
    } else {
      gapDef[0] = m_w[0].x - m_coplan[0];
      for (const auto& wire : m_w) {
        gapDef[0] = std::min(wire.x - m_coplan[0], gapDef[0]);
      }
    }
  }

  if (m_ynplan[1]) {
    if (m_ynplan[0]) {
      gapDef[1] = m_coplan[1] - m_coplan[0];
    } else if (m_nWires <= 0) {
      gapDef[1] = -1.;
    } else {
      gapDef[1] = m_coplan[1] - m_w[0].x;
      for (const auto& wire : m_w) {
        gapDef[1] = std::min(m_coplan[1] - wire.x, gapDef[1]);
      }
    }
  }

  if (m_ynplan[2]) {
    if (m_ynplan[3]) {
      gapDef[2] = m_coplan[3] - m_coplan[2];
    } else if (m_nWires <= 0) {
      gapDef[2] = -1.;
    } else {
      gapDef[2] = m_w[0].y - m_coplan[2];
      for (const auto& wire : m_w) {
        gapDef[2] = std::min(wire.y - m_coplan[2], gapDef[2]);
      }
    }
  }

  if (m_ynplan[3]) {
    if (m_ynplan[2]) {
      gapDef[3] = m_coplan[3] - m_coplan[2];
    } else if (m_nWires <= 0) {
      gapDef[3] = -1.;
    } else {
      gapDef[3] = m_coplan[3] - m_w[0].y;
      for (const auto& wire : m_w) {
        gapDef[3] = std::min(m_coplan[3] - wire.y, gapDef[3]);
      }
    }
  }

  // Assign.
  for (unsigned int i = 0; i < 4; ++i) {
    for (auto& strip : m_planes[i].strips1) {
      if (strip.gap < 0. && gapDef[i] < 0.) {
        std::cerr << m_className << "::PrepareStrips:\n"
                  << "    Not able to set a default anode-cathode gap\n";
        if (m_polar) {
          std::cerr << "    for r/phi-strips of plane " << i << ".\n";
        } else {
          std::cerr << "    for x/y-strips of plane " << i << ".\n";
        }
        return false;
      }
      if (strip.gap < 0.) {
        strip.gap = gapDef[i];
      } else if (m_polar && i < 2) {
        if (i == 0) {
          strip.gap = log1p(strip.gap / exp(m_coplan[i]));
        } else {
          strip.gap = -log1p(-strip.gap / exp(m_coplan[i]));
        }
      }
    }
    for (auto& strip : m_planes[i].strips2) {
      if (strip.gap < 0. && gapDef[i] < 0.) {
        std::cerr << m_className << "::PrepareStrips:\n"
                  << "    Not able to set a default anode-cathode gap\n"
                  << "    for z-strips of plane " << i << ".\n";
        return false;
      }
      if (strip.gap < 0.) {
        strip.gap = gapDef[i];
      } else if (m_polar && i < 2) {
        if (i == 0) {
          strip.gap = log1p(strip.gap / exp(m_coplan[i]));
        } else {
          strip.gap = -log1p(-strip.gap / exp(m_coplan[i]));
        }
      }
    }
    for (auto& pixel : m_planes[i].pixels) {
      if (pixel.gap < 0. && gapDef[i] < 0.) {
        std::cerr << m_className << "::PrepareStrips:\n"
                  << "    Not able to set a default anode-cathode gap\n"
                  << "    for pixels on plane " << i << ".\n";
        return false;
      }
      if (pixel.gap < 0.) {
        pixel.gap = gapDef[i];
      } else if (m_polar && i < 2) {
        if (i == 0) {
          pixel.gap = log1p(pixel.gap / exp(m_coplan[i]));
        } else {
          pixel.gap = -log1p(-pixel.gap / exp(m_coplan[i]));
        }
      }
    }
  }

  return true;
}

void ComponentAnalyticField::AddReadout(const std::string& label) {
  // Check if this readout group already exists.
  if (std::find(m_readout.begin(), m_readout.end(), label) != m_readout.end()) {
    std::cout << m_className << "::AddReadout:\n";
    std::cout << "    Readout group " << label << " already exists.\n";
    return;
  }
  m_readout.push_back(label);

  unsigned int nWiresFound = 0;
  for (const auto& wire : m_w) {
    if (wire.type == label) ++nWiresFound;
  }

  unsigned int nPlanesFound = 0;
  unsigned int nStripsFound = 0;
  unsigned int nPixelsFound = 0;
  for (const auto& plane : m_planes) {
    if (plane.type == label) ++nPlanesFound;
    for (const auto& strip : plane.strips1) {
      if (strip.type == label) ++nStripsFound;
    }
    for (const auto& strip : plane.strips2) {
      if (strip.type == label) ++nStripsFound;
    }
    for (const auto& pixel : plane.pixels) {
      if (pixel.type == label) ++nPixelsFound;
    }
  }

  if (nWiresFound == 0 && nPlanesFound == 0 && nStripsFound == 0 &&
      nPixelsFound == 0) {
    std::cerr << m_className << "::AddReadout:\n";
    std::cerr << "    At present there are no wires, planes or strips\n";
    std::cerr << "    associated to readout group " << label << ".\n";
  } else {
    std::cout << m_className << "::AddReadout:\n";
    std::cout << "    Readout group " << label << " comprises:\n";
    if (nWiresFound > 1) {
      std::cout << "      " << nWiresFound << " wires\n";
    } else if (nWiresFound == 1) {
      std::cout << "      1 wire\n";
    }
    if (nPlanesFound > 1) {
      std::cout << "      " << nPlanesFound << " planes\n";
    } else if (nPlanesFound == 1) {
      std::cout << "      1 plane\n";
    }
    if (nStripsFound > 1) {
      std::cout << "      " << nStripsFound << " strips\n";
    } else if (nStripsFound == 1) {
      std::cout << "      1 strip\n";
    }
    if (nPixelsFound > 1) {
      std::cout << "      " << nPixelsFound << " pixels\n";
    } else if (nPixelsFound == 1) {
      std::cout << "      1 pixel\n";
    }
  }
  m_sigset = false;
}

void ComponentAnalyticField::SetNumberOfCellCopies(const unsigned int nf) {
  if (nf > 0) {
    if ((nf & (nf - 1)) != 0) {
      std::cerr << m_className << "::SetNumberOfCellCopies:\n"
                << "    Argument must be an integral power of 2.\n";
      return;
    }
  }
  m_nFourier = nf;
  m_sigset = false;
}

bool ComponentAnalyticField::Setup() {
  //-----------------------------------------------------------------------
  //     SETUP  - Routine calling the appropriate setup routine.
  //     (Last changed on 19/ 9/07.)
  //-----------------------------------------------------------------------

  // Set a separate set of plane variables to avoid repeated loops.
  if (m_ynplan[0]) {
    m_coplax = m_coplan[0];
    m_ynplax = true;
  } else if (m_ynplan[1]) {
    m_coplax = m_coplan[1];
    m_ynplax = true;
  } else {
    m_ynplax = false;
  }

  if (m_ynplan[2]) {
    m_coplay = m_coplan[2];
    m_ynplay = true;
  } else if (m_ynplan[3]) {
    m_coplay = m_coplan[3];
    m_ynplay = true;
  } else {
    m_ynplay = false;
  }

  // Set the correction parameters for the planes.
  if (m_tube) {
    m_corvta = 0.;
    m_corvtb = 0.;
    m_corvtc = m_vttube;
  } else if ((m_ynplan[0] && m_ynplan[1]) && !(m_ynplan[2] || m_ynplan[3])) {
    m_corvta = (m_vtplan[0] - m_vtplan[1]) / (m_coplan[0] - m_coplan[1]);
    m_corvtb = 0.;
    m_corvtc = (m_vtplan[1] * m_coplan[0] - m_vtplan[0] * m_coplan[1]) /
               (m_coplan[0] - m_coplan[1]);
  } else if ((m_ynplan[2] && m_ynplan[3]) && !(m_ynplan[0] || m_ynplan[1])) {
    m_corvta = 0.;
    m_corvtb = (m_vtplan[2] - m_vtplan[3]) / (m_coplan[2] - m_coplan[3]);
    m_corvtc = (m_vtplan[3] * m_coplan[2] - m_vtplan[2] * m_coplan[3]) /
               (m_coplan[2] - m_coplan[3]);
  } else {
    m_corvta = m_corvtb = m_corvtc = 0.;
    if (m_ynplan[0]) m_corvtc = m_vtplan[0];
    if (m_ynplan[1]) m_corvtc = m_vtplan[1];
    if (m_ynplan[2]) m_corvtc = m_vtplan[2];
    if (m_ynplan[3]) m_corvtc = m_vtplan[3];
  }

  // Skip wire calculations if there aren't any.
  if (m_nWires <= 0) return true;

  bool ok = true;
  // Call the set routine appropriate for the present cell type.
  switch (m_cellType) {
    case A00:
      ok = SetupA00();
      break;
    case B1X:
      ok = SetupB1X();
      break;
    case B1Y:
      ok = SetupB1Y();
      break;
    case B2X:
      ok = SetupB2X();
      break;
    case B2Y:
      ok = SetupB2Y();
      break;
    case C10:
      ok = SetupC10();
      break;
    case C2X:
      ok = SetupC2X();
      break;
    case C2Y:
      ok = SetupC2Y();
      break;
    case C30:
      ok = SetupC30();
      break;
    case D10:
      ok = SetupD10();
      break;
    case D20:
      ok = SetupD20();
      break;
    case D30:
      ok = SetupD30();
      break;
    default:
      std::cerr << m_className << "::Setup: Unknown cell type.\n";
      break;
  }

  if (!ok) {
    std::cerr << m_className << "::Setup:\n"
              << "    Preparing the cell for field calculations"
              << " did not succeed.\n";
    return false;
  }
  return true;
}

bool ComponentAnalyticField::SetupA00() {
  //-----------------------------------------------------------------------
  //   SETA00 - Subroutine preparing the field calculations by calculating
  //            the charges on the wires, for the cell with one charge and
  //            not more than one plane in either x or y.
  //            The potential used is log(r).
  //-----------------------------------------------------------------------
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wire combinations.
  for (size_t i = 0; i < m_nWires; ++i) {
    a[i][i] = m_w[i].r * m_w[i].r;
    // Take care of the equipotential planes.
    const double xi = m_w[i].x;
    const double yi = m_w[i].y;
    if (m_ynplax) a[i][i] /= 4. * pow(xi - m_coplax, 2);
    if (m_ynplay) a[i][i] /= 4. * pow(yi - m_coplay, 2);
    // Take care of combinations of equipotential planes.
    if (m_ynplax && m_ynplay)
      a[i][i] *= 4. * (pow(xi - m_coplax, 2) + pow(yi - m_coplay, 2));
    // Define the final version of a[i][i].
    a[i][i] = -0.5 * log(a[i][i]);
    // Loop over all other wires for the off-diagonal elements.
    for (size_t j = i + 1; j < m_nWires; ++j) {
      const double xj = m_w[j].x;
      const double yj = m_w[j].y;
      a[i][j] = pow(xi - xj, 2) + pow(yi - yj, 2);
      // Take care of equipotential planes.
      if (m_ynplax) {
        a[i][j] = a[i][j] / (pow(xi + xj - 2. * m_coplax, 2) +
                             pow(yi - yj, 2));
      }
      if (m_ynplay) {
        a[i][j] = a[i][j] / (pow(xi - xj, 2) +
                             pow(yi + yj - 2. * m_coplay, 2));
      }
      // Take care of pairs of equipotential planes in different directions.
      if (m_ynplax && m_ynplay) {
        a[i][j] *= pow(xi + xj - 2. * m_coplax, 2) +
                   pow(yi + yj - 2. * m_coplay, 2);
      }
      // Define a final version of a[i][j].
      a[i][j] = -0.5 * log(a[i][j]);
      // Copy this to a[j][i] since the capacitance matrix is symmetric.
      a[j][i] = a[i][j];
    }
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupB1X() {
  //-----------------------------------------------------------------------
  //   SETB1X - Routine preparing the field calculations by filling the
  //            c-matrix, the potential used is re(log(sin Pi/s (z-z0))).
  //   VARIABLES : xx         : Difference in x of two wires * factor.
  //               yy         : Difference in y of two wires * factor.
  //               yymirr     : Difference in y of one wire and the mirror
  //                            image of another * factor.
  //               r2plan     : Periodic length of (xx,yymirr)
  //-----------------------------------------------------------------------

  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wires and calculate the diagonal elements first.
  for (size_t i = 0; i < m_nWires; ++i) {
    a[i][i] = -log(m_w[i].r * Pi / m_sx);
    // Take care of a plane at constant y if it exists.
    if (m_ynplay) {
      const double yy = (Pi / m_sx) * 2. * (m_w[i].y - m_coplay);
      if (fabs(yy) > 20.) {
        a[i][i] += fabs(yy) - CLog2;
      } else {
        a[i][i] += log(fabs(sinh(yy)));
      }
    }
    // Loop over all other wires to obtain off-diagonal elements.
    for (size_t j = i + 1; j < m_nWires; ++j) {
      const double xx = (Pi / m_sx) * (m_w[i].x - m_w[j].x);
      const double yy = (Pi / m_sx) * (m_w[i].y - m_w[j].y);
      if (fabs(yy) > 20.) {
        a[i][j] = -fabs(yy) + CLog2;
      } else {
        const double sinhy = sinh(yy);
        const double sinx = sin(xx);
        a[i][j] = -0.5 * log(sinhy * sinhy + sinx * sinx);
      }
      // Take equipotential planes into account if they exist.
      if (m_ynplay) {
        double r2plan = 0.;
        const double yymirr = (Pi / m_sx) * (m_w[i].y + m_w[j].y - 2. * m_coplay);
        if (fabs(yymirr) > 20.) {
          r2plan = fabs(yymirr) - CLog2;
        } else {
          const double sinhy = sinh(yymirr);
          const double sinx = sin(xx);
          r2plan = 0.5 * log(sinhy * sinhy + sinx * sinx);
        }
        a[i][j] += r2plan;
      }
      // Copy a[i][j] to a[j][i], the capactance matrix is symmetric.
      a[j][i] = a[i][j];
    }
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupB1Y() {
  //-----------------------------------------------------------------------
  //   SETB1Y - Routine preparing the field calculations by setting the
  //            charges. The potential used is Re log(sinh Pi/sy(z-z0)).
  //   VARIABLES : yy         : Difference in y of two wires * factor.
  //               xxmirr     : Difference in x of one wire and the mirror
  //                            image of another * factor.
  //               r2plan     : Periodic length of (xxmirr,yy).
  //-----------------------------------------------------------------------

  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wires and calculate the diagonal elements first.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    a[i][i] = -log(m_w[i].r * Pi / m_sy);
    // Take care of a plane at constant x if it exists.
    if (m_ynplax) {
      const double xx = (Pi / m_sy) * 2. * (m_w[i].x - m_coplax);
      if (fabs(xx) > 20.) {
        a[i][i] += fabs(xx) - CLog2;
      } else {
        a[i][i] += log(fabs(sinh(xx)));
      }
    }
    // Loop over all other wires to obtain off-diagonal elements.
    for (unsigned int j = i + 1; j < m_nWires; ++j) {
      const double xx = (Pi / m_sy) * (m_w[i].x - m_w[j].x);
      const double yy = (Pi / m_sy) * (m_w[i].y - m_w[j].y);
      if (fabs(xx) > 20.) {
        a[i][j] = -fabs(xx) + CLog2;
      } else {
        const double sinhx = sinh(xx);
        const double siny = sin(yy);
        a[i][j] = -0.5 * log(sinhx * sinhx + siny * siny);
      }
      // Take care of a plane at constant x.
      if (m_ynplax) {
        const double xxmirr = (Pi / m_sy) * (m_w[i].x + m_w[j].x - 2. * m_coplax);
        double r2plan = 0.;
        if (fabs(xxmirr) > 20.) {
          r2plan = fabs(xxmirr) - CLog2;
        } else {
          const double sinhx = sinh(xxmirr);
          const double siny = sin(yy);
          r2plan = 0.5 * log(sinhx * sinhx + siny * siny);
        }
        a[i][j] += r2plan;
      }
      // Copy a[i][j] to a[j][i], the capacitance matrix is symmetric.
      a[j][i] = a[i][j];
    }
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupB2X() {
  //-----------------------------------------------------------------------
  //   SETB2X - Routine preparing the field calculations by setting the
  //            charges.
  //   VARIABLES : xx         : Difference in x of two wires * factor.
  //               yy         : Difference in y of two wires * factor.
  //               xxneg      : Difference in x of one wire and the mirror
  //                            image in period direction of another * fac.
  //               yymirr     : Difference in y of one wire and the mirror
  //                            image of another * factor.
  //-----------------------------------------------------------------------

  const double tx = Pi / m_sx;  
  m_b2sin.resize(m_nWires);
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wires and calculate the diagonal elements first.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = tx * (m_w[i].x - m_coplax);
    a[i][i] = 0.5 * tx * m_w[i].r / sin(xx);
    // Take care of a plane at constant y if it exists.
    if (m_ynplay) {
      const double yymirr = tx * (m_w[i].y - m_coplay);
      if (fabs(yymirr) <= 20.) {
        const double sinhy = sinh(yymirr);
        const double sinx = sin(xx);
        a[i][i] *= sqrt(sinhy * sinhy + sinx * sinx) / sinhy;
      }
    }
    // Store the true value of a[i][i].
    a[i][i] = -log(fabs(a[i][i]));
  }

  for (size_t i = 0; i < m_nWires; ++i) {
    // Loop over all other wires to obtain off-diagonal elements.
    for (size_t j = i + 1; j < m_nWires; ++j) {
      const double xx = 0.5 * tx * (m_w[i].x - m_w[j].x);
      const double yy = 0.5 * tx * (m_w[i].y - m_w[j].y);
      const double xxneg = 0.5 * tx * (m_w[i].x + m_w[j].x - 2. * m_coplax);
      if (fabs(yy) <= 20.) {
        const double sinhy = sinh(yy);
        const double sinxx = sin(xx);
        const double sinxxneg = sin(xxneg);
        a[i][j] = (sinhy * sinhy + sinxx * sinxx) /
                  (sinhy * sinhy + sinxxneg * sinxxneg);
      } else {
        a[i][j] = 1.0;
      }
      // Take an equipotential plane at constant y into account.
      if (m_ynplay) {
        const double yymirr = 0.5 * tx * (m_w[i].y + m_w[j].y - 2. * m_coplay);
        if (fabs(yymirr) <= 20.) {
          const double sinhy = sinh(yymirr);
          const double sinxx = sin(xx);
          const double sinxxneg = sin(xxneg);
          a[i][j] *= (sinhy * sinhy + sinxxneg * sinxxneg) /
                     (sinhy * sinhy + sinxx * sinxx);
        }
      }
      // Store the true value of a[i][j] in both a[i][j] and a[j][i].
      a[i][j] = -0.5 * log(a[i][j]);
      a[j][i] = a[i][j];
    }
    // Set the b2sin vector.
    m_b2sin[i] = sin(tx * (m_coplax - m_w[i].x));
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupB2Y() {
  //-----------------------------------------------------------------------
  //   SETB2Y - Routine preparing the field calculations by setting the
  //            charges.
  //   VARIABLES : xx         : Difference in x of two wires * factor.
  //               yy         : Difference in y of two wires * factor.
  //               xxmirr     : Difference in x of one wire and the mirror
  //                            image of another * factor.
  //               yyneg      : Difference in y of one wire and the mirror
  //                            image in period direction of another * fac.
  //-----------------------------------------------------------------------

  const double ty = Pi / m_sy;
  m_b2sin.resize(m_nWires);
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wires and calculate the diagonal elements first.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double yy = ty * (m_w[i].y - m_coplay);
    a[i][i] = 0.5 * ty * m_w[i].r / sin(yy);
    // Take care of a plane at constant x if present.
    if (m_ynplax) {
      const double xxmirr = ty * (m_w[i].x - m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const double sinhx = sinh(xxmirr);
        const double sinyy = sin(yy);
        a[i][i] *= sqrt(sinhx * sinhx + sinyy * sinyy) / sinhx;
      }
    }
    // Store the true value of a[i][i].
    a[i][i] = -log(fabs(a[i][i]));
  }
  for (size_t i = 0; i < m_nWires; ++i) {
    // Loop over all other wires to obtain off-diagonal elements.
    for (size_t j = i + 1; j < m_nWires; j++) {
      const double xx = 0.5 * ty * (m_w[i].x - m_w[j].x);
      const double yy = 0.5 * ty * (m_w[i].y - m_w[j].y);
      const double yyneg = 0.5 * ty * (m_w[i].y + m_w[j].y - 2. * m_coplay);
      if (fabs(xx) <= 20.) {
        const double sinhx = sinh(xx);
        const double sinyy = sin(yy);
        const double sinyyneg = sin(yyneg);
        a[i][j] = (sinhx * sinhx + sinyy * sinyy) /
                  (sinhx * sinhx + sinyyneg * sinyyneg);
      } else {
        a[i][j] = 1.0;
      }
      // Take an equipotential plane at constant x into account.
      if (m_ynplax) {
        const double xxmirr = 0.5 * ty * (m_w[i].x + m_w[j].x - 2. * m_coplax);
        if (fabs(xxmirr) <= 20.) {
          const double sinhx = sinh(xxmirr);
          const double sinyy = sin(yy);
          const double sinyyneg = sin(yyneg);
          a[i][j] *= (sinhx * sinhx + sinyyneg * sinyyneg) /
                     (sinhx * sinhx + sinyy * sinyy);
        }
      }
      // Store the true value of a[i][j] in both a[i][j] and a[j][i].
      a[i][j] = -0.5 * log(a[i][j]);
      a[j][i] = a[i][j];
    }
    // Set the b2sin vector.
    m_b2sin[i] = sin(ty * (m_coplay - m_w[i].y));
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupC10() {
  //-----------------------------------------------------------------------
  //   SETC10 - This initialising routine computes the wire charges E and
  //            sets certain constants in common. The wire are located at
  //            (x[j],y[j])+(LX*SX,LY*SY), J=1(1)NWIRE,
  //            LX=-infinity(1)infinity, LY=-infinity(1)infinity.
  //            Use is made of the function PH2.
  //
  //  (Written by G.A.Erskine/DD, 14.8.1984 modified to some extent)
  //-----------------------------------------------------------------------

  // Initialise the constants.
  double p = 0.;
  m_p1 = m_p2 = 0.;

  m_mode = 0;
  if (m_sx <= m_sy) {
    m_mode = 1;
    if (m_sy / m_sx < 8.) p = exp(-Pi * m_sy / m_sx);
    m_zmult = std::complex<double>(Pi / m_sx, 0.);
  } else {
    m_mode = 0;
    if (m_sx / m_sy < 8.) p = exp(-Pi * m_sx / m_sy);
    m_zmult = std::complex<double>(0., Pi / m_sy);
  }
  m_p1 = p * p;
  if (m_p1 > 1.e-10) m_p2 = pow(p, 6);

  if (m_debug) {
    std::cout << m_className << "::SetupC10:\n"
              << "    p, p1, p2 = " << p << ", " << m_p1 << ", " << m_p2
              << "\n"
              << "    zmult = " << m_zmult << "\n"
              << "    mode = " << m_mode << "\n";
  }

  // Store the capacitance matrix.
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  for (size_t i = 0; i < m_nWires; ++i) {
    for (size_t j = 0; j < m_nWires; ++j) {
      const double xyi = m_mode == 0 ? m_w[i].x : m_w[i].y;
      const double xyj = m_mode == 0 ? m_w[j].x : m_w[j].y;
      const double temp = xyi * xyj * TwoPi / (m_sx * m_sy);
      if (i == j) {
        a[i][j] = Ph2Lim(m_w[i].r) - temp;
      } else {
        a[i][j] = Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) - temp;
      }
    }
  }
  // Invert the capacitance matrix to find the charges.
  if (!Charge(a)) return false;
  // Calculate the non-logarithmic term in the potential.
  double s = 0.;
  for (size_t j = 0; j < m_nWires; ++j) {
    const double xyj = m_mode == 0 ? m_w[j].x : m_w[j].y;
    s += m_w[j].e * xyj;
  }
  m_c1 = -s * 2. * Pi / (m_sx * m_sy);
  return true;
}

bool ComponentAnalyticField::SetupC2X() {
  //-----------------------------------------------------------------------
  //   SETC2X - This initializing subroutine stores the capacitance matrix
  //            for the configuration:
  //            wires at zw(j)+cmplx(lx*2*sx,ly*sy),
  //            j=1(1)n, lx=-infinity(1)infinity, ly=-infinity(1)infinity.
  //            but the signs of the charges alternate in the x-direction
  //-----------------------------------------------------------------------

  // Initialise the constants.
  double p = 0.;
  m_p1 = m_p2 = 0.;

  m_mode = 0;
  if (2. * m_sx <= m_sy) {
    m_mode = 1;
    if (m_sy / m_sx < 25.) p = exp(-HalfPi * m_sy / m_sx);
    m_zmult = std::complex<double>(HalfPi / m_sx, 0.);
  } else {
    m_mode = 0;
    if (m_sx / m_sy < 6.) p = exp(-2. * Pi * m_sx / m_sy);
    m_zmult = std::complex<double>(0., Pi / m_sy);
  }
  m_p1 = p * p;
  if (m_p1 > 1.e-10) m_p2 = pow(p, 6);

  if (m_debug) {
    std::cout << m_className << "::SetupC2X:\n"
              << "    p, p1, p2 = " << p << ", " << m_p1 << ", " << m_p2
              << "\n"
              << "    zmult = " << m_zmult << "\n"
              << "    mode = " << m_mode << "\n";
  }

  // Fill the capacitance matrix.
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  for (size_t i = 0; i < m_nWires; ++i) {
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    for (size_t j = 0; j < m_nWires; ++j) {
      double temp = 0.;
      if (m_mode == 0) {
        temp = (m_w[i].x - cx) * (m_w[j].x - cx) * TwoPi / (m_sx * m_sy);
      }
      if (i == j) {
        a[i][i] = Ph2Lim(m_w[i].r) - Ph2(2. * (m_w[i].x - cx), 0.) - temp;
      } else {
        a[i][j] = Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) -
                  Ph2(m_w[i].x + m_w[j].x - 2. * cx, m_w[i].y - m_w[j].y) -
                  temp;
      }
    }
  }
  // Invert the capacitance matrix to find the charges.
  if (!Charge(a)) return false;
  // Determine the non-logarithmic part of the potential (0 if MODE=1).
  m_c1 = 0.;
  if (m_mode == 0) {
    double s = 0.;
    for (unsigned int i = 0; i < m_nWires; ++i) {
      const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
      s += m_w[i].e * (m_w[i].x - cx);
    }
    m_c1 = -s * TwoPi / (m_sx * m_sy);
  }
  return true;
}

bool ComponentAnalyticField::SetupC2Y() {
  //-----------------------------------------------------------------------
  //   SETC2Y - This initializing subroutine stores the capacitance matrix
  //            for the configuration:
  //            wires at zw(j)+cmplx(lx*sx,ly*2*sy),
  //            j=1(1)n, lx=-infinity(1)infinity, ly=-infinity(1)infinity.
  //            but the signs of the charges alternate in the y-direction
  //-----------------------------------------------------------------------

  // Initialise the constants.
  double p = 0.;
  m_p1 = m_p2 = 0.;

  m_mode = 0;
  if (m_sx <= 2. * m_sy) {
    m_mode = 1;
    if (m_sy / m_sx <= 6.) p = exp(-2. * Pi * m_sy / m_sx);
    m_zmult = std::complex<double>(Pi / m_sx, 0.);
  } else {
    m_mode = 0;
    if (m_sx / m_sy <= 25.) p = exp(-HalfPi * m_sx / m_sy);
    m_zmult = std::complex<double>(0., HalfPi / m_sy);
  }
  m_p1 = p * p;
  if (m_p1 > 1.e-10) m_p2 = pow(p, 6);

  if (m_debug) {
    std::cout << m_className << "::SetupC2Y:\n"
              << "    p, p1, p2 = " << p << ", " << m_p1 << ", " << m_p2
              << "\n"
              << "    zmult = " << m_zmult << "\n"
              << "    mode = " << m_mode << "\n";
  }

  // Fill the capacitance matrix.
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  for (size_t i = 0; i < m_nWires; ++i) {
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    for (size_t j = 0; j < m_nWires; ++j) {
      double temp = 0.;
      if (m_mode == 1) {
        temp = (m_w[i].y - cy) * (m_w[j].y - cy) * TwoPi / (m_sx * m_sy);
      }
      if (i == j) {
        a[i][i] = Ph2Lim(m_w[i].r) - Ph2(0., 2. * (m_w[j].y - cy)) - temp;
      } else {
        a[i][j] = Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) -
                  Ph2(m_w[i].x - m_w[j].x, m_w[i].y + m_w[j].y - 2. * cy) -
                  temp;
      }
    }
  }
  // Invert the capacitance matrix to find the charges.
  if (!Charge(a)) return false;
  // The non-logarithmic part of the potential is zero if MODE=0.
  m_c1 = 0.;
  if (m_mode == 1) {
    double s = 0.;
    for (unsigned int i = 0; i < m_nWires; ++i) {
      const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
      s += m_w[i].e * (m_w[i].y - cy);
    }
    m_c1 = -s * TwoPi / (m_sx * m_sy);
  }
  return true;
}

bool ComponentAnalyticField::SetupC30() {
  //-----------------------------------------------------------------------
  //   SETC30 - This initializing subroutine stores the capacitance matrix
  //            for a configuration with
  //            wires at zw(j)+cmplx(lx*2*sx,ly*2*sy),
  //            j=1(1)n, lx=-infinity(1)infinity, ly=-infinity(1)infinity.
  //            but the signs of the charges alternate in both directions.
  //-----------------------------------------------------------------------

  // Initialise the constants.
  double p = 0.;
  m_p1 = m_p2 = 0.;

  m_mode = 0;
  if (m_sx <= m_sy) {
    m_mode = 1;
    if (m_sy / m_sx <= 13.) p = exp(-Pi * m_sy / m_sx);
    m_zmult = std::complex<double>(HalfPi / m_sx, 0.);
  } else {
    m_mode = 0;
    if (m_sx / m_sy <= 13.) p = exp(-Pi * m_sx / m_sy);
    m_zmult = std::complex<double>(0., HalfPi / m_sy);
  }
  m_p1 = p * p;
  if (m_p1 > 1.e-10) m_p2 = pow(p, 6);

  if (m_debug) {
    std::cout << m_className << "::SetupC30:\n"
              << "    p, p1, p2 = " << p << ", " << m_p1 << ", " << m_p2
              << "\n"
              << "    zmult = " << m_zmult << "\n"
              << "    mode = " << m_mode << "\n";
  }

  // Fill the capacitance matrix.
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  for (size_t i = 0; i < m_nWires; ++i) {
    double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    for (size_t j = 0; j < m_nWires; ++j) {
      if (i == j) {
        a[i][i] = Ph2Lim(m_w[i].r) - Ph2(0., 2. * (m_w[i].y - cy)) -
                  Ph2(2. * (m_w[i].x - cx), 0.) +
                  Ph2(2. * (m_w[i].x - cx), 2. * (m_w[i].y - cy));
      } else {
        a[i][j] =
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) -
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y + m_w[j].y - 2. * cy) -
            Ph2(m_w[i].x + m_w[j].x - 2. * cx, m_w[i].y - m_w[j].y) +
            Ph2(m_w[i].x + m_w[j].x - 2. * cx, m_w[i].y + m_w[j].y - 2. * cy);
      }
    }
  }
  // Invert the capacitance matrix to find the charges.
  if (!Charge(a)) return false;
  // The non-logarithmic part of the potential is zero in this case.
  m_c1 = 0.;
  return true;
}

bool ComponentAnalyticField::SetupD10() {
  //-----------------------------------------------------------------------
  //   SETD10 - Subroutine preparing the field calculations by calculating
  //            the charges on the wires, for cells with a tube.
  //
  //   (Last changed on  4/ 9/95.)
  //-----------------------------------------------------------------------

  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set the diagonal terms.
    a[i][i] = -log(m_w[i].r * m_cotube /
                   (m_cotube2 - (m_w[i].x * m_w[i].x + m_w[i].y * m_w[i].y)));
    // Set a complex wire-coordinate to make things a little easier.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    // Loop over all other wires for the off-diagonal elements.
    for (size_t j = i + 1; j < m_nWires; ++j) {
      // Set a complex wire-coordinate to make things a little easier.
      const std::complex<double> zj(m_w[j].x, m_w[j].y);
      a[i][j] = -log(
        std::abs(m_cotube * (zi - zj) / (m_cotube2 - conj(zi) * zj)));
      // Copy this to a[j][i] since the capacitance matrix is symmetric.
      a[j][i] = a[i][j];
    }
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupD20() {
  //-----------------------------------------------------------------------
  //   SETD20 - Subroutine preparing the field calculations by calculating
  //            the charges on the wires, for cells with a tube and a
  //            phi periodicity. Assymetric capacitance matrix !
  //
  //   (Last changed on 18/ 2/93.)
  //-----------------------------------------------------------------------

  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set a complex wire-coordinate to make things a little easier.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    if (std::abs(zi) < m_w[i].r) {
      // Case of a wire near the centre.
      // Inner loop over the wires.
      for (size_t j = 0; j < m_nWires; ++j) {
        if (i == j) {
          // Set the diagonal terms.
          a[i][i] = -log(m_w[i].r / 
              (m_cotube - (m_w[i].x * m_w[i].x + m_w[i].y * m_w[i].y) / m_cotube));
        } else {
          // Off-diagonal terms.
          const std::complex<double> zj(m_w[j].x, m_w[j].y);
          a[j][i] = -log(std::abs((1. / m_cotube) * (zi - zj) /
                                  (1. - conj(zi) * zj / m_cotube2)));
        }
      }
    } else {
      // Normal case.
      // Inner wire loop.
      for (size_t j = 0; j < m_nWires; ++j) {
        if (i == j) {
          // Diagonal elements.
          a[i][i] = -log(
            std::abs(m_w[i].r * m_mtube * pow(zi, m_mtube - 1) /
                     (pow(m_cotube, m_mtube) *
                      (1. - pow((std::abs(zi) / m_cotube), 2 * m_mtube)))));
        } else {
          // Off-diagonal terms.
          const std::complex<double> zj(m_w[j].x, m_w[j].y);
          a[j][i] = -log(
            std::abs((1. / pow(m_cotube, m_mtube)) *
                     (pow(zj, m_mtube) - pow(zi, m_mtube)) /
                     (1. - pow(zj * conj(zi) / m_cotube2, m_mtube))));
        }
      }
    }
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::SetupD30() {
  //-----------------------------------------------------------------------
  //   SETD30 - Subroutine preparing the field calculations by calculating
  //            the charges on the wires, for cells with wires inside a
  //            polygon.
  //
  //   (Last changed on 21/ 2/94.)
  //-----------------------------------------------------------------------

  m_zw.assign(m_nWires, std::complex<double>(0., 0.));

  std::complex<double> wd = std::complex<double>(0., 0.);

  // Evaluate kappa, a constant needed by ConformalMap.
  m_kappa = tgamma((m_ntube + 1.) / m_ntube) *
            tgamma((m_ntube - 2.) / m_ntube) / tgamma((m_ntube - 1.) / m_ntube);
  DMatrix a(m_nWires, std::vector<double>(m_nWires, 0.));
  // Loop over all wire combinations.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Compute wire mappings only once.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    ConformalMap(zi / m_cotube, m_zw[i], wd);
    // Diagonal elements.
    a[i][i] = -log(
        std::abs((m_w[i].r / m_cotube) * wd / (1. - std::norm(m_zw[i]))));

    // Loop over all other wires for the off-diagonal elements.
    for (size_t j = 0; j < i; ++j) {
      a[i][j] = -log(
          std::abs((m_zw[i] - m_zw[j]) / (1. - conj(m_zw[i]) * m_zw[j])));
      // Copy this to a[j][i] since the capacitance matrix is symmetric.
      a[j][i] = a[i][j];
    }
  }
  // Invert the capacitance matrix and calculate the charges.
  return Charge(a);
}

bool ComponentAnalyticField::Charge(DMatrix& a) {
  //-----------------------------------------------------------------------
  //   CHARGE - Routine actually inverting the capacitance matrix filled in
  //            the SET... routines thereby providing the charges.
  //   (Last changed on 30/ 1/93.)
  //-----------------------------------------------------------------------

  // Transfer the voltages to rhs vector,
  // correcting for the equipotential planes.
  std::vector<double> b(m_nWires, 0.);
  for (size_t i = 0; i < m_nWires; ++i) {
    b[i] = m_w[i].v - (m_corvta * m_w[i].x + m_corvtb * m_w[i].y + m_corvtc);
  }

  bool ok = true;
  // Force sum charges = 0 in case of absence of equipotential planes.
  if (!(m_ynplan[0] || m_ynplan[1] || m_ynplan[2] || m_ynplan[3] || m_tube)) {
    // Add extra elements to A, acting as constraints.
    b.push_back(0.);
    a.resize(m_nWires + 1);
    a[m_nWires].clear();
    for (size_t i = 0; i < m_nWires; ++i) {
      a[i].push_back(1.);
      a[m_nWires].push_back(1.);
    }
    a[m_nWires].push_back(0.);
    // Solve equations to yield charges.
    if (Numerics::CERNLIB::deqinv(m_nWires + 1, a, b) != 0) {
      std::cerr << m_className << "::Charge: Matrix inversion failed.\n";
      return false;
    }
    // Modify A to give true inverse of capacitance matrix.
    if (a[m_nWires][m_nWires] != 0.) {
      const double t = 1. / a[m_nWires][m_nWires];
      for (size_t i = 0; i < m_nWires; ++i) {
        for (size_t j = 0; j < m_nWires; ++j) {
          a[i][j] -= t * a[i][m_nWires] * a[m_nWires][j];
        }
      }
    } else {
      std::cerr << m_className << "::Charge:\n"
                << "    True inverse of the capacitance matrix"
                << " could not be calculated.\n";
      ok = false;
    }
    // Store reference potential.
    m_v0 = b[m_nWires];
  } else {
    // Handle the case when the sum of the charges is zero automatically.
    if (Numerics::CERNLIB::deqinv(m_nWires, a, b) != 0) {
      std::cerr << m_className << "::Charge: Matrix inversion failed.\n";
      return false;
    }
    // Reference potential chosen to be zero.
    m_v0 = 0.;
  }

  // Check the error condition flag.
  if (!ok) {
    std::cerr << m_className << "::Charge:\n"
              << "    Failure to solve the capacitance equations.\n"
              << "    No charges are available.\n";
    return false;
  }

  // Copy the charges to E.
  for (size_t i = 0; i < m_nWires; ++i) m_w[i].e = b[i];

  // If debugging is on, print the capacitance matrix.
  if (m_debug) {
    std::cout << m_className << "::Charge:\n"
              << "    Dump of the capacitance matrix after inversion:\n";
    for (size_t i = 0; i < m_nWires; i += 10) {
      for (size_t j = 0; j < m_nWires; j += 10) {
        std::cout << "    (Block " << i / 10 << ", " << j / 10 << ")\n";
        for (unsigned int ii = 0; ii < 10; ++ii) {
          if (i + ii >= m_nWires) break;
          for (unsigned int jj = 0; jj < 10; ++jj) {
            if (j + jj >= m_nWires) break;
            std::cout << std::setw(6) << a[i + ii][j + jj] << " ";
          }
          std::cout << "\n";
        }
        std::cout << "\n";
      }
    }
    std::cout << "    End of the inverted capacitance matrix.\n";
  }

  // And also check the quality of the matrix inversion.
  if (m_chargeCheck) {
    std::cout << m_className << "::Charge:\n"
              << "    Quality check of the charge calculation.\n"
              << "    Wire       E as obtained        E reconstructed\n";
    for (size_t i = 0; i < m_nWires; ++i) {
      b[i] = 0.;
      for (size_t j = 0; j < m_nWires; ++j) {
        b[i] += a[i][j] *
                (m_w[j].v - m_v0 -
                 (m_corvta * m_w[j].x + m_corvtb * m_w[j].y + m_corvtc));
      }
      std::cout << "    " << i << "      " << m_w[i].e << "    " << b[i]
                << "\n";
    }
  }
  return true;
}

double ComponentAnalyticField::Ph2(const double xpos, const double ypos) const {
  //-----------------------------------------------------------------------
  //   PH2    - Logarithmic contribution to real single-wire potential,
  //            for a doubly priodic wire array.
  //   PH2LIM - Entry, PH2LIM(r) corresponds to z on the surface of a wire
  //            of (small) radius r.
  //
  //            Clenshaw's algorithm is used for the evaluation of the sum
  //            ZTERM = SIN(zeta) - P1*SIN(3*zeta) + P2*SIN(5*zeta).
  //
  //  (G.A.Erskine/DD, 14.8.1984; some minor modifications (i) common block
  //   /EV2COM/ incorporated in /CELDAT/ (ii) large imag(zeta) corrected)
  //-----------------------------------------------------------------------

  // Start of the main subroutine, off diagonal elements.
  std::complex<double> zeta = m_zmult * std::complex<double>(xpos, ypos);
  if (fabs(imag(zeta)) < 10.) {
    std::complex<double> zsin = sin(zeta);
    std::complex<double> zcof = 4. * zsin * zsin - 2.;
    std::complex<double> zu = -m_p1 - zcof * m_p2;
    std::complex<double> zunew = 1. - zcof * zu - m_p2;
    std::complex<double> zterm = (zunew + zu) * zsin;
    return -log(std::abs(zterm));
  }

  return -fabs(imag(zeta)) + CLog2;
}

void ComponentAnalyticField::ConformalMap(const std::complex<double>& z,
                                          std::complex<double>& ww,
                                          std::complex<double>& wd) const {
  //-----------------------------------------------------------------------
  //   EFCMAP - Maps a the interior part of a regular in the unit circle.
  //   Variables: Z     - point to be mapped
  //              W     - the image of Z
  //              WD    - derivative of the mapping at Z
  //              CC1   - coefficients for expansion around centre
  //              CC2   - coefficients for expansion around corner
  //   (Last changed on 19/ 2/94.)
  //-----------------------------------------------------------------------

  // Coefficients for centre expansion in triangles, squares, pentagons,
  // hexagons, heptagons, octogons.
  constexpr std::array<std::array<double, 16>, 6> cc1 = {
      {{{0.1000000000e+01, -.1666666865e+00, 0.3174602985e-01, -.5731921643e-02,
         0.1040112227e-02, -.1886279933e-03, 0.3421107249e-04, -.6204730198e-05,
         0.1125329618e-05, -.2040969207e-06, 0.3701631357e-07, -.6713513301e-08,
         0.1217605794e-08, -.2208327132e-09, 0.4005162868e-10,
         -.7264017512e-11}},
       {{0.1000000000e+01, -.1000000238e+00, 0.8333332837e-02, -.7051283028e-03,
         0.5967194738e-04, -.5049648280e-05, 0.4273189802e-06, -.3616123934e-07,
         0.3060091514e-08, -.2589557457e-09, 0.2191374859e-10, -.1854418528e-11,
         0.1569274224e-12, -.1327975205e-13, 0.1123779363e-14,
         -.9509817570e-16}},
       {{0.1000000000e+01, -.6666666269e-01, 0.1212121220e-02, -.2626262140e-03,
         -.3322110570e-04, -.9413293810e-05, -.2570029210e-05, -.7695705904e-06,
         -.2422486887e-06, -.7945993730e-07, -.2691839640e-07, -.9361642128e-08,
         -.3327319087e-08, -.1204430555e-08, -.4428404310e-09,
         -.1650302672e-09}},
       {{0.1000000000e+01, -.4761904851e-01, -.1221001148e-02, -.3753788769e-03,
         -.9415557724e-04, -.2862767724e-04, -.9587882232e-05, -.3441659828e-05,
         -.1299798896e-05, -.5103651119e-06, -.2066504408e-06, -.8578405186e-07,
         -.3635090096e-07, -.1567239494e-07, -.6857355572e-08,
         -.3038770346e-08}},
       {{0.1000000000e+01, -.3571428731e-01, -.2040816238e-02, -.4936389159e-03,
         -.1446709794e-03, -.4963850370e-04, -.1877940667e-04, -.7600909157e-05,
         -.3232265954e-05, -.1427365532e-05, -.6493634714e-06, -.3026190711e-06,
         -.1438593245e-06, -.6953911225e-07, -.3409525462e-07,
         -.1692310647e-07}},
       {{0.1000000000e+01, -.2777777612e-01, -.2246732125e-02, -.5571441725e-03,
         -.1790652314e-03, -.6708275760e-04, -.2766949183e-04, -.1219387286e-04,
         -.5640039490e-05, -.2706697160e-05, -.1337270078e-05, -.6763995657e-06,
         -.3488264610e-06, -.1828456675e-06, -.9718036154e-07,
         -.5227070332e-07}}}};
  // Coefficients for corner expansion.
  constexpr std::array<std::array<double, 16>, 6> cc2 = {
      {{{0.3333333135e+00, -.5555555597e-01, 0.1014109328e-01, -.1837154618e-02,
         0.3332451452e-03, -.6043842586e-04, 0.1096152027e-04, -.1988050826e-05,
         0.3605655365e-06, -.6539443120e-07, 0.1186035448e-07, -.2151069323e-08,
         0.3901317047e-09, -.7075676156e-10, 0.1283289534e-10,
         -.2327455936e-11}},
       {{0.1000000000e+01, -.5000000000e+00, 0.3000000119e+00, -.1750000119e+00,
         0.1016666889e+00, -.5916666612e-01, 0.3442307562e-01, -.2002724260e-01,
         0.1165192947e-01, -.6779119372e-02, 0.3944106400e-02, -.2294691978e-02,
         0.1335057430e-02, -.7767395582e-03, 0.4519091453e-03,
         -.2629216760e-03}},
       {{0.1248050690e+01, -.7788147926e+00, 0.6355384588e+00, -.4899077415e+00,
         0.3713272810e+00, -.2838423252e+00, 0.2174729109e+00, -.1663445234e+00,
         0.1271933913e+00, -.9728997946e-01, 0.7442557812e-01, -.5692918226e-01,
         0.4354400188e-01, -.3330700099e-01, 0.2547712997e-01,
         -.1948769018e-01}},
       {{0.1333333015e+01, -.8888888955e+00, 0.8395061493e+00, -.7242798209e+00,
         0.6016069055e+00, -.5107235312e+00, 0.4393203855e+00, -.3745460510e+00,
         0.3175755739e+00, -.2703750730e+00, 0.2308617830e+00, -.1966916919e+00,
         0.1672732830e+00, -.1424439549e+00, 0.1214511395e+00,
         -.1034612656e+00}},
       {{0.1359752655e+01, -.9244638681e+00, 0.9593217969e+00, -.8771237731e+00,
         0.7490229011e+00, -.6677658558e+00, 0.6196745634e+00, -.5591596961e+00,
         0.4905325770e+00, -.4393517375e+00, 0.4029803872e+00, -.3631100059e+00,
         0.3199430704e+00, -.2866140604e+00, 0.2627358437e+00,
         -.2368256450e+00}},
       {{0.1362840652e+01, -.9286670089e+00, 0.1035511017e+01, -.9800255299e+00,
         0.8315343261e+00, -.7592730522e+00, 0.7612683773e+00, -.7132136226e+00,
         0.6074471474e+00, -.5554352999e+00, 0.5699443221e+00, -.5357525349e+00,
         0.4329345822e+00, -.3916820884e+00, 0.4401986003e+00,
         -.4197303057e+00}}}};

  constexpr int nterm = 15;
  if (z == 0.) {
    // Z coincides with the centre. Results are trivial.
    ww = 0;
    wd = m_kappa;
  } else if (std::abs(z) < 0.75) {
    // Z is close to the centre. Series expansion.
    std::complex<double> zterm = pow(m_kappa * z, m_ntube);
    std::complex<double> wdsum = 0.;
    std::complex<double> wsum = cc1[m_ntube - 3][nterm];
    for (int i = nterm; i--;) {
      wdsum = wsum + zterm * wdsum;
      wsum = cc1[m_ntube - 3][i] + zterm * wsum;
    }
    // Return the results.
    ww = m_kappa * z * wsum;
    wd = m_kappa * (wsum + double(m_ntube) * zterm * wdsum);
  } else {
    // Z is close to the edge.
    // First rotate Z nearest to 1.
    double arot = -TwoPi * round(std::arg(z) * m_ntube / TwoPi) / m_ntube;
    const std::complex<double> zz =
        z * std::complex<double>(cos(arot), sin(arot));
    // Expand in a series.
    std::complex<double> zterm =
        pow(m_kappa * (1. - zz), m_ntube / (m_ntube - 2.));
    std::complex<double> wdsum = 0.;
    std::complex<double> wsum = cc2[m_ntube - 3][nterm];
    for (int i = nterm; i--;) {
      wdsum = wsum + zterm * wdsum;
      wsum = cc2[m_ntube - 3][i] + zterm * wsum;
    }
    // And return the results.
    ww = std::complex<double>(cos(arot), -sin(arot)) * (1. - zterm * wsum);
    wd = m_ntube * m_kappa * pow(m_kappa * (1. - zz), 2. / (m_ntube - 2.)) *
         (wsum + zterm * wdsum) / (m_ntube - 2.);
  }
}

void ComponentAnalyticField::E2Sum(const double xpos, const double ypos,
                                   double& ex, double& ey) const {
  //-----------------------------------------------------------------------
  //   E2SUM  - Components of the elecrostatic field intensity in a doubly
  //            periodic wire array.
  //            Clenshaw's algorithm is used for the evaluation of the sums
  //                  ZTERM1 = SIN(ZETA) - P1*SIN(3*ZETA) + P2*SIN(5*ZETA),
  //                  ZTERM2 = COS(ZETA)- 3 P1*COS(3*ZETA)+ 5P2*COS(5*ZETA)
  //   VARIABLES : (XPOS,YPOS): Position in the basic cell at which the
  //                            field is to be computed.
  //  (Essentially by G.A.Erskine/DD, 14.8.1984)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);

  std::complex<double> wsum = 0.;
  for (const auto& wire : m_w) {
    const auto zeta =
        m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum += wire.e * (zterm.second / zterm.first);
    }
  }
  ex = -real(-m_zmult * wsum);
  ey = imag(-m_zmult * wsum);
}

void ComponentAnalyticField::FieldA00(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCA00 - Subroutine performing the actual field calculations in case
  //            only one charge and not more than 1 mirror-charge in either
  //            x or y is present.
  //            The potential used is 1/2*pi*eps0  log(r).
  //   VARIABLES : R2         : Potential before taking -log(sqrt(...))
  //               EX, EY     : x,y-component of the electric field.
  //               ETOT       : Magnitude of electric field.
  //               VOLT       : Potential.
  //               EXHELP etc : One term in the series to be summed.
  //               (XPOS,YPOS): The position where the field is calculated.
  //   (Last changed on 25/ 1/96.)
  //-----------------------------------------------------------------------

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  double xxmirr = 0., yymirr = 0.;
  // Loop over all wires.
  for (const auto& wire : m_w) {
    const double xx = xpos - wire.x;
    const double yy = ypos - wire.y;
    double r2 = xx * xx + yy * yy;
    // Calculate the field in case there are no planes.
    double exhelp = xx / r2;
    double eyhelp = yy / r2;
    // Take care of a plane at constant x.
    if (m_ynplax) {
      xxmirr = wire.x + (xpos - 2. * m_coplax);
      const double r2plan = xxmirr * xxmirr + yy * yy;
      exhelp -= xxmirr / r2plan;
      eyhelp -= yy / r2plan;
      r2 /= r2plan;
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      yymirr = wire.y + (ypos - 2. * m_coplay);
      const double r2plan = xx * xx + yymirr * yymirr;
      exhelp -= xx / r2plan;
      eyhelp -= yymirr / r2plan;
      r2 /= r2plan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double r2plan = xxmirr * xxmirr + yymirr * yymirr;
      exhelp += xxmirr / r2plan;
      eyhelp += yymirr / r2plan;
      r2 *= r2plan;
    }
    // Calculate the electric field and potential.
    if (opt) volt -= 0.5 * wire.e * log(r2);
    ex += wire.e * exhelp;
    ey += wire.e * eyhelp;
  }
}

void ComponentAnalyticField::FieldB1X(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCB1X - Routine calculating the potential for a row of positive
  //            charges. The potential used is Re(Log(sin pi/s (z-z0))).
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z,ZZMIRR   : X + I*Y , XXMIRR + I*YYMIRR ; I**2=-1
  //               ECOMPL     : EX + I*EY                   ; I**2=-1
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);

  std::complex<double> ecompl;

  double r2 = 0.;

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  const double tx = Pi / m_sx;
  // Loop over all wires.
  for (const auto& wire : m_w) {
    const double xx = tx * (xpos - wire.x);
    const double yy = tx * (ypos - wire.y);
    // Calculate the field in case there are no equipotential planes.
    if (yy > 20.) {
      ecompl = -icons;
    } else if (yy < -20.) {
      ecompl = icons;
    } else {
      const std::complex<double> zz(xx, yy);
      const auto expzz = exp(2. * icons * zz);
      ecompl = icons * (expzz + 1.) / (expzz - 1.);
    }

    if (opt) {
      if (fabs(yy) > 20.) r2 = -fabs(yy) + CLog2;
      if (fabs(yy) <= 20.) r2 = -0.5 * log(pow(sinh(yy), 2) + pow(sin(xx), 2));
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + wire.y - 2. * m_coplay);
      if (yymirr > 20.) {
        ecompl += icons;
      } else if (yymirr < -20.) {
        ecompl += -icons;
      } else {
        const std::complex<double> zzmirr(xx, yymirr);
        const auto expzzmirr = exp(2. * icons * zzmirr);
        ecompl += -icons * (expzzmirr + 1.) / (expzzmirr - 1.);
      }
      if (opt && fabs(yymirr) > 20.) r2 += fabs(yymirr) - CLog2;
      if (opt && fabs(yymirr) <= 20.)
        r2 += 0.5 * log(pow(sinh(yymirr), 2) + pow(sin(xx), 2));
    }
    // Calculate the electric field and potential.
    ex += wire.e * real(ecompl);
    ey -= wire.e * imag(ecompl);
    if (opt) volt += wire.e * r2;
  }
  ex *= tx;
  ey *= tx;
}

void ComponentAnalyticField::FieldB1Y(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCB1Y - Routine calculating the potential for a row of positive
  //            charges. The potential used is Re(Log(sinh pi/sy(z-z0)).
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z,ZZMIRR   : X + I*Y , XXMIRR + I*YYMIRR ; I**2=-1
  //               ECOMPL     : EX + I*EY                   ; I**2=-1
  //-----------------------------------------------------------------------

  std::complex<double> ecompl;

  double r2 = 0.;

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  const double ty = Pi / m_sy;
  // Loop over all wires.
  for (const auto& wire : m_w) {
    const double xx = ty * (xpos - wire.x);
    const double yy = ty * (ypos - wire.y);
    // Calculate the field in case there are no equipotential planes.
    if (xx > 20.) {
      ecompl = 1.;
    } else if (xx < -20.) {
      ecompl = -1.;
    } else {
      const std::complex<double> zz(xx, yy);
      const auto expzz = exp(2. * zz);
      ecompl = (expzz + 1.) / (expzz - 1.);
    }
    if (opt) {
      if (fabs(xx) > 20.) r2 = -fabs(xx) + CLog2;
      if (fabs(xx) <= 20.) r2 = -0.5 * log(pow(sinh(xx), 2) + pow(sin(yy), 2));
    }
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + wire.x - 2. * m_coplax);
      if (xxmirr > 20.) {
        ecompl -= 1.;
      } else if (xxmirr < -20.) {
        ecompl += 1.;
      } else {
        const std::complex<double> zzmirr(xxmirr, yy);
        const auto expzzmirr = exp(2. * zzmirr);
        ecompl -= (expzzmirr + 1.) / (expzzmirr - 1.);
      }
      if (opt && fabs(xxmirr) > 20.) r2 += fabs(xxmirr) - CLog2;
      if (opt && fabs(xxmirr) <= 20.)
        r2 += 0.5 * log(pow(sinh(xxmirr), 2) + pow(sin(yy), 2));
    }
    // Calculate the electric field and potential.
    ex += wire.e * real(ecompl);
    ey -= wire.e * imag(ecompl);
    if (opt) volt += wire.e * r2;
  }
  ex *= ty;
  ey *= ty;
}

void ComponentAnalyticField::FieldB2X(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCB2X - Routine calculating the potential for a row of alternating
  //            + - charges. The potential used is re log(sin pi/sx (z-z0))
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z, ZZMRR   : X + i*Y , XXMIRR + i*YYMIRR ; i**2=-1
  //               ECOMPL     : EX + i*EY                   ; i**2=-1
  //   (Cray vectorisable)
  //-----------------------------------------------------------------------

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  const double tx = HalfPi / m_sx;
  // Loop over all wires.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const double xx = tx * (xpos - m_w[i].x);
    const double yy = tx * (ypos - m_w[i].y);
    const double xxneg = tx * (xpos + m_w[i].x - 2. * m_coplax);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    double r2 = 1.;
    if (fabs(yy) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xxneg, yy);
      ecompl = -m_b2sin[i] / (sin(zz) * sin(zzneg));
      if (opt) {
        const double sinhy = sinh(yy);
        const double sinxx = sin(xx);
        const double sinxxneg = sin(xxneg);
        r2 = (sinhy * sinhy + sinxx * sinxx) /
             (sinhy * sinhy + sinxxneg * sinxxneg);
      }
    }
    // Take care of a planes at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + m_w[i].y - 2 * m_coplay);
      if (fabs(yymirr) <= 20.) {
        const std::complex<double> zzmirr(xx, yymirr);
        const std::complex<double> zznmirr(xxneg, yymirr);
        ecompl += m_b2sin[i] / (sin(zzmirr) * sin(zznmirr));
        if (opt) {
          const double sinhy = sinh(yymirr);
          const double sinxx = sin(xx);
          const double sinxxneg = sin(xxneg);
          const double r2plan = (sinhy * sinhy + sinxx * sinxx) /
                                (sinhy * sinhy + sinxxneg * sinxxneg);
          r2 /= r2plan;
        }
      }
    }
    // Calculate the electric field and potential.
    ex += m_w[i].e * real(ecompl);
    ey -= m_w[i].e * imag(ecompl);
    if (opt) volt -= 0.5 * m_w[i].e * log(r2);
  }
  ex *= tx;
  ey *= tx;
}

void ComponentAnalyticField::FieldB2Y(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCB2Y - Routine calculating the potential for a row of alternating
  //            + - charges. The potential used is re log(sin pi/sx (z-z0))
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z, ZMIRR   : X + i*Y , XXMIRR + i*YYMIRR ; i**2=-1
  //               ECOMPL     : EX + i*EY                   ; i**2=-1
  //   (Cray vectorisable)
  //-----------------------------------------------------------------------

  const std::complex<double> icons(0., 1.);

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  const double ty = HalfPi / m_sy;
  // Loop over all wires.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const double xx = ty * (xpos - m_w[i].x);
    const double yy = ty * (ypos - m_w[i].y);
    const double yyneg = ty * (ypos + m_w[i].y - 2 * m_coplay);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    double r2 = 1.;
    if (fabs(xx) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xx, yyneg);
      ecompl = icons * m_b2sin[i] / (sin(icons * zz) * sin(icons * zzneg));
      if (opt) {
        const double sinhx = sinh(xx);
        const double sinyy = sin(yy);
        const double sinyyneg = sin(yyneg);
        r2 = (sinhx * sinhx + sinyy * sinyy) /
             (sinhx * sinhx + sinyyneg * sinyyneg);
      }
    }
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + m_w[i].x - 2 * m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const std::complex<double> zzmirr(xxmirr, yy);
        const std::complex<double> zznmirr(xxmirr, yyneg);
        ecompl -=
            icons * m_b2sin[i] / (sin(icons * zzmirr) * sin(icons * zznmirr));
        if (opt) {
          const double sinhx = sinh(xxmirr);
          const double sinyy = sin(yy);
          const double sinyyneg = sin(yyneg);
          const double r2plan = (sinhx * sinhx + sinyy * sinyy) /
                                (sinhx * sinhx + sinyyneg * sinyyneg);
          r2 /= r2plan;
        }
      }
    }
    // Calculate the electric field and potential.
    ex += m_w[i].e * real(ecompl);
    ey -= m_w[i].e * imag(ecompl);
    if (opt) volt -= 0.5 * m_w[i].e * log(r2);
  }
  ex *= ty;
  ey *= ty;
}

void ComponentAnalyticField::FieldC10(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCC10 - Routine returning the potential and electric field. It
  //            calls the routines PH2 and E2SUM written by G.A.Erskine.
  //   VARIABLES : No local variables.
  //-----------------------------------------------------------------------

  // Calculate voltage first, if needed.
  if (opt) {
    if (m_mode == 0) volt = m_v0 + m_c1 * xpos;
    if (m_mode == 1) volt = m_v0 + m_c1 * ypos;
    for (const auto& wire : m_w) {
      volt += wire.e * Ph2(xpos - wire.x, ypos - wire.y);
    }
  }

  // And finally the electric field.
  E2Sum(xpos, ypos, ex, ey);
  if (m_mode == 0) ex -= m_c1;
  if (m_mode == 1) ey -= m_c1;
}

void ComponentAnalyticField::FieldC2X(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCC2X - Routine returning the potential and electric field in a
  //            configuration with 2 x planes and y periodicity.
  //   VARIABLES : see the writeup
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);

  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  volt = 0.;

  // Wire loop.
  for (const auto& wire : m_w) {
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum1 -= wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum1 += wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += wire.e * (zterm.second / zterm.first);
      if (opt) volt -= wire.e * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - wire.x) / m_sx);
    // Mirror contribution.
    zeta =
        m_zmult * std::complex<double>(2. * cx - xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum2 -= wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum2 += wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += wire.e * (zterm.second / zterm.first);
      if (opt) volt += wire.e * log(std::abs(zterm.first));
    }
    // Correct the voltage, if needed (MODE).
    if (opt && m_mode == 0) {
      volt -= TwoPi * wire.e * (xpos - cx) * (wire.x - cx) / (m_sx * m_sy);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2));
  ey = -imag(m_zmult * (wsum1 - wsum2));
  // Constant correction terms.
  if (m_mode == 0) ex -= m_c1;
}

void ComponentAnalyticField::FieldC2Y(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCC2Y - Routine returning the potential and electric field in a
  //            configuration with 2 y planes and x periodicity.
  //   VARIABLES : see the writeup
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);

  // Initial values.
  volt = 0.;
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;

  // Wire loop.
  for (const auto& wire : m_w) {
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum1 -= wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum1 += wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += wire.e * (zterm.second / zterm.first);
      if (opt) volt -= wire.e * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - wire.y) / m_sy);
    // Mirror contribution from the y plane.
    zeta =
        m_zmult * std::complex<double>(xpos - wire.x, 2 * cy - ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum2 -= wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum2 += wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += wire.e * (zterm.second / zterm.first);
      if (opt) volt += wire.e * log(std::abs(zterm.first));
    }
    // Correct the voltage, if needed (MODE).
    if (opt && m_mode == 1) {
      volt -= TwoPi * wire.e * (ypos - cy) * (wire.y - cy) / (m_sx * m_sy);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 - wsum2));
  ey = -imag(m_zmult * (wsum1 + wsum2));
  // Constant correction terms.
  if (m_mode == 1) ey -= m_c1;
}

void ComponentAnalyticField::FieldC30(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCC30 - Routine returning the potential and electric field in a
  //            configuration with 2 y and 2 x planes.
  //   VARIABLES : see the writeup
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);

  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  std::complex<double> wsum3 = 0.;
  std::complex<double> wsum4 = 0.;
  volt = 0.;

  // Wire loop.
  for (const auto& wire : m_w) {
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum1 -= wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum1 += wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += wire.e * (zterm.second / zterm.first);
      if (opt) volt -= wire.e * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - wire.x) / m_sx);
    // Mirror contribution from the x plane.
    zeta =
        m_zmult * std::complex<double>(2. * cx - xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum2 -= wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum2 += wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += wire.e * (zterm.second / zterm.first);
      if (opt) volt += wire.e * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - wire.y) / m_sy);
    // Mirror contribution from the x plane.
    zeta =
        m_zmult * std::complex<double>(xpos - wire.x, 2. * cy - ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum3 -= wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum3 += wire.e * icons;
      if (opt) volt += wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum3 += wire.e * (zterm.second / zterm.first);
      if (opt) volt += wire.e * log(std::abs(zterm.first));
    }
    // Mirror contribution from both the x and the y plane.
    zeta = m_zmult * std::complex<double>(2. * cx - xpos - wire.x,
                                          2. * cy - ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum4 -= wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else if (imag(zeta) < -15.) {
      wsum4 += wire.e * icons;
      if (opt) volt -= wire.e * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum4 += wire.e * (zterm.second / zterm.first);
      if (opt) volt -= wire.e * log(std::abs(zterm.first));
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2 - wsum3 - wsum4));
  ey = -imag(m_zmult * (wsum1 - wsum2 + wsum3 - wsum4));
}

void ComponentAnalyticField::FieldD10(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCD10 - Subroutine performing the actual field calculations for a
  //            cell which has one circular plane and some wires.
  //   VARIABLES : EX, EY, VOLT:Electric field and potential.
  //               ETOT, VOLT : Magnitude of electric field, potential.
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on  4/ 9/95.)
  //-----------------------------------------------------------------------

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (const auto& wire : m_w) {
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(wire.x, wire.y);
    // Compute the contribution to the potential, if needed.
    if (opt) {
      volt -= wire.e * log(
        std::abs(m_cotube * (zpos - zi) / (m_cotube2 - zpos * conj(zi))));
    }
    // Compute the contribution to the electric field, always.
    const auto wi = 1. / conj(zpos - zi) + zi / (m_cotube2 - conj(zpos) * zi);
    ex += wire.e * real(wi);
    ey += wire.e * imag(wi);
  }
}

void ComponentAnalyticField::FieldD20(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCD20 - Subroutine performing the actual field calculations for a
  //            cell which has a tube and phi periodicity.
  //   VARIABLES : EX, EY, VOLT:Electric field and potential.
  //               ETOT, VOLT : Magnitude of electric field, potential.
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on 10/ 2/93.)
  //-----------------------------------------------------------------------

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (const auto& wire : m_w) {
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(wire.x, wire.y);
    // Case of the wire which is not in the centre.
    if (std::abs(zi) > wire.r) {
      // Compute the contribution to the potential, if needed.
      if (opt) {
        volt -= wire.e * log(
          std::abs((1. / pow(m_cotube, m_mtube)) *
                   (pow(zpos, m_mtube) - pow(zi, m_mtube)) /
                   (1. - pow(zpos * conj(zi) / m_cotube2, m_mtube))));
      }
      // Compute the contribution to the electric field, always.
      const auto wi = double(m_mtube) * pow(conj(zpos), m_mtube - 1) *
          (1. / conj(pow(zpos, m_mtube) - pow(zi, m_mtube)) +
           pow(zi, m_mtube) /
               (pow(m_cotube, 2 * m_mtube) - pow(conj(zpos) * zi, m_mtube)));
      ex += wire.e * real(wi);
      ey += wire.e * imag(wi);
    } else {
      // Case of the central wire.
      if (opt) {
        volt -= wire.e * log(
          std::abs((1. / m_cotube) * (zpos - zi) /
                   (1. - zpos * conj(zi) / m_cotube2)));
      }
      const auto wi = 1. / conj(zpos - zi) + zi / (m_cotube2 - conj(zpos) * zi);
      // Compute the contribution to the electric field, always.
      ex += wire.e * real(wi);
      ey += wire.e * imag(wi);
    }
  }
}

void ComponentAnalyticField::FieldD30(const double xpos, const double ypos,
                                      double& ex, double& ey, double& volt,
                                      const bool opt) const {
  //-----------------------------------------------------------------------
  //   EFCD30 - Subroutine performing the actual field calculations for a
  //            cell which has a polygon as tube and some wires.
  //   VARIABLES : EX, EY, VOLT:Electric field and potential.
  //               ETOT, VOLT : Magnitude of electric field, potential.
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on 19/ 2/94.)
  //-----------------------------------------------------------------------

  // Initialise the electric field and potential.
  ex = ey = 0.;
  volt = m_v0;

  // Get the mapping of the position.
  std::complex<double> wpos, wdpos;
  ConformalMap(std::complex<double>(xpos, ypos) / m_cotube, wpos, wdpos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Compute the contribution to the potential, if needed.
    if (opt) {
      volt -= m_w[i].e * log(
        std::abs((wpos - m_zw[i]) / (1. - wpos * conj(m_zw[i]))));
    }
    const auto whelp = wdpos * (1. - std::norm(m_zw[i])) /
            ((wpos - m_zw[i]) * (1. - conj(m_zw[i]) * wpos));
    // Compute the contribution to the electric field, always.
    ex += m_w[i].e * real(whelp);
    ey -= m_w[i].e * imag(whelp);
  }
  ex /= m_cotube;
  ey /= m_cotube;
}

bool ComponentAnalyticField::InTube(const double x0, const double y0,
                                    const double a, const int n) {
  //-----------------------------------------------------------------------
  //   INTUBE - Determines whether a point is located inside a polygon.
  //            ILOC is set to +1 if outside, 0 if inside and -1 if the
  //            arguments are not valid.
  //   (Last changed on 18/ 3/01.)
  //-----------------------------------------------------------------------

  // Special case: x = y = 0
  if (x0 == 0. && y0 == 0.) return true;

  // Special case: round tube.
  if (n == 0) {
    if (x0 * x0 + y0 * y0 > a * a) return false;
    return true;
  }

  if (n < 0 || n == 1 || n == 2) {
    std::cerr << "ComponentAnalyticField::InTube:\n"
              << "    Invalid number of edges (n = " << n << ")\n";
    return false;
  }

  // Truly polygonal tubes.
  // Reduce angle to the first sector.
  double phi = atan2(y0, x0);
  if (phi < 0.) phi += TwoPi;
  phi -= TwoPi * int(0.5 * n * phi / Pi) / n;
  // Compare the length to the local radius.
  if ((x0 * x0 + y0 * y0) * pow(cos(Pi / n - phi), 2) >
      a * a * pow(cos(Pi / n), 2))
    return false;

  return true;
}

void ComponentAnalyticField::Field3dA00(
    const double xpos, const double ypos, const double zpos, 
    double& ex, double& ey, double& ez, double& volt) const {
  //-----------------------------------------------------------------------
  //   E3DA00 - Subroutine adding 3-dimensional charges for A cells.
  //            The potential used is 1/2*pi*eps0  1/r
  //   VARIABLES : EX, EY     : x,y-component of the electric field.
  //               ETOT       : Magnitude of electric field.
  //               VOLT       : Potential.
  //               EXHELP etc : One term in the series to be summed.
  //               (XPOS,YPOS): The position where the field is calculated.
  //   (Last changed on  5/12/94.)
  //-----------------------------------------------------------------------
  // Initialise the electric field and potential.
  ex = ey = ez = volt = 0.;

  // Loop over all charges.
  for (const auto& charge : m_ch3d) {
    // Calculate the field in case there are no planes.
    const double dx = xpos - charge.x;
    const double dy = ypos - charge.y;
    const double dz = zpos - charge.z;
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (fabs(r) < Small) continue;
    const double f = pow(r, -3);
    double exhelp = -dx * f;
    double eyhelp = -dy * f;
    double ezhelp = -dz * f;
    double vhelp = 1. / r;
    // Take care of a plane at constant x.
    double dxm = 0., dym = 0.;
    if (m_ynplax) {
      dxm = charge.x + xpos - 2 * m_coplax;
      const double rplan = sqrt(dxm * dxm + dy * dy);
      if (fabs(rplan) < Small) continue;
      const double fplan = pow(rplan, -3);
      exhelp += dxm * fplan;
      eyhelp += dy * fplan;
      ezhelp += dz * fplan;
      vhelp -= 1. / rplan;
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      dym = charge.y + ypos - 2. * m_coplay;
      const double rplan = sqrt(dx * dx + dym * dym);
      if (fabs(rplan) < Small) continue;
      const double fplan = pow(rplan, -3);
      exhelp += dx * fplan;
      eyhelp += dym * fplan;
      ezhelp += dz * fplan;
      vhelp -= 1. / rplan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double rplan = sqrt(dxm * dxm + dym * dym);
      if (fabs(rplan) < Small) continue;
      const double fplan = pow(rplan, -3);
      exhelp -= dxm * fplan;
      eyhelp -= dym * fplan;
      ezhelp -= dz * fplan;
      vhelp += 1. / rplan;
    }
    // Add the terms to the electric field and the potential.
    ex -= charge.e * exhelp;
    ey -= charge.e * eyhelp;
    ez -= charge.e * ezhelp;
    volt += charge.e * vhelp;
  }
}

void ComponentAnalyticField::Field3dB2X(
    const double xpos, const double ypos, const double zpos, 
    double& ex, double& ey, double& ez, double& volt) const {
  //-----------------------------------------------------------------------
  //   E3DB2X - Routine calculating the potential for a 3 dimensional point
  //            charge between two plates at constant x.
  //            The series expansions for the modified Bessel functions
  //            have been taken from Abramowitz and Stegun.
  //   VARIABLES : See routine E3DA00 for most of the variables.
  //   (Last changed on  5/12/94.)
  //-----------------------------------------------------------------------

  constexpr double rcut = 1.;

  ex = ey = ez = volt = 0.;

  // Shorthand.
  const double a = Pi / m_sx;

  // Loop over all charges.
  for (const auto& charge : m_ch3d) {
    // Skip coordinates that are on the charge.
    if (xpos == charge.x && ypos == charge.y && zpos == charge.z) continue;
    const double dx = xpos - charge.x;
    const double dy = ypos - charge.y;
    const double dz = zpos - charge.z;
    const double dxm = xpos + charge.x - 2 * m_coplax;
    double rho2 = dy * dy + dz * dz;
    double rho = sqrt(rho2);
    double ersum = 0., exsum = 0., vsum = 0.;
    // In the far away zone, sum the modified Bessel function series.
    if (rho > rcut * 2 * m_sx) {
      // Loop over the terms in the series.
      for (unsigned int j = 1; j <= m_nTermBessel; ++j) {
        // Obtain reduced coordinates.
        const double rr = a * j * rho;
        const double zp = a * j * dx;
        const double zn = a * j * dxm;
        // Evaluate the Bessel functions for this R.
        const double k0 = rr < 2. ? Numerics::BesselK0S(rr) :
                                    Numerics::BesselK0L(rr);
        const double k1 = rr < 2. ? Numerics::BesselK1S(rr) :
                                    Numerics::BesselK1L(rr);
        // Get the field components.
        const double cp = cos(zp);
        const double cn = cos(zn);
        vsum += k0 * (cp - cn);
        ersum += j * k1 * (cp - cn);
        exsum += j * k0 * (sin(zp) - sin(zn));
      }
      vsum /= m_sx;
      exsum *= TwoPi / (m_sx * m_sx);
      ersum *= TwoPi / (m_sx * m_sx);
      ersum /= rho;
    } else {
      // Direct polynomial summing, obtain reduced coordinates.
      // Loop over the terms.
      for (unsigned int j = 0; j <= m_nTermPoly; ++j) {
        // Simplify the references to the distances.
        const double sx = j * 2 * m_sx;
        const double rr1 = sqrt(rho2 + pow(dx + sx, 2));
        const double rm1 = sqrt(rho2 + pow(dxm - sx, 2));
        const double fr1 = pow(rr1, -3);
        const double fm1 = pow(rm1, -3);
        // Initialisation of the sum: only a charge and a mirror charge.
        if (j == 0) {
          vsum = 1. / rr1 - 1. / rm1;
          exsum = dx * fr1 - dxm * fm1;
          ersum = fr1 - fm1;
          continue;
        }
        const double rr2 = sqrt(rho2 + pow(dx - sx, 2));
        const double rm2 = sqrt(rho2 + pow(dxm + sx, 2));
        const double fr2 = pow(rr2, -3);
        const double fm2 = pow(rm2, -3);
        // Further terms in the series: 2 charges and 2 mirror charges.
        vsum += 1. / rr1 + 1. / rr2 - 1. / rm1 - 1. / rm2;
        exsum += (dx + sx) * fr1 - (dxm - sx) * fm1 + 
                 (dx - sx) * fr2 - (dxm + sx) * fm2;
        ersum += fr1 + fr2 - fm1 - fm2;
      }
    }
    ex += charge.e * exsum;
    ey += charge.e * ersum * dy;
    ez += charge.e * ersum * dz;
    volt += charge.e * vsum;
    // Take care of a plane at constant y.
    if (!m_ynplay) continue;
    ersum = exsum = vsum = 0.;
    const double dym = ypos + charge.y - 2. * m_coplay;
    rho2 = dym * dym + dz * dz;
    rho = sqrt(rho2);
    if (rho > rcut * 2 * m_sx) {
      // Bessel function series.
      // Loop over the terms in the series.
      for (unsigned int j = 1; j <= m_nTermBessel; ++j) {
        // Obtain reduced coordinates.
        const double rr = a * j * rho;
        const double zp = a * j * dx;
        const double zn = a * j * dxm;
        // Evaluate the Bessel functions for this R.
        const double k0 = rr < 2. ? Numerics::BesselK0S(rr) :
                                    Numerics::BesselK0L(rr);
        const double k1 = rr < 2. ? Numerics::BesselK1S(rr) :
                                    Numerics::BesselK1L(rr);
        // Get the field components.
        const double cp = cos(zp);
        const double cn = cos(zn);
        vsum += k0 * (cp - cn);
        ersum += k1 * (cp - cn);
        exsum += k0 * (sin(zp) - sin(zn));
      }
      vsum /= m_sx;
      exsum *= TwoPi / (m_sx * m_sx);
      ersum *= TwoPi / (m_sx * m_sx);
      ersum /= rho;
    } else {
      // Polynomial sum.
      // Loop over the terms.
      for (unsigned int j = 0; j <= m_nTermPoly; ++j) {
        // Simplify the references to the distances.
        const double sx = j * 2 * m_sx;
        const double rr1 = sqrt(rho2 + pow(dx + sx, 2));
        const double rm1 = sqrt(rho2 + pow(dxm - sx, 2));
        const double fr1 = pow(rr1, -3);
        const double fm1 = pow(rm1, -3);
        // Initialisation of the sum: only a charge and a mirror charge.
        if (j == 0) {
          vsum -= 1. / rr1 - 1. / rm1;
          exsum -= dx * fr1 - dxm * fm1;
          ersum -= fr1 - fm1;
          continue;
        }
        const double rr2 = sqrt(rho2 + pow(dx - sx, 2));
        const double rm2 = sqrt(rho2 + pow(dxm + sx, 2));
        const double fr2 = pow(rr2, -3);
        const double fm2 = pow(rm2, -3);
        // Further terms in the series: 2 charges and 2 mirror charges.
        vsum -= 1. / rr1 + 1. / rr2 - 1. / rm1 - 1. / rm2;
        exsum -= (dx + sx) * fr1 - (dxm - sx) * fm1 + 
                 (dx - sx) * fr2 - (dxm + sx) * fm2;
        ersum -= fr1 + fr2 - fm1 - fm2;
      }
    }
    ex += charge.e * exsum;
    ey += charge.e * ersum * dym;
    ez += charge.e * ersum * dz;
    volt += charge.e * vsum;
  }
}

void ComponentAnalyticField::Field3dB2Y(
    const double xpos, const double ypos, const double zpos, 
    double& ex, double& ey, double& ez, double& volt) const {
  //-----------------------------------------------------------------------
  //   E3DB2Y - Routine calculating the potential for a 3 dimensional point
  //            charge between two plates at constant y.
  //            The series expansions for the modified Bessel functions
  //            have been taken from Abramowitz and Stegun.
  //   VARIABLES : See routine E3DA00 for most of the variables.
  //   (Last changed on  5/12/94.)
  //-----------------------------------------------------------------------

  constexpr double rcut = 1.;

  ex = ey = ez = volt = 0.;

  // Shorthand.
  const double a = Pi / m_sy;

  // Loop over all charges.
  for (const auto& charge : m_ch3d) {
    // Skip wires that are on the charge.
    if (xpos == charge.x && ypos == charge.y && zpos == charge.z) continue;
    const double dx = xpos - charge.x;
    const double dy = ypos - charge.y;
    const double dz = zpos - charge.z;
    const double dym = ypos + charge.y - 2 * m_coplay;
    double rho2 = dx * dx + dz * dz;
    double rho = sqrt(rho2);
    double ersum = 0., eysum = 0., vsum = 0.;
    // In the far away zone, sum the modified Bessel function series.
    if (rho > rcut * 2 * m_sy) {
      // Loop over the terms in the series.
      for (unsigned int j = 1; j <= m_nTermBessel; ++j) {
        // Obtain reduced coordinates.
        const double rr = a * j * rho;
        const double zp = a * j * dy;
        const double zn = a * j * dym;
        // Evaluate the Bessel functions for this R.
        const double k0 = rr < 2. ? Numerics::BesselK0S(rr) :
                                    Numerics::BesselK0L(rr);
        const double k1 = rr < 2. ? Numerics::BesselK1S(rr) :
                                    Numerics::BesselK1L(rr);
        // Get the field components.
        const double cp = cos(zp);
        const double cn = cos(zn);
        vsum += k0 * (cp - cn);
        ersum += j * k1 * (cp - cn);
        eysum += j * k0 * (sin(zp) - sin(zn));
      }
      vsum /= m_sy;
      eysum *= TwoPi / (m_sy * m_sy);
      ersum *= TwoPi / (m_sy * m_sy);
      ersum /= rho;
    } else {
      // Direct polynomial summing, obtain reduced coordinates.
      // Loop over the terms.
      for (unsigned int j = 0; j <= m_nTermPoly; ++j) {
        const double sy = j * 2 * m_sy;
        // Simplify the references to the distances.
        const double rr1 = sqrt(rho2 + pow(dy + sy, 2));
        const double rm1 = sqrt(rho2 + pow(dym - sy, 2));
        const double fr1 = pow(rr1, -3);
        const double fm1 = pow(rm1, -3);
        // Initialisation of the sum: only a charge and a mirror charge.
        if (j == 0) {
          vsum = 1. / rr1 - 1. / rm1;
          ersum = fr1 - fm1;
          eysum = dy * fr1 - dym * fm1;
          continue;
        }
        // Further terms in the series: 2 charges and 2 mirror charges.
        const double rr2 = sqrt(rho2 + pow(dy - sy, 2));
        const double rm2 = sqrt(rho2 + pow(dym + sy, 2));
        const double fr2 = pow(rr2, -3);
        const double fm2 = pow(rm2, -3);
        vsum += 1. / rr1 + 1. / rr2 - 1. / rm1 - 1. / rm2;
        ersum += fr1 + fr2 - fm1 - fm2;
        eysum += (dy + sy) * fr1 - (dym - sy) * fm1 + 
                 (dy - sy) * fr2 - (dym + sy) * fm2;
      }
    }
    ex += charge.e * ersum * dx;
    ey += charge.e * eysum;
    ez += charge.e * ersum * dz;
    volt += charge.e * vsum;
    // Take care of a plane at constant x.
    if (!m_ynplax) continue;
    ersum = eysum = vsum = 0.;
    const double dxm = xpos + charge.x - 2. * m_coplax;
    rho2 = dxm * dxm + dz * dz;
    rho = sqrt(rho2); 
    if (rho > rcut * 2 * m_sy) {
      // Bessel function series.
      // Loop over the terms in the series.
      for (unsigned int j = 1; j <= m_nTermBessel; ++j) {
        // Obtain reduced coordinates.
        const double rr = a * j * rho;
        const double zp = a * j * dy;
        const double zn = a * j * dym;
        // Evaluate the Bessel functions for this R.
        const double k0 = rr < 2. ? Numerics::BesselK0S(rr) :
                                    Numerics::BesselK0L(rr);
        const double k1 = rr < 2. ? Numerics::BesselK1S(rr) :
                                    Numerics::BesselK1L(rr);
        // Get the field components.
        const double cp = cos(zp);
        const double cn = cos(zn);
        vsum += k0 * (cp - cn);
        ersum += k1 * (cp - cn);
        eysum += k0 * (sin(zp) - sin(zn));
      }
      vsum /= m_sy;
      eysum *= TwoPi / (m_sy * m_sy);
      ersum *= TwoPi / (m_sy * m_sy);
      ersum /= rho;
    } else {
      // Polynomial sum.
      // Loop over the terms.
      for (unsigned int j = 0; j <= m_nTermPoly; ++j) {
        // Simplify the references to the distances.
        const double sy = j * 2 * m_sy;
        const double rr1 = sqrt(rho2 + pow(dy + sy, 2));
        const double rm1 = sqrt(rho2 + pow(dym - sy, 2));
        const double fr1 = pow(rr1, -3);
        const double fm1 = pow(rm1, -3);
        // Initialisation of the sum: only a charge and a mirror charge.
        if (j == 0) {
          vsum -= 1. / rr1 - 1. / rm1;
          ersum -= fr1 - fm1;
          eysum -= dy * fr1 - dym * fm1;
          continue;
        }
        const double rr2 = sqrt(rho2 + pow(dy - sy, 2));
        const double rm2 = sqrt(rho2 + pow(dym + sy, 2));
        const double fr2 = pow(rr2, -3);
        const double fm2 = pow(rm2, -3);
        // Further terms in the series: 2 charges and 2 mirror charges.
        vsum -= 1. / rr1 + 1. / rr2 - 1. / rm1 - 1. / rm2;
        ersum -= fr1 + fr2 - fm1 - fm2;
        eysum -= (dy + sy) * fr1 - (dym - sy) * fm1 + 
                 (dy - sy) * fr2 - (dym + sy) * fm2;
      }
    }
    ex += charge.e * ersum * dxm;
    ey += charge.e * eysum;
    ez += charge.e * ersum * dz;
    volt += charge.e * vsum;
  }
}

void ComponentAnalyticField::Field3dD10(const double xin, const double yin,
                                        const double zin, double& fx,
                                        double& fy, double& fz,
                                        double& volt) const {
  //-----------------------------------------------------------------------
  //   E3DD10 - Subroutine adding 3-dimensional charges to tubes with one
  //            wire running down the centre.
  //            The series expansions for the modified Bessel functions
  //            have been taken from Abramowitz and Stegun.
  //   VARIABLES : See routine E3DA00 for most of the variables.
  //   (Last changed on 25/11/95.)
  //-----------------------------------------------------------------------

  constexpr double rcut = 1.;

  // Initialise the sums for the field components.
  fx = fy = fz = volt = 0.;
  double ex = 0., ey = 0., ez = 0.;

  // Ensure that the routine can actually work.
  if (m_nWires < 1) {
    std::cerr << m_className << "::Field3dD10: No wires present.\n";
    return;
  }

  // Define a periodicity and one plane in the mapped frame.
  const double ssx = log(m_cotube / m_w[0].r);
  const double a = Pi / ssx;
  const double cpl = log(m_w[0].r);

  // Transform the coordinates to the mapped frame.
  const double xpos = 0.5 * log(xin * xin + yin * yin);
  const double ypos = atan2(yin, xin);
  const double zpos = zin;

  // Loop over all point charges.
  for (const auto& charge : m_ch3d) {
    const double x3d = 0.5 * log(charge.x * charge.x + charge.y * charge.y);
    const double z3d = charge.z;
    for (int ii = -1; ii <= 1; ++ii) {
      const double y3d = atan2(charge.y, charge.x + ii * TwoPi);
      const double dx = xpos - x3d;
      const double dy = ypos - y3d;
      const double dz = zpos - z3d;
      const double dxm = xpos + x3d - 2 * cpl;
      // Skip wires that are on the charge.
      if (xpos == x3d && ypos == y3d && zpos == z3d) continue;
      const double rho2 = dy * dy + dz * dz;
      const double rho = sqrt(rho2);
      double exsum = 0., ersum = 0., vsum = 0.;
      // In the far away zone, sum the modified Bessel function series.
      if (rho > rcut * 2 * ssx) {
        // Loop over the terms in the series.
        for (unsigned int j = 1; j <= m_nTermBessel; ++j) {
          // Obtain reduced coordinates.
          const double rr = a * j * rho;
          const double zp = a * j * dx;
          const double zn = a * j * dxm;
          // Evaluate the Bessel functions for this R.
          const double k0 = rr < 2. ? Numerics::BesselK0S(rr) :
                                      Numerics::BesselK0L(rr);
          const double k1 = rr < 2. ? Numerics::BesselK1S(rr) :
                                      Numerics::BesselK1L(rr);
          // Get the field components.
          const double cp = cos(zp);
          const double cn = cos(zn);
          vsum += k0 * (cp - cn);
          ersum += j * k1 * (cp - cn);
          exsum += j * k0 * (sin(zp) - sin(zn));
        }
        vsum /= ssx;
        exsum *= TwoPi / (ssx * ssx);
        ersum *= TwoPi / (ssx * ssx);
        ersum /= rho; 
      } else {
        // Direct polynomial summing, obtain reduced coordinates.
        // Loop over the terms.
        for (unsigned int j = 0; j < m_nTermPoly; ++j) {
          // Simplify the references to the distances.
          const double sx = j * 2 * ssx;
          const double rr1 = sqrt(rho2 + pow(dx + sx, 2));
          const double rm1 = sqrt(rho2 + pow(dxm - sx, 2));
          const double fr1 = pow(rr1, -3);
          const double fm1 = pow(rm1, -3);
          // Initialisation of the sum: only a charge and a mirror charge.
          if (j == 0) {
            vsum = 1. / rr1 - 1. / rm1;
            exsum = dxm * fr1 - dxm * fm1;
            ersum = fr1 - fm1;
            continue;
          }
          const double rr2 = sqrt(rho2 + pow(dx - sx, 2));
          const double rm2 = sqrt(rho2 + pow(dxm + sx, 2));
          const double fr2 = pow(rr2, -3);
          const double fm2 = pow(rm2, -3);
          // Further terms in the series: 2 charges and 2 mirror charges.
          vsum += 1. / rr1 + 1. / rr2 - 1. / rm1 - 1. / rm2;
          exsum += (dx + sx) * fr1 - (dxm - sx) * fm1 + 
                   (dx - sx) * fr2 - (dxm + sx) * fm2;
          ersum += fr1 + fr2 - fm1 - fm2;
        }
      }
      ex += charge.e * exsum;
      ey += charge.e * ersum * dy;
      ez += charge.e * ersum * dz;
      volt += charge.e * vsum;
    }
  }

  // Transform the field vectors back to Cartesian coordinates.
  const double ctheta = cos(ypos);
  const double stheta = sin(ypos); 
  fx = exp(-xpos) * (ex * ctheta - ey * stheta);
  fy = exp(-ypos) * (ex * stheta + ey * ctheta);
  fz = ez;
}

size_t ComponentAnalyticField::SignalLayer(const int mx, const int my) const {

  if (m_nFourier == 0) return 0;
  const int m0 = m_nFourier / 2 - 1;
  if (m_fperx && m_fpery) {
    return (mx + m0) * m_nFourier + my + m0;
  } else if (m_fperx) {
    return mx + m0;
  } else if (m_fpery) {
    return my + m0;
  }
  return 0;
}

bool ComponentAnalyticField::PrepareSignals() {
  //-----------------------------------------------------------------------
  //   SIGINI - Initialises signal calculations.
  //   (Last changed on 11/10/06.)
  //-----------------------------------------------------------------------

  if (m_readout.empty()) {
    std::cerr << m_className << "::PrepareSignals:\n"
              << "    There are no readout groups defined.\n"
              << "    Calculation of weighting fields makes no sense.\n";
    return false;
  }

  if (!m_cellset && !Prepare()) {
    std::cerr << m_className << "::PrepareSignals: Cell not set up.\n";
    return false;
  }

  std::lock_guard<std::mutex> guard(m_mutex);

  // If using natural periodicity, copy the cell type.
  // Otherwise, eliminate true periodicities.
  if (m_nFourier == 0) {
    m_cellTypeFourier = m_cellType;
  } else if (m_cellType == A00 || m_cellType == B1X || m_cellType == B1Y ||
             m_cellType == C10) {
    m_cellTypeFourier = A00;
  } else if (m_cellType == B2X || m_cellType == C2X) {
    m_cellTypeFourier = B2X;
  } else if (m_cellType == B2Y || m_cellType == C2Y) {
    m_cellTypeFourier = B2Y;
  } else if (m_cellType == C30) {
    m_cellTypeFourier = C30;
  } else if (m_cellType == D10) {
    m_cellTypeFourier = D10;
  } else if (m_cellType == D30) {
    m_cellTypeFourier = D30;
  } else {
    // Other cases.
    std::cerr << m_className << "::PrepareSignals:\n"
              << "    No potentials available to handle cell type "
              << GetCellType(m_cellType) << "\n.";
    return false;
  }

  // Establish the directions in which convolutions occur.
  m_fperx = m_fpery = false;
  if (m_nFourier == 0) {
    m_mfexp = 0;
  } else {
    if (m_cellType == B1X || m_cellType == C10 || m_cellType == C2Y) {
      m_fperx = true;
    }
    if (m_cellType == B1Y || m_cellType == C10 || m_cellType == C2X) {
      m_fpery = true;
    }
    m_mfexp = int(0.1 + std::log2(m_nFourier));
    if (m_mfexp == 0) {
      m_fperx = m_fpery = false;
    }
  }
  // Set maximum and minimum Fourier terms.
  m_mxmin = m_mymin = m_mxmax = m_mymax = 0;
  if (m_fperx) {
    m_mxmin = std::min(0, 1 - m_nFourier / 2);
    m_mxmax = m_nFourier / 2;
  }
  if (m_fpery) {
    m_mymin = std::min(0, 1 - m_nFourier / 2);
    m_mymax = m_nFourier / 2;
  }

  // Print some debugging output if requested.
  if (m_debug) {
    std::cout << m_className << "::PrepareSignals:\n"
              << "    Cell type:           " << GetCellType(m_cellType) << "\n"
              << "    Fourier cell type:   " << GetCellType(m_cellTypeFourier)
              << "\n";
    std::cout << "    x convolutions:      " << m_fperx << "\n"
              << "    y convolutions:      " << m_fpery << "\n";
    std::cout << "    No of Fourier terms: " << m_nFourier << " (= 2**"
              << m_mfexp << ")\n";
  }

  // Prepare the signal matrices.
  if (!SetupWireSignals()) {
    std::cerr << m_className << "::PrepareSignals:\n"
              << "    Preparing wire signal capacitance matrices failed.\n";
    m_qwire.clear();
    return false;
  }
  if (!SetupPlaneSignals()) {
    std::cerr << m_className << "::PrepareSignals:\n"
              << "    Preparing plane charges failed.\n";
    m_qwire.clear();
    m_qplane.clear();
    return false;
  }

  // Associate wires, planes and strips with readout groups
  const unsigned int nReadout = m_readout.size();
  for (unsigned int i = 0; i < nReadout; ++i) {
    for (auto& wire : m_w) {
      if (wire.type == m_readout[i]) wire.ind = i;
    }
    for (size_t j = 0; j < 5; ++j) {
      if (m_planes[j].type == m_readout[i]) m_planes[j].ind = i;
      for (auto& strip : m_planes[j].strips1) {
        if (strip.type == m_readout[i]) strip.ind = i;
      }
      for (auto& strip : m_planes[j].strips2) {
        if (strip.type == m_readout[i]) strip.ind = i;
      }
      for (auto& pixel : m_planes[j].pixels) {
        if (pixel.type == m_readout[i]) pixel.ind = i;
      }
    }
  }

  // Seems to have worked.
  m_sigset = true;
  return true;
}

bool ComponentAnalyticField::SetupWireSignals() {
  //-----------------------------------------------------------------------
  //   SIGIPR - Prepares the ion tail calculation by filling the signal
  //            matrices (ie non-periodic capacitance matrices),
  //            Fourier transforming them if necessary, inverting them and
  //            Fourier transforming them back. Because of the large number
  //            of terms involved, a (scratch) external file on unit 13 is
  //            used to store the intermediate and final results. This file
  //            is handled by the routines IONBGN and IONIO.
  //   VARIABLES : FFTMAT      : Matrix used for Fourier transforms.
  //   (Last changed on  4/10/06.)
  //-----------------------------------------------------------------------

  std::vector<std::vector<std::complex<double> > > fftmat;
  fftmat.resize(m_nWires);
  for (size_t i = 0; i < m_nWires; ++i) {
    fftmat[i].assign(std::max(m_nFourier, 1), 0.);
  }

  size_t nLayers = 1;
  if (m_fperx && m_fpery) {
    nLayers = m_nFourier * m_nFourier;
  } else if (m_fperx || m_fpery) {
    nLayers = m_nFourier;
  }
  std::vector<std::vector<std::vector<std::complex<double> > > > sigmat;
  sigmat.resize(nLayers);
  for (size_t i = 0; i < nLayers; ++i) {
    sigmat[i].assign(m_nWires, std::vector<std::complex<double> >(m_nWires));
  }

  if (m_debug) {
    std::cout << m_className << "::SetupWireSignals:\n"
              << "    Calculating signal matrices.\n";
  }
  // Have the matrix/matrices filled (and stored).
  for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
    for (int my = m_mymin; my <= m_mymax; ++my) {
      // Select layer to be produced.
      const size_t k = SignalLayer(mx, my);
      if (m_cellTypeFourier == A00) {
        IprA00(mx, my, sigmat[k]);
      } else if (m_cellTypeFourier == B2X) {
        IprB2X(my, sigmat[k]);
      } else if (m_cellTypeFourier == B2Y) {
        IprB2Y(mx, sigmat[k]);
      } else if (m_cellTypeFourier == C2X) {
        IprC2X(sigmat[k]);
      } else if (m_cellTypeFourier == C2Y) {
        IprC2Y(sigmat[k]);
      } else if (m_cellTypeFourier == C30) {
        IprC30(sigmat[k]);
      } else if (m_cellTypeFourier == D10) {
        IprD10(sigmat[k]);
      } else if (m_cellTypeFourier == D30) {
        IprD30(sigmat[k]);
      } else {
        std::cerr << m_className << "::SetupWireSignals:\n"
                  << "    Unknown signal cell type " << m_cellTypeFourier
                  << "\n";
        return false;
      }

      // Dump the signal matrix before inversion, if DEBUG is requested.
      if (m_debug) {
        std::cout << "    Dump of signal matrix (" << mx << ", " << my
                  << ") before inversion:\n";
        for (size_t i = 0; i < m_nWires; i += 10) {
          for (size_t j = 0; j < m_nWires; j += 10) {
            std::cout << "    (Re-Block " << i / 10 << ", " << j / 10 << ")\n";
            for (unsigned int ii = 0; ii < 10; ++ii) {
              if (i + ii >= m_nWires) break;
              for (unsigned int jj = 0; jj < 10; ++jj) {
                if (j + jj >= m_nWires) break;
                std::cout << real(sigmat[k][i + ii][j + jj]) << "  ";
              }
              std::cout << "\n";
            }
            std::cout << "\n";
            std::cout << "    (Im-Block " << i / 10 << ", " << j / 10 << ")\n";
            for (unsigned int ii = 0; ii < 10; ++ii) {
              if (i + ii >= m_nWires) break;
              for (unsigned int jj = 0; jj < 10; ++jj) {
                if (j + jj >= m_nWires) break;
                std::cout << imag(sigmat[k][i + ii][j + jj]) << "  ";
              }
              std::cout << "\n";
            }
            std::cout << "\n";
          }
        }
        std::cout << "    End of the uninverted capacitance matrix dump.\n";
      }
      // Next layer.
    }
  }

  // Have them Fourier transformed (singly periodic case).
  if (m_fperx != m_fpery) {
    if (m_debug) std::cout << "    Calculating Fourier transforms.\n";
    const int mmin = -m_nFourier / 2 + 1;
    const int mmax =  m_nFourier / 2;
    for (size_t i = 0; i < m_nWires; ++i) {
      for (int m = mmin; m <= mmax; ++m) {
        const size_t k = SignalLayer(m, m);
        for (size_t j = 0; j < m_nWires; ++j) {
          fftmat[j][k] = sigmat[k][i][j];
        }
      }
      for (size_t j = 0; j < m_nWires; ++j) {
        Numerics::CERNLIB::cfft(fftmat[j], m_mfexp);
      }
      for (int m = mmin; m <= mmax; ++m) {
        const size_t k = SignalLayer(m, m);
        for (size_t j = 0; j < m_nWires; ++j) {
          sigmat[k][i][j] = fftmat[j][k];
        }
      }
    }
  }
  // Have them Fourier transformed (doubly periodic case).
  if (m_fperx && m_fpery) {
    if (m_debug) std::cout << "    Calculating Fourier transforms.\n";
    const int m0 = m_nFourier / 2 - 1;
    for (size_t i = 0; i < m_nWires; ++i) {
      for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
        for (int my = m_mymin; my <= m_mymax; ++my) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            fftmat[j][my + m0] = sigmat[k][i][j];
          }
        }
        for (size_t j = 0; j < m_nWires; ++j) {
          Numerics::CERNLIB::cfft(fftmat[j], m_mfexp);
        }
        for (int my = m_mymin; my <= m_mymax; ++my) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            sigmat[k][i][j] = fftmat[j][my + m0];
          }
        }
      }
      for (int my = m_mymin; my <= m_mymax; ++my) {
        for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            fftmat[j][mx + m0] = sigmat[k][i][j];
          }
        }
        for (size_t j = 0; j < m_nWires; ++j) {
          Numerics::CERNLIB::cfft(fftmat[j], m_mfexp);
        }
        for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            sigmat[k][i][j] = fftmat[j][mx + m0];
          }
        }
      }
    }
  }

  // Invert the matrices.
  if (m_debug) std::cout << "    Inverting the matrices.\n";
  for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
    for (int my = m_mymin; my <= m_mymax; ++my) {
      // Retrieve the layer.
      const size_t k = SignalLayer(mx, my);
      // Invert.
      if (m_nWires >= 1) {
        if (Numerics::CERNLIB::cinv(m_nWires, sigmat[k]) != 0) {
          std::cerr << m_className << "::PrepareWireSignals:\n"
                    << "    Inversion of signal matrix (" << mx << ", " << my
                    << ") failed.\n"
                    << "    No reliable results.\n"
                    << "    Preparation of weighting fields is abandoned.\n";
          return false;
        }
      }
    }
  }

  // And transform the matrices back to the original domain.
  if (m_fperx != m_fpery) {
    if (m_debug) std::cout << "    Calculating inverse FFT.\n";
    const int mmin = -m_nFourier / 2 + 1;
    const int mmax = m_nFourier / 2;
    for (size_t i = 0; i < m_nWires; ++i) {
      for (int m = mmin; m <= mmax; ++m) {
        const size_t k = SignalLayer(m, m);
        for (size_t j = 0; j < m_nWires; ++j) {
          fftmat[j][k] = sigmat[k][i][j];
        }
      }
      for (size_t j = 0; j < m_nWires; ++j) {
        Numerics::CERNLIB::cfft(fftmat[j], -m_mfexp);
      }
      for (int m = mmin; m <= mmax; ++m) {
        const size_t k = SignalLayer(m, m);
        for (size_t j = 0; j < m_nWires; ++j) {
          sigmat[k][i][j] = fftmat[j][k] / double(m_nFourier);
        }
      }
    }
  }
  // Have them transformed to the original domain (doubly periodic).
  if (m_fperx && m_fpery) {
    if (m_debug) std::cout << "    Calculating inverse FFT.\n";
    const int m0 = m_nFourier / 2 - 1;
    for (size_t i = 0; i < m_nWires; ++i) {
      for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
        for (int my = m_mymin; my <= m_mymax; ++my) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            fftmat[j][my + m0] = sigmat[k][i][j];
          }
        }
        for (size_t j = 0; j < m_nWires; ++j) {
          Numerics::CERNLIB::cfft(fftmat[j], -m_mfexp);
        }
        for (int my = m_mymin; my <= m_mymax; ++my) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            sigmat[k][i][j] = fftmat[j][my + m0] / double(m_nFourier);
          }
        }
      }
      for (int my = m_mymin; my <= m_mymax; ++my) {
        for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            fftmat[j][mx + m0] = sigmat[k][i][j];
          }
        }
        for (size_t j = 0; j < m_nWires; ++j) {
          Numerics::CERNLIB::cfft(fftmat[j], -m_mfexp);
        }
        for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
          const size_t k = SignalLayer(mx, my);
          for (size_t j = 0; j < m_nWires; ++j) {
            sigmat[k][i][j] = fftmat[j][mx + m0] / double(m_nFourier);
          }
        }
      }
    }
  }

  m_qwire.resize(nLayers);
  for (size_t k = 0; k < nLayers; ++k) {
    m_qwire[k].assign(m_nWires, std::vector<double>(m_nWires));
    for (unsigned int i = 0; i < m_nWires; ++i) {
      for (size_t j = 0; j < m_nWires; ++j) {
        m_qwire[k][i][j] = real(sigmat[k][i][j]);
      }
    }
  }

  // Dump the signal matrix after inversion, if DEBUG is requested.
  if (m_debug) {
    for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
      for (int my = m_mymin; my <= m_mymax; ++my) {
        std::cout << m_className << "::SetupWireSignals:\n"
                  << "    Dump of signal matrix (" << mx << ", " << my
                  << ") after inversion:\n";
        const size_t k = SignalLayer(mx, my);
        for (unsigned int i = 0; i < m_nWires; i += 10) {
          for (size_t j = 0; j < m_nWires; j += 10) {
            std::cout << "    (Re-Block " << i / 10 << ", " << j / 10 << ")\n";
            for (unsigned int ii = 0; ii < 10; ++ii) {
              if (i + ii >= m_nWires) break;
              for (unsigned int jj = 0; jj < 10; ++jj) {
                if (j + jj >= m_nWires) break;
                std::cout << real(sigmat[k][i + ii][j + jj]) << "  ";
              }
              std::cout << "\n";
            }
            std::cout << "\n";
            std::cout << "    (Im-Block " << i / 10 << ", " << j / 10 << ")\n";
            for (int ii = 0; ii < 10; ++ii) {
              if (i + ii >= m_nWires) break;
              for (int jj = 0; jj < 10; ++jj) {
                if (j + jj >= m_nWires) break;
                std::cout << imag(sigmat[k][i + ii][j + jj]) << "  ";
              }
              std::cout << "\n";
            }
            std::cout << "\n";
          }
        }
        std::cout << m_className << "::SetupWireSignals:\n"
                  << "    End of the inverted capacitance matrix dump.\n";
      }
    }
  }
  return true;
}

bool ComponentAnalyticField::SetupPlaneSignals() {
  //-----------------------------------------------------------------------
  //   SIGPLP - Computes the weighting field charges for the planes and
  //            the tube.
  //   (Last changed on 14/10/99.)
  //-----------------------------------------------------------------------

  constexpr size_t nPlanes = 5;
  const size_t nLayers = m_qwire.size();
  m_qplane.resize(nLayers);
  for (size_t k = 0; k < nLayers; ++k) {
    m_qplane[k].assign(nPlanes, std::vector<double>(m_nWires, 0.));
  }

  double vw;
  if (m_debug) std::cout << m_className << "::SetupPlaneSignals:\n";
  // Loop over the signal layers.
  for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
    for (int my = m_mymin; my <= m_mymax; ++my) {
      // Load the layers of the signal matrices.
      // CALL IONIO(MX,MY,2,0,IFAIL1)
      const size_t k = SignalLayer(mx, my);
      const auto& qw = m_qwire[k];
      // Charges for plane 1, if present.
      if (m_ynplan[0]) {
        // Set the weighting field voltages.
        for (unsigned int i = 0; i < m_nWires; ++i) {
          if (m_ynplan[1]) {
            vw = -(m_coplan[1] - m_w[i].x) / (m_coplan[1] - m_coplan[0]);
          } else if (m_perx) {
            vw = -(m_coplan[0] + m_sx - m_w[i].x) / m_sx;
          } else {
            vw = -1;
          }
          // Multiply with the matrix.
          for (size_t j = 0; j < m_nWires; ++j) {
            m_qplane[k][0][j] += qw[i][j] * vw;
          }
        }
      }
      // Charges for plane 2, if present.
      if (m_ynplan[1]) {
        // Set the weighting field voltages.
        for (unsigned int i = 0; i < m_nWires; ++i) {
          if (m_ynplan[0]) {
            vw = -(m_coplan[0] - m_w[i].x) / (m_coplan[0] - m_coplan[1]);
          } else if (m_perx) {
            vw = -(m_w[i].x - m_coplan[1] + m_sx) / m_sx;
          } else {
            vw = -1.;
          }
          // Multiply with the matrix.
          for (size_t j = 0; j < m_nWires; ++j) {
            m_qplane[k][1][j] += qw[i][j] * vw;
          }
        }
      }
      // Charges for plane 3, if present.
      if (m_ynplan[2]) {
        // Set the weighting field voltages.
        for (unsigned int i = 0; i < m_nWires; ++i) {
          if (m_ynplan[3]) {
            vw = -(m_coplan[3] - m_w[i].y) / (m_coplan[3] - m_coplan[2]);
          } else if (m_pery) {
            vw = -(m_coplan[2] + m_sy - m_w[i].y) / m_sy;
          } else {
            vw = -1.;
          }
          // Multiply with the matrix.
          for (size_t j = 0; j < m_nWires; ++j) {
            m_qplane[k][2][i] += qw[i][j] * vw;
          }
        }
      }
      // Charges for plane 4, if present.
      if (m_ynplan[3]) {
        // Set the weighting field voltages.
        for (unsigned int i = 0; i < m_nWires; ++i) {
          if (m_ynplan[2]) {
            vw = -(m_coplan[2] - m_w[i].y) / (m_coplan[2] - m_coplan[3]);
          } else if (m_pery) {
            vw = -(m_w[i].y - m_coplan[3] + m_sy) / m_sy;
          } else {
            vw = -1.;
          }
          // Multiply with the matrix.
          for (size_t j = 0; j < m_nWires; ++j) {
            m_qplane[k][3][i] += qw[i][j] * vw;
          }
        }
      }
      // Charges for the tube, if present.
      if (m_tube) {
        for (unsigned int i = 0; i < m_nWires; ++i) {
          for (size_t j = 0; j < m_nWires; ++j) {
            m_qplane[k][4][i] -= qw[i][j];
          }
        }
      }
      if (!m_debug || m_nWires == 0) continue;
      std::cout << "    Matrix (" << mx << ", " << my << "):\n";  
      if (m_tube) {
        std::cout << "    Charges for currents induced in the tube:\n"
                  << "    Wire\n";
        for (unsigned int i = 0; i < m_nWires; ++i) {
          std::printf("   %5d  %15.8f\n", i, m_qplane[k][4][i]);
        }
      } else {
        std::cout << "    Charges for currents induced in the planes:\n"
                  << "    Wire        x-Plane 1        x-Plane 2"
                  << "        y-Plane 1        y-Plane 2\n";
        for (unsigned int i = 0; i < m_nWires; ++i) {
          std::printf("   %5d  %15.8f  %15.8f  %15.8f  %15.8f\n", i,
                      m_qplane[k][0][i], m_qplane[k][1][i], 
                      m_qplane[k][2][i], m_qplane[k][3][i]);
        }
      }
    }
  }
  // Compute the background weighting fields.
  for (size_t i = 0; i < 5; ++i) {
    m_planes[i].ewxcor = 0.;
    m_planes[i].ewycor = 0.;
  }
  // First in x.
  if (m_ynplan[0] && m_ynplan[1]) {
    m_planes[0].ewxcor = 1. / (m_coplan[1] - m_coplan[0]);
    m_planes[1].ewxcor = 1. / (m_coplan[0] - m_coplan[1]);
  } else if (m_ynplan[0] && m_perx) {
    m_planes[0].ewxcor = 1. / m_sx;
    m_planes[1].ewxcor = 0.;
  } else if (m_ynplan[1] && m_perx) {
    m_planes[0].ewxcor = 0.;
    m_planes[1].ewxcor = -1. / m_sx;
  }
  // Next also in y.
  if (m_ynplan[2] && m_ynplan[3]) {
    m_planes[2].ewycor = 1. / (m_coplan[3] - m_coplan[2]);
    m_planes[3].ewycor = 1. / (m_coplan[2] - m_coplan[3]);
  } else if (m_ynplan[2] && m_pery) {
    m_planes[2].ewycor = 1. / m_sy;
    m_planes[3].ewycor = 0.;
  } else if (m_ynplan[3] && m_pery) {
    m_planes[2].ewycor = 0.;
    m_planes[3].ewycor = -1. / m_sy;
  }

  // Debugging output.
  if (m_debug) {
    std::cout << "    Bias fields:\n"
              << "    Plane    x-Bias [1/cm]    y-Bias [1/cm]\n";
    for (int i = 0; i < 4; ++i) {
      std::printf("   %5d  %15.8f  %15.8f\n", i,
                  m_planes[i].ewxcor, m_planes[i].ewycor);
    }
  }

  return true;
}

bool ComponentAnalyticField::IprA00(const int mx, const int my,
                                    CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRA00 - Routine filling the (MX,MY) th layer of the signal matrix
  //            for cells with non-periodic type A (see SIGIPR).
  //-----------------------------------------------------------------------

  const double dx = mx * m_sx;
  const double dy = my * m_sy;
  double aa = 0.;

  for (size_t i = 0; i < m_nWires; ++i) {
    // Diagonal terms.
    if (dx != 0. || dy != 0.) {
      aa = dx * dx + dy * dy;
    } else {
      aa = m_w[i].r * m_w[i].r;
    }
    // Take care of single equipotential planes.
    if (m_ynplax) aa /= 2. * pow(m_w[i].x - m_coplax, 2) + dy * dy;
    if (m_ynplay) aa /= 2. * pow(m_w[i].y - m_coplay, 2) + dx * dx;
    // Take care of pairs of equipotential planes.
    if (m_ynplax && m_ynplay)
      aa *= 4. * (pow(m_w[i].x - m_coplax, 2) + pow(m_w[i].y - m_coplay, 2));
    // Define the final version of a[i][i].
    mat[i][i] = -0.5 * log(aa);
    for (unsigned int j = i + 1; j < m_nWires; ++j) {
      aa = pow(m_w[i].x + dx - m_w[j].x, 2) + pow(m_w[i].y + dy - m_w[j].y, 2);
      // Take care of single planes.
      if (m_ynplax)
        aa /= pow(2. * m_coplax - m_w[i].x - dx - m_w[j].x, 2) +
              pow(m_w[i].y + dy - m_w[j].y, 2);
      if (m_ynplay)
        aa /= pow(m_w[i].x + dx - m_w[j].x, 2) +
              pow(2. * m_coplay - m_w[i].y - dy - m_w[j].y, 2);
      // Take care of pairs of planes.
      if (m_ynplax && m_ynplay) {
        aa *= pow(2. * m_coplax - m_w[i].x - dx - m_w[j].x, 2) +
              pow(2. * m_coplay - m_w[i].y - dy - m_w[j].y, 2);
      }
      // Store the true versions after taking LOGs and SQRT's.
      mat[i][j] = -0.5 * log(aa);
      mat[j][i] = mat[i][j];
    }
  }
  return true;
}

bool ComponentAnalyticField::IprB2X(const int my, CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRB2X - Routine filling the MY th layer of the signal matrix
  //            for cells with non-periodic type B2X (see SIGIPR).
  //   (Last changed on 26/ 4/92.)
  //-----------------------------------------------------------------------

  m_b2sin.resize(m_nWires);

  const double dy = my * m_sy;
  double aa = 0.;
  double xxneg;

  // Loop over all wires and calculate the diagonal elements first.
  for (size_t i = 0; i < m_nWires; ++i) {
    double xx = (Pi / m_sx) * (m_w[i].x - m_coplan[0]);
    if (dy != 0.) {
      aa = pow(sinh(Pi * dy / m_sx) / sin(xx), 2);
    } else {
      aa = pow((0.5 * m_w[i].r * Pi / m_sx) / sin(xx), 2);
    }
    // Take care of a planes at constant y (no dy in this case).
    if (m_ynplay) {
      const double yymirr = (Pi / m_sx) * (m_w[i].y - m_coplay);
      if (fabs(yymirr) <= 20.) {
        const double sinhy = sinh(yymirr);
        const double sinxx = sin(xx);
        aa *= (sinhy * sinhy + sinxx * sinxx) / (sinhy * sinhy);
      }
    }
    // Store the true value of A[i][i].
    mat[i][i] = -0.5 * log(aa);
    // Loop over all other wires to obtain off-diagonal elements.
    for (unsigned int j = i + 1; j < m_nWires; ++j) {
      const double yy = HalfPi * (m_w[i].y + dy - m_w[j].y) / m_sx;
      xx = HalfPi * (m_w[i].x - m_w[j].x) / m_sx;
      xxneg = HalfPi * (m_w[i].x + m_w[j].x - 2. * m_coplan[0]) / m_sx;
      if (fabs(yy) <= 20.) {
        const double sinhy = sinh(yy);
        const double sinxx = sin(xx);
        const double sinxxneg = sin(xxneg);
        aa = (sinhy * sinhy + sinxx * sinxx) /
             (sinhy * sinhy + sinxxneg * sinxxneg);
      } else {
        aa = 1.;
      }
      // Take equipotential planes into account (no dy anyhow).
      if (m_ynplay) {
        const double yymirr =
            HalfPi * (m_w[i].y + m_w[j].y - 2. * m_coplay) / m_sx;
        if (fabs(yymirr) <= 20.) {
          const double sinhy = sinh(yymirr);
          const double sinxx = sin(xx);
          const double sinxxneg = sin(xxneg);
          aa *= (sinhy * sinhy + sinxxneg * sinxxneg) /
                (sinhy * sinhy + sinxx * sinxx);
        }
      }
      // Store the true value of a[i][j] in both a[i][j] and a[j][i].
      mat[i][j] = -0.5 * log(aa);
      mat[j][i] = mat[i][j];
    }
    // Fill the B2SIN vector.
    m_b2sin[i] = sin(Pi * (m_coplan[0] - m_w[i].x) / m_sx);
  }

  return true;
}

bool ComponentAnalyticField::IprB2Y(const int mx, CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRB2Y - Routine filling the MX th layer of the signal matrix
  //            for cells with non-periodic type B2Y (see SIGIPR).
  //   (Last changed on 26/ 4/92.)
  //-----------------------------------------------------------------------

  m_b2sin.resize(m_nWires);

  const double dx = mx * m_sx;
  double aa = 0.;
  double xx, yy, xxmirr, yyneg;

  // Loop over all wires and calculate the diagonal elements first.
  for (size_t i = 0; i < m_nWires; ++i) {
    yy = (Pi / m_sy) * (m_w[i].y - m_coplan[2]);
    if (dx != 0.) {
      aa = pow(sinh(Pi * dx / m_sy) / sin(yy), 2);
    } else {
      aa = pow((0.5 * m_w[i].r * Pi / m_sy) / sin(yy), 2);
    }
    // Take care of a planes at constant x (no dx in this case).
    if (m_ynplax) {
      xxmirr = (Pi / m_sy) * (m_w[i].x - m_coplax);
      if (fabs(xxmirr) <= 20.) {
        aa *= (pow(sinh(xxmirr), 2) + pow(sin(yy), 2)) / pow(sinh(xxmirr), 2);
      }
    }
    // Store the true value of A[i][i].
    mat[i][i] = -0.5 * log(aa);
    // Loop over all other wires to obtain off-diagonal elements.
    for (unsigned int j = i + 1; j < m_nWires; ++j) {
      xx = HalfPi * (m_w[i].x + dx - m_w[j].x) / m_sy;
      yy = HalfPi * (m_w[i].y - m_w[j].y) / m_sy;
      yyneg = HalfPi * (m_w[i].y + m_w[j].y - 2. * m_coplan[2]) / m_sy;
      if (fabs(xx) <= 20.) {
        aa = (pow(sinh(xx), 2) + pow(sin(yy), 2)) /
             (pow(sinh(xx), 2) + pow(sin(yyneg), 2));
      } else {
        aa = 1.;
      }
      // Take equipotential planes into account (no dx anyhow).
      if (m_ynplax) {
        xxmirr = HalfPi * (m_w[i].x + m_w[j].x - 2. * m_coplax) / m_sy;
        if (fabs(xxmirr) <= 20.) {
          aa *= (pow(sinh(xxmirr), 2) + pow(sin(yyneg), 2)) /
                (pow(sinh(xxmirr), 2) + pow(sin(yy), 2));
        }
      }
      // Store the true value of a[i][j] in both a[i][j] and a[j][i].
      mat[i][j] = -0.5 * log(aa);
      mat[j][i] = mat[i][j];
    }
    // Fill the B2SIN vector.
    m_b2sin[i] = sin(Pi * (m_coplan[2] - m_w[i].y) / m_sy);
  }
  return true;
}

bool ComponentAnalyticField::IprC2X(CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRC2X - This initializing subroutine stores the capacitance matrix
  //            for the configuration:
  //            wires at zw(j)+cmplx(lx*2*sx,ly*sy),
  //            j=1(1)n, lx=-infinity(1)infinity, ly=-infinity(1)infinity.
  //            but the signs of the charges alternate in the x-direction
  //   (Last changed on  4/10/06.)
  //-----------------------------------------------------------------------

  // Fill the capacitance matrix.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    for (unsigned int j = 0; j < m_nWires; ++j) {
      double temp = 0.;
      if (m_mode == 0) {
        temp = (m_w[i].x - cx) * (m_w[j].x - cx) * TwoPi / (m_sx * m_sy);
      }
      if (i == j) {
        mat[i][j] = Ph2Lim(m_w[i].r) - Ph2(2. * (m_w[j].x - cx), 0.) - temp;
      } else {
        mat[i][j] =
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) -
            Ph2(m_w[i].x + m_w[j].x - 2. * cx, m_w[i].y - m_w[j].y) - temp;
      }
    }
  }
  return true;
}

bool ComponentAnalyticField::IprC2Y(CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRC2Y - This initializing subroutine stores the capacitance matrix
  //            for the configuration:
  //            wires at zw(j)+cmplx(lx*sx,ly*2*sy),
  //            j=1(1)n, lx=-infinity(1)infinity, ly=-infinity(1)infinity.
  //            but the signs of the charges alternate in the y-direction
  //   (Last changed on  4/10/06.)
  //-----------------------------------------------------------------------

  // Fill the capacitance matrix.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    for (unsigned int j = 0; j < m_nWires; ++j) {
      double temp = 0.;
      if (m_mode == 1) {
        temp = (m_w[i].y - cy) * (m_w[j].y - cy) * TwoPi / (m_sx * m_sy);
      }
      if (i == j) {
        mat[i][j] = Ph2Lim(m_w[i].r) - Ph2(0., 2. * (m_w[j].y - cy)) - temp;
      } else {
        mat[i][j] =
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) -
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y + m_w[j].y - 2. * cy) - temp;
      }
    }
  }
  return true;
}

bool ComponentAnalyticField::IprC30(CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRC30 - Routine filling the signal matrix for cells of type C30.
  //            Since the signal matrix equals the capacitance matrix for
  //            this potential, the routine is identical to SETC30 except
  //            for the C and P parameters.
  //   (Last changed on 11/11/97.)
  //-----------------------------------------------------------------------

  // Fill the capacitance matrix.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    for (unsigned int j = 0; j < m_nWires; ++j) {
      if (i == j) {
        mat[i][i] = Ph2Lim(m_w[i].r) -
                    Ph2(0., 2. * (m_w[i].y - cy)) -
                    Ph2(2. * (m_w[i].x - cx), 0.) +
                    Ph2(2. * (m_w[i].x - cx), 2. * (m_w[i].y - cy));
      } else {
        mat[i][j] =
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y - m_w[j].y) -
            Ph2(m_w[i].x - m_w[j].x, m_w[i].y + m_w[j].y - 2. * cy) -
            Ph2(m_w[i].x + m_w[j].x - 2. * cx, m_w[i].y - m_w[j].y) +
            Ph2(m_w[i].x + m_w[j].x - 2. * cx, m_w[i].y + m_w[j].y - 2. * cy);
      }
    }
  }
  return true;
}

bool ComponentAnalyticField::IprD10(CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRD10 - Signal matrix preparation for D1 cells.
  //   VARIABLES :
  //   (Last changed on  2/ 2/93.)
  //-----------------------------------------------------------------------

  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set the diagonal terms.
    mat[i][i] = -log(
        m_w[i].r /
        (m_cotube - (m_w[i].x * m_w[i].x + m_w[i].y * m_w[i].y) / m_cotube));
    // Set a complex wire-coordinate to make things a little easier.
    std::complex<double> zi(m_w[i].x, m_w[i].y);
    // Loop over all other wires for the off-diagonal elements.
    for (size_t j = i + 1; j < m_nWires; ++j) {
      // Set a complex wire-coordinate to make things a little easier.
      std::complex<double> zj(m_w[j].x, m_w[j].y);
      mat[i][j] = -log(
          std::abs((1. / m_cotube) * (zi - zj) / (1. - conj(zi) * zj / m_cotube2)));
      // Copy this to a[j][i] since the capacitance matrix is symmetric.
      mat[j][i] = mat[i][j];
    }
  }
  return true;
}

bool ComponentAnalyticField::IprD30(CMatrix& mat) {
  //-----------------------------------------------------------------------
  //   IPRD30 - Signal matrix preparation for polygonal cells (type D3).
  //   Variables :
  //   (Last changed on 19/ 6/97.)
  //-----------------------------------------------------------------------

  m_zw.resize(m_nWires);

  std::complex<double> wd;

  // Loop over all wire combinations.
  for (int i = 0; i < int(m_nWires); ++i) {
    // We need to compute the wire mapping again to obtain WD.
    ConformalMap(std::complex<double>(m_w[i].x, m_w[i].y) / m_cotube, m_zw[i],
                 wd);
    // Diagonal elements.
    mat[i][i] = -log(
        std::abs((m_w[i].r / m_cotube) * wd / (1. - std::norm(m_zw[i]))));
    // Loop over all other wires for the off-diagonal elements.
    for (int j = 0; j < i - 1; ++j) {
      mat[i][j] = -log(
          std::abs((m_zw[i] - m_zw[j]) / (1. - conj(m_zw[i]) * m_zw[j])));
      // Copy this to a[j][i] since the capacitance matrix is symmetric.
      mat[j][i] = mat[i][j];
    }
  }
  return true;
}

bool ComponentAnalyticField::Wfield(const double xin, const double yin,
                                    const double zpos, double& exsum,
                                    double& eysum, double& ezsum,
                                    const std::string& label) const {
  //-----------------------------------------------------------------------
  //   SIGFLS - Sums the weighting field components at (XPOS,YPOS,ZPOS).
  //   (Last changed on 11/10/06.)
  //-----------------------------------------------------------------------

  // Initialise the sums.
  exsum = eysum = ezsum = 0.;
  double ex = 0., ey = 0., ez = 0.;

  double xpos = xin, ypos = yin;
  if (m_polar) Cartesian2Internal(xin, yin, xpos, ypos);
  // Stop here if there are no weighting fields defined.
  if (m_readout.empty()) return false;
  if (!m_sigset) {
    std::cerr << m_className << "::Wfield: No weighting fields available.\n";
    return false;
  }
  
  // In case (xpos, ypos) is located behind a plane there is no field.
  if (m_tube) {
    if (!InTube(xpos, ypos, m_cotube, m_ntube)) return false;
  } else {
    if (!m_perx) {
      if ((m_ynplan[0] && xpos < m_coplan[0]) ||
          (m_ynplan[1] && xpos > m_coplan[1])) {
        return false;
      }
    }
    if (!m_pery) {
      if ((m_ynplan[2] && ypos < m_coplan[2]) ||
          (m_ynplan[3] && ypos > m_coplan[3])) {
        return false;
      }
    }
  }

  if (label.empty()) return false;
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end()) return false;
  const auto isw = it - m_readout.begin();

  // Loop over the signal layers.
  for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
    for (int my = m_mymin; my <= m_mymax; ++my) {
      // Get the layers of the wire matrices.
      // CALL IONIO(MX,MY,2,0,IFAIL)
      const size_t k = SignalLayer(mx, my);
      // Loop over all wires.
      for (size_t iw = 0; iw < m_nWires; ++iw) {
        // Pick out those wires that are part of this read out group.
        if (m_w[iw].ind != isw) continue;
        const auto& qw = m_qwire[k][iw];
        ex = ey = ez = 0.;
        if (m_cellTypeFourier == A00) {
          WfieldWireA00(xpos, ypos, ex, ey, mx, my, qw);
        } else if (m_cellTypeFourier == B2X) {
          WfieldWireB2X(xpos, ypos, ex, ey, my, qw);
        } else if (m_cellTypeFourier == B2Y) {
          WfieldWireB2Y(xpos, ypos, ex, ey, mx, qw);
        } else if (m_cellTypeFourier == C2X) {
          WfieldWireC2X(xpos, ypos, ex, ey, qw);
        } else if (m_cellTypeFourier == C2Y) {
          WfieldWireC2Y(xpos, ypos, ex, ey, qw);
        } else if (m_cellTypeFourier == C30) {
          WfieldWireC30(xpos, ypos, ex, ey, qw);
        } else if (m_cellTypeFourier == D10) {
          WfieldWireD10(xpos, ypos, ex, ey, qw);
        } else if (m_cellTypeFourier == D30) {
          WfieldWireD30(xpos, ypos, ex, ey, qw);
        } else {
          std::cerr << m_className << "::Wfield:\n"
                    << "    Unknown signal field type " << m_cellTypeFourier
                    << " received. Program error!\n"
                    << "    Encountered for wire " << iw
                    << ", readout group = " << m_w[iw].ind << "\n";
          exsum = eysum = ezsum = 0.;
          return false;
        }
        exsum += ex;
        eysum += ey;
        ezsum += ez;
      }
      // Loop over all planes.
      for (size_t ip = 0; ip < 5; ++ip) {
        // Pick out those that are part of this read out group.
        if (m_planes[ip].ind != isw) continue;
        ex = ey = ez = 0.;
        const auto& qp = m_qplane[k][ip];
        if (m_cellTypeFourier == A00) {
          WfieldPlaneA00(xpos, ypos, ex, ey, mx, my, qp);
        } else if (m_cellTypeFourier == B2X) {
          WfieldPlaneB2X(xpos, ypos, ex, ey, my, qp);
        } else if (m_cellTypeFourier == B2Y) {
          WfieldPlaneB2Y(xpos, ypos, ex, ey, mx, qp);
        } else if (m_cellTypeFourier == C2X) {
          WfieldPlaneC2X(xpos, ypos, ex, ey, qp);
        } else if (m_cellTypeFourier == C2Y) {
          WfieldPlaneC2Y(xpos, ypos, ex, ey, qp);
        } else if (m_cellTypeFourier == C30) {
          WfieldPlaneC30(xpos, ypos, ex, ey, qp);
        } else if (m_cellTypeFourier == D10) {
          WfieldPlaneD10(xpos, ypos, ex, ey, qp);
        } else if (m_cellTypeFourier == D30) {
          WfieldPlaneD30(xpos, ypos, ex, ey, qp);
        } else {
          std::cerr << m_className << "::Wfield:\n"
                    << "    Unknown field type " << m_cellTypeFourier
                    << " received. Program error!\n"
                    << "    Encountered for plane " << ip
                    << ", readout group = " << m_planes[ip].ind << "\n";
          exsum = eysum = ezsum = 0.;
          return false;
        }
        exsum += ex;
        eysum += ey;
        ezsum += ez;
      }
      // Next signal layer.
    }
  }

  for (size_t ip = 0; ip < 5; ++ip) {
    // Add the field due to the plane itself.
    if (m_planes[ip].ind == isw) {
      exsum += m_planes[ip].ewxcor;
      eysum += m_planes[ip].ewycor;
    }
    // Add strips and pixels, if there are any.
    for (const auto& strip : m_planes[ip].strips1) {
      if (strip.ind != isw) continue;
      WfieldStripXy(xpos, ypos, zpos, ex, ey, ez, ip, strip);
      exsum += ex;
      eysum += ey;
      ezsum += ez;
    }
    for (const auto& strip : m_planes[ip].strips2) {
      if (strip.ind != isw) continue;
      WfieldStripZ(xpos, ypos, ex, ey, ip, strip);
      exsum += ex;
      eysum += ey;
    }
    for (const auto& pixel : m_planes[ip].pixels) {
      if (pixel.ind != isw) continue;
      WfieldPixel(xpos, ypos, zpos, ex, ey, ez, ip, pixel);
      exsum += ex;
      eysum += ey;
      ezsum += ez;
    }
  }
  if (m_polar) {
    const double r = exp(xpos);
    const double er = exsum / r;
    const double ep = eysum / r;
    const double theta = atan2(yin, xin);
    const double ct = cos(theta);
    const double st = sin(theta);
    exsum = +ct * er - st * ep;
    eysum = +st * er + ct * ep; 
  }
  return true;
}

double ComponentAnalyticField::Wpot(const double xin, const double yin,
                                    const double zpos,
                                    const std::string& label) const {

  double vsum = 0.;
  double xpos = xin, ypos = yin;
  if (m_polar) Cartesian2Internal(xin, yin, xpos, ypos);
  // Stop here if there are no weighting fields defined.
  if (m_readout.empty()) return -1.;
  if (!m_sigset) {
    std::cerr << m_className << "::Wpot: No weighting potentials available.\n";
    return -1.;
  }

  // In case (xpos, ypos) is located behind a plane there is no field.
  if (m_tube) {
    if (!InTube(xpos, ypos, m_cotube, m_ntube)) return -1.;
  } else {
    if (!m_perx) {
      if ((m_ynplan[0] && xpos < m_coplan[0]) ||
          (m_ynplan[1] && xpos > m_coplan[1])) {
        return -1.;
      }
    }
    if (!m_pery) {
      if ((m_ynplan[2] && ypos < m_coplan[2]) ||
          (m_ynplan[3] && ypos > m_coplan[3])) {
        return -1.;
      }
    }
  }
  if (label.empty()) return -1.;
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it == m_readout.end()) return -1.;
  const auto isw = it - m_readout.begin();

  // Loop over the signal layers.
  for (int mx = m_mxmin; mx <= m_mxmax; ++mx) {
    for (int my = m_mymin; my <= m_mymax; ++my) {
      const size_t k = SignalLayer(mx, my);
      // Loop over all wires.
      for (size_t iw = 0; iw < m_nWires; ++iw) {
        // Pick out those wires that are part of this read out group.
        if (m_w[iw].ind != isw) continue;
        const auto& qw = m_qwire[k][iw];
        if (m_cellTypeFourier == A00) {
          vsum += WpotWireA00(xpos, ypos, mx, my, qw);
        } else if (m_cellTypeFourier == B2X) {
          vsum += WpotWireB2X(xpos, ypos, my, qw);
        } else if (m_cellTypeFourier == B2Y) {
          vsum += WpotWireB2Y(xpos, ypos, mx, qw);
        } else if (m_cellTypeFourier == C2X) {
          vsum += WpotWireC2X(xpos, ypos, qw);
        } else if (m_cellTypeFourier == C2Y) {
          vsum += WpotWireC2Y(xpos, ypos, qw);
        } else if (m_cellTypeFourier == C30) {
          vsum += WpotWireC30(xpos, ypos, qw);
        } else if (m_cellTypeFourier == D10) {
          vsum += WpotWireD10(xpos, ypos, qw);
        } else if (m_cellTypeFourier == D30) {
          vsum += WpotWireD30(xpos, ypos, qw);
        } else {
          std::cerr << m_className << "::Wpot:\n"
                    << "    Unknown signal field type " << m_cellTypeFourier
                    << " received. Program error!\n"
                    << "    Encountered for wire " << iw
                    << ", readout group = " << m_w[iw].ind << "\n";
          return -1.;
        }
      }
      // Loop over all planes.
      for (size_t ip = 0; ip < 5; ++ip) {
        // Pick out those that are part of this read out group.
        if (m_planes[ip].ind != isw) continue;
        const auto& qp = m_qplane[k][ip];
        if (m_cellTypeFourier == A00) {
          vsum += WpotPlaneA00(xpos, ypos, mx, my, qp);
        } else if (m_cellTypeFourier == B2X) {
          vsum += WpotPlaneB2X(xpos, ypos, my, qp);
        } else if (m_cellTypeFourier == B2Y) {
          vsum += WpotPlaneB2Y(xpos, ypos, mx, qp);
        } else if (m_cellTypeFourier == C2X) {
          vsum += WpotPlaneC2X(xpos, ypos, qp);
        } else if (m_cellTypeFourier == C2Y) {
          vsum += WpotPlaneC2Y(xpos, ypos, qp);
        } else if (m_cellTypeFourier == C30) {
          vsum += WpotPlaneC30(xpos, ypos, qp);
        } else if (m_cellTypeFourier == D10) {
          vsum += WpotPlaneD10(xpos, ypos, qp);
        } else if (m_cellTypeFourier == D30) {
          vsum += WpotPlaneD30(xpos, ypos, qp);
        } else {
          std::cerr << m_className << "::Wpot:\n"
                    << "    Unknown field type " << m_cellTypeFourier
                    << " received. Program error!\n"
                    << "    Encountered for plane " << ip
                    << ", readout group = " << m_planes[ip].ind << "\n";
          return -1.;
        }
      }
      // Next signal layer.
    }
  }
  // Add the potential due to the planes themselves.
  for (size_t ip = 0; ip < 5; ++ip) {
    if (m_planes[ip].ind != isw) continue;
    if (ip == 0 || ip == 1) {
      double xx = xpos;
      if (m_perx) {
        xx -= m_sx * round(xpos / m_sx);
        if (m_ynplan[0] && xx <= m_coplan[0]) xx += m_sx;
        if (m_ynplan[1] && xx >= m_coplan[1]) xx -= m_sx;
      }
      vsum += 1. - m_planes[ip].ewxcor * (xx - m_coplan[ip]);
    } else if (ip == 2 || ip == 3) {
      double yy = ypos;
      if (m_pery) {
        yy -= m_sy * round(ypos / m_sy);
        if (m_ynplan[2] && yy <= m_coplan[2]) yy += m_sy;
        if (m_ynplan[3] && yy >= m_coplan[3]) yy -= m_sy;
      }
      vsum += 1. - m_planes[ip].ewycor * (yy - m_coplan[ip]);
    }
  }

  // Add strips and pixels, if there are any.
  for (size_t ip = 0; ip < 5; ++ip) {
    for (const auto& strip : m_planes[ip].strips1) {
      if (strip.ind != isw) continue;
      vsum += WpotStripXy(xpos, ypos, zpos, ip, strip);
    }
    for (const auto& strip : m_planes[ip].strips2) {
      if (strip.ind != isw) continue;
      vsum += WpotStripZ(xpos, ypos, ip, strip);
    }
    for (const auto& pixel : m_planes[ip].pixels) {
      if (pixel.ind != isw) continue;
      vsum += WpotPixel(xpos, ypos, zpos, ip, pixel);
    }
  }
  return vsum;
}

void ComponentAnalyticField::WfieldWireA00(const double xpos, const double ypos,
    double& ex, double& ey, const int mx, const int my,
    const std::vector<double>& qw) const {
  //-----------------------------------------------------------------------
  //   IONA00 - Routine returning the A I,J [MX,MY] * E terms for A cells.
  //   VARIABLES : R2         : Potential before taking -Log(Sqrt(...))
  //               EX,EY      : x,y-Component of the electric field.
  //               ETOT       : Magnitude of the electric field.
  //               VOLT       : Potential.
  //               EXHELP ETC : One term in the summing series.
  //               (XPOS,YPOS): Position where the field is needed.
  //   (Last changed on 14/ 8/98.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Define a few reduced variables.
    const double xx = xpos - m_w[i].x - mx * m_sx;
    const double yy = ypos - m_w[i].y - my * m_sy;
    // Calculate the field in case there are no planes.
    const double r2 = xx * xx + yy * yy;
    const double s2 = r2 > 0. ? 1. / r2 : 0.;
    double exhelp = xx * s2;
    double eyhelp = yy * s2;
    // Take care of a plane at constant x.
    const double xxmirr = m_ynplax ? xpos + m_w[i].x - 2. * m_coplax : 0.;
    if (m_ynplax) {
      const double r2plan = xxmirr * xxmirr + yy * yy;
      const double s2plan = r2plan > 0. ? 1. / r2plan : 0.;
      exhelp -= xxmirr * s2plan;
      eyhelp -= yy * s2plan;
    }
    // Take care of a plane at constant y.
    const double yymirr = m_ynplay ? ypos + m_w[i].y - 2. * m_coplay : 0.;
    if (m_ynplay) {
      const double r2plan = xx * xx + yymirr * yymirr;
      const double s2plan = r2plan > 0. ? 1. / r2plan : 0.;
      exhelp -= xx * s2plan;
      eyhelp -= yymirr * s2plan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double r2plan = xxmirr * xxmirr + yymirr * yymirr;
      const double s2plan = r2plan > 0. ? 1. / r2plan : 0.;
      exhelp += xxmirr * s2plan;
      eyhelp += yymirr * s2plan;
    }
    // Calculate the electric field.
    ex += qw[i] * exhelp;
    ey += qw[i] * eyhelp;
  }
}

double ComponentAnalyticField::WpotWireA00(const double xpos, const double ypos,
    const int mx, const int my, const std::vector<double>& qw) const { 

  double volt = 0.;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Define a few reduced variables.
    const double xx = xpos - m_w[i].x - mx * m_sx;
    const double yy = ypos - m_w[i].y - my * m_sy;
    // Calculate the field in case there are no planes.
    double r2 = xx * xx + yy * yy;
    if (r2 <= 0.) continue;
    // Take care of a plane at constant x.
    const double xxmirr = m_ynplax ? xpos + m_w[i].x - 2. * m_coplax : 0.;
    if (m_ynplax) {
      const double r2plan = xxmirr * xxmirr + yy * yy;
      if (r2plan <= 0.) continue;
      r2 /= r2plan;
    }
    // Take care of a plane at constant y.
    const double yymirr = m_ynplay ? ypos + m_w[i].y - 2. * m_coplay : 0.;
    if (m_ynplay) {
      const double r2plan = xx * xx + yymirr * yymirr;
      if (r2plan <= 0.) continue;
      r2 /= r2plan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double r2plan = xxmirr * xxmirr + yymirr * yymirr;
      if (r2plan <= 0.) continue;
      r2 *= r2plan;
    }
    // Calculate the potential.
    volt -= qw[i] * log(r2);
  }
  return 0.5 * volt;
}

void ComponentAnalyticField::WfieldWireB2X(const double xpos, const double ypos,
    double& ex, double& ey, const int my, 
    const std::vector<double>& qw) const {
  //-----------------------------------------------------------------------
  //   IONB2X - Routine calculating the MY contribution to the signal on
  //            wire ISW due to a charge at (XPOS,YPOS) for F-B2Y cells.
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z,ZZMIRR   : X + I*Y , XXMIRR + I*YYMIRR ; I**2=-1
  //               ECOMPL     : EX + I*EY                   ; I**2=-1
  //   (Last changed on 20/ 2/90.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  const double tx = HalfPi / m_sx;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = tx * (xpos - m_w[i].x);
    const double yy = tx * (ypos - m_w[i].y - my * m_sy);
    const double xxneg = tx * (xpos + m_w[i].x - 2 * m_coplan[0]);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    if (fabs(yy) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xxneg, yy);
      ecompl = -m_b2sin[i] / (sin(zz) * sin(zzneg));
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + m_w[i].y - 2. * m_coplay);
      if (fabs(yymirr) <= 20.) {
        const std::complex<double> zzmirr(xx, yymirr);
        const std::complex<double> zznmirr(xxneg, yymirr);
        ecompl += m_b2sin[i] / (sin(zzmirr) * sin(zznmirr));
      }
    }
    // Calculate the electric field.
    ex += qw[i] * real(ecompl);
    ey -= qw[i] * imag(ecompl);
  }
  ex *= tx;
  ey *= tx;
}

double ComponentAnalyticField::WpotWireB2X(const double xpos, const double ypos,
    const int my, const std::vector<double>& qw) const {
  double volt = 0.;
  const double tx = HalfPi / m_sx;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = tx * (xpos - m_w[i].x);
    const double yy = tx * (ypos - m_w[i].y - my * m_sy);
    const double xxneg = tx * (xpos + m_w[i].x - 2 * m_coplan[0]);
    // Calculate the potential in case there are no equipotential planes.
    double r2 = 1.;
    if (fabs(yy) <= 20.) {
      const double sinhy = sinh(yy);
      const double sinxx = sin(xx);
      const double sinxxneg = sin(xxneg);
      r2 = (sinhy * sinhy + sinxx * sinxx) /
           (sinhy * sinhy + sinxxneg * sinxxneg);
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + m_w[i].y - 2. * m_coplay);
      if (fabs(yymirr) <= 20.) {
        const double sinhy = sinh(yymirr);
        const double sinxx = sin(xx);
        const double sinxxneg = sin(xxneg);
        const double r2plan = (sinhy * sinhy + sinxx * sinxx) /
                              (sinhy * sinhy + sinxxneg * sinxxneg);
        r2 /= r2plan;
      }
    }
    // Calculate the potential.
    volt -= qw[i] * log(r2);
  }
  return 0.5 * volt;
}

void ComponentAnalyticField::WfieldWireB2Y(const double xpos, const double ypos,
    double& ex, double& ey, const int mx, 
    const std::vector<double>& qw) const {
  //-----------------------------------------------------------------------
  //   IONB2Y - Routine calculating the MX contribution to the signal on
  //            wire ISW due to a charge at (XPOS,YPOS) for F-B2X cells.
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z,ZZMIRR   : X + I*Y , XXMIRR + I*YYMIRR ; I**2=-1
  //               ECOMPL     : EX + I*EY                   ; I**2=-1
  //   (Last changed on 20/ 2/90.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  ex = ey = 0.;
  const double ty = HalfPi / m_sy;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = ty * (xpos - m_w[i].x - mx * m_sx);
    const double yy = ty * (ypos - m_w[i].y);
    const double yyneg = ty * (ypos + m_w[i].y - 2. * m_coplan[2]);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    if (fabs(xx) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xx, yyneg);
      ecompl = icons * m_b2sin[i] / (sin(icons * zz) * sin(icons * zzneg));
    }
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + m_w[i].x - 2 * m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const std::complex<double> zzmirr(xxmirr, yy);
        const std::complex<double> zznmirr(xxmirr, yyneg);
        ecompl -=
            icons * m_b2sin[i] / (sin(icons * zzmirr) * sin(icons * zznmirr));
      }
    }
    // Calculate the electric field.
    ex += qw[i] * real(ecompl);
    ey -= qw[i] * imag(ecompl);
  }
  ex *= ty;
  ey *= ty;
}

double ComponentAnalyticField::WpotWireB2Y(const double xpos, const double ypos,
    const int mx, const std::vector<double>& qw) const {
  double volt = 0.;
  const double ty = HalfPi / m_sy;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = ty * (xpos - m_w[i].x - mx * m_sx);
    const double yy = ty * (ypos - m_w[i].y);
    const double yyneg = ty * (ypos + m_w[i].y - 2. * m_coplan[2]);
    // Calculate the field in case there are no equipotential planes.
    double r2 = 1.;
    if (fabs(xx) <= 20.) {
      const double sinhx = sinh(xx);
      const double sinyy = sin(yy);
      const double sinyyneg = sin(yyneg);
      r2 = (sinhx * sinhx + sinyy * sinyy) /
           (sinhx * sinhx + sinyyneg * sinyyneg);
    }
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + m_w[i].x - 2 * m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const double sinhx = sinh(xxmirr);
        const double sinyy = sin(yy);
        const double sinyyneg = sin(yyneg);
        const double r2plan = (sinhx * sinhx + sinyy * sinyy) /
                              (sinhx * sinhx + sinyyneg * sinyyneg);
        r2 /= r2plan;
      }
    }
    volt -= qw[i] * log(r2);
  }
  return 0.5 * volt;
}

void ComponentAnalyticField::WfieldWireC2X(const double xpos, const double ypos,
    double& ex, double& ey, const std::vector<double>& qw) const {
  //-----------------------------------------------------------------------
  //   IONC2X - Routine returning the potential and electric field in a
  //            configuration with 2 x planes and y periodicity.
  //   VARIABLES : see the writeup
  //   (Last changed on 12/10/06.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  double s = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (imag(zeta) > 15.) {
      wsum1 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum1 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += qw[i] * (zterm.second / zterm.first);
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    // Constant terms sum
    s += qw[i] * (m_w[i].x - cx);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(2. * cx - xpos - m_w[i].x, yy);
    if (imag(zeta) > +15.) {
      wsum2 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += qw[i] * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2));
  ey = -imag(m_zmult * (wsum1 - wsum2));
  // Constant correction terms.
  if (m_mode == 0) ex += s * TwoPi / (m_sx * m_sy);
}

double ComponentAnalyticField::WpotWireC2X(const double xpos, const double ypos,
    const std::vector<double>& qw) const {
  double volt = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qw[i] * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(2. * cx - xpos - m_w[i].x, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt += qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qw[i] * log(std::abs(zterm.first));
    }
    // Correct the voltage, if needed (MODE).
    if (m_mode == 0) {
      volt -= TwoPi * qw[i] * (xpos - cx) * (m_w[i].x - cx) / (m_sx * m_sy);
    }
  }
  return volt;
}

void ComponentAnalyticField::WfieldWireC2Y(const double xpos, const double ypos,
    double& ex, double& ey, const std::vector<double>& qw) const { 
  //-----------------------------------------------------------------------
  //   IONC2Y - Routine returning the potential and electric field in a
  //            configuration with 2 y planes and x periodicity.
  //   VARIABLES : see the writeup
  //   (Last changed on 12/10/06.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  double s = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (imag(zeta) > +15.) {
      wsum1 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum1 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += qw[i] * (zterm.second / zterm.first);
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    // Constant terms sum
    s += qw[i] * (m_w[i].y - cy);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(xx, 2. * cy - ypos - m_w[i].y);
    if (imag(zeta) > +15.) {
      wsum2 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += qw[i] * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 - wsum2));
  ey = -imag(m_zmult * (wsum1 + wsum2));
  // Constant correction terms.
  if (m_mode == 1) ey += s * TwoPi / (m_sx * m_sy);
}

double ComponentAnalyticField::WpotWireC2Y(const double xpos, const double ypos,
    const std::vector<double>& qw) const {

  double volt = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qw[i] * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(xx, 2. * cy - ypos - m_w[i].y);
    if (fabs(imag(zeta)) > 15.) {
      volt += qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qw[i] * log(std::abs(zterm.first));
    }
    // Correct the voltage, if needed (MODE).
    if (m_mode == 1) {
      volt -= TwoPi * qw[i] * (ypos - cy) * (m_w[i].y - cy) / (m_sx * m_sy);
    }
  }
  return volt;
}

void ComponentAnalyticField::WfieldWireC30(const double xpos, const double ypos,
    double& ex, double& ey, const std::vector<double>& qw) const { 
  //-----------------------------------------------------------------------
  //   IONC30 - Routine returning the weighting field field in a
  //            configuration with 2 y and 2 x planes. This routine is
  //            basically the same as EFCC30.
  //   (Last changed on 11/11/97.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  std::complex<double> wsum3 = 0.;
  std::complex<double> wsum4 = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (imag(zeta) > +15.) {
      wsum1 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum1 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += qw[i] * (zterm.second / zterm.first);
    }
    // Mirror contribution from the x plane.
    const double xxmirr = MirrorCoordinate(xpos, m_coplax, m_w[i].x, m_sx);
    zeta = m_zmult * std::complex<double>(xxmirr, yy);
    if (imag(zeta) > +15.) {
      wsum2 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += qw[i] * (zterm.second / zterm.first);
    }
    // Mirror contribution from the y plane.
    const double yymirr = MirrorCoordinate(ypos, m_coplay, m_w[i].y, m_sy);
    zeta = m_zmult * std::complex<double>(xx, yymirr);
    if (imag(zeta) > +15.) {
      wsum3 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum3 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum3 += qw[i] * (zterm.second / zterm.first);
    }
    // Mirror contribution from both the x and the y plane.
    zeta = m_zmult * std::complex<double>(xxmirr, yymirr);
    if (imag(zeta) > +15.) {
      wsum4 -= qw[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum4 += qw[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum4 += qw[i] * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2 - wsum3 - wsum4));
  ey = -imag(m_zmult * (wsum1 - wsum2 + wsum3 - wsum4));
}

double ComponentAnalyticField::WpotWireC30(const double xpos, const double ypos,
    const std::vector<double>& qw) const {
  double volt = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qw[i] * log(std::abs(zterm.first));
    }
    // Mirror contribution from the x plane.
    const double xxmirr = MirrorCoordinate(xpos, m_coplax, m_w[i].x, m_sx);
    zeta = m_zmult * std::complex<double>(xxmirr, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt += qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qw[i] * log(std::abs(zterm.first));
    }
    // Mirror contribution from the y plane.
    const double yymirr = MirrorCoordinate(ypos, m_coplay, m_w[i].y, m_sy);
    zeta = m_zmult * std::complex<double>(xx, yymirr);
    if (fabs(imag(zeta)) > 15.) {
      volt += qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qw[i] * log(std::abs(zterm.first));
    }
    // Mirror contribution from both the x and the y plane.
    zeta = m_zmult * std::complex<double>(xxmirr, yymirr);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qw[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qw[i] * log(std::abs(zterm.first));
    }
  }
  return volt;
}

void ComponentAnalyticField::WfieldWireD10(const double xpos, const double ypos,
    double& ex, double& ey, const std::vector<double>& qw) const {
  //-----------------------------------------------------------------------
  //   IOND10 - Subroutine computing the signal on wire ISW due to a charge
  //            at (XPOS,YPOS). This is effectively routine EFCD10.
  //   VARIABLES : EX, EY, VOLT:Electric field and potential.
  //               ETOT, VOLT : Magnitude of electric field, potential.
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on  2/ 2/93.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    // Compute the contribution to the electric field.
    const auto wi = 1. / conj(zpos - zi) + zi / (m_cotube2 - conj(zpos) * zi);
    ex += qw[i] * real(wi);
    ey += qw[i] * imag(wi);
  }
}

double ComponentAnalyticField::WpotWireD10(const double xpos, const double ypos,
    const std::vector<double>& qw) const {
  double volt = 0.;
  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    volt -= qw[i] * log(std::abs(m_cotube * (zpos - zi) / 
                                 (m_cotube2 - zpos * conj(zi))));
  }
  return volt;
}

void ComponentAnalyticField::WfieldWireD30(const double xpos, const double ypos,
    double& ex, double& ey, const std::vector<double>& qw) const { 
  //-----------------------------------------------------------------------
  //   IOND30 - Subroutine computing the weighting field for a polygonal
  //            cells without periodicities, type D3.
  //   VARIABLES : EX, EY     :Electric field
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on 19/ 6/97.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  // Get the mapping of the position.
  std::complex<double> wpos, wdpos;
  ConformalMap(std::complex<double>(xpos, ypos) / m_cotube, wpos, wdpos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Compute the contribution to the electric field.
    const auto whelp = wdpos * (1. - std::norm(m_zw[i])) /
            ((wpos - m_zw[i]) * (1. - conj(m_zw[i]) * wpos));
    ex += qw[i] * real(whelp);
    ey -= qw[i] * imag(whelp);
  }
  ex /= m_cotube;
  ey /= m_cotube;
}

double ComponentAnalyticField::WpotWireD30(const double xpos, const double ypos,
    const std::vector<double>& qw) const {
  double volt = 0.;
  // Get the mapping of the position.
  std::complex<double> wpos, wdpos;
  ConformalMap(std::complex<double>(xpos, ypos) / m_cotube, wpos, wdpos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Compute the contribution to the potential.
    volt -= qw[i] * log(std::abs((wpos - m_zw[i]) / (1. - wpos * conj(m_zw[i]))));
  }
  return volt;
}

void ComponentAnalyticField::WfieldPlaneA00(
    const double xpos, const double ypos, double& ex, double& ey,
    const int mx, const int my, const std::vector<double>& qp) const {
  //-----------------------------------------------------------------------
  //   IPLA00 - Routine returning the A I,J [MX,MY] * E terms for A cells.
  //   VARIABLES : R2         : Potential before taking -Log(Sqrt(...))
  //               EX,EY      : x,y-Component of the electric field.
  //               EXHELP ETC : One term in the summing series.
  //               (XPOS,YPOS): Position where the field is needed.
  //   (Last changed on  9/11/98.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Define a few reduced variables.
    const double xx = xpos - m_w[i].x - mx * m_sx;
    const double yy = ypos - m_w[i].y - my * m_sy;
    // Calculate the field in case there are no planes.
    const double r2 = xx * xx + yy * yy;
    if (r2 <= 0.) continue;
    const double s2 = 1. / r2;
    double exhelp = xx * s2;
    double eyhelp = yy * s2;
    // Take care of a plane at constant x.
    const double xxmirr = m_ynplax ? xpos + m_w[i].x - 2 * m_coplax : 0.;
    if (m_ynplax) {
      const double r2plan = xxmirr * xxmirr + yy * yy;
      if (r2plan <= 0.) continue;
      const double s2plan = 1. / r2plan;
      exhelp -= xxmirr * s2plan;
      eyhelp -= yy * s2plan;
    }
    // Take care of a plane at constant y.
    const double yymirr = m_ynplay ? ypos + m_w[i].y - 2 * m_coplay : 0.;
    if (m_ynplay) {
      const double r2plan = xx * xx + yymirr * yymirr;
      if (r2plan <= 0.) continue;
      const double s2plan = 1. / r2plan;
      exhelp -= xx * s2plan;
      eyhelp -= yymirr * s2plan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double r2plan = xxmirr * xxmirr + yymirr * yymirr;
      if (r2plan <= 0.) continue;
      const double s2plan = 1. / r2plan;
      exhelp += xxmirr * s2plan;
      eyhelp += yymirr * s2plan;
    }
    ex += qp[i] * exhelp;
    ey += qp[i] * eyhelp;
  }
}

double ComponentAnalyticField::WpotPlaneA00(
    const double xpos, const double ypos, 
    const int mx, const int my, const std::vector<double>& qp) const {

  double volt = 0.;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Define a few reduced variables.
    const double xx = xpos - m_w[i].x - mx * m_sx;
    const double yy = ypos - m_w[i].y - my * m_sy;
    // Calculate the direct contribution.
    double r2 = xx * xx + yy * yy;
    if (r2 <= 0.) continue;
    // Take care of a plane at constant x.
    const double xxmirr = m_ynplax ? xpos + m_w[i].x - 2 * m_coplax : 0.;
    if (m_ynplax) {
      const double r2plan = xxmirr * xxmirr + yy * yy;
      if (r2plan <= 0.) continue;
      r2 /= r2plan;
    }
    // Take care of a plane at constant y.
    const double yymirr = m_ynplay ? ypos + m_w[i].y - 2 * m_coplay : 0.;
    if (m_ynplay) {
      const double r2plan = xx * xx + yymirr * yymirr;
      if (r2plan <= 0.) continue;
      r2 /= r2plan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double r2plan = xxmirr * xxmirr + yymirr * yymirr;
      if (r2plan <= 0.) continue;
      r2 *= r2plan;
    }
    // Calculate the potential.
    volt -= qp[i] * log(r2);
  }
  return 0.5 * volt;
}

void ComponentAnalyticField::WfieldPlaneB2X(
    const double xpos, const double ypos, double& ex, double& ey, 
    const int my, const std::vector<double>& qp) const { 
  //-----------------------------------------------------------------------
  //   IPLB2X - Routine calculating the MY contribution to the signal on
  //            wire IPLANE due to a charge at (XPOS,YPOS) for F-B2Y cells.
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z,ZZMIRR   : X + I*Y , XXMIRR + I*YYMIRR ; I**2=-1
  //               ECOMPL     : EX + I*EY                   ; I**2=-1
  //   (Last changed on 12/11/98.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  const double tx = HalfPi / m_sx;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = tx * (xpos - m_w[i].x);
    const double yy = tx * (ypos - m_w[i].y - my * m_sy);
    const double xxneg = tx * (xpos + m_w[i].x - 2 * m_coplan[0]);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    if (fabs(yy) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xxneg, yy);
      ecompl = -m_b2sin[i] / (sin(zz) * sin(zzneg));
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + m_w[i].y - 2 * m_coplay);
      if (fabs(yymirr) <= 20.) {
        const std::complex<double> zzmirr(yy, yymirr);
        const std::complex<double> zznmirr(xxneg, yymirr);
        ecompl += m_b2sin[i] / (sin(zzmirr) * sin(zznmirr));
      }
    }
    // Calculate the electric field.
    ex += qp[i] * real(ecompl);
    ey -= qp[i] * imag(ecompl);
  }
  ex *= tx;
  ey *= tx;
}

double ComponentAnalyticField::WpotPlaneB2X(
    const double xpos, const double ypos,
    const int my, const std::vector<double>& qp) const {

  double volt = 0.;
  const double tx = HalfPi / m_sx;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = tx * (xpos - m_w[i].x);
    const double yy = tx * (ypos - m_w[i].y - my * m_sy);
    const double xxneg = tx * (xpos + m_w[i].x - 2 * m_coplan[0]);
    // Calculate the direct contribution.
    double r2 = 1.;
    if (fabs(yy) <= 20.) {
      const double sinhy = sinh(yy);
      const double sinxx = sin(xx);
      const double sinxxneg = sin(xxneg);
      r2 = (sinhy * sinhy + sinxx * sinxx) /
           (sinhy * sinhy + sinxxneg * sinxxneg);
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + m_w[i].y - 2 * m_coplay);
      if (fabs(yymirr) <= 20.) {
        const double sinhy = sinh(yymirr);
        const double sinxx = sin(xx);
        const double sinxxneg = sin(xxneg);
        const double r2plan = (sinhy * sinhy + sinxx * sinxx) /
                              (sinhy * sinhy + sinxxneg * sinxxneg);
        r2 /= r2plan;
      }
    }
    volt -= qp[i] * log(r2);
  }
  return 0.5 * volt;
}

void ComponentAnalyticField::WfieldPlaneB2Y(
    const double xpos, const double ypos, double& ex, double& ey,
    const int mx, const std::vector<double>& qp) const { 
  //-----------------------------------------------------------------------
  //   IPLB2Y - Routine calculating the MX contribution to the signal on
  //            wire IPLANE due to a charge at (XPOS,YPOS) for F-B2X cells.
  //   VARIABLES : See routine EFCA00 for most of the variables.
  //               Z,ZZMIRR   : X + I*Y , XXMIRR + I*YYMIRR ; I**2=-1
  //               ECOMPL     : EX + I*EY                   ; I**2=-1
  //   (Last changed on 12/11/98.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  ex = ey = 0.;
  const double ty = HalfPi / m_sy;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = ty * (xpos - m_w[i].x - mx * m_sx);
    const double yy = ty * (ypos - m_w[i].y);
    const double yyneg = ty * (ypos + m_w[i].y - 2 * m_coplan[2]);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    if (fabs(xx) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xx, yyneg);
      ecompl = icons * m_b2sin[i] / (sin(icons * zz) * sin(icons * zzneg));
    }
    // Take care of a plane at constant y.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + m_w[i].x - 2 * m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const std::complex<double> zzmirr(xxmirr, yy);
        const std::complex<double> zznmirr(xxmirr, yyneg);
        ecompl -= m_b2sin[i] / (sin(icons * zzmirr) * sin(icons * zznmirr));
      }
    }
    ex += qp[i] * real(ecompl);
    ey -= qp[i] * imag(ecompl);
  }
  ex *= ty;
  ey *= ty;
}

double ComponentAnalyticField::WpotPlaneB2Y(
    const double xpos, const double ypos,
    const int mx, const std::vector<double>& qp) const {
  double volt = 0.;
  const double ty = HalfPi / m_sy;
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = ty * (xpos - m_w[i].x - mx * m_sx);
    const double yy = ty * (ypos - m_w[i].y);
    const double yyneg = ty * (ypos + m_w[i].y - 2 * m_coplan[2]);
    // Calculate the direct contribution.
    double r2 = 1.;
    if (fabs(xx) <= 20.) {
      const double sinhx = sinh(xx);
      const double sinyy = sin(yy);
      const double sinyyneg = sin(yyneg);
      r2 = (sinhx * sinhx + sinyy * sinyy) /
           (sinhx * sinhx + sinyyneg * sinyyneg);
    }
    // Take care of a plane at constant y.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + m_w[i].x - 2 * m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const double sinhx = sinh(xxmirr);
        const double sinyy = sin(yy);
        const double sinyyneg = sin(yyneg);
        const double r2plan = (sinhx * sinhx + sinyy * sinyy) /
                              (sinhx * sinhx + sinyyneg * sinyyneg);
        r2 /= r2plan;
      }
    }
    volt -= qp[i] * log(r2);
  }
  return 0.5 * volt;
}

void ComponentAnalyticField::WfieldPlaneC2X(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<double>& qp) const {
  //-----------------------------------------------------------------------
  //   IPLC2X - Routine returning the potential and electric field in a
  //            configuration with 2 x planes and y periodicity.
  //   VARIABLES : see the writeup
  //   (Last changed on 12/10/06.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  double s = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (imag(zeta) > +15.) {
      wsum1 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum1 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += qp[i] * (zterm.second / zterm.first);
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    // Constant terms sum
    s += qp[i] * (m_w[i].x - cx);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(2. * cx - xpos - m_w[i].x, yy);
    if (imag(zeta) > 15.) {
      wsum2 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += qp[i] * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2));
  ey = -imag(m_zmult * (wsum1 - wsum2));
  // Constant correction terms.
  if (m_mode == 0) ex += s * TwoPi / (m_sx * m_sy);
}

double ComponentAnalyticField::WpotPlaneC2X(
    const double xpos, const double ypos, 
    const std::vector<double>& qp) const {
  double volt = 0.;
  const double c0 = m_mode == 0 ? TwoPi / (m_sx * m_sy) : 0.;  
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qp[i] * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - m_w[i].x) / m_sx);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(2. * cx - xpos - m_w[i].x, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt += qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qp[i] * log(std::abs(zterm.first));
    }
    if (m_mode == 0) {
      volt -= c0 * qp[i] * (xpos - cx) * (m_w[i].x - cx);
    }
  }
  return volt;
}

void ComponentAnalyticField::WfieldPlaneC2Y(
    const double xpos, const double ypos, double& ex, double& ey, 
    const std::vector<double>& qp) const {
  //-----------------------------------------------------------------------
  //   IPLC2Y - Routine returning the potential and electric field in a
  //            configuration with 2 y planes and x periodicity.
  //   VARIABLES : see the writeup
  //   (Last changed on 12/10/06.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  double s = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (imag(zeta) > +15.) {
      wsum1 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum1 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += qp[i] * (zterm.second / zterm.first);
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    // Constant terms sum
    s += qp[i] * (m_w[i].y - cy);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(xx, 2. * cy - ypos - m_w[i].y);
    if (imag(zeta) > 15.) {
      wsum2 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += qp[i] * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 - wsum2));
  ey = -imag(m_zmult * (wsum1 + wsum2));
  // Constant correction terms.
  if (m_mode == 1) ey += s * TwoPi / (m_sx * m_sy);
}

double ComponentAnalyticField::WpotPlaneC2Y(
    const double xpos, const double ypos, 
    const std::vector<double>& qp) const {
  double volt = 0.;
  const double c1 = m_mode == 1 ? TwoPi / (m_sx * m_sy) : 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qp[i] * log(std::abs(zterm.first));
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - m_w[i].y) / m_sy);
    // Mirror contribution.
    zeta = m_zmult * std::complex<double>(xx, 2. * cy - ypos - m_w[i].y);
    if (fabs(imag(zeta)) > 15.) {
      volt += qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qp[i] * log(std::abs(zterm.first));
    }
    // Correct the voltage, if needed (MODE).
    if (m_mode == 1) {
      volt -= c1 * qp[i] * (ypos - cy) * (m_w[i].y - cy);
    }
  }
  return volt;
}

void ComponentAnalyticField::WfieldPlaneC30(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<double>& qp) const {
  //-----------------------------------------------------------------------
  //   IPLC30 - Routine returning the weighting field field in a
  //            configuration with 2 y and 2 x planes. This routine is
  //            basically the same as EFCC30.
  //   (Last changed on  9/11/98.)
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  std::complex<double> wsum3 = 0.;
  std::complex<double> wsum4 = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y;
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (imag(zeta) > +15.) {
      wsum1 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum1 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum1 += qp[i] * zterm.second / zterm.first;
    }
    // Mirror contribution from the x plane.
    const double xxmirr = MirrorCoordinate(xpos, m_coplax, m_w[i].x, m_sx);
    zeta = m_zmult * std::complex<double>(xxmirr, yy);
    if (imag(zeta) > 15.) {
      wsum2 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += qp[i] * zterm.second / zterm.first;
    }
    // Mirror contribution from the y plane.
    const double yymirr = MirrorCoordinate(ypos, m_coplay, m_w[i].y, m_sy);
    zeta = m_zmult * std::complex<double>(xx, yymirr);
    if (imag(zeta) > 15.) {
      wsum3 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum3 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum3 += qp[i] * zterm.second / zterm.first;
    }
    // Mirror contribution from both the x and the y plane.
    zeta = m_zmult * std::complex<double>(xxmirr, yymirr);
    if (imag(zeta) > 15.) {
      wsum4 -= qp[i] * icons;
    } else if (imag(zeta) < -15.) {
      wsum4 += qp[i] * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum4 += qp[i] * zterm.second / zterm.first;
    }
  }
  ex = real(m_zmult * (wsum1 + wsum2 - wsum3 - wsum4));
  ey = -imag(m_zmult * (wsum1 - wsum2 + wsum3 - wsum4));
}

double ComponentAnalyticField::WpotPlaneC30(
    const double xpos, const double ypos,
    const std::vector<double>& qp) const {
  double volt = 0.;
  // Wire loop.
  for (size_t i = 0; i < m_nWires; ++i) {
    const double xx = xpos - m_w[i].x;
    const double yy = ypos - m_w[i].y; 
    // Compute the direct contribution.
    auto zeta = m_zmult * std::complex<double>(xx, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qp[i] * log(std::abs(zterm.first));
    }
    // Mirror contribution from the x plane.
    const double xxmirr = MirrorCoordinate(xpos, m_coplax, m_w[i].x, m_sx);
    zeta = m_zmult * std::complex<double>(xxmirr, yy);
    if (fabs(imag(zeta)) > 15.) {
      volt += qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qp[i] * log(std::abs(zterm.first));
    }
    // Mirror contribution from the y plane.
    const double yymirr = MirrorCoordinate(ypos, m_coplay, m_w[i].y, m_sy);
    zeta = m_zmult * std::complex<double>(xx, yymirr);
    if (fabs(imag(zeta)) > 15.) {
      volt += qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt += qp[i] * log(std::abs(zterm.first));
    }
    // Mirror contribution from both the x and the y plane.
    zeta = m_zmult * std::complex<double>(xxmirr, yymirr);
    if (fabs(imag(zeta)) > 15.) {
      volt -= qp[i] * (fabs(imag(zeta)) - CLog2);
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      volt -= qp[i] * log(std::abs(zterm.first));
    }
  }
  return volt;
}

void ComponentAnalyticField::WfieldPlaneD10(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<double>& qp) const {
  //-----------------------------------------------------------------------
  //   IPLD10 - Subroutine computing the signal on wire IPLANE due to a
  //            charge at (XPOS,YPOS). This is effectively routine EFCD10.
  //   VARIABLES : EX, EY     : Electric field.
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on  9/11/98.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    // Compute the contribution to the electric field.
    const auto wi = 1. / conj(zpos - zi) + zi / (m_cotube2 - conj(zpos) * zi);
    ex += qp[i] * real(wi);
    ey += qp[i] * imag(wi);
  }
}

double ComponentAnalyticField::WpotPlaneD10(
    const double xpos, const double ypos,
    const std::vector<double>& qp) const {
  double volt = 0.;
  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(m_w[i].x, m_w[i].y);
    volt -= qp[i] * log(
        std::abs(m_cotube * (zpos - zi) / (m_cotube2 - zpos * conj(zi))));
  }
  return volt;
}

void ComponentAnalyticField::WfieldPlaneD30(
    const double xpos, const double ypos, double& ex, double& ey, 
    const std::vector<double>& qp) const {
  //-----------------------------------------------------------------------
  //   IPLD30 - Subroutine computing the weighting field for a polygonal
  //            cells without periodicities, type D3.
  //   VARIABLES : EX, EY     : Electric field
  //               (XPOS,YPOS): The position where the field is calculated.
  //               ZI, ZPOS   : Shorthand complex notations.
  //   (Last changed on  9/11/98.)
  //-----------------------------------------------------------------------
  ex = ey = 0.;
  // Get the mapping of the position.
  std::complex<double> wpos, wdpos;
  ConformalMap(std::complex<double>(xpos, ypos) / m_cotube, wpos, wdpos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    // Compute the contribution to the electric field.
    const auto whelp = wdpos * (1. - std::norm(m_zw[i])) /
            ((wpos - m_zw[i]) * (1. - conj(m_zw[i]) * wpos));
    ex += qp[i] * real(whelp);
    ey -= qp[i] * imag(whelp);
  }
  ex /= m_cotube;
  ey /= m_cotube;
}

double ComponentAnalyticField::WpotPlaneD30(
    const double xpos, const double ypos,
    const std::vector<double>& qp) const {
  double volt = 0.;
  // Get the mapping of the position.
  std::complex<double> wpos, wdpos;
  ConformalMap(std::complex<double>(xpos, ypos) / m_cotube, wpos, wdpos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    volt -= qp[i] * log(
        std::abs((wpos - m_zw[i]) / (1. - wpos * conj(m_zw[i]))));
  }
  return volt;
}

void ComponentAnalyticField::WfieldStrip(const double x, const double y, 
                                         const double g, const double w, 
                                         double& fx, double& fy) const {

  // Define shorthand notations.
  const double invg = 1. / g;
  const double s = sin(Pi * y * invg);
  const double c = cos(Pi * y * invg);
  const double s2 = s * s;
  // Evaluate the field.
  const double a1 = Pi * (w - x) * invg;
  if (a1 < 500.) {
    const double e1 = exp(a1);
    // Check for singularities.
    if (c == e1) return;
    const double d1 = c - e1;
    const double t1 = 1. / (s2 + d1 * d1);
    fx += e1 * t1;
    fy -= (1. - c * e1) * t1;
  }
  const double a2 = -Pi * (w + x) * invg;
  if (a2 < 500.) {
    const double e2 = exp(a2);
    // Check for singularities.
    if (c == e2) return;
    const double d2 = c - e2;
    const double t2 = 1. / (s2 + d2 * d2);
    fx -= e2 * t2;
    fy += (1. - c * e2) * t2;
  }
  fx *= s * invg;
  fy *= invg;

}

void ComponentAnalyticField::WfieldStripZ(
    const double xpos, const double ypos, double& ex, double& ey, 
    const int ip, const Strip& strip) const {
  //-----------------------------------------------------------------------
  //   IONEST - Weighting field for strips.
  //   (Last changed on  6/12/00.)
  //-----------------------------------------------------------------------

  // Initialise the weighting field.
  ex = ey = 0.;

  // Transform to normalised coordinates.
  double xw = 0., yw = 0.;
  switch (ip) {
    case 0:
      xw = -ypos + 0.5 * (strip.smin + strip.smax);
      yw = xpos - m_coplan[ip];
      break;
    case 1:
      xw = ypos - 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - xpos;
      break;
    case 2:
      xw = xpos - 0.5 * (strip.smin + strip.smax);
      yw = ypos - m_coplan[ip];
      break;
    case 3:
      xw = -xpos + 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - ypos;
      break;
    default:
      return;
  }

  // Make sure we are in the fiducial part of the weighting map.
  if (yw <= 0. || yw > strip.gap) return;

  double fx = 0., fy = 0.;
  WfieldStrip(xw, yw, strip.gap, 0.5 * fabs(strip.smax - strip.smin), fx, fy);

  // Rotate the field back to the original coordinates.
  switch (ip) {
    case 0:
      ex = fy;
      ey = -fx;
      break;
    case 1:
      ex = -fy;
      ey = fx;
      break;
    case 2:
      ex = fx;
      ey = fy;
      break;
    case 3:
      ex = -fx;
      ey = -fy;
      break;
  }
}

double ComponentAnalyticField::WpotStripZ(
    const double xpos, const double ypos, 
    const int ip, const Strip& strip) const {

  // Transform to normalised coordinates.
  double xw = 0., yw = 0.;
  switch (ip) {
    case 0:
      xw = -ypos + 0.5 * (strip.smin + strip.smax);
      yw = xpos - m_coplan[ip];
      break;
    case 1:
      xw = ypos - 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - xpos;
      break;
    case 2:
      xw = xpos - 0.5 * (strip.smin + strip.smax);
      yw = ypos - m_coplan[ip];
      break;
    case 3:
      xw = -xpos + 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - ypos;
      break;
    default:
      return 0.;
  }

  // Make sure we are in the fiducial part of the weighting map.
  if (yw <= 0. || yw > strip.gap) return 0.;

  // Define shorthand notations.
  const double a = Pi / strip.gap;
  const double c = cos(a * yw);
  // Strip halfwidth.
  const double w = 0.5 * fabs(strip.smax - strip.smin);
  const double e1 = exp(a * (w - xw));
  const double e2 = exp(-a * (w + xw));
  // Check for singularities.
  if (c == e1 || c == e2) return 0.;
  const double invs = 1. / sin(a * yw);
  constexpr double invPi = 1. / Pi;
  return (atan((c - e2) * invs) - atan((c - e1) * invs)) * invPi;
}

void ComponentAnalyticField::WfieldStripXy(const double xpos, const double ypos,
                                           const double zpos, double& ex,
                                           double& ey, double& ez,
                                           const int ip, 
                                           const Strip& strip) const {
  //-----------------------------------------------------------------------
  //   IONEST - Weighting field for strips.
  //   (Last changed on  6/12/00.)
  //-----------------------------------------------------------------------

  // Initialise the weighting field.
  ex = ey = ez = 0.;

  // Transform to normalised coordinates.
  double xw = 0., yw = 0.;
  switch (ip) {
    case 0:
      xw = -zpos + 0.5 * (strip.smin + strip.smax);
      yw = xpos - m_coplan[ip];
      break;
    case 1:
      xw = zpos - 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - xpos;
      break;
    case 2:
      xw = zpos - 0.5 * (strip.smin + strip.smax);
      yw = ypos - m_coplan[ip];
      break;
    case 3:
      xw = -zpos + 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - ypos;
      break;
    default:
      return;
  }

  // Make sure we are in the fiducial part of the weighting map.
  if (yw <= 0. || yw > strip.gap) return;

  double fx = 0., fy = 0.;
  WfieldStrip(xw, yw, strip.gap, 0.5 * fabs(strip.smax - strip.smin), fx, fy);

  // Rotate the field back to the original coordinates.
  switch (ip) {
    case 0:
      ex = fy;
      ey = 0.;
      ez = -fx;
      break;
    case 1:
      ex = -fy;
      ey = 0.;
      ez = fx;
      break;
    case 2:
      ex = 0.;
      ey = fy;
      ez = fx;
      break;
    case 3:
      ex = 0.;
      ey = -fy;
      ez = -fx;
      break;
  }
}

double ComponentAnalyticField::WpotStripXy(const double xpos, const double ypos,
                                           const double zpos, 
                                           const int ip, const Strip& strip) const {
  // Transform to normalised coordinates.
  double xw = 0., yw = 0.;
  switch (ip) {
    case 0:
      xw = -zpos + 0.5 * (strip.smin + strip.smax);
      yw = xpos - m_coplan[ip];
      break;
    case 1:
      xw = zpos - 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - xpos;
      break;
    case 2:
      xw = zpos - 0.5 * (strip.smin + strip.smax);
      yw = ypos - m_coplan[ip];
      break;
    case 3:
      xw = -zpos + 0.5 * (strip.smin + strip.smax);
      yw = m_coplan[ip] - ypos;
      break;
    default:
      return 0.;
  }

  // Make sure we are in the fiducial part of the weighting map.
  if (yw <= 0. || yw > strip.gap) return 0.;

  // Define shorthand notations.
  const double a = Pi / strip.gap;
  const double c = cos(a * yw);
  // Strip halfwidth.
  const double w = 0.5 * fabs(strip.smax - strip.smin);
  const double e1 = exp(a * (w - xw));
  const double e2 = exp(-a * (w + xw));
  // Check for singularities.
  if (c == e1 || c == e2) return 0.;
  const double invs = 1. / sin(a * yw);
  constexpr double invPi = 1. / Pi;
  return (atan((c - e2) * invs) - atan((c - e1) * invs)) * invPi;
}

void ComponentAnalyticField::WfieldPixel(const double xpos, const double ypos,
                                         const double zpos, 
                                         double& ex, double& ey, double& ez,
                                         const int ip, 
                                         const Pixel& pixel) const {
  //-----------------------------------------------------------------------
  //   Weighting field for pixels.
  //-----------------------------------------------------------------------
  // W. Riegler, G. Aglieri Rinella,
  // Point charge potential and weighting field of a pixel or pad
  // in a plane condenser,
  // Nucl. Instr. Meth. A 767, 2014, 267 - 270
  // http://dx.doi.org/10.1016/j.nima.2014.08.044

  // Initialise the weighting field.
  ex = ey = ez = 0.;

  // Transform to standard coordinates.
  double x = 0., y = 0., z = 0.;

  // Pixel centre and widths.
  const double ps = 0.5 * (pixel.smin + pixel.smax);
  const double pz = 0.5 * (pixel.zmin + pixel.zmax);
  const double wx = pixel.smax - pixel.smin;
  const double wy = pixel.zmax - pixel.zmin;
  switch (ip) {
    case 0:
      x = ypos - ps;
      y = zpos - pz;
      z = xpos - m_coplan[ip];
      break;
    case 1:
      x = ypos - ps;
      y = -zpos + pz;
      z = -xpos + m_coplan[ip];
      break;
    case 2:
      x = xpos - ps;
      y = -zpos + pz;
      z = ypos - m_coplan[ip];
      break;
    case 3:
      x = xpos - ps;
      y = zpos - pz;
      z = -ypos + m_coplan[ip];
      break;
    default:
      return;
  }
  // If needed, rotate into place.
  const bool rot = fabs(pixel.sphi) > 1.e-9;
  if (rot) {
    const double xx = x;
    const double yy = y;
    x = pixel.cphi * xx + pixel.sphi * yy;
    y = -pixel.sphi * xx + pixel.cphi * yy;
  }
  // if (z < 0.) std::cerr << " z = " << z << std::endl;
  // Make sure we are in the fiducial part of the weighting map.
  // Commenting out this lines either breaks the simulation or the plot!
  // TODO!
  // if (z <= 0. || z > d) return;

  // Define shorthand notations and common terms.
  const double x1 = x - 0.5 * wx;
  const double x2 = x + 0.5 * wx;
  const double y1 = y - 0.5 * wy;
  const double y2 = y + 0.5 * wy;
  const double x1s = x1 * x1;
  const double x2s = x2 * x2;
  const double y1s = y1 * y1;
  const double y2s = y2 * y2;

  // Calculate number of terms needed to have sufficiently small error.
  const double maxError = 1.e-5;
  const double d = pixel.gap;
  const double d3 = d * d * d;
  const unsigned int nz = std::ceil(sqrt(wx * wy / (8 * Pi * d3 * maxError)));
  const unsigned int nx = std::ceil(sqrt(wy * z / (4 * Pi * d3 * maxError)));
  const unsigned int ny = std::ceil(sqrt(wx * z / (4 * Pi * d3 * maxError)));
  const unsigned int nn = std::max(ny, std::max(nx, nz));
  for (unsigned int i = 1; i <= nn; ++i) {
    const double u1 = 2 * i * d - z;
    const double u2 = 2 * i * d + z;
    const double u1s = u1 * u1;
    const double u2s = u2 * u2;
    const double u1x1y1 = sqrt(x1s + y1s + u1s);
    const double u1x1y2 = sqrt(x1s + y2s + u1s);
    const double u1x2y1 = sqrt(x2s + y1s + u1s);
    const double u1x2y2 = sqrt(x2s + y2s + u1s);
    const double u2x1y1 = sqrt(x1s + y1s + u2s);
    const double u2x1y2 = sqrt(x1s + y2s + u2s);
    const double u2x2y1 = sqrt(x2s + y1s + u2s);
    const double u2x2y2 = sqrt(x2s + y2s + u2s);

    if (i <= nx) {
      //-df/dx(x,y,2nd-z)
      ex -= u1 * y1 / ((u1s + x2s) * u1x2y1) -
            u1 * y1 / ((u1s + x1s) * u1x1y1) +
            u1 * y2 / ((u1s + x1s) * u1x1y2) - u1 * y2 / ((u1s + x2s) * u1x2y2);

      //-df/dx(x,y,2nd+z)
      ex += u2 * y1 / ((u2s + x2s) * u2x2y1) -
            u2 * y1 / ((u2s + x1s) * u2x1y1) +
            u2 * y2 / ((u2s + x1s) * u2x1y2) - u2 * y2 / ((u2s + x2s) * u2x2y2);
    }
    if (i <= ny) {
      //-df/dy(x,y,2nd-z)
      ey -= u1 * x1 / ((u1s + y2s) * u1x1y2) -
            u1 * x1 / ((u1s + y1s) * u1x1y1) +
            u1 * x2 / ((u1s + y1s) * u1x2y1) - u1 * x2 / ((u1s + y2s) * u1x2y2);

      //-df/dy(x,y,2nd+z)
      ey += u2 * x1 / ((u2s + y2s) * u2x1y2) -
            u2 * x1 / ((u2s + y1s) * u2x1y1) +
            u2 * x2 / ((u2s + y1s) * u2x2y1) - u2 * x2 / ((u2s + y2s) * u2x2y2);
    }
    if (i <= nz) {
      //-df/dz(x,y,2nd-z)
      ez += x1 * y1 * (x1s + y1s + 2 * u1s) /
                ((x1s + u1s) * (y1s + u1s) * u1x1y1) +
            x2 * y2 * (x2s + y2s + 2 * u1s) /
                ((x2s + u1s) * (y2s + u1s) * u1x2y2) -
            x1 * y2 * (x1s + y2s + 2 * u1s) /
                ((x1s + u1s) * (y2s + u1s) * u1x1y2) -
            x2 * y1 * (x2s + y1s + 2 * u1s) /
                ((x2s + u1s) * (y1s + u1s) * u1x2y1);

      //-df/dz(x,y,2nd+z)
      ez += x1 * y1 * (x1s + y1s + 2 * u2s) /
                ((x1s + u2s) * (y1s + u2s) * u2x1y1) +
            x2 * y2 * (x2s + y2s + 2 * u2s) /
                ((x2s + u2s) * (y2s + u2s) * u2x2y2) -
            x1 * y2 * (x1s + y2s + 2 * u2s) /
                ((x1s + u2s) * (y2s + u2s) * u2x1y2) -
            x2 * y1 * (x2s + y1s + 2 * u2s) /
                ((x2s + u2s) * (y1s + u2s) * u2x2y1);
    }
  }

  const double zs = z * z;
  const double x1y1 = sqrt(x1s + y1s + zs);
  const double x1y2 = sqrt(x1s + y2s + zs);
  const double x2y1 = sqrt(x2s + y1s + zs);
  const double x2y2 = sqrt(x2s + y2s + zs);
  //-df/dx(x,y,z)
  ex += z * y1 / ((zs + x2s) * x2y1) - z * y1 / ((zs + x1s) * x1y1) +
        z * y2 / ((zs + x1s) * x1y2) - z * y2 / ((zs + x2s) * x2y2);

  //-df/y(x,y,z)
  ey += z * x1 / ((zs + y2s) * x1y2) - z * x1 / ((zs + y1s) * x1y1) +
        z * x2 / ((zs + y1s) * x2y1) - z * x2 / ((zs + y2s) * x2y2);

  //-df/dz(x,y,z)
  ez += x1 * y1 * (x1s + y1s + 2 * zs) / ((x1s + zs) * (y1s + zs) * x1y1) +
        x2 * y2 * (x2s + y2s + 2 * zs) / ((x2s + zs) * (y2s + zs) * x2y2) -
        x1 * y2 * (x1s + y2s + 2 * zs) / ((x1s + zs) * (y2s + zs) * x1y2) -
        x2 * y1 * (x2s + y1s + 2 * zs) / ((x2s + zs) * (y1s + zs) * x2y1);

  constexpr double invTwoPi = 1. / TwoPi;
  ex *= invTwoPi;
  ey *= invTwoPi;
  ez *= invTwoPi;

  // Rotate the field back to the original coordinates.
  const double fx = rot ? pixel.cphi * ex - pixel.sphi * ey : ex; 
  const double fy = rot ? pixel.sphi * ex + pixel.cphi * ey : ey;
  const double fz = ez;
  switch (ip) {
    case 0:
      ex = fz;
      ey = fx;
      ez = fy;
      break;
    case 1:
      ex = -fz;
      ey = fx;
      ez = -fy;
      break;
    case 2:
      ex = fx;
      ey = fz;
      ez = -fy;
      break;
    case 3:
      ex = fx;
      ey = -fz;
      ez = fy;
      break;
  }
}

double ComponentAnalyticField::WpotPixel(const double xpos, const double ypos,
                                         const double zpos, 
                                         const int ip, const Pixel& pixel) const {
  // Transform to standard coordinates.
  double x = 0., y = 0., z = 0.;

  // Pixel centre and widths.
  const double ps = 0.5 * (pixel.smin + pixel.smax);
  const double pz = 0.5 * (pixel.zmin + pixel.zmax);
  const double wx = pixel.smax - pixel.smin;
  const double wy = pixel.zmax - pixel.zmin;
  switch (ip) {
    case 0:
      x = ypos - ps;
      y = zpos - pz;
      z = xpos - m_coplan[ip];
      break;
    case 1:
      x = ypos - ps;
      y = -zpos + pz;
      z = -xpos + m_coplan[ip];
      break;
    case 2:
      x = xpos - ps;
      y = -zpos + pz;
      z = ypos - m_coplan[ip];
      break;
    case 3:
      x = xpos - ps;
      y = zpos - pz;
      z = -ypos + m_coplan[ip];
      break;
    default:
      return 0.;
  }
  // If needed, rotate into place.
  if (fabs(pixel.sphi) > 1.e-9) {
    const double xx = x;
    const double yy = y;
    x = pixel.cphi * xx + pixel.sphi * yy;
    y = -pixel.sphi * xx + pixel.cphi * yy;
  }

  // Define shorthand notations and common terms.
  const double x1 = x - 0.5 * wx;
  const double x2 = x + 0.5 * wx;
  const double y1 = y - 0.5 * wy;
  const double y2 = y + 0.5 * wy;
  const double x1s = x1 * x1;
  const double x2s = x2 * x2;
  const double y1s = y1 * y1;
  const double y2s = y2 * y2;

  // Calculate number of terms needed to have sufficiently small error.
  const double maxError = 1.e-5;
  const double d = pixel.gap;
  const double d3 = d * d * d;
  const unsigned int nn = std::ceil(sqrt(wx * wy * z / (8 * Pi * d3 * maxError)));
  double volt = 0.;
  for (unsigned int i = 1; i <= nn; ++i) {
    const double u1 = 2 * i * d - z;
    const double u2 = 2 * i * d + z;
    const double u1s = u1 * u1;
    const double u2s = u2 * u2;
    const double u1x1y1 = sqrt(x1s + y1s + u1s);
    const double u1x1y2 = sqrt(x1s + y2s + u1s);
    const double u1x2y1 = sqrt(x2s + y1s + u1s);
    const double u1x2y2 = sqrt(x2s + y2s + u1s);
    const double u2x1y1 = sqrt(x1s + y1s + u2s);
    const double u2x1y2 = sqrt(x1s + y2s + u2s);
    const double u2x2y1 = sqrt(x2s + y1s + u2s);
    const double u2x2y2 = sqrt(x2s + y2s + u2s);

    volt -= atan(x1 * y1 / (u1 * u1x1y1)) + atan(x2 * y2 / (u1 * u1x2y2)) -
            atan(x1 * y2 / (u1 * u1x1y2)) - atan(x2 * y1 / (u1 * u1x2y1));
    volt += atan(x1 * y1 / (u2 * u2x1y1)) + atan(x2 * y2 / (u2 * u2x2y2)) -
            atan(x1 * y2 / (u2 * u2x1y2)) - atan(x2 * y1 / (u2 * u2x2y1));
  }

  const double zs = z * z;
  const double x1y1 = sqrt(x1s + y1s + zs);
  const double x1y2 = sqrt(x1s + y2s + zs);
  const double x2y1 = sqrt(x2s + y1s + zs);
  const double x2y2 = sqrt(x2s + y2s + zs);

  volt += atan(x1 * y1 / (z * x1y1)) + atan(x2 * y2 / (z * x2y2)) -
          atan(x1 * y2 / (z * x1y2)) - atan(x2 * y1 / (z * x2y1));
  constexpr double invTwoPi = 1. / TwoPi;
  return volt * invTwoPi;
}

void ComponentAnalyticField::FieldAtWireA00(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCA00 - Subroutine performing the actual field calculations in case
  //          only one charge and not more than 1 mirror-charge in either
  //          x or y is present.
  //          The potential used is 1/2*pi*eps0  log(r).
  // (Last changed on 27/ 1/96.)
  //-----------------------------------------------------------------------

  ex = ey = 0.;
  // Loop over all wires.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    // Calculate the field in case there are no planes.
    double exhelp = 0.;
    double eyhelp = 0.;
    const double xx = xpos - wire.x;
    const double yy = ypos - wire.y;
    if (cnalso[i]) {
      const double r2 = xx * xx + yy * yy;
      exhelp = xx / r2;
      eyhelp = yy / r2;
    }
    // Take care of a plane at constant x.
    double xxmirr = 0.;
    if (m_ynplax) {
      xxmirr = wire.x + xpos - 2 * m_coplax;
      const double r2plan = xxmirr * xxmirr + yy * yy;
      exhelp -= xxmirr / r2plan;
      eyhelp -= yy / r2plan;
    }
    // Take care of a plane at constant y.
    double yymirr = 0.;
    if (m_ynplay) {
      yymirr = wire.y + ypos - 2 * m_coplay;
      const double r2plan = xx * xx + yymirr * yymirr;
      exhelp -= xx / r2plan;
      eyhelp -= yymirr / r2plan;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double r2plan = xxmirr * xxmirr + yymirr * yymirr;
      exhelp += xxmirr / r2plan;
      eyhelp += yymirr / r2plan;
    }
    ex += wire.e * exhelp;
    ey += wire.e * eyhelp;
  }
}

void ComponentAnalyticField::FieldAtWireB1X(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCB1X - Routine calculating the potential for a row of positive
  //          charges. The potential used is Re(Log(sin pi/s (z-z0))).
  //-----------------------------------------------------------------------

  constexpr std::complex<double> icons(0., 1.);
  std::complex<double> ecompl;
  ex = ey = 0.;
  const double tx = Pi / m_sx;
  if (m_ynplay) {
    // With a y plane.
    for (unsigned int i = 0; i < m_nWires; ++i) {
      const auto& wire = m_w[i];
      const double xx = tx * (xpos - wire.x);
      const double yy = tx * (ypos - wire.y);
      if (!cnalso[i]) {
        ecompl = 0.;
      } else if (yy > 20.) {
        ecompl = -icons;
      } else if (yy < -20.) {
        ecompl = icons;
      } else {
        const std::complex<double> zz(xx, yy);
        const auto expzz = exp(2. * icons * zz);
        ecompl = icons * (expzz + 1.) / (expzz - 1.);
      }
      const double yymirr = tx * (ypos + wire.y - 2. * m_coplay);
      if (yymirr > 20.) {
        ecompl += icons;
      } else if (yymirr < -20.) {
        ecompl += -icons;
      } else {
        const std::complex<double> zzmirr(xx, yymirr);
        const auto expzzmirr = exp(2. * icons * zzmirr);
        ecompl += -icons * (expzzmirr + 1.) / (expzzmirr - 1.);
      }
      // Update the field.
      ex += wire.e * real(ecompl);
      ey -= wire.e * imag(ecompl);
    }
  } else {
    // Without a y plane.
    for (unsigned int i = 0; i < m_nWires; ++i) {
      if (!cnalso[i]) continue;
      const auto& wire = m_w[i];
      const double xx = tx * (xpos - wire.x);
      const double yy = tx * (ypos - wire.y);
      if (yy > 20.) {
        ecompl = -icons;
      } else if (yy < -20.) {
        ecompl = icons;
      } else {
        const std::complex<double> zz(xx, yy);
        const auto expzz = exp(2. * icons * zz);
        ecompl = icons * (expzz + 1.) / (expzz - 1.);
      }
      ex += wire.e * real(ecompl);
      ey -= wire.e * imag(ecompl);
    }
  }
  ex *= tx;
  ey *= tx;
}

void ComponentAnalyticField::FieldAtWireB1Y(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCB1Y - Routine calculating the potential for a row of positive
  //          charges. The potential used is Re(Log(sinh pi/sy(z-z0)).
  //-----------------------------------------------------------------------

  std::complex<double> ecompl;
  ex = ey = 0.;
  const double ty = Pi / m_sy;
  if (m_ynplax) {
    // With an x plane.
    for (unsigned int i = 0; i < m_nWires; ++i) {
      const auto& wire = m_w[i];
      const double xx = ty * (xpos - wire.x);
      const double yy = ty * (ypos - wire.y);
      if (!cnalso[i]) {
        ecompl = 0.;
      } else if (xx > 20.) {
        ecompl = 1.;
      } else if (xx < -20.) {
        ecompl = -1.;
      } else {
        const std::complex<double> zz(xx, yy);
        const auto expzz = exp(2. * zz);
        ecompl = (expzz + 1.) / (expzz - 1.);
      }
      const double xxmirr = ty * (xpos + wire.x - 2. * m_coplax);
      if (xxmirr > 20.) {
        ecompl -= 1.;
      } else if (xxmirr < -20.) {
        ecompl += 1.;
      } else {
        const std::complex<double> zzmirr(xxmirr, yy);
        const auto expzzmirr = exp(2. * zzmirr);
        ecompl -= (expzzmirr + 1.) / (expzzmirr - 1.);
      }
      ex += wire.e * real(ecompl);
      ey -= wire.e * imag(ecompl);
    }
  } else {
    // Without an x plane.
    for (unsigned int i = 0; i < m_nWires; ++i) {
      if (!cnalso[i]) continue;
      const auto& wire = m_w[i];
      const double xx = ty * (xpos - wire.x);
      const double yy = ty * (ypos - wire.y);
      if (xx > 20.) {
        ecompl = 1.;
      } else if (xx < -20.) {
        ecompl = -1.;
      } else {
        const std::complex<double> zz(xx, yy);
        const auto expzz = exp(2. * zz);
        ecompl = (expzz + 1.) / (expzz - 1.);
      }
      ex += wire.e * real(ecompl);
      ey -= wire.e * imag(ecompl);
    }
  }
  ex *= ty;
  ey *= ty;
}

void ComponentAnalyticField::FieldAtWireB2X(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCB2X - Routine calculating the potential for a row of alternating
  //          + - charges. The potential used is re log(sin pi/sx (z-z0))
  //-----------------------------------------------------------------------
  constexpr std::complex<double> icons(0., 1.);
  ex = ey = 0.;
  const double tx = HalfPi / m_sx;
  // Loop over all wires.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const double xx = tx * (xpos - m_w[i].x);
    const double yy = tx * (ypos - m_w[i].y);
    const double xxneg = tx * (xpos + m_w[i].x - 2 * m_coplax);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    if (cnalso[i] && fabs(yy) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xxneg, yy);
      ecompl = -m_b2sin[i] / (sin(zz) * sin(zzneg));
    } else if (fabs(yy) <= 20.) {
      const std::complex<double> zzneg(xxneg, yy);
      const auto expzzneg = exp(2. * icons * zzneg);
      ecompl = -icons * (expzzneg + 1.) / (expzzneg - 1.);
    }
    // Take care of planes at constant y.
    if (m_ynplay) {
      const double yymirr = tx * (ypos + m_w[i].y - 2 * m_coplay);
      if (fabs(yymirr) <= 20.) {
        const std::complex<double> zzmirr(xx, yymirr);
        const std::complex<double> zznmirr(xxneg, yymirr);
        ecompl += m_b2sin[i] / (sin(zzmirr) * sin(zznmirr));
      }
    }
    ex += m_w[i].e * real(ecompl);
    ey -= m_w[i].e * imag(ecompl);
  }
  ex *= tx;
  ey *= tx;
}

void ComponentAnalyticField::FieldAtWireB2Y(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCB2Y - Routine calculating the potential for a row of alternating
  //          + - charges. The potential used is re log(sin pi/sx (z-z0))
  //-----------------------------------------------------------------------
  constexpr std::complex<double> icons(0., 1.);
  ex = ey = 0.;
  const double ty = HalfPi / m_sy;
  // Loop over all wires.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const double xx = ty * (xpos - m_w[i].x);
    const double yy = ty * (ypos - m_w[i].y);
    const double yyneg = ty * (ypos + m_w[i].y - 2 * m_coplay);
    // Calculate the field in case there are no equipotential planes.
    std::complex<double> ecompl(0., 0.);
    if (cnalso[i] && fabs(xx) <= 20.) {
      const std::complex<double> zz(xx, yy);
      const std::complex<double> zzneg(xx, yyneg);
      ecompl = icons * m_b2sin[i] / (sin(icons * zz) * sin(icons * zzneg));
    } else if (fabs(xx) <= 20.) {
      const std::complex<double> zzneg(xx, yyneg);
      const auto expzzneg = exp(2. * zzneg);
      ecompl = -(expzzneg + 1.) / (expzzneg - 1.);
    }
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double xxmirr = ty * (xpos + m_w[i].x - 2 * m_coplax);
      if (fabs(xxmirr) <= 20.) {
        const std::complex<double> zzmirr(xxmirr, yy);
        const std::complex<double> zznmirr(xxmirr, yyneg);
        ecompl -=
            icons * m_b2sin[i] / (sin(icons * zzmirr) * sin(icons * zznmirr));
      }
    }
    ex += m_w[i].e * real(ecompl);
    ey -= m_w[i].e * imag(ecompl);
  }
  ex *= ty;
  ey *= ty;
}

void ComponentAnalyticField::FieldAtWireC10(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCC10 - Routine returning the potential and electric field. It
  //          calls the routines PH2 and E2SUM written by G.A.Erskine.
  //-----------------------------------------------------------------------
  constexpr std::complex<double> icons(0., 1.);
  std::complex<double> wsum(0., 0.);
  // Loop over the wires.
  for (unsigned int j = 0; j < m_nWires; ++j) {
    if (!cnalso[j]) continue;
    const auto& wire = m_w[j];
    auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum += wire.e * (zterm.second / zterm.first);
    }
  }
  ex = -real(-m_zmult * wsum);
  ey = imag(-m_zmult * wsum);
  if (m_mode == 0) ex -= m_c1;
  if (m_mode == 1) ey -= m_c1;
}

void ComponentAnalyticField::FieldAtWireC2X(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCC2X - Routine returning the potential and electric field in a
  //          configuration with 2 x planes and y periodicity.
  //-----------------------------------------------------------------------
  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    if (cnalso[i]) {
      auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
      if (imag(zeta) > 15.) {
        wsum1 -= wire.e * icons;
      } else if (imag(zeta) < -15.) {
        wsum1 += wire.e * icons;
      } else {
        const auto zterm = Th1(zeta, m_p1, m_p2);
        wsum1 += wire.e * (zterm.second / zterm.first);
      }
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - wire.x) / m_sx);
    // Mirror contribution.
    auto zeta =
        m_zmult * std::complex<double>(2. * cx - xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum2 -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += wire.e * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2));
  ey = -imag(m_zmult * (wsum1 - wsum2));
  // Constant correction terms.
  if (m_mode == 0) ex -= m_c1;
}

void ComponentAnalyticField::FieldAtWireC2Y(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCC2Y - Routine returning the potential and electric field in a
  //          configuration with 2 y planes and x periodicity.
  //-----------------------------------------------------------------------
  const std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    if (cnalso[i]) {
      // Compute the direct contribution.
      auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
      if (imag(zeta) > 15.) {
        wsum1 -= wire.e * icons;
      } else if (imag(zeta) < -15.) {
        wsum1 += wire.e * icons;
      } else {
        const auto zterm = Th1(zeta, m_p1, m_p2);
        wsum1 += wire.e * (zterm.second / zterm.first);
      }
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - wire.y) / m_sy);
    // Mirror contribution from the y plane.
    auto zeta =
        m_zmult * std::complex<double>(xpos - wire.x, 2 * cy - ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum2 -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += wire.e * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 - wsum2));
  ey = -imag(m_zmult * (wsum1 + wsum2));
  // Constant correction terms.
  if (m_mode == 1) ey -= m_c1;
}

void ComponentAnalyticField::FieldAtWireC30(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCC30 - Routine returning the potential and electric field in a
  //          configuration with 2 y and 2 x planes.
  //-----------------------------------------------------------------------
  constexpr std::complex<double> icons(0., 1.);
  // Initial values.
  std::complex<double> wsum1 = 0.;
  std::complex<double> wsum2 = 0.;
  std::complex<double> wsum3 = 0.;
  std::complex<double> wsum4 = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    if (cnalso[i]) {
      auto zeta = m_zmult * std::complex<double>(xpos - wire.x, ypos - wire.y);
      if (imag(zeta) > 15.) {
        wsum1 -= wire.e * icons;
      } else if (imag(zeta) < -15.) {
        wsum1 += wire.e * icons;
      } else {
        const auto zterm = Th1(zeta, m_p1, m_p2);
        wsum1 += wire.e * (zterm.second / zterm.first);
      }
    }
    // Find the plane nearest to the wire.
    const double cx = m_coplax - m_sx * round((m_coplax - wire.x) / m_sx);
    // Mirror contribution from the x plane.
    auto zeta =
        m_zmult * std::complex<double>(2. * cx - xpos - wire.x, ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum2 -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum2 += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum2 += wire.e * (zterm.second / zterm.first);
    }
    // Find the plane nearest to the wire.
    const double cy = m_coplay - m_sy * round((m_coplay - wire.y) / m_sy);
    // Mirror contribution from the x plane.
    zeta =
        m_zmult * std::complex<double>(xpos - wire.x, 2. * cy - ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum3 -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum3 += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum3 += wire.e * (zterm.second / zterm.first);
    }
    // Mirror contribution from both the x and the y plane.
    zeta = m_zmult * std::complex<double>(2. * cx - xpos - wire.x,
                                          2. * cy - ypos - wire.y);
    if (imag(zeta) > 15.) {
      wsum4 -= wire.e * icons;
    } else if (imag(zeta) < -15.) {
      wsum4 += wire.e * icons;
    } else {
      const auto zterm = Th1(zeta, m_p1, m_p2);
      wsum4 += wire.e * (zterm.second / zterm.first);
    }
  }
  // Convert the two contributions to a real field.
  ex = real(m_zmult * (wsum1 + wsum2 - wsum3 - wsum4));
  ey = -imag(m_zmult * (wsum1 - wsum2 + wsum3 - wsum4));
}

void ComponentAnalyticField::FieldAtWireD10(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCD10 - Subroutine performing the actual field calculations for a
  //          cell which has a one circular plane and some wires.
  //-----------------------------------------------------------------------
  ex = ey = 0.;
  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  // Loop over all wires.
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(wire.x, wire.y);
    // First the case that the wire has to be taken fully.
    if (cnalso[i]) {
      const std::complex<double> wi =
          1. / conj(zpos - zi) + zi / (m_cotube2 - conj(zpos) * zi);
      ex += wire.e * real(wi);
      ey += wire.e * imag(wi);
    } else {
      const std::complex<double> wi = zi / (m_cotube2 - conj(zpos) * zi);
      ex += wire.e * real(wi);
      ey += wire.e * imag(wi);
    }
  }
}

void ComponentAnalyticField::FieldAtWireD20(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCD20 - Subroutine performing the actual field calculations for a
  //          cell which has a tube and phi periodicity.
  //-----------------------------------------------------------------------
  ex = ey = 0.;
  // Set the complex position coordinates.
  const std::complex<double> zpos(xpos, ypos);
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    // Set the complex version of the wire-coordinate for simplicity.
    const std::complex<double> zi(wire.x, wire.y);
    if (cnalso[i]) {
      // Case of the wire which is not in the centre.
      if (std::abs(zi) > wire.r) {
        const std::complex<double> wi =
            double(m_mtube) * pow(conj(zpos), m_mtube - 1) *
            (1. / conj(pow(zpos, m_mtube) - pow(zi, m_mtube)) +
             pow(zi, m_mtube) /
                 (pow(m_cotube, 2 * m_mtube) - pow(conj(zpos) * zi, m_mtube)));
        ex += wire.e * real(wi);
        ey += wire.e * imag(wi);
      } else {
        const std::complex<double> wi =
            1. / conj(zpos - zi) + zi / (m_cotube2 - conj(zpos) * zi);
        ex += wire.e * real(wi);
        ey += wire.e * imag(wi);
      }
    } else {
      if (std::abs(zi) > wire.r) {
        const std::complex<double> wi =
            double(m_mtube) * pow(conj(zpos), m_mtube - 1) *
            (pow(zi, m_mtube) /
             (pow(m_cotube, 2 * m_mtube) - pow(conj(zpos) * zi, m_mtube)));
        ex += wire.e * real(wi);
        ey += wire.e * imag(wi);
      } else {
        const std::complex<double> wi = zi / (m_cotube2 - conj(zpos) * zi);
        ex += wire.e * real(wi);
        ey += wire.e * imag(wi);
      }
    }
  }
}

void ComponentAnalyticField::FieldAtWireD30(
    const double xpos, const double ypos, double& ex, double& ey,
    const std::vector<bool>& cnalso) const {
  //-----------------------------------------------------------------------
  // FFCD30 - Subroutine performing the actual field calculations for a
  //          cell which has a polygon as tube and some wires.
  //-----------------------------------------------------------------------
  ex = ey = 0.;
  // Get the mapping of the position.
  std::complex<double> wpos, wdpos;
  ConformalMap(std::complex<double>(xpos, ypos) / m_cotube, wpos, wdpos);
  // Loop over all wires.
  for (size_t i = 0; i < m_nWires; ++i) {
    if (cnalso[i]) {
      // Full contribution.
      const std::complex<double> whelp =
          wdpos * (1. - std::norm(m_zw[i])) /
          ((wpos - m_zw[i]) * (1. - conj(m_zw[i]) * wpos));
      ex += m_w[i].e * real(whelp);
      ey -= m_w[i].e * imag(whelp);
    } else {
      // Mirror charges only.
      const std::complex<double> whelp =
          wdpos * conj(m_zw[i]) / (1. - conj(m_zw[i]) * wpos);
      ex += m_w[i].e * real(whelp);
      ey -= m_w[i].e * imag(whelp);
    }
  }
  ex /= m_cotube;
  ey /= m_cotube;
}

bool ComponentAnalyticField::SagDetailed(
    const Wire& wire, const std::vector<double>& xMap,
    const std::vector<double>& yMap,
    const std::vector<std::vector<double> >& fxMap,
    const std::vector<std::vector<double> >& fyMap, std::vector<double>& csag,
    std::vector<double>& xsag, std::vector<double>& ysag) const {
  //-----------------------------------------------------------------------
  //   OPTSAG - Computes the wire sag due to electrostatic and gravitational
  //            forces, using a Runge-Kutta-Nystrom multiple shoot method,
  //            where the intermediate conditions are imposed through a
  //            Broyden rank-1 zero search.
  //-----------------------------------------------------------------------

  csag.clear();
  xsag.clear();
  ysag.clear();
  const unsigned int np = m_nSteps * (m_nShots + 1);
  // Compute the step width.
  const double h = wire.u / np;
  // Compute expected maximum sag, constant-force approximation.
  std::array<double, 2> xst = {0., 0.};
  std::array<double, 2> dxst = {0., 0.};
  double fxmean = 0.;
  double fymean = 0.;
  // Loop over the whole wire.
  for (unsigned int i = 0; i <= np; ++i) {
    const double z = i * h;
    std::array<double, 2> force = {0., 0.};
    if (!GetForceRatio(wire, z, xst, dxst, force, xMap, yMap, fxMap, fyMap)) {
      std::cerr << m_className << "::SagDetailed:\n"
                << "    Wire at nominal position outside scanning area.\n";
      return false;
    }
    fxmean += force[0];
    fymean += force[1];
  }
  const double u2 = wire.u * wire.u;
  // Compute expected sag.
  const double s = u2 / (8. * (1. + np));
  double sagx0 = -fxmean * s;
  double sagy0 = -fymean * s;
  if (m_debug) {
    std::cout << m_className << "::SagDetailed: Parabolic sag.\n";
    std::printf("    dx = %12.5e, dy = %12.5e [cm]\n", sagx0, sagy0);
  }
  // Starting position: parabolic sag.
  std::vector<double> xx(4 * m_nShots + 2);
  // Derivative first point.
  xx[0] = 4 * sagx0 / wire.u;
  xx[1] = 4 * sagy0 / wire.u;
  // Intermediate points, both position and derivative.
  for (unsigned int i = 1; i <= m_nShots; ++i) {
    // Position along the wire.
    const double z = -0.5 * wire.u + i * m_nSteps * h;
    const unsigned int k = 4 * i - 2;
    // Deflection.
    const double f = 1. - 4 * z * z / u2;
    xx[k] = sagx0 * f;
    xx[k + 1] = sagy0 * f;
    // Derivative.
    const double fp = -8 * z / u2;
    xx[k + 2] = sagx0 * fp;
    xx[k + 3] = sagy0 * fp;
  }
  // Search for solution.
  if (!FindZeroes(wire, h, xx, xMap, yMap, fxMap, fyMap)) {
    std::cerr << m_className << "::SagDetailed:\n"
              << "    Failed to solve the differential equation for the sag.\n";
    return false;
  }
  // And return the detailed solution, first the starting point.
  csag.assign(np + 1, 0.);
  xsag.assign(np + 1, 0.);
  ysag.assign(np + 1, 0.);
  csag[0] = -0.5 * wire.u;
  double coor = -0.5 * wire.u;
  for (unsigned int i = 0; i <= m_nShots; ++i) {
    // Set the starting value and starting derivative.
    if (i == 0) {
      xst[0] = 0;
      xst[1] = 0;
      dxst[0] = xx[0];
      dxst[1] = xx[1];
    } else {
      xst[0] = xx[4 * i - 2];
      xst[1] = xx[4 * i - 1];
      dxst[0] = xx[4 * i];
      dxst[1] = xx[4 * i + 1];
    }
    // Store the intermediate values.
    for (unsigned int j = 1; j <= m_nSteps; ++j) {
      StepRKN(wire, h, coor, xst, dxst, xMap, yMap, fxMap, fyMap);
      csag[i * m_nSteps + j] = coor;
      xsag[i * m_nSteps + j] = xst[0];
      ysag[i * m_nSteps + j] = xst[1];
    }
  }
  // Seems to have worked.
  return true;
}

bool ComponentAnalyticField::GetForceRatio(
    const Wire& wire, const double /*coor*/, const std::array<double, 2>& bend,
    const std::array<double, 2>& /*dbend*/, std::array<double, 2>& f,
    const std::vector<double>& xMap, const std::vector<double>& yMap,
    const std::vector<std::vector<double> >& fxMap,
    const std::vector<std::vector<double> >& fyMap) const {
  //-----------------------------------------------------------------------
  //    OPTSTP - Returns the electrostatic and gravitational force divided
  //             by the stretching force acting on a wire at position COOR,
  //             with deflection BEND and bending derivative DBEND.
  //-----------------------------------------------------------------------

  // TODO: COOR and DBEND don't seem to be used?

  // Initialise the forces.
  f.fill(0.);

  const double xw = wire.x + bend[0];
  const double yw = wire.y + bend[1];
  if (m_useElectrostaticForce) {
    // Electrostatic force.
    if (xMap.empty() || yMap.empty() || fxMap.empty() || fyMap.empty()) {
      return false;
    }
    // In case extrapolation is not permitted, check range.
    if (!m_extrapolateForces) {
      if ((xMap.front() - xw) * (xw - xMap.back()) < 0 ||
          (yMap.front() - yw) * (yw - yMap.back()) < 0) {
        return false;
      }
    }
    // Interpolation order.
    constexpr int order = 2;
    // Electrostatic force: interpolate the table, first along the y-lines.
    const unsigned int nX = xMap.size();
    const unsigned int nY = yMap.size();
    std::vector<double> xaux(nX, 0.);
    std::vector<double> yaux(nX, 0.);
    for (unsigned int i = 0; i < nX; ++i) {
      xaux[i] = Numerics::Divdif(fxMap[i], yMap, nY, yw, order);
      yaux[i] = Numerics::Divdif(fyMap[i], yMap, nY, yw, order);
    }
    // Then along the x-lines.
    f[0] += Numerics::Divdif(xaux, xMap, nX, xw, order);
    f[1] += Numerics::Divdif(yaux, xMap, nX, xw, order);
  }
  // Add the gravity term.
  if (m_useGravitationalForce) {
    // Mass per unit length [kg / cm]. 
    const double m = 1.e-3 * wire.density * Pi * wire.r * wire.r;
    f[0] -= m_down[0] * m * GravitationalAcceleration;
    f[1] -= m_down[1] * m * GravitationalAcceleration;
  }
  // Divide by the stretching force.
  const double s = 1000. / (GravitationalAcceleration * wire.tension);
  f[0] *= s;
  f[1] *= s;
  return true;
}

// Subroutine subprograms RRKNYS and DRKNYS advance the solution of the system
// of n >= 1 second-order differential equations
// by a single step of length h in the independent variable x.
bool ComponentAnalyticField::StepRKN(
    const Wire& wire, const double h, double& x, std::array<double, 2>& y,
    std::array<double, 2>& yp, const std::vector<double>& xMap,
    const std::vector<double>& yMap,
    const std::vector<std::vector<double> >& fxMap,
    const std::vector<std::vector<double> >& fyMap) const {
  constexpr double r2 = 1. / 2.;
  constexpr double r6 = 1. / 6.;
  constexpr double r8 = 1. / 8.;
  if (h == 0) return true;
  const double h2 = r2 * h;
  const double h6 = r6 * h;
  const double hh2 = h * h2;
  const double hh6 = h * h6;
  const double hh8 = r8 * h * h;

  constexpr unsigned int n = 2;
  std::array<std::array<double, n>, 6> w;
  if (!GetForceRatio(wire, x, y, yp, w[0], xMap, yMap, fxMap, fyMap)) {
    return false;
  }
  for (unsigned int j = 0; j < n; ++j) {
    w[3][j] = y[j] + h2 * yp[j];
    w[4][j] = w[3][j] + hh8 * w[0][j];
    w[5][j] = yp[j] + h2 * w[0][j];
  }
  const double xh2 = x + h2;
  if (!GetForceRatio(wire, xh2, w[4], w[5], w[1], xMap, yMap, fxMap, fyMap)) {
    return false;
  }
  for (unsigned int j = 0; j < n; ++j) {
    w[5][j] = yp[j] + h2 * w[1][j];
    w[0][j] = w[0][j] + w[1][j];
    w[1][j] = w[0][j] + w[1][j];
  }
  if (!GetForceRatio(wire, xh2, w[4], w[5], w[2], xMap, yMap, fxMap, fyMap)) {
    return false;
  }
  for (unsigned int j = 0; j < n; ++j) {
    w[3][j] = w[3][j] + h2 * yp[j];
    w[4][j] = w[3][j] + hh2 * w[2][j];
    w[5][j] = yp[j] + h * w[2][j];
    w[0][j] = w[0][j] + w[2][j];
    w[1][j] = w[1][j] + 2 * w[2][j];
  }
  const double xh = x + h;
  if (!GetForceRatio(wire, xh, w[4], w[5], w[2], xMap, yMap, fxMap, fyMap)) {
    return false;
  }
  for (unsigned int j = 0; j < n; ++j) {
    y[j] = w[3][j] + hh6 * w[0][j];
    yp[j] += h6 * (w[1][j] + w[2][j]);
  }
  x = xh;
  return true;
}

bool ComponentAnalyticField::FindZeroes(
    const Wire& wire, const double h, std::vector<double>& x,
    const std::vector<double>& xMap, const std::vector<double>& yMap,
    const std::vector<std::vector<double> >& fxMap,
    const std::vector<std::vector<double> >& fyMap) const {
  //-----------------------------------------------------------------------
  //   OPTZRO - Tries to find zeroes of a set of functions F. Uses the
  //            Broyden rank-1 update variant of an n-dimensional Newton-
  //            Raphson zero search in most steps, except every 5th step
  //            and whenever the step length update becomes less than 0.5,
  //            when a new derivative is computed.
  //-----------------------------------------------------------------------

  if (x.empty()) {
    std::cerr << m_className << "::FindZeroes: Empty vector.\n";
    return false;
  }
  const unsigned int n = x.size();
  // Initial deviation.
  std::vector<double> fold(n, 0.);
  if (!Trace(wire, h, x, fold, xMap, yMap, fxMap, fyMap)) {
    std::cerr << m_className << "::FindZeroes: Zero search stopped.\n    "
              << "Initial position outside scanning area.\n";
    return false;
  }
  double fnorml =
      std::inner_product(fold.begin(), fold.end(), fold.begin(), 0.);
  // Debugging output for initial situation.
  constexpr unsigned int nbsmax = 10;
  constexpr unsigned int nitmax = 20;
  constexpr double eps = 1.e-4;
  constexpr double epsx = 1.e-4;
  constexpr double epsf = 1.e-4;
  if (m_debug) {
    std::cout << m_className << "::FindZeroes: Start of zero search.\n"
              << "    Number of parameters:      " << n << "\n"
              << "    Maximum bisections:        " << nbsmax << "\n"
              << "    Maximum iterations:        " << nitmax << "\n"
              << "    Epsilon differentiation:   " << eps << "\n"
              << "    Required location change:  " << epsx << "\n"
              << "    Required function norm:    " << epsf << "\n"
              << "    Initial function norm:     " << sqrt(fnorml) << "\n"
              << " Parameter        Value     Function\n";
    for (unsigned int i = 0; i < n; ++i) {
      std::printf(" %9d %12.5e %12.5e\n", i, x[i], fold[i]);
    }
  }
  // Derivative matrix.
  std::vector<std::vector<double> > b(n, std::vector<double>(n, 0.));
  std::vector<std::vector<double> > bb(n, std::vector<double>(n, 0.));
  // Flag whether the matrix needs to be recomputed.
  bool updateMatrix = true;
  // Count function calls.
  unsigned int nCalls = 0;
  bool converged = false;
  for (unsigned int iter = 0; iter < nitmax; ++iter) {
    // If needed, (re-)compute the derivative matrix.
    if (updateMatrix) {
      std::vector<double> f1(n, 0.);
      std::vector<double> f2(n, 0.);
      for (unsigned int i = 0; i < n; ++i) {
        const double epsdif = eps * (1. + std::abs(x[i]));
        x[i] += 0.5 * epsdif;
        if (!Trace(wire, h, x, f1, xMap, yMap, fxMap, fyMap)) {
          std::cerr << m_className << "::FindZeroes: Zero search stopped.\n    "
                    << "Differential matrix requires a point "
                    << "outside scanning area.\n";
          return false;
        }
        x[i] -= epsdif;
        if (!Trace(wire, h, x, f2, xMap, yMap, fxMap, fyMap)) {
          std::cerr << m_className << "::FindZeroes: Zero search stopped.\n    "
                    << "Differential matrix requires a point "
                    << "outside scanning area.\n";
          return false;
        }
        x[i] += 0.5 * epsdif;
        for (unsigned int j = 0; j < n; ++j) {
          b[j][i] = (f1[j] - f2[j]) / epsdif;
        }
      }
      nCalls += 2 * n;
      updateMatrix = false;
    }
    if (m_debug) {
      std::cout << "    Start of iteration " << iter << "\n";
      for (unsigned int i = 0; i < m_nShots; ++i) {
        const unsigned int k = 4 * i + 2;
        std::printf("     x = %12.5e,  y = %12.5e\n", x[k], x[k + 1]);
      }
    }
    // Find the correction vector to 0th order.
    std::vector<double> dx = fold;
    bb = b;
    if (Numerics::CERNLIB::deqn(n, bb, dx) != 0) {
      std::cerr << m_className << "::FindZeroes: Zero search stopped.\n"
                << "    Solving the update equation failed.\n";
      break;
    }
    if (m_debug) {
      for (unsigned int i = 0; i < m_nShots; ++i) {
        const unsigned int k = 4 * i + 2;
        std::printf("    dx = %12.5e, dy = %12.5e\n", dx[k], dx[k + 1]);
      }
    }
    // Scale the correction vector to improve FNORM, AUX3: f.
    double scale = 1.;
    double fnorm = 2 * fnorml;
    std::vector<double> xnew(n, 0.);
    std::vector<double> fnew(n, 0.);
    for (unsigned int kbs = 0; kbs < nbsmax; ++kbs) {
      for (unsigned int i = 0; i < n; ++i) {
        xnew[i] = x[i] - scale * dx[i];
      }
      if (!Trace(wire, h, xnew, fnew, xMap, yMap, fxMap, fyMap)) {
        std::cerr
            << m_className << "::FindZeroes: Zero search stopped.\n    "
            << "Step update leads to a point outside the scanning area.\n";
        return false;
      }
      ++nCalls;
      fnorm = std::inner_product(fnew.begin(), fnew.end(), fnew.begin(), 0.);
      if (fnorm <= fnorml) {
        if (m_debug) std::cout << "    Scaling factor: " << scale << "\n";
        break;
      }
      scale *= 0.5;
    }
    if (fnorm > fnorml) {
      std::cerr << m_className << "::FindZeroes: Zero search stopped.\n    "
                << "Bisection search for scaling factor did not converge.\n";
      break;
    }
    // Update the estimate.
    std::vector<double> df(n, 0.);
    for (unsigned int i = 0; i < n; ++i) {
      dx[i] = xnew[i] - x[i];
      x[i] = xnew[i];
      df[i] = fnew[i] - fold[i];
      fold[i] = fnew[i];
    }
    double xnorm = std::inner_product(x.begin(), x.end(), x.begin(), 0.);
    double dxnorm = std::inner_product(dx.begin(), dx.end(), dx.begin(), 0.);
    double dfnorm = std::inner_product(df.begin(), df.end(), df.begin(), 0.);
    // Debugging output to show current status.
    if (m_debug) {
      std::cout << "    After this iteration...\n";
      std::printf("      Norm and change of position: %12.5e %12.5e\n",
                  sqrt(xnorm), sqrt(dxnorm));
      std::printf("      Norm and change of function: %12.5e %12.5e\n",
                  sqrt(fnorm), sqrt(dfnorm));
    }
    // See whether convergence has been achieved.
    if (sqrt(dxnorm) < epsx * sqrt(xnorm)) {
      if (m_debug) {
        std::cout << "    Positional convergence criterion is satisfied.\n";
      }
      converged = true;
      break;
    } else if (sqrt(fnorm) < epsf) {
      if (m_debug) {
        std::cout << "    Function value convergence criterion is satisfied.\n";
      }
      converged = true;
      break;
    }
    // Update the difference.
    fnorml = fnorm;
    if (scale > 0.4 && iter != 5 * (iter / 5)) {
      // If the scaling factor is small, then update (rank-1 Broyden).
      if (m_debug) std::cout << "    Performing a Broyden rank-1 update.\n";
      // Compute the "df - B dx" term.
      std::vector<double> corr(n, 0.);
      for (unsigned int i = 0; i < n; ++i) {
        corr[i] = df[i];
        for (unsigned int j = 0; j < n; ++j) {
          corr[i] -= b[i][j] * dx[j];
        }
      }
      // Update the matrix.
      for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
          b[i][j] += corr[i] * dx[j] / dxnorm;
        }
      }
    } else {
      // Otherwise, recompute the differential.
      if (m_debug) std::cout << "    Recomputing the covariance matrix.\n";
      updateMatrix = true;
    }
  }
  if (!converged) {
    std::cerr << m_className << "::FindZeroes: Search did not converge.\n";
  }
  if (m_debug) {
    std::vector<double> f(n, 0.);
    Trace(wire, h, x, f, xMap, yMap, fxMap, fyMap);
    ++nCalls;
    std::cout << "    Final values:\n"
              << " Parameter        Value     Function\n";
    for (unsigned int i = 0; i < n; ++i) {
      std::printf(" %9d %12.5e %12.5e\n", i, x[i], f[i]);
    }
    std::cout << "    Total number of function calls: " << nCalls << "\n";
  }
  return converged;
}

bool ComponentAnalyticField::Trace(
    const Wire& wire, const double h, const std::vector<double>& xx,
    std::vector<double>& delta, const std::vector<double>& xMap,
    const std::vector<double>& yMap,
    const std::vector<std::vector<double> >& fxMap,
    const std::vector<std::vector<double> >& fyMap) const {
  //-----------------------------------------------------------------------
  //   OPTSHT - Auxiliary routine for the wire sag routines which computes
  //            for a given set of positions and derivatives the next set
  //            which is used by OPTZRO to match the sections. Uses a
  //            2nd order Runge-Kutta-Nystrom integration routine (D203).
  //-----------------------------------------------------------------------

  delta.assign(xx.size(), 0.);
  // For the starting set in XX, compute the next round.
  double z = -0.5 * wire.u;
  std::array<double, 2> xst = {0., 0.};
  std::array<double, 2> dxst = {xx[0], xx[1]};
  for (unsigned int i = 0; i <= m_nShots; ++i) {
    const unsigned int k = 4 * i;
    // Set the starting value and starting derivative.
    if (i > 0) {
      xst = {xx[k - 2], xx[k - 1]};
      dxst = {xx[k], xx[k + 1]};
    }
    // Compute the end value and end derivative.
    for (unsigned int j = 0; j < m_nSteps; ++j) {
      if (!StepRKN(wire, h, z, xst, dxst, xMap, yMap, fxMap, fyMap)) {
        return false;
      }
    }
    // Store the differences as function value.
    if (i < m_nShots) {
      delta[k] = xst[0] - xx[k + 2];
      delta[k + 1] = xst[1] - xx[k + 3];
      delta[k + 2] = dxst[0] - xx[k + 4];
      delta[k + 3] = dxst[1] - xx[k + 5];
    } else {
      delta[k] = xst[0];
      delta[k + 1] = xst[1];
    }
  }
  return true;
}

bool ComponentAnalyticField::SetupDipoleTerms() {

  //-----------------------------------------------------------------------
  //   SETDIP - Subroutine computing coefficients for the dipole terms in
  //            cells without periodicities.
  //-----------------------------------------------------------------------

  // Parameters.
  constexpr double epsp = 1.e-3;
  constexpr double epsa = 1.e-3;

  const unsigned int nWires = m_w.size();
  // Initial dipole moments.
  std::vector<double> phi2(nWires, 0.);
  m_cosph2.assign(nWires, 1.);
  m_sinph2.assign(nWires, 0.);
  m_amp2.assign(nWires, 0.);

  // Iterate until the dipole terms have converged.
  std::vector<double> phit2(nWires, 0.);
  std::vector<double> ampt2(nWires, 0.);
  constexpr unsigned int nMaxIter = 10;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    if (m_debug) {
      std::cout << "Iteration " << iter << "/" << nMaxIter << "\n"
                << "  Wire  correction angle [deg]  amplitude\n";
    }
    // Loop over the wires.
    for (unsigned int iw = 0; iw < nWires; ++iw) {
      const double xw = m_w[iw].x;
      const double yw = m_w[iw].y;
      const double rw = m_w[iw].r;
      // Set the radius of the wire to 0.
      m_w[iw].r = 0.;
      // Loop around the wire.
      constexpr unsigned int nAngles = 20;
      std::vector<double> angle(nAngles, 0.);
      std::vector<double> volt(nAngles, 0.);
      constexpr double rmult = 1.;
      const double r = rw * rmult;
      int status = 0;
      for (unsigned int i = 0; i < nAngles; ++i) {
        angle[i] = TwoPi * (i + 1.) / nAngles;
        const double x = xw + r * cos(angle[i]);
        const double y = yw + r * sin(angle[i]);
        double ex = 0., ey = 0., ez = 0.;
        status = Field(x, y, 0., ex, ey, ez, volt[i], true);
        if (status != 0) {
          std::cerr << "Unexpected status code; computation stopped.\n";
          break;
        }
        volt[i] -= m_w[iw].v;
      }
      // Restore wire radius.
      m_w[iw].r = rw;
      if (status != 0) continue;
      // Determine the dipole term.
      double ampdip = 0., phidip = 0.;
      FitDipoleMoment(angle, volt, ampdip, phidip, false);
      // Store the parameters, removing the radial dependence.
      phit2[iw] = phidip;
      ampt2[iw] = ampdip * r;
      if (m_debug) {
        std::printf(" %3d  %10.3f  %12.5e\n", 
                    iw, RadToDegree * phit2[iw], ampt2[iw]);
      }
    }
    // Transfer to the arrays where the dipole moments have impact
    bool reiter = false;
    if (m_debug) std::cout << "  Wire  new angle [deg]  amplitude\n";
    for (unsigned int iw = 0; iw < nWires; ++iw) {
      // See whether we need further refinements.
      bool converged = true;
      if (std::abs(phi2[iw]) > epsp * (1. + std::abs(phi2[iw])) || 
          std::abs(m_amp2[iw]) > epsa * (1. + std::abs(m_amp2[iw]))) {
        reiter = true;
        converged = false;
      }
      // Add the new term to the existing one.
      const double s0 = m_sinph2[iw] * m_amp2[iw] + sin(phit2[iw]) * ampt2[iw];
      const double c0 = m_cosph2[iw] * m_amp2[iw] + cos(phit2[iw]) * ampt2[iw];
      phi2[iw] = atan2(s0, c0);
      m_cosph2[iw] = cos(phi2[iw]);
      m_sinph2[iw] = sin(phi2[iw]);
      const double s1 = m_sinph2[iw] * m_amp2[iw] + sin(phit2[iw]) * ampt2[iw];
      const double c1 = m_cosph2[iw] * m_amp2[iw] + cos(phit2[iw]) * ampt2[iw];
      m_amp2[iw] = sqrt(s1 * s1 + c1 * c1);
      if (m_debug) {
        std::printf(" %3d  %10.3f  %12.5e %s\n",
                    iw, RadToDegree * phi2[iw], m_amp2[iw], 
                    converged ? "CONVERGED" : "");
      }
    }
    // Next iteration?
    if (!reiter) return true;
  }
  // Maximum number of iterations exceeded.
  std::cerr << m_className << "::SetupDipoleTerms:\n"
            << "    Maximum number of dipole iterations exceeded "
            << "without convergence; abandoned.\n";
  return false;
}

void ComponentAnalyticField::DipoleFieldA00(
    const double xpos, const double ypos, 
    double& ex, double& ey, double& volt, const bool opt) const {

  //-----------------------------------------------------------------------
  //   EMCA00 - Subroutine computing dipole terms in a field generated by
  //            wires without periodicities.
  //-----------------------------------------------------------------------

  // Initialise the potential and the electric field.
  ex = 0.;
  ey = 0.;
  volt = 0.;
  // Loop over all wires.
  double v = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    const double dx = xpos - wire.x;
    const double dy = ypos - wire.y;
    const double dxm = xpos + wire.x - 2. * m_coplax;
    const double dym = ypos + wire.y - 2. * m_coplay;
    // Calculate the field in case there are no planes.
    const double a = dx * dx - dy * dy;
    const double b = 2 * dx * dy;
    const double d2 = dx * dx + dy * dy;
    const double d4 = d2 * d2;
    double fx = (a * m_cosph2[i] + b * m_sinph2[i]) / d4;
    double fy = (b * m_cosph2[i] - a * m_sinph2[i]) / d4;
    if (opt) v = (dx * m_cosph2[i] + dy * m_sinph2[i]) / d2;
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double am = dxm * dxm - dy * dy;
      const double bm = 2 * dxm * dy;
      const double d2m = dxm * dxm + dy * dy;
      const double d4m = d2m * d2m;
      fx -= (am * m_cosph2[i] + bm * m_sinph2[i]) / d4m;
      fy -= (bm * m_cosph2[i] - am * m_sinph2[i]) / d4m;
      if (opt) v -= (dxm * m_cosph2[i] + dy * m_sinph2[i]) / d2m;
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double am = dx * dx - dym * dym;
      const double bm = 2 * dx * dym;
      const double d2m = dx * dx + dym * dym;
      const double d4m = d2m * d2m;
      fx -= (am * m_cosph2[i] + bm * m_sinph2[i]) / d4m;
      fy -= (bm * m_cosph2[i] - am * m_sinph2[i]) / d4m;
      if (opt) v -= (dx * m_cosph2[i] + dym * m_sinph2[i]) / d2m;
    }
    // Take care of pairs of planes.
    if (m_ynplax && m_ynplay) {
      const double am = dxm * dxm - dym * dym;
      const double bm = 2 * dxm * dym;
      const double d2m = dxm * dxm + dym * dym;
      const double d4m = d2m * d2m;
      fx += (am * m_cosph2[i] + bm * m_sinph2[i]) / d4m;
      fy += (bm * m_cosph2[i] - am * m_sinph2[i]) / d4m;
      if (opt) v += (dxm * m_cosph2[i] + dym * m_sinph2[i]) / d2m;
    }
    // Normalise.
    volt -= m_amp2[i] * v; 
    ex -= m_amp2[i] * fx;
    ey -= m_amp2[i] * fy;
  }
}

void ComponentAnalyticField::DipoleFieldB1X(
    const double xpos, const double ypos, 
    double& ex, double& ey, double& volt, const bool opt) const {

  //-----------------------------------------------------------------------
  //   EMCB1X - Subroutine computing dipole terms in a field generated by
  //            a row of wires along the x-axis.
  //-----------------------------------------------------------------------

  // Initialise the potential and the electric field.
  ex = 0.;
  ey = 0.;
  volt = 0.;
  // Shorthand.
  const double tx = Pi / m_sx;
  const double tx2 = tx * tx;
  // Loop over all wires.
  double v = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    // Calculate the field in case there are no planes.
    const double dx = tx * (xpos - wire.x);
    const double dy = tx * (ypos - wire.y);
    const double a = 1 - cos(2 * dx) * cosh(2 * dy);
    const double b = sin(2 * dx) * sinh(2 * dy);
    const double sx = sin(dx);
    const double shy = sinh(dy);
    const double d2 = sx * sx + shy * shy;
    const double d4 = d2 * d2;
    double fx = ( m_cosph2[i] * a + m_sinph2[i] * b) / d4;
    double fy = (-m_sinph2[i] * a + m_cosph2[i] * b) / d4;
    if (opt) {
      v = (m_cosph2[i] * sin(2 * dx) + m_sinph2[i] * sinh(2 * dy)) / d2;
    }
    // Take care of a plane at constant y.
    if (m_ynplay) {
      const double dym = tx * (ypos + wire.y - 2. * m_coplay);
      const double am = 1 - cos(2 * dx) * cosh(2 * dym);
      const double bm = sin(2 * dx) * sinh(2 * dym);
      const double shym = sinh(dym);
      const double d2m = sx * sx + shym * shym;
      const double d4m = d2m * d2m;
      fx -= (m_cosph2[i] * am - m_sinph2[i] * bm) / d4m;
      fy -= (m_sinph2[i] * am + m_cosph2[i] * bm) / d4m;
      if (opt) {
        v -= (m_cosph2[i] * sin(2 * dx) - m_sinph2[i] * sinh(2 * dym)) / d2m;
      }
    }
    // Calculate the electric field and the potential.
    ex -= m_amp2[i] * 0.5 * tx2 * fx;
    ey -= m_amp2[i] * 0.5 * tx2 * fy;
    if (opt) volt -= 0.5 * tx * m_amp2[i] * v;
  } 
}

void ComponentAnalyticField::DipoleFieldB1Y(
    const double xpos, const double ypos, 
    double& ex, double& ey, double& volt, const bool opt) const {

  //-----------------------------------------------------------------------
  //   EMCB1Y - Subroutine computing dipole terms in a field generated by
  //            a row of wires along the y-axis.
  //-----------------------------------------------------------------------

  // Initialise the potential and the electric field.
  ex = 0.;
  ey = 0.;
  volt = 0.;
  // Shorthand.
  const double ty = Pi / m_sy;
  const double ty2 = ty * ty;
  // Loop over all wires.
  double v = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    // Calculate the field in case there are no planes.
    const double dx = ty * (xpos - wire.x);
    const double dy = ty * (ypos - wire.y);
    const double a = 1 - cosh(2 * dx) * cos(2 * dy);
    const double b = sinh(2 * dx) * sin(2 * dy);
    const double shx = sinh(dx);
    const double sy = sin(dy);
    const double d2 = shx * shx + sy * sy;
    const double d4 = d2 * d2;
    double fx = (-m_cosph2[i] * a + m_sinph2[i] * b) / d4;
    double fy = ( m_sinph2[i] * a + m_cosph2[i] * b) / d4;
    if (opt) {
      v = (m_cosph2[i] * sinh(2 * dx) + m_sinph2[i] * sin(2 * dy)) / d2;
    }
    // Take care of a plane at constant x.
    if (m_ynplax) {
      const double dxm = ty * (xpos + wire.x - 2. * m_coplax);
      const double am = 1 - cosh(2 * dxm) * cos(2 * dy);
      const double bm = sinh(2 * dxm) * sin(2 * dy);
      const double shxm = sinh(dxm);
      const double d2m = shxm * shxm + sy * sy;
      const double d4m = d2m * d2m;
      fx -= (m_cosph2[i] * am + m_sinph2[i] * bm) / d4m;
      fy -= (m_sinph2[i] * am - m_cosph2[i] * bm) / d4m;
      if (opt) {
        v -= (-m_cosph2[i] * sinh(2 * dxm) + m_sinph2[i] * sin(2 * dy)) / d2m;
      }
    }
    // Calculate the electric field and the potential.
    ex -= m_amp2[i] * 0.5 * ty2 * fx;
    ey -= m_amp2[i] * 0.5 * ty2 * fy;
    if (opt) volt -= 0.5 * ty * m_amp2[i] * v;
  }
}

void ComponentAnalyticField::DipoleFieldB2X(
    const double xpos, const double ypos, 
    double& ex, double& ey, double& volt, const bool opt) const {

  //-----------------------------------------------------------------------
  //   EMCB2X - Routine calculating the dipole terms for a charge between
  //            parallel conducting planes.
  //-----------------------------------------------------------------------

  // Initialise the potential and the electric field.
  ex = 0.;
  ey = 0.;
  volt = 0.;
  // Shorthand.
  const double tx = HalfPi / m_sx;
  const double tx2 = tx * tx;
  // Loop over all wires.
  double v = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    const double dx = tx * (xpos - wire.x);
    const double dy = tx * (ypos - wire.y);
    // Calculate the field in case there are no equipotential planes.
    const double a = 1 - cos(2 * dx) * cosh(2 * dy);
    const double b = sin(2 * dx) * sinh(2 * dy);
    const double sx = sin(dx);
    const double shy = sinh(dy);
    const double d2 = sx * sx + shy * shy;
    const double d4 = d2 * d2;
    double fx = ( m_cosph2[i] * a + m_sinph2[i] * b) / d4;
    double fy = (-m_sinph2[i] * a + m_cosph2[i] * b) / d4;
    const double dxn = tx * (xpos + wire.x - 2. * m_coplax);
    const double an = 1 - cos(2 * dxn) * cosh(2 * dy);
    const double bn = sin(2 * dxn) * sinh(2 * dy);
    const double sxn = sin(dxn);
    const double d2n = sxn * sxn + shy * shy; 
    const double d4n = d2n * d2n;
    fx += (m_cosph2[i] * an - m_sinph2[i] * bn) / d4n;
    fy += (m_sinph2[i] * an + m_cosph2[i] * bn) / d4n;
    if (opt) {
      v = -sin(dx + dxn) * (-2 * m_cosph2[i] * (
          (1 + shy * shy) * sx * sxn + cos(dx) * cos(dxn) * shy * shy) + 
          m_sinph2[i] * m_b2sin[i] * sinh(2 * dy)) / (d2 * d2n);
    }
    // Take care of planes at constant y.
    if (m_ynplay) {
      const double dym = tx * (ypos + wire.y - 2. * m_coplay); 
      const double am = 1 - cos(2 * dx) * cosh(2 * dym);
      const double bm = sin(2 * dx) * sinh(2 * dym);
      const double shym = sinh(dym);
      const double d2m = sx * sx + shym * shym;
      const double d4m = d2m * d2m;
      fx -= (m_cosph2[i] * am - m_sinph2[i] * bm) / d4m;
      fy -= (m_sinph2[i] * am + m_cosph2[i] * bm) / d4m;
      const double amn = 1 - cos(2 * dxn) * cosh(2 * dym);
      const double bmn = sin(2 * dxn) * sinh(2 * dym);
      const double d2mn = sxn * sxn + shym * shym;
      const double d4mn = d2mn * d2mn;
      fx -= ( m_cosph2[i] * amn + m_sinph2[i] * bmn) / d4mn;
      fy -= (-m_sinph2[i] * amn + m_cosph2[i] * bmn) / d4mn;
      if (opt) {
        v += sin(dx + dxn) * (-2 * m_cosph2[i] * (
            (1 + shym * shym) * sx * sxn + cos(dx) * cos(dxn) * shym * shym) -
            m_sinph2[i]  * m_b2sin[i] * sinh(2 * dym)) / (d2m * d2mn);
      }
    }
    // Calculate the electric field and the potential.
    ex -= m_amp2[i] * 0.5 * tx2 * fx;
    ey -= m_amp2[i] * 0.5 * tx2 * fy;
    if (opt) volt -= 0.5 * tx * m_amp2[i] * v;
  }
}

void ComponentAnalyticField::DipoleFieldB2Y(
    const double xpos, const double ypos, 
    double& ex, double& ey, double& volt, const bool opt) const {

  //-----------------------------------------------------------------------
  //   EMCB2Y - Routine calculating the dipole terms for a charge between
  //            parallel conducting planes.
  //-----------------------------------------------------------------------

  // Initialise the potential and the electric field.
  ex = 0.;
  ey = 0.;
  volt = 0.;
  // Shorthand.
  const double ty = HalfPi / m_sy;
  const double ty2 = ty * ty;
  // Loop over all wires.
  double v = 0.;
  for (unsigned int i = 0; i < m_nWires; ++i) {
    const auto& wire = m_w[i];
    const double dx = ty * (xpos - wire.x);
    const double dy = ty * (ypos - wire.y);
    // Calculate the field in case there are no equipotential planes.
    const double a = 1 - cosh(2 * dx) * cos(2 * dy);
    const double b = sinh(2 * dx) * sin(2 * dy);
    const double shx = sinh(dx);
    const double sy = sin(dy);
    const double d2 = shx * shx + sy * sy;
    const double d4 = d2 * d2;
    double fx = (-m_cosph2[i] * a + m_sinph2[i] * b) / d4; 
    double fy = ( m_sinph2[i] * a + m_cosph2[i] * b) / d4; 
    const double dyn = ty * (ypos + wire.y - 2. * m_coplay);
    const double an = 1 - cosh(2 * dx) * sin(2 * dyn);
    const double bn = sinh(2 * dx) * sin(2 * dyn);
    const double syn = sin(dyn);
    const double d2n = shx * shx + syn * syn;
    const double d4n = d2n * d2n; 
    fx += (m_cosph2[i] * an + m_sinph2[i] * bn) / d4n;
    fy += (m_sinph2[i] * an - m_cosph2[i] * bn) / d4n;
    if (opt) {
      // TODO: check!
      v = sin(dy + dyn) * (
          m_sinph2[i] * (-cos(dy + dyn) + cos(dy - dyn) * cosh(2 * dx)) -
          m_cosph2[i] * sinh(2 * dx) * (cos(dyn) * sy + cos(dy) * syn)) /
          (d2 * d2n);
    }
    // Take care of planes at constant x.
    if (m_ynplax) {
      const double dxm = ty * (xpos + wire.x - 2. * m_coplax);
      const double am = 1 - cosh(2 * dxm) * cos(2 * dy);
      const double bm = sinh(2 * dxm) * sin(2 * dy);
      const double shxm = sinh(dxm);
      const double d2m = shxm * shxm + sy * sy;
      const double d4m = d2m * d2m;
      fx -= (m_cosph2[i] * am + m_sinph2[i] * bm) / d4m;
      fy -= (m_sinph2[i] * am - m_cosph2[i] * bm) / d4m;
      const double amn = 1 - cosh(2 * dxm) * cos(2 * dyn);
      const double bmn = sinh(2 * dxm) * sin(2 * dyn);
      const double d2mn = shxm * shxm + syn * syn;
      const double d4mn = d2mn * d2mn;
      fx -= (-m_cosph2[i] * amn + m_sinph2[i] * bmn) / d4mn; 
      fy -= ( m_sinph2[i] * amn + m_cosph2[i] * bmn) / d4mn;
      if (opt) {
        // TODO: check!
        v -= sin(dy + dyn) * (
            m_sinph2[i] * (-cos(dy + dyn) + cos(dy - dyn) * cosh(2 * dxm)) +  
            m_cosph2[i] * sinh(2 * dxm) * (cos(dyn) * sy - cos(dy) * syn)) /
            (d2m * d2mn);

      }
    }
    // Calculate the electric field and the potential.
    ex -= m_amp2[i] * 0.5 * ty2 * fx;
    ey -= m_amp2[i] * 0.5 * ty2 * fy;
    if (opt) volt -= 0.5 * ty * m_amp2[i] * v;
  }
}

bool ComponentAnalyticField::MultipoleMoments(const unsigned int iw, 
    const unsigned int nPoles, const bool print, const bool plot, 
    const double rmult, const double eps,
    const unsigned int nMaxIter) {

  //-----------------------------------------------------------------------
  //   EFMWIR - Computes the dipole moment of a given wire.
  //-----------------------------------------------------------------------
  if (!m_cellset && !Prepare()) return false;
  // Check input parameters.
  if (iw >= m_nWires) {
    std::cerr << m_className << "::MultipoleMoments:\n"
              << "    Wire index out of range.\n";
    return false;
  }
  if (eps <= 0.) {
    std::cerr << m_className << "::MultipoleMoments:\n"
              << "    Epsilon must be positive.\n";
    return false;
  }
  if (nPoles == 0) {
    std::cerr << m_className << "::MultipoleMoments:\n"
              << "    Multipole order out of range.\n";
    return false;
  }
  if (rmult <= 0.) {
    std::cerr << m_className << "::MultipoleMoments:\n"
              << "    Radius multiplication factor out of range.\n";
    return false;
  }

  const double xw = m_w[iw].x;
  const double yw = m_w[iw].y;
  // Set the radius of the wire to 0.
  const double rw = m_w[iw].r;
  m_w[iw].r = 0.;

  // Loop around the wire.
  constexpr unsigned int nPoints = 20000;
  std::vector<double> angle(nPoints, 0.);
  std::vector<double> volt(nPoints, 0.);
  std::vector<double> weight(nPoints, 1.);
  for (unsigned int i = 0; i < nPoints; ++i) {
    // Set angle around wire.
    angle[i] = TwoPi * (i + 1.) / nPoints;
    // Compute E field, make sure the point is in a free region.
    const double x = xw + rmult * rw * cos(angle[i]); 
    const double y = yw + rmult * rw * sin(angle[i]);
    double ex = 0., ey = 0., ez = 0., v = 0.;
    if (Field(x, y, 0., ex, ey, ez, v, true) != 0) { 
      std::cerr << m_className << "::MultipoleMoments:\n"
                << "    Unexpected location code. Computation stopped.\n";
      m_w[iw].r = rw;
      return false;
    }
    // Assign the result to the fitting array.
    volt[i] = v;
  } 
  // Restore the wire diameter.
  m_w[iw].r = rw;

  // Determine the maximum, minimum and average.
  double vmin = *std::min_element(volt.cbegin(), volt.cend());
  double vmax = *std::max_element(volt.cbegin(), volt.cend());
  double vave = std::accumulate(volt.cbegin(), volt.cend(), 0.) / nPoints;
  // Subtract the wire potential to centre the data more or less.
  for (unsigned int i = 0; i < nPoints; ++i) volt[i] -= vave;
  vmax -= vave;
  vmin -= vave;

  // Perform the fit.
  const double vm = 0.5 * fabs(vmin) + fabs(vmax);
  double chi2 = 1.e-6 * nPoints * vm * vm;
  const double dist = 1.e-3 * (1. + vm);
  const unsigned int nPar = 2 * nPoles + 1;
  std::vector<double> pars(nPar, 0.);
  std::vector<double> epar(nPar, 0.);
  pars[0] = 0.5 * (vmax + vmin);
  for (unsigned int i = 1; i <= nPoles; ++i) {
    pars[2 * i - 1] = 0.5 * (vmax - vmin);
    pars[2 * i] = 0.;
  } 

  auto f = [nPoles](const double x, const std::vector<double>& par) {
    // EFMFUN
    // Sum the series, initial value is the monopole term.
    double sum = par[0];
    for (unsigned int k = 1; k <= nPoles; ++k) {
      // Obtain the Legendre polynomial of this order and add to the series.
      const float cphi = cos(x - par[2 * k]);
      sum += par[2 * k - 1] * sqrt(k + 0.5) * Numerics::Legendre(k, cphi);
    }
    return sum;
  };

  if (!Numerics::LeastSquaresFit(f, pars, epar, angle, volt, weight, 
                                 nMaxIter, dist, chi2, eps, m_debug, print)) {
    std::cerr << m_className << "::MultipoleMoments:\n"
              << "    Fitting the multipoles failed; computation stopped.\n";
  }
  // Plot the result of the fit.
  if (plot) {
    const std::string name = ViewBase::FindUnusedCanvasName("cMultipole");
    TCanvas* cfit = new TCanvas(name.c_str(), "Multipoles");
    cfit->SetGridx();
    cfit->SetGridy();
    cfit->DrawFrame(0., vmin, TwoPi, vmax,
                    ";Angle around the wire [rad]; Potential - average [V]");
    TGraph graph;
    graph.SetLineWidth(2);
    graph.SetLineColor(kBlack);
    graph.DrawGraph(angle.size(), angle.data(), volt.data(), "lsame");
    // Sum of contributions.
    constexpr unsigned int nP = 1000;
    std::array<double, nP> xp;
    std::array<double, nP> yp;
    for (unsigned int i = 0; i < nP; ++i) {
      xp[i] = TwoPi * (i + 1.) / nP;
      yp[i] = f(xp[i], pars);
    }
    graph.SetLineColor(kViolet + 3);
    graph.DrawGraph(nP, xp.data(), yp.data(), "lsame");
    // Individual contributions.
    std::vector<double> parres = pars;
    for (unsigned int i = 1; i <= nPoles; ++i) parres[2 * i - 1] = 0.;
    for (unsigned int j = 1; j <= nPoles; ++j) {
      parres[2 * j - 1] = pars[2 * j - 1];
      for (unsigned int i = 0; i < nP; ++i) {
        yp[i] = f(xp[i], parres);
      }
      parres[2 * j - 1] = 0.;
      graph.SetLineColor(kAzure + j);
      graph.DrawGraph(nP, xp.data(), yp.data(), "lsame"); 
    }
    gPad->Update();
  }

  // Print the results.
  std::cout << m_className << "::MultipoleMoments:\n"
            << "    Multipole moments for wire " << iw << ":\n"
            << "  Moment            Value       Angle [degree]\n";
  std::printf("  %6u  %15.8f        Arbitrary\n", 0, vave);
  for (unsigned int i = 1; i <= nPoles; ++i) {
    // Remove radial term from the multipole moment.
    const double val = pow(rmult * rw, i) * pars[2 * i - 1];
    const double phi = RadToDegree * fmod(pars[2 * i], Pi);
    std::printf("  %6u  %15.8f  %15.8f\n", i, val, phi);
  }
  return true;
}

}  // namespace Garfield
