#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>

#include "Garfield/ComponentGrid.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Utilities.hh"

namespace {

void PrintError(const std::string& fcn, const unsigned int line,
                const std::string& par) {
  std::cerr << fcn << ": Error reading line " << line << ".\n"
            << "    Could not read " << par << ".\n";
}

void PrintNotReady(const std::string& fcn) {
  std::cerr << fcn << ": Map not available.\n";
}

void PrintProgress(const double f) {
  if (f < 0.) return;
  constexpr unsigned int width = 70;
  const unsigned int n = static_cast<unsigned int>(std::floor(width * f));
  std::string bar = "[";
  if (n < 1) {
    bar += std::string(width, ' ');
  } else if (n >= width) {
    bar += std::string(width, '=');
  } else {
    bar += std::string(n, '=') + ">" + std::string(width - n - 1, ' ');
  }
  bar += "]";
  std::cout << bar << "\r" << std::flush;
}

bool IsComment(const std::string& line) {
  if (line.empty()) return false;
  if (line[0] == '#') return true;
  if (line.size() > 1 && (line[0] == '/' && line[1] == '/')) return true;
  return false;
}

}  // namespace

namespace Garfield {

ComponentGrid::ComponentGrid() : Component("Grid") {}

void ComponentGrid::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, double& p, Medium*& m,
                                  int& status) {
  m = nullptr;
  status = 0;

  // Make sure the field map has been loaded.
  if (m_efields.empty()) {
    PrintNotReady(m_className + "::ElectricField");
    status = -10;
    return;
  }

  status = 0;
  bool active = true;
  if (!GetField(x, y, z, m_efields, ex, ey, ez, p, active)) {
    status = -11;
    return;
  }
  if (!active) {
    status = -5;
    return;
  }
  m = m_medium;
  if (!m) status = -5;
}

void ComponentGrid::ElectricField(const double x, const double y,
                                  const double z, double& ex, double& ey,
                                  double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentGrid::WeightingField(const double x, const double y,
                                   const double z, double& wx, double& wy,
                                   double& wz, const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (m_wfields.empty()) return;
  const double xx = x - m_wFieldOffset[0];
  const double yy = y - m_wFieldOffset[1];
  const double zz = z - m_wFieldOffset[2];
  double wp = 0.;
  bool active = true;
  GetField(xx, yy, zz, m_wfields, wx, wy, wz, wp, active);
}

double ComponentGrid::WeightingPotential(const double x, const double y,
                                         const double z,
                                         const std::string& /*label*/) {
  if (m_wfields.empty()) return 0.;
  const double xx = x - m_wFieldOffset[0];
  const double yy = y - m_wFieldOffset[1];
  const double zz = z - m_wFieldOffset[2];
  double wx = 0., wy = 0., wz = 0.;
  double wp = 0.;
  bool active = true;
  if (!GetField(xx, yy, zz, m_wfields, wx, wy, wz, wp, active)) return 0.;
  return wp;
}

void ComponentGrid::DelayedWeightingField(const double x, const double y,
                                          const double z, const double t,
                                          double& wx, double& wy, double& wz,
                                          const std::string& /*label*/) {
  wx = wy = wz = 0.;
  if (m_wdtimes.empty()) return;
  // Assume no weighting field for times outside the range of available maps.
  if (t < m_wdtimes.front() || t > m_wdtimes.back()) return;

  const double xx = x - m_wFieldOffset[0];
  const double yy = y - m_wFieldOffset[1];
  const double zz = z - m_wFieldOffset[2];

  const auto it1 = std::upper_bound(m_wdtimes.cbegin(), m_wdtimes.cend(), t);
  const auto it0 = std::prev(it1);

  const double dt = t - *it0;
  double wp = 0.;
  const unsigned int i0 = it0 - m_wdtimes.cbegin();
  double wx0 = 0., wy0 = 0., wz0 = 0.;
  bool active = true;
  if (!GetField(xx, yy, zz, m_wdfields[i0], wx0, wy0, wz0, wp, active)) return;

  if (dt < Small || it1 == m_wdtimes.cend()) {
    wx = wx0;
    wy = wy0;
    wz = wz0;
    return;
  }
  const unsigned int i1 = it1 - m_wdtimes.cbegin();
  double wx1 = 0., wy1 = 0., wz1 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i1], wx1, wy1, wz1, wp, active)) return;

  const double f1 = dt / (*it1 - *it0);
  const double f0 = 1. - f1;
  wx = f0 * wx0 + f1 * wx1;
  wy = f0 * wy0 + f1 * wy1;
  wz = f0 * wz0 + f1 * wz1;
}

double ComponentGrid::DelayedWeightingPotential(
    const double x, const double y, const double z, const double t,
    const std::string& /*label*/) {

  if (m_wdtimes.empty()) return 0.;
  // Outside the range of the available maps?
  if (t < m_wdtimes.front() || t > m_wdtimes.back()) return 0.;

  const double xx = x - m_wFieldOffset[0];
  const double yy = y - m_wFieldOffset[1];
  const double zz = z - m_wFieldOffset[2];

  const auto it1 = std::upper_bound(m_wdtimes.cbegin(), m_wdtimes.cend(), t);
  const auto it0 = std::prev(it1);

  const double dt = t - *it0;
  const unsigned int i0 = it0 - m_wdtimes.cbegin();
  double wp0 = 0., wx0 = 0., wy0 = 0., wz0 = 0.;
  bool active = true;
  if (!GetField(xx, yy, zz, m_wdfields[i0], wx0, wy0, wz0, wp0, active)) return 0.;

  if (dt < Small || it1 == m_wdtimes.cend()) return wp0;

  const unsigned int i1 = it1 - m_wdtimes.cbegin();
  double wp1 = 0., wx1 = 0., wy1 = 0., wz1 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i1], wx1, wy1, wz1, wp1, active)) return 0.;

  const double f1 = dt / (*it1 - *it0);
  const double f0 = 1. - f1;
  return f0 * wp0 + f1 * wp1;
}

void ComponentGrid::SetWeightingFieldOffset(const double x, const double y,
                                            const double z) {
  m_wFieldOffset = {x, y, z};
}

void ComponentGrid::MagneticField(const double x, const double y,
                                  const double z, double& bx, double& by,
                                  double& bz, int& status) {
  status = 0;
  if (m_bfields.empty()) {
    return Component::MagneticField(x, y, z, bx, by, bz, status);
  }

  double p = 0.;
  bool active = true;
  if (!GetField(x, y, z, m_bfields, bx, by, bz, p, active)) {
    status = -11;
  }
}

bool ComponentGrid::HasMagneticField() const {
  return m_bfields.empty() ? Component::HasMagneticField() : true;
}

Medium* ComponentGrid::GetMedium(const double x, const double y,
                                 const double z) {
  // Make sure the field map has been loaded.
  if (m_efields.empty()) {
    PrintNotReady(m_className + "::GetMedium");
    return nullptr;
  }

  std::array<double, 3> xx = {x, y, z};
  if (m_coordinates == Coordinates::Cylindrical) {
    if (fabs(x) > Small || fabs(y) > Small) {
      xx[0] = sqrt(x * x + y * y);
      xx[1] = atan2(y, x);
    }
  } 
  for (size_t i = 0; i < 3; ++i) {
    if (m_periodic[i] || m_mirrorPeriodic[i]) continue;
    if (xx[i] < m_xMin[i] || xx[i] > m_xMax[i]) return nullptr;
  }
  if (m_active.empty()) return m_medium;
  for (size_t i = 0; i < 3; ++i) {
    bool mirrored = false;
    xx[i] = Reduce(xx[i], m_xMin[i], m_xMax[i], 
                   m_periodic[i], m_mirrorPeriodic[i], mirrored);
  }
  // Get the indices.
  const double sx = (xx[0] - m_xMin[0]) * m_sX[0];
  const double sy = (xx[1] - m_xMin[1]) * m_sX[1];
  const double sz = (xx[2] - m_xMin[2]) * m_sX[2];
  const unsigned int i0 = static_cast<unsigned int>(std::floor(sx));
  const unsigned int j0 = static_cast<unsigned int>(std::floor(sy));
  const unsigned int k0 = static_cast<unsigned int>(std::floor(sz));
  const unsigned int i1 = std::min(i0 + 1, m_nX[0] - 1);
  const unsigned int j1 = std::min(j0 + 1, m_nX[1] - 1);
  const unsigned int k1 = std::min(k0 + 1, m_nX[2] - 1);
  if (m_active[i0][j0][k0] && m_active[i0][j0][k1] && m_active[i0][j1][k0] &&
      m_active[i0][j1][k1] && m_active[i1][j0][k0] && m_active[i1][j0][k1] &&
      m_active[i1][j1][k0] && m_active[i1][j1][k1]) {
    return m_medium;
  }
  return nullptr;
}

bool ComponentGrid::SetMesh(const unsigned int nx, const unsigned int ny,
                            const unsigned int nz, const double xmin,
                            const double xmax, const double ymin,
                            const double ymax, const double zmin,
                            const double zmax) {
  Reset();
  if (nx == 0 || ny == 0 || nz == 0) {
    std::cerr << m_className << "::SetMesh:\n"
              << "    Number of mesh elements must be positive.\n";
    return false;
  }
  if (xmin >= xmax) {
    std::cerr << m_className << "::SetMesh: Invalid x range.\n";
    return false;
  } else if (ymin >= ymax) {
    std::cerr << m_className << "::SetMesh: Invalid y range.\n";
    return false;
  } else if (zmin >= zmax) {
    std::cerr << m_className << "::SetMesh: Invalid z range.\n";
    return false;
  }
  if (m_coordinates == Coordinates::Cylindrical) {
    if (xmin < 0.) {
      std::cerr << m_className << "::SetMesh: Invalid range.\n"
                << "    Radial coordinates must be >= 0.\n";
      return false; 
    }
  }
  m_nX[0] = nx;
  m_nX[1] = ny;
  m_nX[2] = nz;
  m_xMin[0] = xmin;
  m_xMin[1] = ymin;
  m_xMin[2] = zmin;
  m_xMax[0] = xmax;
  m_xMax[1] = ymax;
  m_xMax[2] = zmax;
  constexpr double tol = 1.e-10;
  for (size_t i = 0; i < 3; ++i) {
    if (m_xMax[i] - m_xMin[i] > tol) {
      m_sX[i] = std::max(m_nX[i] - 1., 1.) / (m_xMax[i] - m_xMin[i]);
    } else {
      m_sX[i] = 0.;
    }
  }
  if (m_coordinates == Coordinates::Cylindrical) {
    if (fabs(m_xMax[1] - m_xMin[1] - TwoPi) < tol) {
      if (!m_periodic[1]) {
        std::cerr << m_className << "::SetMesh: Enabling theta periodicity.\n";
      }
      m_periodic[1] = true;
      m_mirrorPeriodic[1] = false;
    }
  }
  m_hasMesh = true;
  return true;
}

bool ComponentGrid::GetMesh(unsigned int& nx, unsigned int& ny,
                            unsigned int& nz, double& xmin, double& xmax,
                            double& ymin, double& ymax, double& zmin,
                            double& zmax) const {
  if (!m_hasMesh) return false;
  nx = m_nX[0];
  ny = m_nX[1];
  nz = m_nX[2];
  xmin = m_xMin[0];
  ymin = m_xMin[1];
  zmin = m_xMin[2];
  xmax = m_xMax[0];
  ymax = m_xMax[1];
  zmax = m_xMax[2];
  return true;
}

void ComponentGrid::SetCylindricalCoordinates() {

  if (m_xMin[0] < 0. || m_xMax[0] < 0.) {
    std::cerr << m_className << "::SetCylindricalCoordinates:\n"
              << "    Invalid mesh limits. Radial coordinate must be >= 0.\n"; 
    return;
  }
  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    const double s = m_xMax[1] - m_xMin[1];
    if (std::abs(TwoPi - s * int(TwoPi / s)) > 1.e-4) {
      std::cerr << m_className << "::SetCylindricalCoordinates:\n"
                << "    Angular range does not divide 2 pi.\n"
                << "    Switching off periodicity.\n";
      m_periodic[1] = false;
      m_mirrorPeriodic[1] = false;
    }
  }
  m_coordinates = Coordinates::Cylindrical;
}

bool ComponentGrid::LoadElectricField(const std::string& fname,
                                      const std::string& fmt, const bool withP,
                                      const bool withFlag, const double scaleX,
                                      const double scaleE,
                                      const double scaleP) {
  m_efields.clear();
  m_hasPotential = false;
  m_active.assign(m_nX[0], std::vector<std::vector<bool> >(
                            m_nX[1], std::vector<bool>(m_nX[2], true)));
  // Read the file.
  m_pMin = withP ? +1. : 0.;
  m_pMax = withP ? -1. : 0.;
  if (!LoadData(fname, fmt, withP, withFlag, scaleX, scaleE, scaleP,
                 m_efields)) {
    m_efields.clear();
    return false;
  }
  if (withP) m_hasPotential = true;
  return true;
}

bool ComponentGrid::LoadWeightingField(const std::string& fname,
                                       const std::string& fmt, const bool withP,
                                       const double scaleX, const double scaleE,
                                       const double scaleP) {
  // Read the file.
  if (!LoadData(fname, fmt, withP, false, scaleX, scaleE, scaleP, m_wfields)) {
    m_wfields.clear();
    return false;
  }
  return true;
}

bool ComponentGrid::LoadWeightingField(const std::string& fname,
                                       const std::string& fmt, const double t,
                                       const bool withP, const double scaleX,
                                       const double scaleE,
                                       const double scaleP) {
  std::vector<std::vector<std::vector<Node> > > wfield;
  // Read the file.
  if (!LoadData(fname, fmt, withP, false, scaleX, scaleE, scaleP, wfield)) {
    return false;
  }
  if (m_wdtimes.empty() || t > m_wdtimes.back()) {
    m_wdtimes.push_back(t);
    m_wdfields.push_back(std::move(wfield));
  } else {
    const auto it = std::upper_bound(m_wdtimes.begin(), m_wdtimes.end(), t);
    const auto n = std::distance(m_wdtimes.begin(), it);
    m_wdtimes.insert(it, t);
    m_wdfields.insert( m_wdfields.begin() + n, std::move(wfield));
  }
  return true;
}

bool ComponentGrid::LoadMagneticField(const std::string& fname,
                                      const std::string& fmt,
                                      const double scaleX,
                                      const double scaleB) {
  // Read the file.
  if (!LoadData(fname, fmt, false, false, scaleX, scaleB, 1., m_bfields)) {
    m_bfields.clear();
    return false;
  }
  return true;
}

bool ComponentGrid::SaveElectricField(Component* cmp,
                                      const std::string& filename,
                                      const std::string& format) {
  if (!cmp) {
    std::cerr << m_className << "::SaveElectricField: Null pointer.\n";
    return false;
  }
  if (!m_hasMesh) {
    std::cerr << m_className << "::SaveElectricField: Mesh not set.\n";
    return false;
  }
  const auto fmt = GetFormat(format);
  if (fmt == Format::Unknown) {
    std::cerr << m_className << "::SaveElectricField:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile) {
    std::cerr << m_className << "::SaveElectricField:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }
  std::cout << m_className << "::SaveElectricField:\n"
            << "    Exporting field/potential to " << filename << ".\n"
            << "    Be patient...\n";
  PrintProgress(0.);
  outfile << "# XMIN = " << m_xMin[0] << ", XMAX = " << m_xMax[0] 
          << ", NX = " << m_nX[0] << "\n";
  outfile << "# YMIN = " << m_xMin[1] << ", YMAX = " << m_xMax[1] 
          << ", NY = " << m_nX[1] << "\n";
  outfile << "# ZMIN = " << m_xMin[2] << ", ZMAX = " << m_xMax[2] 
          << ", NZ = " << m_nX[2] << "\n";

  const unsigned int nValues = m_nX[0] * m_nX[1] * m_nX[2];
  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nValues)) - 1, 1.)));
  unsigned int nLines = 0;
  Medium* medium = nullptr;
  int status = 0;
  const double dx = (m_xMax[0] - m_xMin[0]) / std::max(m_nX[0] - 1., 1.);
  const double dy = (m_xMax[1] - m_xMin[1]) / std::max(m_nX[1] - 1., 1.);
  const double dz = (m_xMax[2] - m_xMin[2]) / std::max(m_nX[2] - 1., 1.);
  for (unsigned int i = 0; i < m_nX[0]; ++i) {
    const double x = m_xMin[0] + i * dx;
    for (unsigned int j = 0; j < m_nX[1]; ++j) {
      const double y = m_xMin[1] + j * dy;
      for (unsigned int k = 0; k < m_nX[2]; ++k) {
        const double z = m_xMin[2] + k * dz;
        if (fmt == Format::XY) {
          outfile << x << "  " << y << "  ";
        } else if (fmt == Format::XZ) {
          outfile << x << "  " << z << "  ";
        } else if (fmt == Format::XYZ) {
          outfile << x << "  " << y << "  " << z << "  ";
        } else if (fmt == Format::IJ) {
          outfile << i << "  " << j << "  ";
        } else if (fmt == Format::IK) {
          outfile << i << "  " << k << "  ";
        } else if (fmt == Format::IJK) {
          outfile << i << "  " << j << "  " << k << "  ";
        } else if (fmt == Format::YXZ) {
          outfile << y << "  " << x << "  " << z << "  ";
        }
        if (m_coordinates == Coordinates::Cylindrical) {
          const double ct = cos(y);
          const double st = sin(y);
          double ex = 0., ey = 0., ez = 0., v = 0.;
          cmp->ElectricField(x * ct, x * st, z, ex, ey, ez, v, medium, status);
          const double er = +ex * ct + ey * st;
          const double et = -ex * st + ey * ct;
          outfile << er << "  " << et << "  " << ez << "  " << v << "\n";
        } else { 
          double ex = 0., ey = 0., ez = 0., v = 0.;
          cmp->ElectricField(x, y, z, ex, ey, ez, v, medium, status);
          outfile << ex << "  " << ey << "  " << ez << "  " << v << "\n";
        }
        ++nLines;
        if (nLines % nPrint == 0) PrintProgress(double(nLines) / nValues);
      }
    }
  }
  outfile.close();
  std::cout << std::endl << m_className << "::SaveElectricField: Done.\n";
  return true;
}

bool ComponentGrid::SaveWeightingField(Component* cmp,
                                       const std::string& id,
                                       const std::string& filename,
                                       const std::string& format) {
  if (!cmp) {
    std::cerr << m_className << "::SaveWeightingField: Null pointer.\n";
    return false;
  }
  if (!m_hasMesh) {
    std::cerr << m_className << "::SaveWeightingField: Mesh not set.\n";
    return false;
  }
  const auto fmt = GetFormat(format);
  if (fmt == Format::Unknown) {
    std::cerr << m_className << "::SaveWeightingField:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile) {
    std::cerr << m_className << "::SaveWeightingField:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }
  std::cout << m_className << "::SaveWeightingField:\n"
            << "    Exporting field/potential to " << filename << ".\n"
            << "    Be patient...\n";
  PrintProgress(0.);
  outfile << "# XMIN = " << m_xMin[0] << ", XMAX = " << m_xMax[0] 
          << ", NX = " << m_nX[0] << "\n";
  outfile << "# YMIN = " << m_xMin[1] << ", YMAX = " << m_xMax[1] 
          << ", NY = " << m_nX[1] << "\n";
  outfile << "# ZMIN = " << m_xMin[2] << ", ZMAX = " << m_xMax[2] 
          << ", NZ = " << m_nX[2] << "\n";
  const unsigned int nValues = m_nX[0] * m_nX[1] * m_nX[2];
  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nValues)) - 1, 1.)));
  unsigned int nLines = 0;
  const double dx = (m_xMax[0] - m_xMin[0]) / std::max(m_nX[0] - 1., 1.);
  const double dy = (m_xMax[1] - m_xMin[1]) / std::max(m_nX[1] - 1., 1.);
  const double dz = (m_xMax[2] - m_xMin[2]) / std::max(m_nX[2] - 1., 1.);
  for (unsigned int i = 0; i < m_nX[0]; ++i) {
    const double x = m_xMin[0] + i * dx;
    for (unsigned int j = 0; j < m_nX[1]; ++j) {
      const double y = m_xMin[1] + j * dy;
      for (unsigned int k = 0; k < m_nX[2]; ++k) {
        const double z = m_xMin[2] + k * dz;
        if (fmt == Format::XY) {
          outfile << x << "  " << y << "  ";
        } else if (fmt == Format::XZ) {
          outfile << x << "  " << z << "  ";
        } else if (fmt == Format::XYZ) {
          outfile << x << "  " << y << "  " << z << "  ";
        } else if (fmt == Format::IJ) {
          outfile << i << "  " << j << "  ";
        } else if (fmt == Format::IK) {
          outfile << i << "  " << k << "  ";
        } else if (fmt == Format::IJK) {
          outfile << i << "  " << j << "  " << k << "  ";
        } else if (fmt == Format::YXZ) {
          outfile << y << "  " << x << "  " << z << "  ";
        }
        if (m_coordinates == Coordinates::Cylindrical) {
          const double ct = cos(y);
          const double st = sin(y);
          double wx = 0., wy = 0., wz = 0.;
          cmp->WeightingField(x * ct, x * st, z, wx, wy, wz, id);
          const double v = cmp->WeightingPotential(x * ct, x * st, z, id);
          const double wr = +wx * ct + wy * st;
          const double wt = -wx * st + wy * ct;
          outfile << wr << "  " << wt << "  " << wz << "  " << v << "\n";
        } else { 
          double wx = 0., wy = 0., wz = 0.;
          cmp->WeightingField(x, y, z, wx, wy, wz, id);
          const double v = cmp->WeightingPotential(x, y, z, id);
          outfile << wx << "  " << wy << "  " << wz << "  " << v << "\n";
        }
        ++nLines;
        if (nLines % nPrint == 0) PrintProgress(double(nLines) / nValues);
      }
    }
  }
  outfile.close();
  std::cout << std::endl << m_className << "::SaveWeightingField: Done.\n";
  return true;
}

bool ComponentGrid::LoadMesh(const std::string& filename, std::string format,
                             const double scaleX) {
  const auto fmt = GetFormat(format);
  if (fmt == Format::Unknown) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }

  // Keep track of which mesh parameters we have found.
  std::bitset<9> found;
  found.reset();
  double xmin = 0., ymin = 0., zmin = 0.;
  double xmax = 0., ymax = 0., zmax = 0.;
  unsigned int nx = 0, ny = 0, nz = 0;
  bool cylindrical = (m_coordinates == Coordinates::Cylindrical);
  // Parse the comment lines in the file.
  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }
  std::string line;
  unsigned int nLines = 0;
  // Read the file line by line.
  while (std::getline(infile, line)) {
    ++nLines;
    // Strip white space from the beginning of the line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip lines that are not comments.
    if (!IsComment(line)) continue;
    std::size_t pos0 = 0;
    std::size_t pos1 = line.find("=", pos0);
    while (pos1 != std::string::npos) {
      std::string key = line.substr(pos0, pos1 - pos0);
      std::transform(key.begin(), key.end(), key.begin(), toupper);
      const std::size_t pos2 = line.find_first_of(",;", pos1 + 1);
      std::istringstream val(line.substr(pos1 + 1, pos2 - pos1 - 1));
      if (key.find("XMIN") != std::string::npos) {
        val >> xmin;
        found.set(0);
      } else if (key.find("YMIN") != std::string::npos) {
        val >> ymin;
        found.set(1);
      } else if (key.find("ZMIN") != std::string::npos) {
        val >> zmin;
        found.set(2);
      } else if (key.find("XMAX") != std::string::npos) {
        val >> xmax;
        found.set(3);
      } else if (key.find("YMAX") != std::string::npos) {
        val >> ymax;
        found.set(4);
      } else if (key.find("ZMAX") != std::string::npos) {
        val >> zmax;
        found.set(5);
      } else if (key.find("NX") != std::string::npos) {
        val >> nx;
        found.set(6);
      } else if (key.find("NY") != std::string::npos) {
        val >> ny;
        found.set(7);
      } else if (key.find("NZ") != std::string::npos) {
        val >> nz;
        found.set(8);
      } else if (key.find("CYLINDRICAL") != std::string::npos) {
        cylindrical = true;
      }
      if (pos2 == std::string::npos) break;
      pos0 = pos2 + 1;
      pos1 = line.find("=", pos0);
    }
  }
  infile.close();

  if (fmt == Format::XY || fmt == Format::IJ) {
    // Try to complement missing information on the z-range.
    if (!found[8]) {
      nz = 1;
      found.set(8);
    }
    if (!found[2]) {
      if (found[0] || found[1] || found[3] || found[4] || found[5]) {
        zmin = -std::max(
            {fabs(xmin), fabs(xmax), fabs(ymin), fabs(ymax), fabs(zmax)});
      } else {
        zmin = -100.;
      }
      found.set(2);
    }
    if (!found[5]) {
      zmax = std::max(
          {fabs(xmin), fabs(xmax), fabs(ymin), fabs(ymax), fabs(zmin)});
      found.set(5);
    }
  } else if (fmt == Format::XZ || fmt == Format::IK) {
    // Try to complement missing information on the y/theta-range.
    if (!found[7]) {
      ny = 1;
      found.set(7);
    }
    if (cylindrical) {
      if (!found[1]) ymin = -Pi;
      if (!found[4]) ymax = +Pi;
      found.set(1);
      found.set(4);
    } else {
      if (!found[1]) {
        if (found[0] || found[2] || found[3] || found[4] || found[5]) {
          ymin = -std::max(
              {fabs(xmin), fabs(xmax), fabs(zmin), fabs(zmax), fabs(ymax)});
        } else {
          ymin = -100.;
        }
        found.set(1);
      }
      if (!found[4]) {
        ymax = std::max(
            {fabs(xmin), fabs(xmax), fabs(zmin), fabs(zmax), fabs(ymin)});
        found.set(4);
      }
    }
  }
  if (found.all()) {
    if (cylindrical && xmin < 0.) {
      std::cerr << m_className << "::LoadMesh:\n"
                << "    Radial coordinate must be >= 0.\n";
      return false;
    }
    std::cout << m_className << "::LoadMesh:\n";
    std::printf("%12.6f < x [cm] < %12.6f, %5u points\n", xmin, xmax, nx);
    std::printf("%12.6f < y [cm] < %12.6f, %5u points\n", ymin, ymax, ny);
    std::printf("%12.6f < z [cm] < %12.6f, %5u points\n", zmin, zmax, nz);
    if (cylindrical) {
      m_coordinates = Coordinates::Cylindrical;
    } else {
      m_coordinates = Coordinates::Cartesian;
    }
    return SetMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax);
  }
  if ((fmt == Format::IJ || fmt == Format::IJK || fmt == Format::IK) && 
      !(found[0] && found[3])) {
    std::cerr << m_className << "::LoadMesh: x-limits not found.\n";
    return false;
  } else if ((fmt == Format::IJ || fmt == Format::IJK) && 
             !(found[1] && found[4])) {
    std::cerr << m_className << "::LoadMesh: y-limits not found.\n";
    return false;
  } else if ((fmt == Format::IK || fmt == Format::IJK) && 
             !(found[2] && found[5])) {
    std::cerr << m_className << "::LoadMesh: z-limits not found.\n";
    return false;
  }

  unsigned int nValues = 0;
  infile.open(filename, std::ios::in);
  if (!infile) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  if (!found[0]) xmin = std::numeric_limits<double>::max();
  if (!found[1]) ymin = std::numeric_limits<double>::max();
  if (!found[2]) zmin = std::numeric_limits<double>::max();
  if (!found[3]) xmax = std::numeric_limits<double>::min();
  if (!found[4]) ymax = std::numeric_limits<double>::min();
  if (!found[5]) zmax = std::numeric_limits<double>::min();
  constexpr double tol = 1.e-10;
  auto cmp = [](double x, double y) {
    return x < y - tol * (std::fabs(x) + std::fabs(y));
  };
  std::set<double, decltype(cmp)> xLines(cmp);
  std::set<double, decltype(cmp)> yLines(cmp);
  std::set<double, decltype(cmp)> zLines(cmp);
  nLines = 0;
  bool bad = false;
  // Read the file line by line.
  while (std::getline(infile, line)) {
    ++nLines;
    // Strip white space from the beginning of the line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip comments.
    if (IsComment(line)) continue;
    std::istringstream data(line);
    if (fmt == Format::XY) {
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[1]) ymin = std::min(y, ymin);
      if (!found[3]) xmax = std::max(x, xmax);
      if (!found[4]) ymax = std::max(y, ymax);
      xLines.insert(x);
      yLines.insert(y);
    } else if (fmt == Format::XZ) {
      double x, z;
      data >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      z *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[2]) zmin = std::min(z, zmin);
      if (!found[3]) xmax = std::max(x, xmax);
      if (!found[5]) zmax = std::max(z, zmax);
      xLines.insert(x);
      zLines.insert(z);
    } else if (fmt == Format::XYZ) {
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[1]) ymin = std::min(y, ymin);
      if (!found[2]) zmin = std::min(z, zmin);
      if (!found[3]) xmax = std::max(x, xmax);
      if (!found[4]) ymax = std::max(y, ymax);
      if (!found[5]) zmax = std::max(z, zmax);
      xLines.insert(x);
      yLines.insert(y);
      zLines.insert(z);
    } else if (fmt == Format::IJ) {
      unsigned int i = 0, j = 0;
      data >> i >> j;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "indices");
        bad = true;
        break;
      }
      if (!found[6]) nx = std::max(nx, i);
      if (!found[7]) ny = std::max(ny, j);
    } else if (fmt == Format::IK) {
      unsigned int i = 0, k = 0;
      data >> i >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "indices");
        bad = true;
        break;
      }
      if (!found[6]) nx = std::max(nx, i);
      if (!found[8]) nz = std::max(nz, k);
    } else if (fmt == Format::IJK) {
      unsigned int i = 0, j = 0, k = 0;
      data >> i >> j >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "indices");
        bad = true;
        break;
      }
      if (!found[6]) nx = std::max(nx, i);
      if (!found[7]) ny = std::max(ny, j);
      if (!found[8]) nz = std::max(nz, k);
    } else if (fmt == Format::YXZ) {
      double x, y, z;
      data >> y >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadMesh", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (!found[0]) xmin = std::min(x, xmin);
      if (!found[1]) ymin = std::min(y, ymin);
      if (!found[2]) zmin = std::min(z, zmin);
      if (!found[3]) xmax = std::max(x, xmax);
      if (!found[4]) ymax = std::max(y, ymax);
      if (!found[5]) zmax = std::max(z, zmax);
      xLines.insert(x);
      yLines.insert(y);
      zLines.insert(z);
    }
    ++nValues;
  }
  infile.close();
  if (bad) return false;

  if (fmt == Format::XY || fmt == Format::XYZ || 
      fmt == Format::XZ || fmt == Format::YXZ) {
    if (!found[6]) nx = xLines.size();
    if (!found[7]) ny = yLines.size();
    if (!found[8]) nz = zLines.size();
  }

  if (cylindrical && xmin < 0.) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "    Radial coordinate must be >= 0.\n";
    return false;
  }
  std::cout << m_className << "::LoadMesh:\n";
  if (cylindrical) {
    std::printf("%12.6f < r [cm] < %12.6f, %5u points\n", xmin, xmax, nx);
    std::printf("%12.6f < theta  < %12.6f, %5u points\n", ymin, ymax, ny);
  } else {
    std::printf("%12.6f < x [cm] < %12.6f, %5u points\n", xmin, xmax, nx);
    std::printf("%12.6f < y [cm] < %12.6f, %5u points\n", ymin, ymax, ny);
  }
  std::printf("%12.6f < z [cm] < %12.6f, %5u points\n", zmin, zmax, nz);
  unsigned int nExpected = nx;
  if (fmt == Format::XZ || fmt == Format::IK) {
    nExpected *= nz;
  } else {
    nExpected *= ny;
  }
  if (fmt == Format::XYZ || fmt == Format::IJK || fmt == Format::YXZ) {
    nExpected *= nz;
  }
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadMesh:\n"
              << "   Warning: Expected " << nExpected << " lines, read "
              << nValues << ".\n";
  }
  if (cylindrical) {
    m_coordinates = Coordinates::Cylindrical;
  } else {
    m_coordinates = Coordinates::Cartesian;
  }
  return SetMesh(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax);
}

bool ComponentGrid::LoadData(
    const std::string& filename, std::string format, const bool withPotential,
    const bool withFlag, const double scaleX, const double scaleF,
    const double scaleP,
    std::vector<std::vector<std::vector<Node> > >& fields) {
  if (!m_hasMesh) {
    if (!LoadMesh(filename, format, scaleX)) {
      std::cerr << m_className << "::LoadData: Mesh not set.\n";
      return false;
    }
  }

  const auto fmt = GetFormat(format);
  if (fmt == Format::Unknown) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }

  // Set up the grid.
  Initialise(fields);

  unsigned int nValues = 0;
  // Keep track of which elements have been read.
  std::vector<std::vector<std::vector<bool> > > isSet(
      m_nX[0],
      std::vector<std::vector<bool> >(m_nX[1], std::vector<bool>(m_nX[2], false)));

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  std::string line;
  unsigned int nLines = 0;
  bool bad = false;
  // Read the file line by line.
  while (std::getline(infile, line)) {
    ++nLines;
    // Strip white space from beginning of line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip comments.
    if (IsComment(line)) continue;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    double fx = 0.;
    double fy = 0.;
    double fz = 0.;
    double p = 0.;
    std::istringstream data(line);
    if (fmt == Format::XY) {
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
        i = u < 0. ? 0 : static_cast<unsigned int>(u);
        if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[1] > 1) {
        const double v = std::round((y - m_xMin[1]) * m_sX[1]);
        j = v < 0. ? 0 : static_cast<unsigned int>(v);
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
    } else if (fmt == Format::XZ) {
      double x, z;
      data >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      z *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
        i = u < 0. ? 0 : static_cast<unsigned int>(u);
        if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[2] > 1) {
        const double v = std::round((z - m_xMin[2]) * m_sX[2]);
        k = v < 0. ? 0 : static_cast<unsigned int>(v);
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
    } else if (fmt == Format::XYZ) {
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
        i = u < 0. ? 0 : static_cast<unsigned int>(u);
        if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[1] > 1) {
        const double v = std::round((y - m_xMin[1]) * m_sX[1]);
        j = v < 0. ? 0 : static_cast<unsigned int>(v);
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
      if (m_nX[2] > 1) {
        const double w = std::round((z - m_xMin[2]) * m_sX[2]);
        k = w < 0. ? 0 : static_cast<unsigned int>(w);
        if (k >= m_nX[2]) k = m_nX[2] - 1;
      }
    } else if (fmt == Format::IJ) {
      data >> i >> j;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == Format::IJK) {
      data >> i >> j >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == Format::YXZ) {
      double x, y, z;
      data >> y >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
        i = u < 0. ? 0 : static_cast<unsigned int>(u);
        if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[1] > 1) {
        const double v = std::round((y - m_xMin[1]) * m_sX[1]);
        j = v < 0. ? 0 : static_cast<unsigned int>(v);
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
      if (m_nX[2] > 1) {
        const double w = std::round((z - m_xMin[2]) * m_sX[2]);
        k = w < 0. ? 0 : static_cast<unsigned int>(w);
        if (k >= m_nX[2]) k = m_nX[2] - 1;
      }
    }
    // Check the indices.
    if (i >= m_nX[0] || j >= m_nX[1] || k >= m_nX[2]) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Index (" << i << ", " << j << ", " << k
                << ") out of range.\n";
      continue;
    }
    if (isSet[i][j][k]) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Node (" << i << ", " << j << ", " << k
                << ") has already been set.\n";
      continue;
    }
    // Get the field values.
    if (fmt == Format::XY || fmt == Format::IJ) {
      // Two-dimensional map.
      fz = 0.;
      data >> fx >> fy;
    } else if (fmt == Format::XZ || fmt == Format::IK) {
      // Two-dimensional map.
      fy = 0.;
      data >> fx >> fz;
    } else if (fmt == Format::YXZ) {
      data >> fy >> fx >> fz;
    } else {
      data >> fx >> fy >> fz;
    }
    if (data.fail()) {
      PrintError(m_className + "::LoadData", nLines, "field components");
      bad = true;
      break;
    }
    fx *= scaleF;
    fy *= scaleF;
    fz *= scaleF;
    if (withPotential) {
      data >> p;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "potential");
        bad = true;
        break;
      }
      p *= scaleP;
      if (m_pMin > m_pMax) {
        // First value.
        m_pMin = p;
        m_pMax = p;
      } else {
        if (p < m_pMin) m_pMin = p;
        if (p > m_pMax) m_pMax = p;
      }
    }
    int flag = 0;
    if (withFlag) {
      data >> flag;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "region");
        bad = true;
        break;
      }
    }
    const bool isActive = flag == 0 ? false : true;
    if (fmt == Format::XY || fmt == Format::IJ) {
      // Two-dimensional map.
      for (unsigned int kk = 0; kk < m_nX[2]; ++kk) {
        fields[i][j][kk].fx = fx;
        fields[i][j][kk].fy = fy;
        fields[i][j][kk].fz = fz;
        fields[i][j][kk].v = p;
        if (withFlag) m_active[i][j][kk] = isActive;
        isSet[i][j][kk] = true;
      }
    } else if (fmt == Format::XZ || fmt == Format::IK) {
      // Two-dimensional map.
      for (unsigned int jj = 0; jj < m_nX[1]; ++jj) {
        fields[i][jj][k].fx = fx;
        fields[i][jj][k].fy = fy;
        fields[i][jj][k].fz = fz;
        fields[i][jj][k].v = p;
        if (withFlag) m_active[i][jj][k] = isActive;
        isSet[i][jj][k] = true;
      }
    } else {
      fields[i][j][k].fx = fx;
      fields[i][j][k].fy = fy;
      fields[i][j][k].fz = fz;
      fields[i][j][k].v = p;
      isSet[i][j][k] = true;
    }
    ++nValues;
  }
  infile.close();
  if (bad) return false;
  std::cout << m_className << "::LoadData:\n"
            << "    Read " << nValues << " values from " << filename << ".\n";
  unsigned int nExpected = m_nX[0];
  if (fmt == Format::XY || fmt == Format::IJ) {
    nExpected *= m_nX[1];
  } else if (fmt == Format::XZ || fmt == Format::IK) {
    nExpected *= m_nX[2];
  } else {
    nExpected *= m_nX[1] * m_nX[2];
  }
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadData:\n"
              << "   Expected " << nExpected << " values.\n";
  }
  return true;
}

bool ComponentGrid::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                   double& xmax, double& ymax, double& zmax) {
  if (m_efields.empty() && m_wfields.empty() && m_bfields.empty()) {
    return false;
  }
  if (m_coordinates == Coordinates::Cylindrical) {
    const double rmax = m_xMax[0];
    xmin = -rmax;
    ymin = -rmax;
    xmax = +rmax;
    ymax = +rmax;
    zmin = (m_periodic[2] || m_mirrorPeriodic[2]) ? -INFINITY : m_xMin[2];
    zmax = (m_periodic[2] || m_mirrorPeriodic[2]) ? +INFINITY : m_xMax[2];
    return true;
  }
  if (m_periodic[0] || m_mirrorPeriodic[0]) {
    xmin = -INFINITY;
    xmax = +INFINITY;
  } else {
    xmin = m_xMin[0];
    xmax = m_xMax[0];
  }

  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    ymin = -INFINITY;
    ymax = +INFINITY;
  } else {
    ymin = m_xMin[1];
    ymax = m_xMax[1];
  }

  if (m_periodic[2] || m_mirrorPeriodic[2]) {
    zmin = -INFINITY;
    zmax = +INFINITY;
  } else {
    zmin = m_xMin[2];
    zmax = m_xMax[2];
  }
  return true;
}

bool ComponentGrid::GetElementaryCell(
    double& xmin, double& ymin, double& zmin,
    double& xmax, double& ymax, double& zmax) {

  if (m_efields.empty() && m_wfields.empty() && m_bfields.empty()) {
    return false;
  }
  if (m_coordinates == Coordinates::Cylindrical) {
    const double rmax = m_xMax[0];
    xmin = -rmax;
    ymin = -rmax;
    xmax = +rmax;
    ymax = +rmax;
  } else {
    xmin = m_xMin[0];
    xmax = m_xMax[0];
    ymin = m_xMin[1];
    ymax = m_xMax[1];
  }
  zmin = m_xMin[2];
  zmax = m_xMax[2];
  return true;
}

bool ComponentGrid::GetVoltageRange(double& vmin, double& vmax) {
  if (m_efields.empty()) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

bool ComponentGrid::GetElectricFieldRange(double& exmin, double& exmax,
                                          double& eymin, double& eymax,
                                          double& ezmin, double& ezmax) {
  if (m_efields.empty()) {
    PrintNotReady(m_className + "::GetElectricFieldRange");
    return false;
  }

  exmin = exmax = m_efields[0][0][0].fx;
  eymin = eymax = m_efields[0][0][0].fy;
  ezmin = ezmax = m_efields[0][0][0].fz;
  for (unsigned int i = 0; i < m_nX[0]; ++i) {
    for (unsigned int j = 0; j < m_nX[1]; ++j) {
      for (unsigned int k = 0; k < m_nX[2]; ++k) {
        const Node& node = m_efields[i][j][k];
        if (node.fx < exmin) exmin = node.fx;
        if (node.fx > exmax) exmax = node.fx;
        if (node.fy < eymin) eymin = node.fy;
        if (node.fy > eymax) eymax = node.fy;
        if (node.fz < ezmin) ezmin = node.fz;
        if (node.fz > ezmax) ezmax = node.fz;
      }
    }
  }
  return true;
}

void ComponentGrid::SetMedium(Medium* m) {
  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
  }
  m_medium = m;
}

bool ComponentGrid::GetField(
    const double xi, const double yi, const double zi,
    const std::vector<std::vector<std::vector<Node> > >& field, double& fx,
    double& fy, double& fz, double& p, bool& active) {
  if (!m_hasMesh) {
    std::cerr << m_className << "::GetField: Mesh is not set.\n";
    return false;
  }

  // Reduce the point to the basic cell (in case of periodicity) and
  // check if it is inside the mesh.
  std::array<bool, 3> mirrored = {false, false, false};
  std::array<double, 3> xx = {xi, yi, zi};
  double theta = 0.;
  if (m_coordinates == Coordinates::Cylindrical) {
    if (fabs(xi) > Small || fabs(yi) > Small) {
      theta = atan2(yi, xi);
      xx[0] = sqrt(xi * xi + yi * yi);
      xx[1] = theta;
    } 
  }
  for (size_t i = 0; i < 3; ++i) {
    xx[i] = Reduce(xx[i], m_xMin[i], m_xMax[i], 
                   m_periodic[i], m_mirrorPeriodic[i], mirrored[i]);
    if (xx[i] < m_xMin[i] || xx[i] > m_xMax[i]) return false;
  }
  // Get the indices.
  const double sx = (xx[0] - m_xMin[0]) * m_sX[0];
  const double sy = (xx[1] - m_xMin[1]) * m_sX[1];
  const double sz = (xx[2] - m_xMin[2]) * m_sX[2];
  // Get the indices.
  const unsigned int i0 = static_cast<unsigned int>(std::floor(sx));
  const unsigned int j0 = static_cast<unsigned int>(std::floor(sy));
  const unsigned int k0 = static_cast<unsigned int>(std::floor(sz));
  const double ux = sx - i0;
  const double uy = sy - j0;
  const double uz = sz - k0;
  const unsigned int i1 = std::min(i0 + 1, m_nX[0] - 1);
  const unsigned int j1 = std::min(j0 + 1, m_nX[1] - 1);
  const unsigned int k1 = std::min(k0 + 1, m_nX[2] - 1);
  const double vx = 1. - ux;
  const double vy = 1. - uy;
  const double vz = 1. - uz;
  if (!m_active.empty()) {
    active = m_active[i0][j0][k0] && m_active[i0][j0][k1] &&
             m_active[i0][j1][k0] && m_active[i0][j1][k1] &&
             m_active[i1][j0][k0] && m_active[i1][j0][k1] &&
             m_active[i1][j1][k0] && m_active[i1][j1][k1];
  }
  const Node& n000 = field[i0][j0][k0];
  const Node& n100 = field[i1][j0][k0];
  const Node& n010 = field[i0][j1][k0];
  const Node& n110 = field[i1][j1][k0];
  const Node& n001 = field[i0][j0][k1];
  const Node& n101 = field[i1][j0][k1];
  const Node& n011 = field[i0][j1][k1];
  const Node& n111 = field[i1][j1][k1];

  if (m_debug) {
    std::cout << m_className << "::GetField: Determining field at (" << xi
              << ", " << yi << ", " << zi << ").\n"
              << "    X: " << i0 << " (" << ux << ") - " << i1 << " (" << vx
              << ").\n"
              << "    Y: " << j0 << " (" << uy << ") - " << j1 << " (" << vy
              << ").\n"
              << "    Z: " << k0 << " (" << uz << ") - " << k1 << " (" << vz
              << ").\n";
  }
  fx = ((n000.fx * vx + n100.fx * ux) * vy +
        (n010.fx * vx + n110.fx * ux) * uy) *
           vz +
       ((n001.fx * vx + n101.fx * ux) * vy +
        (n011.fx * vx + n111.fx * ux) * uy) *
           uz;
  fy = ((n000.fy * vx + n100.fy * ux) * vy +
        (n010.fy * vx + n110.fy * ux) * uy) *
           vz +
       ((n001.fy * vx + n101.fy * ux) * vy +
        (n011.fy * vx + n111.fy * ux) * uy) *
           uz;
  fz = ((n000.fz * vx + n100.fz * ux) * vy +
        (n010.fz * vx + n110.fz * ux) * uy) *
           vz +
       ((n001.fz * vx + n101.fz * ux) * vy +
        (n011.fz * vx + n111.fz * ux) * uy) *
           uz;
  p = ((n000.v * vx + n100.v * ux) * vy + (n010.v * vx + n110.v * ux) * uy) *
          vz +
      ((n001.v * vx + n101.v * ux) * vy + (n011.v * vx + n111.v * ux) * uy) *
          uz;
  if (mirrored[0]) fx = -fx;
  if (mirrored[1]) fy = -fy;
  if (mirrored[2]) fz = -fz;
  if (m_coordinates == Coordinates::Cylindrical) {
    const double ct = cos(theta);
    const double st = sin(theta);
    const double fr = fx;
    const double ft = fy;
    fx = ct * fr - st * ft;
    fy = st * fr + ct * ft;
  }
  return true;
}

bool ComponentGrid::GetElectricField(const unsigned int i, const unsigned int j,
                                     const unsigned int k, double& v,
                                     double& ex, double& ey, double& ez) const {
  v = ex = ey = ez = 0.;
  if (m_efields.empty()) {
    if (!m_hasMesh) {
      std::cerr << m_className << "::GetElectricField: Mesh not set.\n";
      return false;
    }
    PrintNotReady(m_className + "::GetElectricField");
    return false;
  }
  if (i >= m_nX[0] || j >= m_nX[1] || k >= m_nX[2]) {
    std::cerr << m_className << "::GetElectricField: Index out of range.\n";
    return false;
  }
  const Node& node = m_efields[i][j][k];
  v = node.v;
  ex = node.fx;
  ey = node.fy;
  ez = node.fz;
  return true;
}

void ComponentGrid::Print() {

  std::cout << m_className << "::Print:\n";
  if (!m_hasMesh) {
    std::cout << "    Mesh not set.\n";
    return;
  }
  std::printf("    %15.8f < x [cm] < %15.8f, %10u nodes\n", 
              m_xMin[0], m_xMax[0], m_nX[0]); 
  std::printf("    %15.8f < y [cm] < %15.8f, %10u nodes\n", 
              m_xMin[1], m_xMax[1], m_nX[1]); 
  std::printf("    %15.8f < z [cm] < %15.8f, %10u nodes\n", 
              m_xMin[2], m_xMax[2], m_nX[2]);
  if (m_efields.empty() && m_bfields.empty() && 
      m_wfields.empty() && m_wdfields.empty() && 
      m_eAttachment.empty() && m_hAttachment.empty() &&
      m_eVelocity.empty() && m_hVelocity.empty()) {
    std::cout << "    Available data: None.\n";
    return;
  }
  std::cout << "    Available data:\n";
  if (!m_efields.empty()) std::cout << "      Electric field.\n"; 
  if (!m_bfields.empty()) std::cout << "      Magnetic field.\n"; 
  if (!m_wfields.empty()) std::cout << "      Weighting field.\n"; 
  if (!m_wdfields.empty()) {
    std::cout << "      Delayed weighting field.\n";
  }
  if (!m_eVelocity.empty()) {
    std::cout << "      Electron drift velocity.\n";
  }
  if (!m_hVelocity.empty()) {
    std::cout << "      Hole drift velocity.\n";
  }
  if (!m_eAttachment.empty()) {
    std::cout << "      Electron attachment coefficient.\n";
  }
  if (!m_hAttachment.empty()) {
    std::cout << "      Hole attachment coefficient.\n";
  }
}

void ComponentGrid::Reset() {
  m_efields.clear();
  m_bfields.clear();
  m_wfields.clear();
  m_eAttachment.clear();
  m_hAttachment.clear();
  m_eVelocity.clear();
  m_hVelocity.clear();

  m_wdfields.clear();
  m_wdtimes.clear();

  m_active.clear();

  m_nX.fill(1);
  m_xMin.fill(0.);
  m_xMax.fill(0.);
  m_sX[0] = m_sX[1] = m_sX[2] = 0.;
  m_pMin = m_pMax = 0.;
  m_medium = nullptr;

  m_hasMesh = false;
  m_hasPotential = false;

  m_wFieldOffset.fill(0.);
}

void ComponentGrid::UpdatePeriodicity() {

  // Check for conflicts.
  for (size_t i = 0; i < 3; ++i) {
    if (m_periodic[i] && m_mirrorPeriodic[i]) {
      std::cerr << m_className << "::UpdatePeriodicity:\n"
                << "    Both simple and mirror periodicity requested. Reset.\n";
      m_periodic[i] = m_mirrorPeriodic[i] = false;
    }
  }

  if (m_axiallyPeriodic[0] || m_axiallyPeriodic[1] || m_axiallyPeriodic[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Axial symmetry is not supported. Reset.\n";
    m_axiallyPeriodic.fill(false);
  }

  if (m_rotationSymmetric[0] || m_rotationSymmetric[1] ||
      m_rotationSymmetric[2]) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Rotation symmetry is not supported. Reset.\n";
    m_rotationSymmetric.fill(false);
  }
}

double ComponentGrid::Reduce(const double xin, const double xmin,
                             const double xmax, const bool simplePeriodic,
                             const bool mirrorPeriodic, bool& mirrored) const {
  // In case of periodicity, reduce the coordinate to the basic cell.
  double x = xin;
  const double lx = xmax - xmin;
  if (simplePeriodic) {
    x = xmin + fmod(x - xmin, lx);
    if (x < xmin) x += lx;
  } else if (mirrorPeriodic) {
    double xNew = xmin + fmod(x - xmin, lx);
    if (xNew < xmin) xNew += lx;
    const int nx = int(floor(0.5 + (xNew - x) / lx));
    if (nx != 2 * (nx / 2)) {
      xNew = xmin + xmax - xNew;
      mirrored = true;
    }
    x = xNew;
  }
  return x;
}

void ComponentGrid::Initialise(
    std::vector<std::vector<std::vector<Node> > >& fields) {
  fields.resize(m_nX[0]);
  for (unsigned int i = 0; i < m_nX[0]; ++i) {
    fields[i].resize(m_nX[1]);
    for (unsigned int j = 0; j < m_nX[1]; ++j) {
      fields[i][j].resize(m_nX[2]);
      for (unsigned int k = 0; k < m_nX[2]; ++k) {
        fields[i][j][k].fx = 0.;
        fields[i][j][k].fy = 0.;
        fields[i][j][k].fz = 0.;
        fields[i][j][k].v = 0.;
      }
    }
  }
}

bool ComponentGrid::LoadElectronVelocity(const std::string& fname, 
                                         const std::string& fmt,
                                         const double scaleX,
                                         const double scaleV) {
  // Read the file.
  if (!LoadData(fname, fmt, false, false, scaleX, scaleV, 1., m_eVelocity)) {
    return false;
  }
  return true;
}

bool ComponentGrid::LoadHoleVelocity(const std::string& fname, 
                                     const std::string& fmt,
                                     const double scaleX,
                                     const double scaleV) {
  // Read the file.
  if (!LoadData(fname, fmt, false, false, scaleX, scaleV, 1., m_hVelocity)) {
    return false;
  }
  return true;
}

bool ComponentGrid::ElectronVelocity(const double x, const double y, 
                                     const double z,
                                     double& vx, double& vy, double& vz) {
  if (m_eVelocity.empty()) {
    PrintNotReady(m_className + "::ElectronVelocity");
    return false;
  }
  double p = 0.;
  bool active = true;
  return GetField(x, y, z, m_eVelocity, vx, vy, vz, p, active);
}

bool ComponentGrid::HoleVelocity(const double x, const double y, 
                                 const double z,
                                 double& vx, double& vy, double& vz) {
  if (m_hVelocity.empty()) {
    PrintNotReady(m_className + "::HoleVelocity");
    return false;
  }
  double p = 0.;
  bool active = true;
  return GetField(x, y, z, m_hVelocity, vx, vy, vz, p, active);
}

bool ComponentGrid::LoadElectronAttachment(const std::string& fname,
                                           const std::string& fmt, 
                                           const unsigned int col,
                                           const double scaleX) {
  // Read the file.
  return LoadData(fname, fmt, scaleX, m_eAttachment, col);
}

bool ComponentGrid::LoadHoleAttachment(const std::string& fname,
                                       const std::string& fmt, 
                                       const unsigned int col,
                                       const double scaleX) {
  // Read the file.
  return LoadData(fname, fmt, scaleX, m_hAttachment, col);
}

bool ComponentGrid::LoadData(
    const std::string& filename, std::string format, const double scaleX,
    std::vector<std::vector<std::vector<double> > >& tab, 
    const unsigned int col) {
  if (!m_hasMesh) {
    if (!LoadMesh(filename, format, scaleX)) {
      std::cerr << m_className << "::LoadData: Mesh not set.\n";
      return false;
    }
  }

  const auto fmt = GetFormat(format);
  if (fmt == Format::Unknown) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Unknown format (" << format << ").\n";
    return false;
  }
  // Check the column index.
  unsigned int offset = 0;
  if (fmt == Format::XY || fmt == Format::IJ) {
    if (col < 2) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Unexpected column index (" << col << ").\n";
      return false; 
    }
    offset = 2;
  } else {
    if (col < 3) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Unexpected column index (" << col << ").\n";
      return false; 
    }
    offset = 3;
  } 

  // Set up the grid.
  tab.assign(
      m_nX[0], 
      std::vector<std::vector<double> >(m_nX[1], std::vector<double>(m_nX[2], 0.)));

  unsigned int nValues = 0;
  // Keep track of which elements have been read.
  std::vector<std::vector<std::vector<bool> > > isSet(
      m_nX[0],
      std::vector<std::vector<bool> >(m_nX[1], std::vector<bool>(m_nX[2], false)));

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  std::string line;
  unsigned int nLines = 0;
  bool bad = false;
  // Read the file line by line.
  while (std::getline(infile, line)) {
    ++nLines;
    // Strip white space from beginning of line.
    ltrim(line);
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip comments.
    if (IsComment(line)) continue;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    double val = 0;
    std::istringstream data(line);
    if (fmt == Format::XY) {
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
         i = u < 0. ? 0 : static_cast<unsigned int>(u);
         if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[1] > 1) {
        const double v = std::round((y - m_xMin[1]) * m_sX[1]);
        j = v < 0. ? 0 : static_cast<unsigned int>(v);
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
    } else if (fmt == Format::XYZ) {
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
        i = u < 0. ? 0 : static_cast<unsigned int>(u);
        if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[1] > 1) {
        const double v = std::round((y - m_xMin[1]) * m_sX[1]);
        j = v < 0. ? 0 : static_cast<unsigned int>(v); 
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
      if (m_nX[2] > 1) {
        const double w = std::round((z - m_xMin[2]) * m_sX[2]);
        k = w < 0. ? 0 : static_cast<unsigned int>(w);
        if (k >= m_nX[2]) k = m_nX[2] - 1;
      }
    } else if (fmt == Format::IJ) {
      data >> i >> j;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == Format::IJK) {
      data >> i >> j >> k;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "indices");
        bad = true;
        break;
      }
    } else if (fmt == Format::YXZ) {
      double x, y, z;
      data >> y >> x >> z;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, "coordinates");
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      if (m_nX[0] > 1) {
        const double u = std::round((x - m_xMin[0]) * m_sX[0]);
        i = u < 0. ? 0 : static_cast<unsigned int>(u);
        if (i >= m_nX[0]) i = m_nX[0] - 1;
      }
      if (m_nX[1] > 1) {
        const double v = std::round((y - m_xMin[1]) * m_sX[1]);
        j = v < 0. ? 0 : static_cast<unsigned int>(v);
        if (j >= m_nX[1]) j = m_nX[1] - 1;
      }
      if (m_nX[2] > 1) {
        const double w = std::round((z - m_xMin[2]) * m_sX[2]);
        k = w < 0. ? 0 : static_cast<unsigned int>(w);
        if (k >= m_nX[2]) k = m_nX[2] - 1;
      }
    }
    // Check the indices.
    if (i >= m_nX[0] || j >= m_nX[1] || k >= m_nX[2]) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Index (" << i << ", " << j << ", " << k
                << ") out of range.\n";
      continue;
    }
    if (isSet[i][j][k]) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Node (" << i << ", " << j << ", " << k
                << ") has already been set.\n";
      continue;
    }

    // Skip to the requested column.
    for (unsigned int ii = 0; ii < col - offset; ++ii) {
      double dummy = 0.;
      data >> dummy;
      if (data.fail()) {
        PrintError(m_className + "::LoadData", nLines, 
                   "column " + std::to_string(offset + ii));
        break;
      }
    }
    if (data.fail()) {
      bad = true;
      break;
    }
    data >> val;

    if (data.fail()) {
      PrintError(m_className + "::LoadData", nLines, 
                 "column " + std::to_string(col));
      bad = true;
      break;
    }

    if (fmt == Format::XY || fmt == Format::IJ) {
      // Two-dimensional map
      for (unsigned int kk = 0; kk < m_nX[2]; ++kk) {
        tab[i][j][kk] = val;
        isSet[i][j][kk] = true;
      }
    } else {
      tab[i][j][k] = val;
      isSet[i][j][k] = true;
    }
    ++nValues;
  }
  infile.close();
  if (bad) return false;
  std::cout << m_className << "::LoadData:\n"
            << "    Read " << nValues << " values from " << filename << ".\n";
  unsigned int nExpected = m_nX[0] * m_nX[1];
  if (fmt == Format::XYZ || fmt == Format::IJK || fmt == Format::YXZ) {
    nExpected *= m_nX[2];
  }
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadData:\n"
              << "   Expected " << nExpected << " values.\n";
  }
  return true;
}

bool ComponentGrid::GetData(
    const double xi, const double yi, const double zi,
    const std::vector<std::vector<std::vector<double> > >& tab, double& val) {
  if (!m_hasMesh) {
    std::cerr << m_className << "::GetData: Mesh is not set.\n";
    return false;
  }

  // Reduce the point to the basic cell (in case of periodicity) and
  // check if it is inside the mesh.
  bool xMirrored = false;
  const double x =
      Reduce(xi, m_xMin[0], m_xMax[0], m_periodic[0], m_mirrorPeriodic[0], xMirrored);
  if (x < m_xMin[0] || x > m_xMax[0]) return false;
  bool yMirrored = false;
  const double y =
      Reduce(yi, m_xMin[1], m_xMax[1], m_periodic[1], m_mirrorPeriodic[1], yMirrored);
  if (y < m_xMin[1] || y > m_xMax[1]) return false;
  bool zMirrored = false;
  const double z =
      Reduce(zi, m_xMin[2], m_xMax[2], m_periodic[2], m_mirrorPeriodic[2], zMirrored);
  if (z < m_xMin[2] || z > m_xMax[2]) return false;

  // Get the indices.
  const double sx = (x - m_xMin[0]) * m_sX[0];
  const double sy = (y - m_xMin[1]) * m_sX[1];
  const double sz = (z - m_xMin[2]) * m_sX[2];
  const unsigned int i0 = static_cast<unsigned int>(std::floor(sx));
  const unsigned int j0 = static_cast<unsigned int>(std::floor(sy));
  const unsigned int k0 = static_cast<unsigned int>(std::floor(sz));
  const double ux = sx - i0;
  const double uy = sy - j0;
  const double uz = sz - k0;
  const unsigned int i1 = std::min(i0 + 1, m_nX[0] - 1);
  const unsigned int j1 = std::min(j0 + 1, m_nX[1] - 1);
  const unsigned int k1 = std::min(k0 + 1, m_nX[2] - 1);
  const double vx = 1. - ux;
  const double vy = 1. - uy;
  const double vz = 1. - uz;
  const double n000 = tab[i0][j0][k0];
  const double n100 = tab[i1][j0][k0];
  const double n010 = tab[i0][j1][k0];
  const double n110 = tab[i1][j1][k0];
  const double n001 = tab[i0][j0][k1];
  const double n101 = tab[i1][j0][k1];
  const double n011 = tab[i0][j1][k1];
  const double n111 = tab[i1][j1][k1];

  if (m_debug) {
    std::cout << m_className << "::GetData: Interpolating at (" << xi
              << ", " << yi << ", " << zi << ").\n"
              << "    X: " << i0 << " (" << ux << ") - " << i1 << " (" << vx
              << ").\n"
              << "    Y: " << j0 << " (" << uy << ") - " << j1 << " (" << vy
              << ").\n"
              << "    Z: " << k0 << " (" << uz << ") - " << k1 << " (" << vz
              << ").\n";
  }
  val = ((n000 * vx + n100 * ux) * vy + (n010 * vx + n110 * ux) * uy) * vz +
        ((n001 * vx + n101 * ux) * vy + (n011 * vx + n111 * ux) * uy) * uz;

  return true;
}

bool ComponentGrid::ElectronAttachment(const double x, const double y,
                                       const double z, double& att) {
  // Make sure the map has been loaded.
  if (m_eAttachment.empty()) {
    PrintNotReady(m_className + "::ElectronAttachment");
    return false;
  }
  return GetData(x, y, z, m_eAttachment, att);
}

bool ComponentGrid::HoleAttachment(const double x, const double y,
                                   const double z, double& att) {
  // Make sure the map has been loaded.
  if (m_hAttachment.empty()) {
    PrintNotReady(m_className + "::HoleAttachment");
    return false;
  }
  return GetData(x, y, z, m_hAttachment, att);
}

ComponentGrid::Format ComponentGrid::GetFormat(std::string format) {
  std::transform(format.begin(), format.end(), format.begin(), toupper);
  if (format == "XY") {
    return Format::XY;
  } else if (format == "XZ") {
    return Format::XZ;
  } else if (format == "XYZ") {
    return Format::XYZ;
  } else if (format == "IJ") {
    return Format::IJ;
  } else if (format == "IK") {
    return Format::IK;
  } else if (format == "IJK") {
    return Format::IJK;
  } else if (format == "YXZ") {
    return Format::YXZ;
  }
  return Format::Unknown;
}

}  // namespace Garfield
