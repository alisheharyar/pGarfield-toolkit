#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Garfield/ComponentVoxel.hh"
#include "Garfield/Utilities.hh"
#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

ComponentVoxel::ComponentVoxel() : Component("Voxel") {}

void ComponentVoxel::ElectricField(
    const double x, const double y, const double z, 
    double& ex, double& ey, double& ez, double& p, Medium*& m, int& status) {
  m = nullptr;
  status = 0;

  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::ElectricField:\n"
              << "    Field map is not available for interpolation.\n";
    status = -10;
    return;
  }

  status = 0;
  int region = -1;
  if (!GetField(x, y, z, m_efields, ex, ey, ez, p, region)) {
    status = -6;
    return;
  }

  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  m = m_media[region];
  if (!m) status = -5;
}

void ComponentVoxel::ElectricField(
    const double x, const double y, const double z, 
    double& ex, double& ey, double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentVoxel::WeightingField(
    const double x, const double y, const double z, 
    double& wx, double& wy, double& wz, const std::string& /*label*/) {

  wx = wy = wz = 0.;
  if (!m_hasWfield) return;
  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;
  double wp = 0.;
  int region = 0;
  GetField(xx, yy, zz, m_wfields, wx, wy, wz, wp, region);
}

double ComponentVoxel::WeightingPotential(
    const double x, const double y, const double z,
    const std::string& /*label*/) {

  if (!m_hasWfield) return 0.;
  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;
  double wx = 0., wy = 0., wz = 0.;
  double wp = 0.;
  int region = 0;
  if (!GetField(xx, yy, zz, m_wfields, wx, wy, wz, wp, region)) return 0.;
  return wp;
}

void ComponentVoxel::DelayedWeightingField(
    const double x, const double y, const double z, const double t,
    double& wx, double& wy, double& wz, const std::string& /*label*/) {

  wx = wy = wz = 0.;
  if (m_wdtimes.empty()) return;
  // Assume no weighting field for times outside the range of available maps.
  if (t < m_wdtimes.front() || t > m_wdtimes.back()) return;

  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;

  const auto it1 = std::upper_bound(m_wdtimes.cbegin(), m_wdtimes.cend(), t);
  const auto it0 = std::prev(it1);
 
  const double dt = t - *it0; 
  double wp = 0.;
  int region = 0;
  const unsigned int i0 = it0 - m_wdtimes.cbegin();
  double wx0 = 0., wy0 = 0., wz0 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i0], wx0, wy0, wz0, wp, region)) {
    return;
  } 
  if (dt < Small || it1 == m_wdtimes.cend()) {
    wx = wx0;
    wy = wy0;
    wz = wz0;
    return; 
  }
  const unsigned int i1 = it1 - m_wdtimes.cbegin();
  double wx1 = 0., wy1 = 0., wz1 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i1], wx1, wy1, wz1, wp, region)) {
    return;
  } 
  const double f1 = dt / (*it1 - *it0);
  const double f0 = 1. - f1;
  wx = f0 * wx0 + f1 * wx1;
  wy = f0 * wy0 + f1 * wy1;
  wz = f0 * wz0 + f1 * wz1;
}

double ComponentVoxel::DelayedWeightingPotential(
    const double x, const double y, const double z, const double t,
    const std::string& /*label*/) {

  if (m_wdtimes.empty()) return 0.;
  // Outside the range of the available maps?
  if (t < m_wdtimes.front() || t > m_wdtimes.back()) return 0.;

  const double xx = x - m_wField_xOffset;
  const double yy = y - m_wField_yOffset;
  const double zz = z - m_wField_zOffset;

  const auto it1 = std::upper_bound(m_wdtimes.cbegin(), m_wdtimes.cend(), t);
  const auto it0 = std::prev(it1);
 
  const double dt = t - *it0; 
  int region = 0;
  const unsigned int i0 = it0 - m_wdtimes.cbegin();
  double wp0 = 0., wx0 = 0., wy0 = 0., wz0 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i0], wx0, wy0, wz0, wp0, region)) {
    return 0.;
  } 
  if (dt < Small || it1 == m_wdtimes.cend()) return 0.;

  const unsigned int i1 = it1 - m_wdtimes.cbegin();
  double wp1 = 0., wx1 = 0., wy1 = 0., wz1 = 0.;
  if (!GetField(xx, yy, zz, m_wdfields[i1], wx1, wy1, wz1, wp1, region)) {
    return 0.;
  } 
  const double f1 = dt / (*it1 - *it0);
  const double f0 = 1. - f1;
  return f0 * wp0 + f1 * wp1;
}

void ComponentVoxel::SetWeightingFieldOffset(const double x, const double y,
                                             const double z) {
  m_wField_xOffset = x;
  m_wField_yOffset = y;
  m_wField_zOffset = z;
}

void ComponentVoxel::MagneticField(const double x, const double y,
                                   const double z, double& bx, double& by,
                                   double& bz, int& status) {
  status = 0;
  if (!m_hasBfield) {
    return Component::MagneticField(x, y, z, bx, by, bz, status);
  }

  int region = -1;
  double p = 0.;
  if (!GetField(x, y, z, m_bfields, bx, by, bz, p, region)) {
    status = -6;
  }
}

bool ComponentVoxel::HasMagneticField() const {
  return m_hasBfield ? true : Component::HasMagneticField();
}

Medium* ComponentVoxel::GetMedium(const double x, const double y,
                                  const double z) {
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::GetMedium:\n"
              << "    Field map is not available for interpolation.\n";
    return nullptr;
  }

  unsigned int i, j, k;
  bool xMirrored, yMirrored, zMirrored;
  if (!GetElement(x, y, z, i, j, k, xMirrored, yMirrored, zMirrored)) {
    return nullptr;
  }
  const int region = m_regions[i][j][k];
  if (region < 0 || region > (int)m_media.size()) return nullptr;
  return m_media[region];
}

void ComponentVoxel::SetMesh(const unsigned int nx, const unsigned int ny,
                             const unsigned int nz, const double xmin,
                             const double xmax, const double ymin,
                             const double ymax, const double zmin,
                             const double zmax) {
  Reset();
  if (nx == 0 || ny == 0 || nz == 0) {
    std::cerr << m_className << "::SetMesh:\n"
              << "    Number of mesh elements must be positive.\n";
    return;
  }
  if (xmin >= xmax) {
    std::cerr << m_className << "::SetMesh: Invalid x range.\n";
    return;
  } else if (ymin >= ymax) {
    std::cerr << m_className << "::SetMesh: Invalid y range.\n";
    return;
  } else if (zmin >= zmax) {
    std::cerr << m_className << "::SetMesh: Invalid z range.\n";
    return;
  }
  m_nX = nx;
  m_nY = ny;
  m_nZ = nz;
  m_xMin = xmin;
  m_yMin = ymin;
  m_zMin = zmin;
  m_xMax = xmax;
  m_yMax = ymax;
  m_zMax = zmax;
  m_dx = (m_xMax - m_xMin) / m_nX;
  m_dy = (m_yMax - m_yMin) / m_nY;
  m_dz = (m_zMax - m_zMin) / m_nZ;
  m_hasMesh = true;
}

bool ComponentVoxel::LoadElectricField(const std::string& fname,
                                       const std::string& fmt,
                                       const bool withP, const bool withR,
                                       const double scaleX, const double scaleE,
                                       const double scaleP) {
  m_ready = false;
  m_efields.clear();
  m_hasPotential = m_hasEfield = false;
  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadElectricField:\n"
              << "    Mesh is not set. Call SetMesh first.\n";
    return false;
  }

  // Set up the grid.
  Initialise(m_efields);
  InitialiseRegions();

  m_pMin = m_pMax = 0.;
  if (withP) {
    m_pMin = 1.;
    m_pMax = -1.;
  }
  if (!LoadData(fname, fmt, withP, withR, scaleX, scaleE, scaleP, m_efields)) {
    return false;
  }
  m_hasEfield = true;
  m_ready = true;
  if (withP) m_hasPotential = true;
  return true;
}

bool ComponentVoxel::LoadWeightingField(const std::string& fname, 
                                        const std::string& fmt,
                                        const bool withP, 
                                        const double scaleX, 
                                        const double scaleE,
                                        const double scaleP) {
  m_hasWfield = false;
  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadWeightingField:\n"
              << "    Mesh is not set. Call SetMesh first.\n";
    return false;
  }

  // Set up the grid.
  Initialise(m_wfields);
  if (m_regions.empty()) InitialiseRegions();

  // Read the file.
  if (!LoadData(fname, fmt, withP, false, scaleX, scaleE, scaleP, m_wfields)) {
    return false;
  }
  m_hasWfield = true;
  return true;
}
bool ComponentVoxel::LoadWeightingField(const std::string& fname, 
                                        const std::string& fmt,
                                        const double t, const bool withP, 
                                        const double scaleX, const double scaleE,
                                        const double scaleP) {

  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadWeightingField:\n"
              << "    Mesh is not set. Call SetMesh first.\n";
    return false;
  }

  std::vector<std::vector<std::vector<Element> > > wfield;
  Initialise(wfield);
  if (m_regions.empty()) InitialiseRegions();
 
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
    m_wdfields.insert(m_wdfields.begin() + n, std::move(wfield));
  }
  return true;
}

bool ComponentVoxel::LoadMagneticField(const std::string& fname,
                                       const std::string& fmt,
                                       const double scaleX,
                                       const double scaleB) {
  m_hasBfield = false;
  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadMagneticField:\n"
              << "    Mesh is not set. Call SetMesh first.\n";
    return false;
  }

  // Set up the grid.
  Initialise(m_bfields);
  InitialiseRegions();

  // Read the file.
  if (!LoadData(fname, fmt, false, false, scaleX, scaleB, 1., m_bfields)) {
    return false;
  }
  m_hasBfield = true;
  return true;
}

bool ComponentVoxel::LoadData(const std::string& filename, std::string format,
    const bool withPotential, const bool withRegion,
    const double scaleX, const double scaleF, const double scaleP,
    std::vector<std::vector<std::vector<Element> > >& fields) {

  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadData: Mesh has not been set.\n";
    return false;
  }

  unsigned int nValues = 0;
  // Keep track of which elements have been read.
  std::vector<std::vector<std::vector<bool> > > isSet(
      m_nX,
      std::vector<std::vector<bool> >(m_nY, std::vector<bool>(m_nZ, false)));

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << m_className << "::LoadData:\n"
              << "    Could not open file " << filename << ".\n";
    return false;
  }

  std::transform(format.begin(), format.end(), format.begin(), toupper);
  unsigned int fmt = 0;
  if (format == "XY") {
    fmt = 1;
  } else if (format == "XYZ") {
    fmt = 2;
  } else if (format == "IJ") {
    fmt = 3;
  } else if (format == "IJK") {
    fmt = 4;
  } else if (format == "YXZ") {
    fmt = 5;
  } else {
    std::cerr << m_className << "::LoadData:\n"
              << "    Unknown format (" << format << ").\n";
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
    if (line[0] == '#') continue;
    if (line[0] == '/' && line[1] == '/') continue;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    double fx = 0.;
    double fy = 0.;
    double fz = 0.;
    double v = 0.;
    int region = 0;
    std::istringstream data(line);
    if (fmt == 1) {
      // "XY"
      double x, y;
      data >> x >> y;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot retrieve element coordinates.\n";
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      const double z = 0.5 * (m_zMin + m_zMax);
      bool xMirrored, yMirrored, zMirrored;
      if (!GetElement(x, y, z, i, j, k, xMirrored, yMirrored, zMirrored)) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Point is outside mesh.\n";
        bad = true;
        break;
      }
    } else if (fmt == 2) {
      // "XYZ"
      double x, y, z;
      data >> x >> y >> z;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot retrieve element coordinates.\n";
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      bool xMirrored, yMirrored, zMirrored;
      if (!GetElement(x, y, z, i, j, k, xMirrored, yMirrored, zMirrored)) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Point is outside mesh.\n";
        bad = true;
        break;
      }
    } else if (fmt == 3) {
      // "IJ"
      k = 0;
      data >> i >> j;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot retrieve element index.\n";
        bad = true;
        break;
      }
    } else if (fmt == 4) {
      // "IJK"
      data >> i >> j >> k;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot retrieve element index.\n";
        bad = true;
        break;
      }
    } else if (fmt == 5) {
      // "YXZ"
      double x, y, z, temp;
      data >> y >> x >> temp;
      z = temp;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot retrieve element coordinates.\n";
        bad = true;
        break;
      }
      x *= scaleX;
      y *= scaleX;
      z *= scaleX;
      bool xMirrored, yMirrored, zMirrored;
      if (!GetElement(x, y, z, i, j, k, xMirrored, yMirrored, zMirrored)) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Point is outside mesh.\n";
        bad = true;
        break;
      }
    }
    // Check the indices.
    if (i >= m_nX || j >= m_nY || k >= m_nZ) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Index (" << i << ", " << j << ", " << k
                << ") out of range.\n";
      continue;
    }
    if (isSet[i][j][k]) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Mesh element (" << i << ", " << j << ", " << k
                << ") has already been set.\n";
      continue;
    }
    // Get the field values.
    if (fmt == 1 || fmt == 3) {
      // Two-dimensional field-map
      fz = 0.;
      data >> fx >> fy;
    } else if (fmt == 5) {
      double temp;
      data >> fy >> fx >> temp;
      fz = temp;
    } else {
      data >> fx >> fy >> fz;
    }
    if (data.fail()) {
      std::cerr << m_className << "::LoadData:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot read field values.\n";
      bad = true;
      break;
    }
    fx *= scaleF;
    fy *= scaleF;
    fz *= scaleF;
    if (withPotential) {
      data >> v;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot read potential.\n";
        bad = true;
        break;
      }
      v *= scaleP;
      if (m_pMin > m_pMax) {
        // First value.
        m_pMin = v;
        m_pMax = v;
      } else {
        if (v < m_pMin) m_pMin = v;
        if (v > m_pMax) m_pMax = v;
      }
    }
    if (withRegion) {
      data >> region;
      if (data.fail()) {
        std::cerr << m_className << "::LoadData:\n"
                  << "    Error reading line " << nLines << ".\n"
                  << "    Cannot read region.\n";
        bad = true;
        break;
      }
    }
    if (fmt == 1 || fmt == 3) {
      // Two-dimensional field-map
      for (unsigned int kk = 0; kk < m_nZ; ++kk) {
        fields[i][j][kk].fx = fx;
        fields[i][j][kk].fy = fy;
        fields[i][j][kk].fz = fz;
        fields[i][j][kk].v = v;
        if (withRegion) m_regions[i][j][kk] = region;
        isSet[i][j][kk] = true;
      }
    } else {
      fields[i][j][k].fx = fx;
      fields[i][j][k].fy = fy;
      fields[i][j][k].fz = fz;
      fields[i][j][k].v = v;
      if (withRegion) m_regions[i][j][k] = region;
      isSet[i][j][k] = true;
    }
    ++nValues;
  }
  if (bad) return false;
  std::cout << m_className << "::LoadData:\n"
            << "    Read " << nValues << " values from " << filename << ".\n";
  unsigned int nExpected = m_nX * m_nY;
  if (fmt == 2 || fmt == 4 || fmt == 5) nExpected *= m_nZ;
  if (nExpected != nValues) {
    std::cerr << m_className << "::LoadData:\n"
              << "   Expected " << nExpected << " values.\n";
  }
  return true;
}

bool ComponentVoxel::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                                    double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  if (m_periodic[0] || m_mirrorPeriodic[0]) {
    xmin = -INFINITY;
    xmax = +INFINITY;
  } else {
    xmin = m_xMin;
    xmax = m_xMax;
  }

  if (m_periodic[1] || m_mirrorPeriodic[1]) {
    ymin = -INFINITY;
    ymax = +INFINITY;
  } else {
    ymin = m_yMin;
    ymax = m_yMax;
  }

  if (m_periodic[2] || m_mirrorPeriodic[2]) {
    zmin = -INFINITY;
    zmax = +INFINITY;
  } else {
    zmin = m_zMin;
    zmax = m_zMax;
  }
  return true;
}

bool ComponentVoxel::GetElementaryCell(
    double& xmin, double& ymin, double& zmin,
    double& xmax, double& ymax, double& zmax) {
  if (!m_ready) return false;
  xmin = m_xMin;
  xmax = m_xMax;
  ymin = m_yMin;
  ymax = m_yMax;
  zmin = m_zMin;
  zmax = m_zMax;
  return true;
}

bool ComponentVoxel::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_ready) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

bool ComponentVoxel::GetElectricFieldRange(double& exmin, double& exmax,
                                           double& eymin, double& eymax,
                                           double& ezmin, double& ezmax) {
  if (!m_ready) {
    std::cerr << m_className << "::GetElectricFieldRange:\n"
              << "    Field map is not ready for interpolation.\n";
    return false;
  }

  exmin = exmax = m_efields[0][0][0].fx;
  eymin = eymax = m_efields[0][0][0].fy;
  ezmin = ezmax = m_efields[0][0][0].fz;
  for (unsigned int i = 0; i < m_nX; ++i) {
    for (unsigned int j = 0; j < m_nY; ++j) {
      for (unsigned int k = 0; k < m_nZ; ++k) {
        const Element& element = m_efields[i][j][k];
        if (element.fx < exmin) exmin = element.fx;
        if (element.fx > exmax) exmax = element.fx;
        if (element.fy < eymin) eymin = element.fy;
        if (element.fy > eymax) eymax = element.fy;
        if (element.fz < ezmin) ezmin = element.fz;
        if (element.fz > ezmax) ezmax = element.fz;
      }
    }
  }
  return true;
}

void ComponentVoxel::PrintRegions() const {
  // Do not proceed if not properly initialised.
  if (!m_ready) {
    std::cerr << m_className << "::PrintRegions:\n"
              << "    Field map not yet initialised.\n";
    return;
  }

  if (m_media.empty()) {
    std::cerr << m_className << "::PrintRegions: No regions defined.\n";
    return;
  }

  std::cout << m_className << "::PrintRegions:\n";
  std::cout << "      Index     Medium\n";
  const unsigned int nMedia = m_media.size();
  for (unsigned int i = 0; i < nMedia; ++i) {
    const std::string name = m_media[i] ? m_media[i]->GetName() : "none";
    std::cout << "      " << i << "            " << name << "\n";
  }
}

void ComponentVoxel::SetMedium(const unsigned int i, Medium* m) {
  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
    if (m_media.empty()) return;
  }
  if (i >= m_media.size()) m_media.resize(i + 1, nullptr);
  m_media[i] = m;
}

Medium* ComponentVoxel::GetMedium(const unsigned int i) const {
  if (i >= m_media.size()) {
    std::cerr << m_className << "::GetMedium: Index out of range.\n";
    return nullptr;
  }
  return m_media[i];
}

bool ComponentVoxel::GetField(
    const double xi, const double yi, const double zi,
    const std::vector<std::vector<std::vector<Element> > >& field, double& fx,
    double& fy, double& fz, double& p, int& region) {
  if (!m_hasMesh) {
    std::cerr << m_className << "::GetField: Mesh is not set.\n";
    return false;
  }

  // Reduce the point to the basic cell (in case of periodicity) and
  // check if it is inside the mesh.
  bool xMirrored = false;
  const double x =
      Reduce(xi, m_xMin, m_xMax, m_periodic[0], m_mirrorPeriodic[0], xMirrored);
  if (x < m_xMin || x > m_xMax) return false;
  bool yMirrored = false;
  const double y =
      Reduce(yi, m_yMin, m_yMax, m_periodic[1], m_mirrorPeriodic[1], yMirrored);
  if (y < m_yMin || y > m_yMax) return false;
  bool zMirrored = false;
  const double z =
      Reduce(zi, m_zMin, m_zMax, m_periodic[2], m_mirrorPeriodic[2], zMirrored);
  if (z < m_zMin || z > m_zMax) return false;

  // Get the indices.
  const double sx = (x - m_xMin) / m_dx;
  const double sy = (y - m_yMin) / m_dy;
  const double sz = (z - m_zMin) / m_dz;
  unsigned int i = static_cast<unsigned int>(sx);
  unsigned int j = static_cast<unsigned int>(sy);
  unsigned int k = static_cast<unsigned int>(sz);
  if (i >= m_nX) i = m_nX - 1;
  if (j >= m_nY) j = m_nY - 1;
  if (k >= m_nZ) k = m_nZ - 1;
  region = m_regions[i][j][k];

  // Get the field and potential.
  if (m_interpolate) {
    // Get the "nodes" (voxel centres) surrounding the point.
    const double tx = sx - 0.5;
    const double ty = sy - 0.5;
    const double tz = sz - 0.5;
    int i0 = static_cast<int>(std::floor(tx));
    int j0 = static_cast<int>(std::floor(ty));
    int k0 = static_cast<int>(std::floor(tz));
    double vx = tx - i0;
    double vy = ty - j0;
    double vz = tz - k0;
    unsigned int i1 = i0 + 1;
    unsigned int j1 = j0 + 1;
    unsigned int k1 = k0 + 1;
    const bool perx = m_periodic[0] || m_mirrorPeriodic[0];
    const bool pery = m_periodic[1] || m_mirrorPeriodic[1];
    const bool perz = m_periodic[2] || m_mirrorPeriodic[2];
    if (i0 < 0) {
      if (perx) {
        i0 = m_nX - 1;
      } else {
        i0 = 0;
        vx = 0.;
      }
    } 
    if (j0 < 0) {
      if (pery) {
        j0 = m_nY - 1;
      } else {
        j0 = 0;
        vy = 0.;
      }
    } 
    if (k0 < 0) {
      if (perz) {
        k0 = m_nZ - 1;
      } else {
        k0 = 0;
        vz = 0.;
      }
    } 
    if (i1 >= m_nX) i1 = perx ? 0 : m_nX - 1;
    if (j1 >= m_nY) j1 = pery ? 0 : m_nY - 1;
    if (k1 >= m_nZ) k1 = perz ? 0 : m_nZ - 1;
    const Element& n000 = field[i0][j0][k0];
    const Element& n100 = field[i1][j0][k0];
    const Element& n010 = field[i0][j1][k0];
    const Element& n110 = field[i1][j1][k0];
    const Element& n001 = field[i0][j0][k1];
    const Element& n101 = field[i1][j0][k1];
    const Element& n011 = field[i0][j1][k1];
    const Element& n111 = field[i1][j1][k1];

    const double ux = 1. - vx;
    const double uy = 1. - vy;
    const double uz = 1. - vz;
    if (m_debug) {
      std::cout << m_className << "::GetField:\n    Determining field at ("
                << xi << ", " << yi << ", " << zi << ").\n"
                << "    X: " << i0 << " (" << ux << ") - " 
                             << i1 << " (" << vx << ").\n"
                << "    Y: " << j0 << " (" << uy << ") - " 
                             << j1 << " (" << vy << ").\n"
                << "    Z: " << k0 << " (" << uz << ") - " 
                             << k1 << " (" << vz << ").\n";
    } 
    fx = ((n000.fx * ux + n100.fx * vx) * uy +
          (n010.fx * ux + n110.fx * vx) * vy) *
             uz +
         ((n001.fx * ux + n101.fx * vx) * uy +
          (n011.fx * ux + n111.fx * vx) * vy) *
             vz;
    fy = ((n000.fy * ux + n100.fy * vx) * uy +
          (n010.fy * ux + n110.fy * vx) * vy) *
             uz +
         ((n001.fy * ux + n101.fy * vx) * uy +
          (n011.fy * ux + n111.fy * vx) * vy) *
             vz;
    fz = ((n000.fz * ux + n100.fz * vx) * uy +
          (n010.fz * ux + n110.fz * vx) * vy) *
             uz +
         ((n001.fz * ux + n101.fz * vx) * uy +
          (n011.fz * ux + n111.fz * vx) * vy) *
             vz;
    p = ((n000.v * ux + n100.v * vx) * uy + (n010.v * ux + n110.v * vx) * vy) *
            uz +
        ((n001.v * ux + n101.v * vx) * uy + (n011.v * ux + n111.v * vx) * vy) *
            vz;
  } else {
    const Element& element = field[i][j][k];
    fx = element.fx;
    fy = element.fy;
    fz = element.fz;
    p = element.v;
  }
  if (xMirrored) fx = -fx;
  if (yMirrored) fy = -fy;
  if (zMirrored) fz = -fz;
  return true;
}

bool ComponentVoxel::GetElement(const double xi, const double yi,
                                const double zi, unsigned int& i,
                                unsigned int& j, unsigned int& k,
                                bool& xMirrored, bool& yMirrored,
                                bool& zMirrored) const {
  if (!m_hasMesh) {
    std::cerr << m_className << "::GetElement: Mesh is not set.\n";
    return false;
  }

  // Reduce the point to the basic cell (in case of periodicity) and
  // check if it is inside the mesh.
  const double x =
      Reduce(xi, m_xMin, m_xMax, m_periodic[0], m_mirrorPeriodic[0], xMirrored);
  if (x < m_xMin || x > m_xMax) return false;
  const double y =
      Reduce(yi, m_yMin, m_yMax, m_periodic[1], m_mirrorPeriodic[1], yMirrored);
  if (y < m_yMin || y > m_yMax) return false;
  const double z =
      Reduce(zi, m_zMin, m_zMax, m_periodic[2], m_mirrorPeriodic[2], zMirrored);
  if (z < m_zMin || z > m_zMax) return false;

  // Get the indices.
  i = (unsigned int)((x - m_xMin) / m_dx);
  j = (unsigned int)((y - m_yMin) / m_dy);
  k = (unsigned int)((z - m_zMin) / m_dz);
  if (i >= m_nX) i = m_nX - 1;
  if (j >= m_nY) j = m_nY - 1;
  if (k >= m_nZ) k = m_nZ - 1;
  return true;
}

bool ComponentVoxel::GetElement(const unsigned int i, const unsigned int j,
                                const unsigned int k, double& v, double& ex,
                                double& ey, double& ez) const {
  v = ex = ey = ez = 0.;
  if (!m_ready) {
    if (!m_hasMesh) {
      std::cerr << m_className << "::GetElement: Mesh not set.\n";
      return false;
    }
    std::cerr << m_className << "::GetElement: Field map not set.\n";
    return false;
  }
  if (i >= m_nX || j >= m_nY || k >= m_nZ) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }
  const Element& element = m_efields[i][j][k];
  v = element.v;
  ex = element.fx;
  ey = element.fy;
  ez = element.fz;
  return true;
}

void ComponentVoxel::Reset() {
  m_regions.clear();
  m_efields.clear();
  m_bfields.clear();
  m_wfields.clear();

  m_wdfields.clear();
  m_wdtimes.clear();

  m_nX = m_nY = m_nZ = 0;
  m_xMin = m_yMin = m_zMin = 0.;
  m_xMax = m_yMax = m_zMax = 0.;
  m_pMin = m_pMax = 0.;
  m_media.clear();

  m_hasMesh = false;
  m_hasPotential = false;
  m_hasEfield = false;
  m_hasBfield = false;
  m_hasWfield = false;
  m_ready = false;
  
  m_wField_xOffset = 0.;
  m_wField_yOffset = 0.;
  m_wField_zOffset = 0.;
}

void ComponentVoxel::UpdatePeriodicity() {
  if (!m_ready) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Field map not available.\n";
    return;
  }

  // Check for conflicts.
  for (unsigned int i = 0; i < 3; ++i) {
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

double ComponentVoxel::Reduce(const double xin, const double xmin,
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

void ComponentVoxel::Initialise(
    std::vector<std::vector<std::vector<Element> > >& fields) {

  fields.resize(m_nX);
  for (unsigned int i = 0; i < m_nX; ++i) {
    fields[i].resize(m_nY);
    for (unsigned int j = 0; j < m_nY; ++j) {
      fields[i][j].resize(m_nZ);
      for (unsigned int k = 0; k < m_nZ; ++k) {
        fields[i][j][k].fx = 0.;
        fields[i][j][k].fy = 0.;
        fields[i][j][k].fz = 0.;
        fields[i][j][k].v = 0.;
      }
    }
  }
}

void ComponentVoxel::InitialiseRegions() {
  if (!m_hasMesh) return; 
  m_regions.resize(m_nX);
  for (unsigned int i = 0; i < m_nX; ++i) {
    m_regions[i].resize(m_nY);
    for (unsigned int j = 0; j < m_nY; ++j) {
      m_regions[i][j].assign(m_nZ, 0);
    }
  }
}
}
