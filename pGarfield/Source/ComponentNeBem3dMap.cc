#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Garfield/ComponentNeBem3dMap.hh"
#include "Garfield/Utilities.hh"

namespace Garfield {

ComponentNeBem3dMap::ComponentNeBem3dMap() : Component("NeBem3dMap") {}

void ComponentNeBem3dMap::ElectricField(const double x, const double y,
                                        const double z, double& ex, double& ey,
                                        double& ez, double& p, Medium*& m,
                                        int& status) {
  if (m_debug) std::cout << m_className << ": In ElectricField\n";

  m = nullptr;
  // Make sure the field map has been loaded.
  if (!m_ready) {
    std::cerr << m_className << "::ElectricField:\n"
              << "    Field map is not available for interpolation.\n";
    status = -10;
    return;
  }

  // Get the mesh element - values of the lowest corner, i,j,k are returned
  // For trilinear interpolation, we need to use till the i+1, j+1, k+1 corner
  unsigned int i = 0, j = 0, k = 0;
  bool xMirrored = false, yMirrored = false, zMirrored = false;
  if (!GetElement(x, y, z, i, j, k, xMirrored, yMirrored, zMirrored)) {
    status = -11;
    return;
  }
  status = 0;
  if (m_debug) {
    std::cout << "x, y, z: " << x << ", " << y << ", " << z << "\n"
              << "i, j, k: " << i << ", " << j << ", " << k << std::endl;
  }
  {  // adjustment, if needed. int adj not needed outside this block
    int adj = 0;
    if (i >= m_nX - 1) {
      i = m_nX - 1;  // CHECK! arbitrary adjustment to avoid index overflow
      adj = 1;
    }
    if (j >= m_nY - 1) {
      j = m_nY - 1;  // This should not happen
      adj = 1;
    }
    if (k >= m_nZ - 1) {
      k = m_nZ - 1;
      adj = 1;
    }
    if (adj) {
      std::cout << "x, y, z: " << x << ", " << y << ", " << z << "\n"
                << "adjusted indices:\n "
                << "i, j, k: " << i << ", " << j << ", " << k << std::endl;
    }
  }  // adjustment block

  // Get the electric field and potential using values of six corners.
  Element& element = m_efields[i][j][k];
  double ex000 = element.fx;
  double ey000 = element.fy;
  double ez000 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p000 = element.v;
  // Get the medium.
  int region = m_regions[i][j][k];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m000 = m_media[region];
  if (!m000) status = -5;

  element = m_efields[i + 1][j][k];
  double ex100 = element.fx;
  double ey100 = element.fy;
  double ez100 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p100 = element.v;
  // Get the medium.
  region = m_regions[i + 1][j][k];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m100 = m_media[region];
  if (!m100) status = -5;

  element = m_efields[i][j + 1][k];
  double ex010 = element.fx;
  double ey010 = element.fy;
  double ez010 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p010 = element.v;
  // Get the medium.
  region = m_regions[i][j + 1][k];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m010 = m_media[region];
  if (!m010) status = -5;

  element = m_efields[i][j][k + 1];
  double ex001 = element.fx;
  double ey001 = element.fy;
  double ez001 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p001 = element.v;
  // Get the medium.
  region = m_regions[i][j][k + 1];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m001 = m_media[region];
  if (!m001) status = -5;

  element = m_efields[i + 1][j + 1][k];
  double ex110 = element.fx;
  double ey110 = element.fy;
  double ez110 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p110 = element.v;
  // Get the medium.
  region = m_regions[i + 1][j + 1][k];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m110 = m_media[region];
  if (!m110) status = -5;

  element = m_efields[i + 1][j][k + 1];
  double ex101 = element.fx;
  double ey101 = element.fy;
  double ez101 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p101 = element.v;
  // Get the medium.
  region = m_regions[i + 1][j][k + 1];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m101 = m_media[region];
  if (!m101) status = -5;

  element = m_efields[i][j + 1][k + 1];
  double ex011 = element.fx;
  double ey011 = element.fy;
  double ez011 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p011 = element.v;
  // Get the medium.
  region = m_regions[i][j + 1][k + 1];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m011 = m_media[region];
  if (!m011) status = -5;

  element = m_efields[i + 1][j + 1][k + 1];
  double ex111 = element.fx;
  double ey111 = element.fy;
  double ez111 = element.fz;
  if (xMirrored) ex = -ex;
  if (yMirrored) ey = -ey;
  if (zMirrored) ez = -ez;
  double p111 = element.v;
  // Get the medium.
  region = m_regions[i + 1][j + 1][k + 1];
  if (region < 0 || region > (int)m_media.size()) {
    m = nullptr;
    status = -5;
    return;
  }
  Medium* m111 = m_media[region];
  if (!m111) status = -5;

  double delx = (m_xMax - m_xMin) / double(m_nX - 1);
  double x0 = m_xMin + double(i) * delx;
  double x1 = m_xMin + double(i + 1) * delx;
  double dely = (m_yMax - m_yMin) / double(m_nY - 1);
  double y0 = m_yMin + double(j) * dely;
  double y1 = m_yMin + double(j + 1) * dely;
  double delz = (m_zMax - m_zMin) / double(m_nZ - 1);
  double z0 = m_zMin + double(k) * delz;
  double z1 = m_zMin + double(k + 1) * delz;
  double dx0 = (x - x0);
  double dx1 = (x1 - x);
  double dy0 = (y - y0);
  double dy1 = (y1 - y);
  double dz0 = (z - z0);
  double dz1 = (z1 - z);
  double xd = (x - x0) / (x1 - x0);
  if (xd < 0.0) xd = 0.0;
  if (xd > 1.0) xd = 1.0;
  double yd = (y - y0) / (y1 - y0);
  if (yd < 0.0) yd = 0.0;
  if (yd > 1.0) yd = 1.0;
  double zd = (z - z0) / (z1 - z0);
  if (zd < 0.0) zd = 0.0;
  if (zd > 1.0) zd = 1.0;
  ex = TriLinInt(xd, yd, zd, ex000, ex100, ex010, ex001, ex110, ex101, ex011,
                 ex111);
  ey = TriLinInt(xd, yd, zd, ey000, ey100, ey010, ey001, ey110, ey101, ey011,
                 ey111);
  ez = TriLinInt(xd, yd, zd, ez000, ez100, ez010, ez001, ez110, ez101, ez011,
                 ez111);
  p = TriLinInt(xd, yd, zd, p000, p100, p010, p001, p110, p101, p011, p111);
  /*
  m = either this material or that!
  Material value is that of the nearest node
  If there are equidistant nodes, m is the lower value to give drift a chance.
  */
  double d000 = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  double d100 = sqrt(dx1 * dx1 + dy0 * dy0 + dz0 * dz0);
  double d010 = sqrt(dx0 * dx0 + dy1 * dy1 + dz0 * dz0);
  double d001 = sqrt(dx0 * dx0 + dy0 * dy0 + dz1 * dz1);
  double d110 = sqrt(dx1 * dx1 + dy1 * dy1 + dz0 * dz0);
  double d101 = sqrt(dx1 * dx1 + dy0 * dy0 + dz1 * dz1);
  double d011 = sqrt(dx0 * dx0 + dy1 * dy1 + dz1 * dz1);
  double d111 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

  // The lower value notion is not implemented yet
  // At present the algo works benefiting the last match for a condition.
  // int mindx = 0;
  double mindst = d000;
  m = m000;
  if (d100 <= mindst) {
    // mindx = 1;
    mindst = d100;
    m = m100;
  } else if (d010 <= mindst) {
    // mindx = 2;
    mindst = d010;
    m = m010;
  } else if (d001 <= mindst) {
    // mindx = 3;
    mindst = d001;
    m = m001;
  } else if (d110 <= mindst) {
    // mindx = 4;
    mindst = d110;
    m = m110;
  } else if (d101 <= mindst) {
    // mindx = 5;
    mindst = d101;
    m = m101;
  } else if (d011 <= mindst) {
    // mindx = 6;
    mindst = d011;
    m = m011;
  } else if (d111 <= mindst) {
    // mindx = 7;
    mindst = d111;
    m = m111;
  }

  if (m_debug) {
    std::cout << "x, y, z: " << x << ", " << y << ", " << z << "\n"
              << "i, j, k: " << i << ", " << j << ", " << k << "\n"
              << "000=> ex, ey, ez, p, m: " << ex000 << ", " << ey000 << ", "
              << ez000 << ", " << p000 << ", " << m000 << "\n"
              << "100=> ex, ey, ez, p, m: " << ex100 << ", " << ey100 << ", "
              << ez100 << ", " << p100 << ", " << m100 << std::endl
              << "010=> ex, ey, ez, p, m: " << ex010 << ", " << ey010 << ", "
              << ez010 << ", " << p010 << ", " << m010 << std::endl
              << "001=> ex, ey, ez, p, m: " << ex001 << ", " << ey001 << ", "
              << ez001 << ", " << p001 << ", " << m001 << std::endl
              << "110=> ex, ey, ez, p, m: " << ex110 << ", " << ey110 << ", "
              << ez110 << ", " << p110 << ", " << m110 << std::endl
              << "101=> ex, ey, ez, p, m: " << ex101 << ", " << ey101 << ", "
              << ez101 << ", " << p101 << ", " << m101 << std::endl
              << "011=> ex, ey, ez, p, m: " << ex011 << ", " << ey011 << ", "
              << ez011 << ", " << p011 << ", " << m011 << std::endl
              << "111=> ex, ey, ez, p, m: " << ex111 << ", " << ey111 << ", "
              << ez111 << ", " << p111 << ", " << m111 << std::endl
              << "delx, x, x0, x1, dx0, dx1, xd: " << delx << ", " << x << ", "
              << x0 << ", " << x1 << ", " << dx0 << ", " << dx1 << ", " << xd
              << std::endl
              << "dely, y, y0, y1, dy0, dy1, yd: " << dely << ", " << y << ", "
              << y0 << ", " << y1 << ", " << dy0 << ", " << dy1 << ", " << yd
              << std::endl
              << "delz, z, z0, z1, dz0, dz1, zd: " << delz << ", " << z << ", "
              << z0 << ", " << z1 << ", " << dz0 << ", " << dz1 << ", " << zd
              << std::endl
              << "Values after LinInt=> ex, ey, ez, p, m: " << ex << ", " << ey
              << ", " << ez << ", " << p << ", " << m << std::endl;
  }
  /*
  ex = ex000;	// manual override interpolation
  ey = ey000;
  ez = ez000;
  p = p000;
  m = m000;	// override the nearest node material assignment
  */
}

void ComponentNeBem3dMap::ElectricField(const double x, const double y,
                                        const double z, double& ex, double& ey,
                                        double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

void ComponentNeBem3dMap::WeightingField(const double x, const double y,
                                         const double z, double& wx, double& wy,
                                         double& wz,
                                         const std::string& /*label*/) {
  int status = 0;
  Medium* med = nullptr;
  double v = 0.;
  const double x1 = x - m_wField_xOffset;
  const double y1 = y - m_wField_yOffset;
  const double z1 = z - m_wField_zOffset;
  ElectricField(x1, y1, z1, wx, wy, wz, v, med, status);
}

double ComponentNeBem3dMap::WeightingPotential(const double x, const double y,
                                               const double z,
                                               const std::string& /*label*/) {
  int status = 0;
  Medium* med = nullptr;
  double v = 0.;
  const double x1 = x - m_wField_xOffset;
  const double y1 = y - m_wField_yOffset;
  const double z1 = z - m_wField_zOffset;
  double wx = 0., wy = 0., wz = 0.;
  ElectricField(x1, y1, z1, wx, wy, wz, v, med, status);
  return v;
}

void ComponentNeBem3dMap::SetWeightingFieldOffset(const double x,
                                                  const double y,
                                                  const double z) {
  m_wField_xOffset = x;
  m_wField_yOffset = y;
  m_wField_zOffset = z;
}

void ComponentNeBem3dMap::MagneticField(const double x, const double y,
                                        const double z, double& bx, double& by,
                                        double& bz, int& status) {
  if (!m_hasBfield) {
    return Component::MagneticField(x, y, z, bx, by, bz, status);
  }

  // Get the mesh element.
  unsigned int i = 0, j = 0, k = 0;
  bool xMirrored = false, yMirrored = false, zMirrored = false;
  if (!GetElement(x, y, z, i, j, k, xMirrored, yMirrored, zMirrored)) {
    status = -11;
    return;
  }
  status = 0;
  // Get the field.
  const Element& element = m_bfields[i][j][k];
  bx = element.fx;
  by = element.fy;
  bz = element.fz;
  if (xMirrored) bx = -bx;
  if (yMirrored) by = -by;
  if (zMirrored) bz = -bz;
}

Medium* ComponentNeBem3dMap::GetMedium(const double x, const double y,
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

/// Read map information from file.
bool ComponentNeBem3dMap::LoadMapInfo(
    const std::string& MapInfoFile, std::string& MapVersion, int& OptMap,
    int& OptStaggerMap, unsigned int& NbOfXCells, unsigned int& NbOfYCells,
    unsigned int& NbOfZCells, double& Xmin, double& Xmax, double& Ymin,
    double& Ymax, double& Zmin, double& Zmax, double& XStagger,
    double& YStagger, double& ZStagger, std::string& MapDataFile) {
  std::ifstream infile;
  infile.open(MapInfoFile.c_str(), std::ios::in);
  if (!infile) {
    std::cerr << m_className << "::LoadMapInfo:\n"
              << "    Could not open file " << MapInfoFile << ".\n";
    return false;
  }

  std::string line;
  unsigned int nLines = 0;

  // read version string
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> MapVersion;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve MapVersion.\n";
      return false;
    }
  }  // version string

  // read OptMap
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> OptMap;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve OptMap.\n";
      return false;
    }
  }  // OptMap

  // read OptStaggerMap
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> OptStaggerMap;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve OptStaggerMap.\n";
      return false;
    }
  }  // OptStaggerMap

  // read NbOfXCells
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> NbOfXCells;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve NbOfXCells.\n";
      return false;
    }
  }  // NbOfXCells

  // read NbOfYCells
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> NbOfYCells;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve NbOfYCells.\n";
      return false;
    }
  }  // NbOfYCells

  // read NbOfZCells
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> NbOfZCells;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve NbOfZCells.\n";
      return false;
    }
  }  // NbOfZCells

  // read Xmin, Xmax
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> Xmin >> Xmax;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve Xmin, Xmax.\n";
      return false;
    }
  }  // Xmin, Xmax

  // read Ymin, Ymax
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> Ymin >> Ymax;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve Ymin, Ymax.\n";
      return false;
    }
  }  // Ymin, Ymax

  // read Zmin, Zmax
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> Zmin >> Zmax;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve Zmin, Zmax.\n";
      return false;
    }
  }  // Zmin, Zmax

  // read XStagger
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> XStagger;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve XStagger.\n";
      return false;
    }
  }  // XStagger

  // read YStagger
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> YStagger;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve YStagger.\n";
      return false;
    }
  }  // YStagger

  // read ZStagger
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> ZStagger;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve ZStagger.\n";
      return false;
    }
  }  // ZStagger

  // read MapDataFile
  if (!infile.fail()) {
    // Read one line.
    std::getline(infile, line);
    ++nLines;

    std::istringstream data;
    data.str(line);
    data >> MapDataFile;

    if (data.fail()) {
      std::cerr << m_className << "::LoadMapInfo:\n"
                << "    Error reading line " << nLines << ".\n"
                << "    Cannot retrieve MapDataFile.\n";
      return false;
    }
  }  // MapDataFile

  return true;
}  // end of LoadMapInfo

void ComponentNeBem3dMap::SetMesh(const unsigned int nx, const unsigned int ny,
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
  m_nX = nx;  // Note that these are the number of nodes
  m_nY = ny;  // Number of elements is one less than the number of nodes
  m_nZ = nz;
  m_xMin = xmin;
  m_yMin = ymin;
  m_zMin = zmin;
  m_xMax = xmax;
  m_yMax = ymax;
  m_zMax = zmax;
  m_hasMesh = true;
}

bool ComponentNeBem3dMap::LoadElectricField(
    const std::string& filename, const std::string& format,
    const bool withPotential, const bool withRegion, const double scaleX,
    const double scaleE, const double scaleP) {
  m_ready = false;
  m_hasPotential = m_hasEfield = false;
  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadElectricField:\n"
              << "    Mesh is not set. Call SetMesh first.\n";
    return false;
  }

  // Set up the grid.
  m_efields.resize(m_nX);
  m_regions.resize(m_nX);
  for (unsigned int i = 0; i < m_nX; ++i) {
    m_efields[i].resize(m_nY);
    m_regions[i].resize(m_nY);
    for (unsigned int j = 0; j < m_nY; ++j) {
      m_efields[i][j].resize(m_nZ);
      m_regions[i][j].resize(m_nZ);
      for (unsigned int k = 0; k < m_nZ; ++k) {
        m_efields[i][j][k].fx = 0.;
        m_efields[i][j][k].fy = 0.;
        m_efields[i][j][k].fz = 0.;
        m_efields[i][j][k].v = 0.;
        m_regions[i][j][k] = 0;
      }
    }
  }

  m_pMin = m_pMax = 0.;
  if (withPotential) {
    m_pMin = 1.;
    m_pMax = -1.;
  }
  return LoadData(filename, format, withPotential, withRegion, scaleX, scaleE,
                  scaleP, 'e');
}

bool ComponentNeBem3dMap::LoadMagneticField(const std::string& filename,
                                            const std::string& format,
                                            const double scaleX,
                                            const double scaleB) {
  m_hasBfield = false;
  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadMagneticField:\n"
              << "    Mesh is not set. Call SetMesh first.\n";
    return false;
  }

  // Set up the grid.
  m_bfields.resize(m_nX);
  for (unsigned int i = 0; i < m_nX; ++i) {
    m_bfields[i].resize(m_nY);
    for (unsigned int j = 0; j < m_nY; ++j) {
      m_bfields[i][j].resize(m_nZ);
      for (unsigned int k = 0; k < m_nZ; ++k) {
        m_bfields[i][j][k].fx = 0.;
        m_bfields[i][j][k].fy = 0.;
        m_bfields[i][j][k].fz = 0.;
        m_bfields[i][j][k].v = 0.;
      }
    }
  }

  return LoadData(filename, format, false, false, scaleX, scaleB, 1., 'b');
}

bool ComponentNeBem3dMap::LoadData(const std::string& filename,
                                   std::string format, const bool withPotential,
                                   const bool withRegion, const double scaleX,
                                   const double scaleF, const double scaleP,
                                   const char field) {
  if (!m_hasMesh) {
    std::cerr << m_className << "::LoadData: Mesh has not been set.\n";
    return false;
  }

  unsigned int nValues = 0;
  // Keep track of which elements have been read.
  std::vector<std::vector<std::vector<bool> > > isSet(
      m_nX,
      std::vector<std::vector<bool> >(m_nY, std::vector<bool>(m_nZ, false)));

  std::ifstream infile;
  infile.open(filename.c_str(), std::ios::in);
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
              << "    Unkown format (" << format << ").\n";
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
    std::istringstream data;
    data.str(line);
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
        if (field == 'e') {
          m_efields[i][j][kk].fx = fx;
          m_efields[i][j][kk].fy = fy;
          m_efields[i][j][kk].fz = fz;
          m_efields[i][j][kk].v = v;
          m_regions[i][j][kk] = region;
        } else if (field == 'b') {
          m_bfields[i][j][kk].fx = fx;
          m_bfields[i][j][kk].fy = fy;
          m_bfields[i][j][kk].fz = fz;
        }
        isSet[i][j][kk] = true;
      }
    } else {
      if (field == 'e') {
        m_efields[i][j][k].fx = fx;
        m_efields[i][j][k].fy = fy;
        m_efields[i][j][k].fz = fz;
        m_efields[i][j][k].v = v;
        m_regions[i][j][k] = region;
      } else if (field == 'b') {
        m_bfields[i][j][k].fx = fx;
        m_bfields[i][j][k].fy = fy;
        m_bfields[i][j][k].fz = fz;
      }
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
  if (field == 'e') {
    m_hasEfield = true;
    m_ready = true;
    if (withPotential) m_hasPotential = true;
  } else if (field == 'b') {
    m_hasBfield = true;
  }
  return true;
}

bool ComponentNeBem3dMap::GetBoundingBox(double& xmin, double& ymin,
                                         double& zmin, double& xmax,
                                         double& ymax, double& zmax) {
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

bool ComponentNeBem3dMap::GetVoltageRange(double& vmin, double& vmax) {
  if (!m_ready) return false;
  vmin = m_pMin;
  vmax = m_pMax;
  return true;
}

bool ComponentNeBem3dMap::GetElectricFieldRange(double& exmin, double& exmax,
                                                double& eymin, double& eymax,
                                                double& ezmin, double& ezmax) {
  if (!m_ready) {
    std::cerr << m_className << "::GetElectricFieldRange:\n";
    std::cerr << "    Field map not available.\n";
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

void ComponentNeBem3dMap::PrintRegions() const {
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

void ComponentNeBem3dMap::SetMedium(const unsigned int i, Medium* m) {
  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
    if (m_media.empty()) return;
  }
  if (i >= m_media.size()) m_media.resize(i + 1, nullptr);
  m_media[i] = m;
}

Medium* ComponentNeBem3dMap::GetMedium(const unsigned int i) const {
  if (i >= m_media.size()) {
    std::cerr << m_className << "::GetMedium: Index out of range.\n";
    return nullptr;
  }
  return m_media[i];
}

bool ComponentNeBem3dMap::GetElement(const double xi, const double yi,
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
  const double dx =
      (m_xMax - m_xMin) /
      (m_nX - 1);  // m_nX is number of nodes, (m_nX-1) number of elements
  const double dy = (m_yMax - m_yMin) / (m_nY - 1);
  const double dz = (m_zMax - m_zMin) / (m_nZ - 1);
  i = (unsigned int)((x - m_xMin) / dx);
  j = (unsigned int)((y - m_yMin) / dy);
  k = (unsigned int)((z - m_zMin) / dz);
  if (i >= m_nX) i = m_nX - 1;
  if (j >= m_nY) j = m_nY - 1;
  if (k >= m_nZ) k = m_nZ - 1;
  if (m_debug) {
    std::cout << m_className << ":In GetElement\n"
              << "x, y, z: " << x << ", " << y << ", " << z << std::endl
              << "m_xMax, m_yMax, m_zMax: " << m_xMax << ", " << m_yMax << ", "
              << m_zMax << std::endl
              << "m_xMin, m_yMin, m_zMin: " << m_xMin << ", " << m_yMin << ", "
              << m_zMin << std::endl
              << "m_nX, m_nY, m_nZ: " << m_nX << ", " << m_nY << ", " << m_nZ
              << std::endl
              << "dx, dy, dz: " << dx << ", " << dy << ", " << dz << std::endl
              << "x-m_xMin, y-m_yMin, z-m_zMin: " << x - m_xMin << ", "
              << y - m_yMin << ", " << z - m_zMin << std::endl
              << "i, j, k: " << i << ", " << j << ", " << k << std::endl
              << "End GetElement" << std::endl;
  }
  return true;
}

bool ComponentNeBem3dMap::GetElement(const unsigned int i, const unsigned int j,
                                     const unsigned int k, double& v,
                                     double& ex, double& ey, double& ez) const {
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

void ComponentNeBem3dMap::Reset() {
  m_efields.clear();
  m_bfields.clear();
  m_regions.clear();
  m_nX = m_nY = m_nZ = 0;
  m_xMin = m_yMin = m_zMin = 0.;
  m_xMax = m_yMax = m_zMax = 0.;
  m_pMin = m_pMax = 0.;
  m_media.clear();

  m_hasMesh = false;
  m_hasPotential = false;
  m_hasEfield = false;
  m_hasBfield = false;
  m_ready = false;
}

void ComponentNeBem3dMap::UpdatePeriodicity() {
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

double ComponentNeBem3dMap::Reduce(const double xin, const double xmin,
                                   const double xmax, const bool simplePeriodic,
                                   const bool mirrorPeriodic,
                                   bool& mirrored) const {
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

double ComponentNeBem3dMap::TriLinInt(const double xd, const double yd,
                                      const double zd, const double c000,
                                      const double c100, const double c010,
                                      const double c001, const double c110,
                                      const double c101, const double c011,
                                      const double c111) {
  double c00 = c000 * (1.0 - xd) + c100 * xd;
  double c10 = c010 * (1.0 - xd) + c110 * xd;
  double c01 = c001 * (1.0 - xd) + c101 * xd;
  double c11 = c011 * (1.0 - xd) + c111 * xd;
  double c0 = c00 * (1.0 - yd) + c10 * yd;
  double c1 = c01 * (1.0 - yd) + c11 * yd;
  return (c0 * (1.0 - zd) + c1 * zd);
}
}
