#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "Garfield/ComponentCST.hh"

namespace {

bool ReadHeader(FILE* f, const int fileSize, const bool debug, 
                int& nX, int& nY, int& nZ, int& nNS,
                int& nES, int& nEM, int& nMaterials) {

  if (!f) return false;
  // Size of the header in binary files used in the CST export
  static constexpr int headerSize = 1000;
  if (fileSize < headerSize) {
    std::cerr << "ComponentCST::ReadHeader:\n"
              << "     Error. The file is extremely short and does not seem to "
              << "contain a header or data." << std::endl;
    return false;
  }
  char header[headerSize];
  size_t result = fread(header, sizeof(char), headerSize, f);
  if (result != headerSize) {
    std::cerr << "ComponentCST::ReadHeader: Could not read the header.\n";
    return false;
  }

  int nMeshX = 0, nMeshY = 0, nMeshZ = 0;
  int nNx = 0, nNy = 0, nNz = 0;
  int nEx = 0, nEy = 0, nEz = 0;

  std::string fmt = "mesh_nx=%d mesh_ny=%d mesh_nz=%d\n";
  fmt += "mesh_xlines=%d mesh_ylines=%d mesh_zlines=%d\n";
  fmt += "nodes_scalar=%d ";
  fmt += "nodes_vector_x=%d nodes_vector_y=%d nodes_vector_z=%d\n";
  fmt += "elements_scalar=%d elements_vector_x=%d ";
  fmt += "elements_vector_y=%d elements_vector_z=%d\n";
  fmt += "elements_material=%d\n";
  fmt += "n_materials=%d\n";
  int filled = std::sscanf(header, fmt.c_str(),
      &nMeshX, &nMeshY, &nMeshZ, &nX, &nY, &nZ, 
      &nNS, &nNx, &nNy, &nNz, &nES, &nEx, &nEy, &nEz, &nEM, &nMaterials);
  if (filled != 16) {
    std::cerr << "ComponentCST::ReadHeader: File header is broken.\n";
    return false;
  }
  if (fileSize < 1000 + (nX + nY + nZ) * 8 +
                 (nNS + nNx + nNy + nNz + nES + nEx + nEy + nEz) * 4 +
                  nEM * 1 + nMaterials * 20) {
    std::cerr << "ComponentCST::ReadHeader: Unexpected file size.\n";
    return false;
  }
  if (debug) {
    std::cout << "ComponentCST::ReadHeader:\n"
              << "  Mesh (nx): " << nMeshX << "\t Mesh (ny): " << nMeshY
              << "\t Mesh (nz): " << nMeshZ << std::endl
              << "  Mesh (x_lines): " << nX << "\t Mesh (y_lines): " << nY 
              << "\t Mesh (z_lines): " << nZ << std::endl
              << "  Nodes (scalar): " << nNS << "\t Nodes (x): " << nNx
              << "\t Nodes (y): " << nNy << "\t Nodes (z): " << nNz << "\n"
              << "  Field (scalar): " << nES << "\t Field (x): " << nEx
              << "\t Field (y): " << nEy << "\t Field (z): " << nEz << "\n"
              << "  Elements: " << nEM << "\t Materials: " << nMaterials
              << std::endl;
  }
  return true;
}

}
namespace Garfield {

ComponentCST::ComponentCST() : ComponentFieldMap("CST") {
  // Default bounding box
  m_minBoundingBox[2] = -50.;
  m_maxBoundingBox[2] = 50.;
  m_deleteBackground = false;
}

bool ComponentCST::Initialise(std::string elist, std::string nlist,
                              std::string mplist, std::string prnsol,
                              std::string unit) {
  Reset();
  // Keep track of the success
  bool ok = true;

  // Buffer for reading
  const int size = 200;
  char line[size];
  // Open the material list
  std::ifstream fmplist(mplist);
  if (!fmplist) {
    PrintCouldNotOpen("Initialise", mplist);
    return false;
  }

  // Read the material list
  int il = 0;
  bool readerror = false;
  while (fmplist.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13)
      continue;
    // Read number of materials,
    // ensure it does not exceed the maximum and initialize the list
    if (strcmp(token, "Materials") == 0) {
      token = strtok(NULL, " ");
      const int nMaterials = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Error reading file " << mplist << " (line " << il
                  << ")." << std::endl;
        fmplist.close();
        return false;
      }
      m_materials.resize(nMaterials);
      for (auto& material : m_materials) {
        material.ohm = -1;
        material.eps = -1;
        material.medium = nullptr;
      }
      if (m_debug) {
        std::cout << m_className << "::Initialise:" << std::endl;
        std::cout << "    Number of materials: " << nMaterials << "\n";
      }
    } else if (strcmp(token, "Material") == 0) {
      token = strtok(NULL, " ");
      int imat = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:" << std::endl;
        std::cerr << "     Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      } else if (imat < 1 || imat > (int)m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Found out-of-range material index " << imat << "in\n"
                  << "    material properties file " << mplist << ".\n";
        ok = false;
      } else {
        token = strtok(NULL, " ");
        int itype = 0;
        if (strcmp(token, "PERX") == 0) {
          itype = 1;
        } else if (strcmp(token, "RSVX") == 0) {
          itype = 2;
        } else {
          std::cerr << m_className << "::Initialise:\n"
                    << "    Unknown material property flag " << token << "\n"
                    << "    in material properties file " << mplist 
                    << " (line " << il << ").\n";
          ok = false;
        }
        token = strtok(NULL, " ");
        if (itype == 1) {
          m_materials[imat - 1].eps = ReadDouble(token, -1, readerror);
        } else if (itype == 2) {
          m_materials[imat - 1].ohm = ReadDouble(token, -1, readerror);
          token = strtok(NULL, " ");
          if (strcmp(token, "PERX") != 0) {
            std::cerr << m_className << "::Initialise:\n"
                      << "   Unknown material property flag " << token << "\n"
                      << "   in material file " << mplist << " (material "
                      << imat << ").\n";
            ok = false;
          } else {
            token = strtok(NULL, " ");
            m_materials[imat - 1].eps = ReadDouble(token, -1, readerror);
          }
        }
        if (readerror) {
          std::cerr << m_className << "::Initialise:\n"
                    << "     Error reading file " << mplist 
                    << " (line " << il << ")." << std::endl;
          fmplist.close();
          return false;
        }
        if (m_debug) {
          std::cout << m_className << "::Initialise:" << std::endl;
          std::cout << "    Read material properties for material "
                    << (imat - 1) << "" << std::endl;
          if (itype == 2) {
            std::cout << "    eps = " << m_materials[imat - 1].eps
                      << " ohm = " << m_materials[imat - 1].ohm << ""
                      << std::endl;
          } else {
            std::cout << "    eps = " << m_materials[imat - 1].eps << ""
                      << std::endl;
          }
        }
      }
    }
  }
  // Close the file
  fmplist.close();

  // Find lowest epsilon, check for eps = 0, set default drift medium.
  if (!SetDefaultDriftMedium()) ok = false;

  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read properties of " << m_materials.size() << " materials\n"
            << "    from file " << mplist << "." << std::endl;
  if (m_debug) PrintMaterials();

  // Check the value of the unit
  double funit = ScalingFactor(unit);
  if (funit <= 0.) {
    std::cerr << m_className << "::Initialise:\n" 
              << "    Unknown length unit " << unit << ".\n";
    ok = false;
    funit = 1.0;
  }
  if (m_debug) {
    std::cout << m_className << "::Initialise: Unit scaling factor = " 
              << funit << ".\n";
  }

  // Open the node list
  std::ifstream fnlist(nlist);
  if (!fnlist) {
    PrintCouldNotOpen("Initialise", nlist);
    return false;
  }
  // Read the node list
  m_nNodes = 0;
  il = 0;
  int xlines = 0, ylines = 0, zlines = 0;
  int lines_type = -1;
  double line_tmp;
  while (fnlist.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13)
      continue;
    // Read max sizes
    if (strcmp(token, "xmax") == 0) {
      token = strtok(NULL, " ");
      xlines = ReadInteger(token, -1, readerror);
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      ylines = ReadInteger(token, -1, readerror);
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      zlines = ReadInteger(token, -1, readerror);
      if (readerror) break;
      continue;
    }
    if (strcmp(token, "x-lines\n") == 0 || strcmp(token, "x-lines") == 0) {
      lines_type = 1;
      if (m_debug) {
        std::cout << m_className << "::Initialise:\n"
                  << "    Reading x-lines from file  " << nlist << ".\n";
      }
      continue;
    }
    if (strcmp(token, "y-lines\n") == 0 || strcmp(token, "y-lines") == 0) {
      lines_type = 2;
      if (m_debug) {
        std::cout << m_className << "::Initialise:\n"
                  << "    Reading y-lines from file  " << nlist << ".\n";
      }
      continue;
    }
    if (strcmp(token, "z-lines\n") == 0 || strcmp(token, "z-lines") == 0) {
      lines_type = 3;
      if (m_debug) {
        std::cout << m_className << "::Initialise:\n"
                  << "    Reading z-lines from file  " << nlist << ".\n";
      }
      continue;
    }
    line_tmp = ReadDouble(token, -1, readerror);
    if (lines_type == 1)
      m_xlines.push_back(line_tmp * funit);
    else if (lines_type == 2)
      m_ylines.push_back(line_tmp * funit);
    else if (lines_type == 3)
      m_zlines.push_back(line_tmp * funit);
    else {
      std::cerr << m_className << "::Initialise:" << std::endl;
      std::cerr << "    Line type was not set in  " << nlist << " (line " << il
                << ", token = " << token << ")." << std::endl;
      std::cerr << "    Maybe it is in the wrong format" << std::endl;
      std::cerr << "    e.g. missing tailing space after x-lines." << std::endl;
      ok = false;
      break;
    }
    if (readerror) break;
  }
  // Check syntax
  if (readerror) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Error reading file " << nlist 
              << " (line " << il << ").\n";
    fnlist.close();
    return false;
  }
  // Close the file
  fnlist.close();

  if ((unsigned)xlines == m_xlines.size() &&
      (unsigned)ylines == m_ylines.size() &&
      (unsigned)zlines == m_zlines.size()) {
    std::cout << m_className << "::Initialise:" << std::endl;
    std::cout << "    Found in file " << nlist << "\n    " << xlines
              << " x-lines\n    " << ylines << " y-lines\n    " << zlines
              << " z-lines" << std::endl;
  } else {
    std::cerr << m_className << "::Initialise:" << std::endl;
    std::cerr << "    There should be " << xlines << " x-lines, " << ylines
              << " y-lines and " << zlines << " z-lines in file " << nlist
              << " but I found :\n    " << m_xlines.size() << " x-lines, "
              << m_ylines.size() << " x-lines, " << m_zlines.size()
              << " z-lines." << std::endl;
  }
  m_nx = m_xlines.size();
  m_ny = m_ylines.size();
  m_nz = m_zlines.size();
  m_nNodes = m_nx * m_ny * m_nz;
  m_nElements = (m_nx - 1) * (m_ny - 1) * (m_nz - 1);

  // Tell how many lines read
  std::cout << m_className << "::Initialise:" << std::endl;
  std::cout << "    Read " << m_nNodes << " nodes from file " << nlist << "."
            << std::endl;

  // Open the element list
  std::ifstream felist(elist);
  if (!felist) {
    PrintCouldNotOpen("Initialise", elist);
    return false;
  }
  // Read the element list
  m_elementMaterial.resize(m_nElements);
  il = 0;
  while (felist.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "ELEM") == 0)
      continue;
    // Read the element
    int ielem = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    unsigned char imat = atoi(token);
    // construct node numbers
    std::vector<int> node_nb;
    try {
      // Read element material - the number of the material is stored (1, 2,
      // ...) but we need the index (0, 1, ...)
      m_elementMaterial.at(ielem) = (imat - 1);
    } catch (...) {
      std::cerr << m_className << "::Initialise:" << std::endl;
      std::cerr << "    Error reading file " << elist << " (line " << il << ")."
                << std::endl;
      std::cerr << "    The element index (" << ielem
                << ") is not in the expected range: 0 - " << m_nElements
                << std::endl;
      ok = false;
    }
    // Check the material number and ensure that epsilon is non-negative
    //    int check_mat = imat;
    if (imat < 1 || imat > m_materials.size()) {
      std::cerr << m_className << "::Initialise:" << std::endl;
      std::cerr << "   Out-of-range material number on file " << elist
                << " (line " << il << ")." << std::endl;
      std::cerr << "    Element: " << ielem << ", material: " << imat
                << std::endl;
      ok = false;
    }
    if (m_materials[imat - 1].eps < 0) {
      std::cerr << m_className << "::Initialise:" << std::endl;
      std::cerr << "    Element " << ielem << " in element list " << elist
                << " uses material " << imat << " which" << std::endl;
      std::cerr << "    has not been assigned a positive permittivity"
                << std::endl;
      std::cerr << "    in material list " << mplist << "." << std::endl;
      ok = false;
    }
  }
  // Close the file
  felist.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:" << std::endl;
  std::cout << "    Read " << m_nElements << " elements.\n";

  // Open the voltage list
  m_potential.resize(m_nNodes);
  std::ifstream fprnsol(prnsol);
  if (!fprnsol) {
    PrintCouldNotOpen("Initialise", prnsol);
    return false;
  }
  // Read the voltage list
  il = 0;
  int nread = 0;
  readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 || strcmp(token, "Max") == 0)
      continue;
    // Read node potential (in prnsol node id starts with 1 and here we will use
    // 0 as first node...)
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    double volt = ReadDouble(token, -1, readerror);

    try {
      m_potential.at(inode - 1) = volt;
      nread++;
    } catch (...) {
      std::cerr << m_className << "::Initialise:" << std::endl;
      std::cerr << "    Error reading file " << prnsol << " (line " << il
                << ")." << std::endl;
      std::cerr << "    The node index (" << inode - 1
                << ") is not in the expected range: 0 - " << m_nNodes
                << std::endl;
      ok = false;
    }
  }
  // Close the file
  fprnsol.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:" << std::endl;
  std::cout << "    Read " << nread << "/" << m_nNodes
            << " (expected) potentials from file " << prnsol << "."
            << std::endl;
  // Check number of nodes
  if (nread != (int)m_nNodes) {
    std::cerr << m_className << "::Initialise:" << std::endl;
    std::cerr << "    Number of nodes read (" << nread << ") on potential file "
              << prnsol << " does not" << std::endl;
    std::cerr << "    match the node list (" << m_nNodes << ").\n";
    ok = false;
  }
  // Set the ready flag
  if (ok) {
    m_ready = true;
  } else {
    std::cerr << m_className << "::Initialise:" << std::endl;
    std::cerr << "    Field map could not be read and cannot be interpolated."
              << std::endl;
    return false;
  }
  Prepare();
  return true;
}

bool ComponentCST::Initialise(std::string dataFile, std::string unit) {
  m_ready = false;

  // Keep track of the success
  bool ok = true;
  // Check the value of the unit
  double funit = ScalingFactor(unit);
  if (funit <= 0.) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Unknown length unit " << unit << ".\n";
    ok = false;
    funit = 1.0;
  }
  if (m_debug) {
    std::cout << m_className << "::Initialise: Unit scaling factor = " 
              << funit << ".\n";
  }
  FILE* f = fopen(dataFile.c_str(), "rb");
  if (f == nullptr) {
    PrintCouldNotOpen("Initialise", dataFile);
    return false;
  }

  struct stat fileStatus;
  stat(dataFile.c_str(), &fileStatus);
  int fileSize = fileStatus.st_size;

  int nLinesX = 0, nLinesY = 0, nLinesZ = 0;
  int nNS = 0, nES = 0, nEM = 0;
  int nMaterials = 0;
  if (!ReadHeader(f, fileSize, m_debug, nLinesX, nLinesY, nLinesZ,
                  nNS, nES, nEM, nMaterials)) {
    if (f) fclose(f);
    return false;
  } 
  m_nx = nLinesX;
  m_ny = nLinesY;
  m_nz = nLinesZ;
  m_nNodes = m_nx * m_ny * m_nz;
  m_nElements = (m_nx - 1) * (m_ny - 1) * (m_nz - 1);

  m_xlines.resize(nLinesX);
  m_ylines.resize(nLinesY);
  m_zlines.resize(nLinesZ);
  m_potential.resize(nNS);
  m_elementMaterial.resize(nEM);
  m_materials.resize(nMaterials);
  auto result = fread(m_xlines.data(), sizeof(double), m_xlines.size(), f);
  if (result != m_xlines.size()) {
    fputs("Reading error while reading xlines.", stderr);
    exit(3);
  } else if (result == 0) {
    fputs("No xlines are stored in the data file.", stderr);
    exit(3);
  }
  result = fread(m_ylines.data(), sizeof(double), m_ylines.size(), f);
  if (result != m_ylines.size()) {
    fputs("Reading error while reading ylines", stderr);
    exit(3);
  } else if (result == 0) {
    fputs("No ylines are stored in the data file.", stderr);
    exit(3);
  }
  result = fread(m_zlines.data(), sizeof(double), m_zlines.size(), f);
  if (result != m_zlines.size()) {
    fputs("Reading error while reading zlines", stderr);
    exit(3);
  } else if (result == 0) {
    fputs("No zlines are stored in the data file.", stderr);
    exit(3);
  }
  result = fread(m_potential.data(), sizeof(float), m_potential.size(), f);
  if (result != m_potential.size()) {
    fputs("Reading error while reading nodes.", stderr);
    exit(3);
  } else if (result == 0) {
    fputs("No potentials are stored in the data file.", stderr);
    exit(3);
  }
  fseek(f, nES * sizeof(float), SEEK_CUR);
  // not needed in principle - thus it is ok if nothing is read
  result = fread(m_elementMaterial.data(), sizeof(unsigned char),
                 m_elementMaterial.size(), f);
  if (result != m_elementMaterial.size()) {
    fputs("Reading error while reading element material", stderr);
    exit(3);
  }
  std::stringstream st;
  st << m_className << "::Initialise:" << std::endl;
  /*
   *  The material vector is filled according to the material id!
   *  Thus material.at(0) is material with id 0.
   */
  for (unsigned int i = 0; i < m_materials.size(); i++) {
    float id;
    result = fread(&(id), sizeof(float), 1, f);
    if (result != 1) {
      fputs("Input error while reading material id.", stderr);
      exit(3);
    }
    // const unsigned int index = id;
    const unsigned int index = i;
    unsigned int description_size = 0;
    result = fread(&(description_size), sizeof(int), 1, f);
    if (result != 1) {
      fputs("Input error while reading material description size.", stderr);
      exit(3);
    }
    char* c = new char[description_size];
    result = fread(c, sizeof(char), description_size, f);
    if (result != description_size) {
      fputs("Input error while reading material description.", stderr);
      exit(3);
    }
    std::string name = c;
    st << "  Read material: " << name;
    if (name.compare("gas") == 0) {
      st << " (considered as drift medium)";
      m_materials.at(index).driftmedium = true;
    } else {
      m_materials.at(index).driftmedium = false;
    }
    delete[] c;
    float eps;
    result = fread(&(eps), sizeof(float), 1, f);
    m_materials.at(index).eps = eps;
    if (result != 1) {
      fputs("Reading error while reading eps.", stderr);
      exit(3);
    }
    st << "; eps is: " << m_materials.at(index).eps;
    // float mue;
    // result = fread(&(mue), sizeof(float), 1, f);
    // if (result != 1) {
    //   fputs ("Reading error while reading mue.", stderr);
    //   exit (3);
    // }
    // st << "\t mue is: " << mue;
    // float rho;
    // result = fread(&(rho), sizeof(float), 1, f);
    // if (result != 1) {
    //   fputs ("Reading error while reading rho.", stderr);
    //   exit (3);
    // }
    // st << "\t rho is: " << rho;
    st << "\t id is: " << id << std::endl;
    // Skip mue and rho
    fseek(f, 2 * sizeof(float), SEEK_CUR);
    // ToDo: Check if rho should be used to decide, which material is driftable
  }
  // To be sure that they are sorted (should be already be the case)
  std::sort(m_xlines.begin(), m_xlines.end());
  std::sort(m_ylines.begin(), m_ylines.end());
  std::sort(m_zlines.begin(), m_zlines.end());
  if (funit != 1) {
    std::transform(m_xlines.begin(), m_xlines.end(), m_xlines.begin(), [funit](double x) { return x * funit;});
    std::transform(m_ylines.begin(), m_ylines.end(), m_ylines.begin(), [funit](double x) { return x * funit;});
    std::transform(m_zlines.begin(), m_zlines.end(), m_zlines.begin(), [funit](double x) { return x * funit;});
  }

  std::cout << m_className << "::Initialise" << std::endl;
  std::cout << "    x range: " << *(m_xlines.begin()) << " - "
            << *(m_xlines.end() - 1) << std::endl;
  std::cout << "    y range: " << *(m_ylines.begin()) << " - "
            << *(m_ylines.end() - 1) << std::endl;
  std::cout << "    z range: " << *(m_zlines.begin()) << " - "
            << *(m_zlines.end() - 1) << std::endl;
  fclose(f);
  // Set the ready flag
  if (ok) {
    m_ready = true;
  } else {
    std::cerr << m_className << "::Initialise:" << std::endl;
    std::cerr << "    Field map could not be read and cannot be interpolated."
              << std::endl;
    return false;
  }

  SetRange();
  UpdatePeriodicity();
  return true;
}

bool ComponentCST::SetWeightingField(std::string prnsol, std::string label,
                                     bool isBinary) {
  std::vector<float> potentials(m_nNodes);
  if (!m_ready) {
    std::cerr << m_className << "::SetWeightingField:" << std::endl;
    std::cerr << "    No valid field map is present." << std::endl;
    std::cerr << "    Weighting field cannot be added." << std::endl;
    return false;
  }

  // Open the voltage list
  std::ifstream fprnsol(prnsol);
  if (!fprnsol) {
    PrintCouldNotOpen("SetWeightingField", prnsol);
    return false;
  }
  // Check if a weighting field with the same label already exists.
  auto it = m_weightingFields.find(label);
  if (it != m_weightingFields.end()) {
    std::cout << m_className << "::SetWeightingField:" << std::endl;
    std::cout << "    Replacing existing weighting field " << label << "."
              << std::endl;
  } else {
    m_wfields.push_back(label);
    m_wfieldsOk.push_back(false);
  }

  if (std::distance(m_weightingFields.begin(), it) !=
      std::distance(m_wfields.begin(),
                    find(m_wfields.begin(), m_wfields.end(), label))) {
    std::cerr << m_className << "::SetWeightingField:" << std::endl;
    std::cerr << "    Indices of the weighting fields and the weighting field "
                 "counter are not equal!"
              << std::endl;
    return false;
  }
  unsigned int iField = std::distance(m_weightingFields.begin(), it);
  int nread = 0;
  bool ok = true;

  if (isBinary) {
    std::cout << m_className << "::SetWeightingField:" << std::endl;
    std::cout << "    Reading weighting field from binary file:"
              << prnsol << std::endl;
    FILE* f = fopen(prnsol.c_str(), "rb");
    if (f == nullptr) {
      PrintCouldNotOpen("SetWeightingField", prnsol);
      return false;
    }

    struct stat fileStatus;
    stat(prnsol.c_str(), &fileStatus);
    int fileSize = fileStatus.st_size;

    int nLinesX = 0, nLinesY = 0, nLinesZ = 0;
    int nES = 0, nEM = 0;
    int nMaterials = 0;
    if (!ReadHeader(f, fileSize, m_debug, nLinesX, nLinesY, nLinesZ,
                    nread, nES, nEM, nMaterials)) {
      if (f) fclose(f);
      return false;
    } 
    // Skip everything, but the potential
    fseek(f, nLinesX * sizeof(double), SEEK_CUR);
    fseek(f, nLinesY * sizeof(double), SEEK_CUR);
    fseek(f, nLinesZ * sizeof(double), SEEK_CUR);
    auto result = fread(potentials.data(), sizeof(float), potentials.size(), f);
    if (result != potentials.size()) {
      fputs("Reading error while reading nodes.", stderr);
      exit(3);
    } else if (result == 0) {
      fputs("No weighting potentials are stored in the data file.", stderr);
      exit(3);
    }
    fprnsol.close();
  } else {
    std::cout << m_className << "::SetWeightingField:" << std::endl;
    std::cout << "    Reading weighting field from text file:" << prnsol
              << std::endl;
    // Buffer for reading
    const int size = 100;
    char line[size];

    // Read the voltage list
    int il = 0;

    bool readerror = false;
    while (fprnsol.getline(line, size, '\n')) {
      il++;
      // Split the line in tokens
      char* token = NULL;
      token = strtok(line, " ");
      // Skip blank lines and headers
      if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
          int(token[0]) == 10 || int(token[0]) == 13 ||
          strcmp(token, "PRINT") == 0 || strcmp(token, "*****") == 0 ||
          strcmp(token, "LOAD") == 0 || strcmp(token, "TIME=") == 0 ||
          strcmp(token, "MAXIMUM") == 0 || strcmp(token, "VALUE") == 0 ||
          strcmp(token, "NODE") == 0)
        continue;
      // Read the element
      int inode = ReadInteger(token, -1, readerror);
      token = strtok(NULL, " ");
      double volt = ReadDouble(token, -1, readerror);
      try {
        potentials.at(inode - 1) = volt;
        nread++;
      } catch (...) {
        std::cerr << m_className << "::SetWeightingField:" << std::endl;
        std::cerr << "    Node number " << inode << " out of range."
                  << std::endl;
        std::cerr << "    on potential file " << prnsol << " (line " << il
                  << ")." << std::endl;
        std::cerr << "    Size of the potential vector is: "
                  << potentials.size() << std::endl;
        ok = false;
      }
    }
    // Close the file
    fprnsol.close();
  }
  // Tell how many lines read
  std::cout << m_className << "::SetWeightingField:" << std::endl;
  std::cout << "    Read " << nread << "/" << m_nNodes
            << " (expected) potentials from file " << prnsol << "."
            << std::endl;
  // Check number of nodes
  if (nread != (int)m_nNodes) {
    std::cerr << m_className << "::SetWeightingField:" << std::endl;
    std::cerr << "    Number of nodes read (" << nread << ")"
              << " on potential file (" << prnsol << ")" << std::endl;
    std::cerr << "     does not match the node list (" << m_nNodes << ")."
              << std::endl;
    ok = false;
  }
  if (!ok) {
    std::cerr << m_className << "::SetWeightingField:" << std::endl;
    std::cerr << "    Field map could not be read "
              << "and cannot be interpolated." << std::endl;
    return false;
  }

  m_weightingFields[label] = potentials;

  // Set the ready flag.
  m_wfieldsOk[iField] = ok;
  return true;
}

void ComponentCST::ShiftComponent(const double xShift, const double yShift,
                                  const double zShift) {
  std::transform(m_xlines.begin(), m_xlines.end(), m_xlines.begin(), [xShift](double x) { return x + xShift;});
  std::transform(m_ylines.begin(), m_ylines.end(), m_ylines.begin(), [yShift](double x) { return x + yShift;});
  std::transform(m_zlines.begin(), m_zlines.end(), m_zlines.begin(), [zShift](double x) { return x + zShift;});
  SetRange();
  UpdatePeriodicity();

  std::cout << m_className << "::ShiftComponent:" << std::endl;
  std::cout << "    Shifted component in x-direction: " << xShift
            << "\t y-direction: " << yShift << "\t z-direction: " << zShift
            << std::endl;
}

void ComponentCST::ElectricField(const double xin, const double yin,
                                 const double zin, double& ex, double& ey,
                                 double& ez, Medium*& m, int& status) {
  double volt;
  ElectricFieldBinary(xin, yin, zin, ex, ey, ez, volt, m, status);
}

void ComponentCST::ElectricField(const double xin, const double yin,
                                 const double zin, double& ex, double& ey,
                                 double& ez, double& volt, Medium*& m,
                                 int& status) {
  ElectricFieldBinary(xin, yin, zin, ex, ey, ez, volt, m, status, true);
}

void ComponentCST::WeightingField(const double xin, const double yin,
                                  const double zin, double& wx, double& wy,
                                  double& wz, const std::string& label) {
  // Initial values
  wx = wy = wz = 0;

  // Do not proceed if not properly initialised.
  if (!m_ready) return;

  // Look for the label.
  auto it = m_weightingFields.find(label);
  if (it == m_weightingFields.end()) {
    // Do not proceed if the requested weighting field does not exist.
    std::cerr << "No weighting field named " << label << " found!\n";
    return;
  }

  // Check if the weighting field is properly initialised.
  if (!m_wfieldsOk[std::distance(m_weightingFields.begin(), it)]) return;

  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates and get indexes
  bool mirrored[3];
  unsigned int i, j, k;
  double pos[3] = {0., 0., 0.};
  if (!Coordinate2Index(x, y, z, i, j, k, pos, mirrored)) {
    return;
  }

  double rx = (pos[0] - m_xlines.at(i)) / (m_xlines.at(i + 1) - m_xlines.at(i));
  double ry = (pos[1] - m_ylines.at(j)) / (m_ylines.at(j + 1) - m_ylines.at(j));
  double rz = (pos[2] - m_zlines.at(k)) / (m_zlines.at(k + 1) - m_zlines.at(k));

  float fwx = 0., fwy = 0., fwz = 0.;
  if (!disableFieldComponent[0])
    fwx = GetFieldComponent(i, j, k, rx, ry, rz, 'x', (*it).second);
  if (!disableFieldComponent[1])
    fwy = GetFieldComponent(i, j, k, rx, ry, rz, 'y', (*it).second);
  if (!disableFieldComponent[2])
    fwz = GetFieldComponent(i, j, k, rx, ry, rz, 'z', (*it).second);

  if (m_elementMaterial.size() > 0 && doShaping) {
    ShapeField(fwx, fwy, fwz, rx, ry, rz, i, j, k, (*it).second);
  }
  if (mirrored[0]) fwx *= -1.f;
  if (mirrored[1]) fwy *= -1.f;
  if (mirrored[2]) fwz *= -1.f;
  if (m_warning) PrintWarning("WeightingField");
  if (m_materials.at(m_elementMaterial.at(Index2Element(i, j, k))).driftmedium) {
    if (!disableFieldComponent[0]) wx = fwx;
    if (!disableFieldComponent[1]) wy = fwy;
    if (!disableFieldComponent[2]) wz = fwz;
  }
}

double ComponentCST::WeightingPotential(const double xin, const double yin,
                                        const double zin,
                                        const std::string& label) {
  // Do not proceed if not properly initialised.
  if (!m_ready) return 0.;

  // Look for the label.
  auto it = m_weightingFields.find(label);
  if (it == m_weightingFields.end()) {
    // Do not proceed if the requested weighting field does not exist.
    std::cerr << "No weighting field named " << label << " found!\n";
    return 0.;
  }

  // Check if the weighting field is properly initialised.
  if (!m_wfieldsOk[std::distance(m_weightingFields.begin(), it)]) return 0.;

  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool mirrored[3];
  unsigned int i, j, k;
  double pos[3] = {0., 0., 0.};
  if (!Coordinate2Index(x, y, z, i, j, k, pos, mirrored)) {
    return 0.;
  }
  double rx = (pos[0] - m_xlines.at(i)) /
              (m_xlines.at(i + 1) - m_xlines.at(i));
  double ry = (pos[1] - m_ylines.at(j)) /
              (m_ylines.at(j + 1) - m_ylines.at(j));
  double rz = (pos[2] - m_zlines.at(k)) /
              (m_zlines.at(k + 1) - m_zlines.at(k));

  double potential = GetPotential(i, j, k, rx, ry, rz, (*it).second);

  if (m_debug) {
    std::cout << m_className << "::WeightingPotential:" << std::endl;
    std::cout << "    Global: (" << x << "," << y << "," << z << "),"
              << std::endl;
    std::cout << "    Local: (" << rx << "," << ry << "," << rz
              << ") in element with indexes: i=" << i << ", j=" << j
              << ", k=" << k << std::endl;
    std::cout << "  Node xyzV:" << std::endl;
    std::cout << "Node 0 position: " << Index2Node(i + 1, j, k)
              << "\t potential: " << ((*it).second).at(Index2Node(i + 1, j, k))
              << "Node 1 position: " << Index2Node(i + 1, j + 1, k)
              << "\t potential: "
              << ((*it).second).at(Index2Node(i + 1, j + 1, k))
              << "Node 2 position: " << Index2Node(i, j + 1, k)
              << "\t potential: " << ((*it).second).at(Index2Node(i, j + 1, k))
              << "Node 3 position: " << Index2Node(i, j, k)
              << "\t potential: " << ((*it).second).at(Index2Node(i, j, k))
              << "Node 4 position: " << Index2Node(i + 1, j, k + 1)
              << "\t potential: "
              << ((*it).second).at(Index2Node(i + 1, j, k + 1))
              << "Node 5 position: " << Index2Node(i + 1, j + 1, k + 1)
              << "\t potential: "
              << ((*it).second).at(Index2Node(i + 1, j + 1, k + 1))
              << "Node 6 position: " << Index2Node(i, j + 1, k + 1)
              << "\t potential: "
              << ((*it).second).at(Index2Node(i, j + 1, k + 1))
              << "Node 7 position: " << Index2Node(i, j, k + 1)
              << "\t potential: " << ((*it).second).at(Index2Node(i, j, k))
              << std::endl;
  }
  return potential;
}

void ComponentCST::GetNumberOfMeshLines(unsigned int& n_x, unsigned int& n_y,
                                        unsigned int& n_z) const {
  n_x = m_xlines.size();
  n_y = m_ylines.size();
  n_z = m_zlines.size();
}

bool ComponentCST::GetElement(const size_t element, size_t& mat, bool& drift,
                              std::vector<size_t>& nodes) const {
  if (element >= m_nElements || element >= m_elementMaterial.size()) {
    std::cerr << m_className << "::GetElement: Index out of range.\n";
    return false;
  }
  mat = m_elementMaterial[element];
  drift = m_materials[mat].driftmedium;
  nodes.clear(); 
  unsigned int i0 = 0, j0 = 0, k0 = 0;
  Element2Index(element, i0, j0, k0);
  const auto i1 = i0 + 1;
  const auto j1 = j0 + 1;
  const auto k1 = k0 + 1; 
  nodes.push_back(Index2Node(i0, j0, k0));
  nodes.push_back(Index2Node(i1, j0, k0));
  nodes.push_back(Index2Node(i0, j1, k0));
  nodes.push_back(Index2Node(i1, j1, k0));
  nodes.push_back(Index2Node(i0, j0, k1));
  nodes.push_back(Index2Node(i1, j0, k1));
  nodes.push_back(Index2Node(i0, j1, k1));
  nodes.push_back(Index2Node(i1, j1, k1));
  return true;
}

bool ComponentCST::GetNode(const size_t node, 
                           double& x, double& y, double& z) const {
  if (node >= m_nNodes) {
    std::cerr << m_className << "::GetNode: Index out of range.\n";
    return false;
  }
  unsigned int i = 0, j = 0, k = 0;
  Node2Index(node, i, j, k);
  x = m_xlines[i];
  y = m_ylines[j];
  z = m_zlines[k];
  return true; 
}

void ComponentCST::GetElementBoundaries(unsigned int element, double& xmin,
                                        double& xmax, double& ymin,
                                        double& ymax, double& zmin,
                                        double& zmax) const {
  unsigned int i, j, k;
  Element2Index(element, i, j, k);
  xmin = m_xlines.at(i);
  xmax = m_xlines.at(i + 1);
  ymin = m_ylines.at(j);
  ymax = m_ylines.at(j + 1);
  zmin = m_zlines.at(k);
  zmax = m_zlines.at(k + 1);
}

Medium* ComponentCST::GetMedium(const double x, const double y,
                                const double z) {
  unsigned int i, j, k;
  Coordinate2Index(x, y, z, i, j, k);
  if (m_debug) {
    std::cout << m_className << "::GetMedium:\n"
              << "    Position (" << x << ", " << y << ", " << z << "):\n"
              << "    Indices are: x: " << i << "/" << m_xlines.size()
              << "\t y: " << j << "/" << m_ylines.size() 
              << "\t z: " << k << "/" << m_zlines.size() << std::endl;
    const auto element = Index2Element(i, j, k);
    std::cout << "    Element index: " << element << std::endl
              << "    Material index: "
              << (int)m_elementMaterial.at(element) << std::endl;
  }
  return m_materials.at(m_elementMaterial.at(Index2Element(i, j, k))).medium;
}

void ComponentCST::SetRange() {
  // Establish the ranges
  m_mapmin[0] = m_xlines.front();
  m_mapmax[0] = m_xlines.back();
  m_mapmin[1] = m_ylines.front();
  m_mapmax[1] = m_ylines.back();
  m_mapmin[2] = m_zlines.front();
  m_mapmax[2] = m_zlines.back();
  m_mapvmin = *std::min_element(m_potential.begin(), m_potential.end());
  m_mapvmax = *std::max_element(m_potential.begin(), m_potential.end());

  // Set provisional cell dimensions.
  m_minBoundingBox[0] = m_mapmin[0];
  m_maxBoundingBox[0] = m_mapmax[0];
  m_minBoundingBox[1] = m_mapmin[1];
  m_maxBoundingBox[1] = m_mapmax[1];
  if (m_is3d) {
    m_minBoundingBox[2] = m_mapmin[2];
    m_maxBoundingBox[2] = m_mapmax[2];
  } else {
    m_mapmin[2] = m_minBoundingBox[2];
    m_mapmax[2] = m_maxBoundingBox[2];
  }
  m_hasBoundingBox = true;
}

void ComponentCST::SetRangeZ(const double zmin, const double zmax) {
  if (fabs(zmax - zmin) <= 0.) {
    std::cerr << m_className << "::SetRangeZ:" << std::endl;
    std::cerr << "    Zero range is not permitted." << std::endl;
    return;
  }
  m_minBoundingBox[2] = std::min(zmin, zmax);
  m_maxBoundingBox[2] = std::max(zmin, zmax);
}

bool ComponentCST::Coordinate2Index(const double x, const double y,
                                    const double z, unsigned int& i,
                                    unsigned int& j, unsigned int& k) const {
  bool mirrored[3] = {false, false, false};
  double pos[3] = {0., 0., 0.};
  return Coordinate2Index(x, y, z, i, j, k, pos, mirrored);
}


bool ComponentCST::Coordinate2Index(const double xin, const double yin,
                                    const double zin, unsigned int& i,
                                    unsigned int& j, unsigned int& k,
                                    double* pos, bool* mirrored) const {
  // Map the coordinates onto field map coordinates
  pos[0] = xin;
  pos[1] = yin;
  pos[2] = zin;
  double rcoordinate = 0.;
  double rotation = 0.;
  MapCoordinates(pos[0], pos[1], pos[2],
                 mirrored[0], mirrored[1], mirrored[2], rcoordinate, rotation);

  auto it_x =
      std::lower_bound(m_xlines.begin(), m_xlines.end(), pos[0]);
  auto it_y =
      std::lower_bound(m_ylines.begin(), m_ylines.end(), pos[1]);
  auto it_z =
      std::lower_bound(m_zlines.begin(), m_zlines.end(), pos[2]);
  if (it_x == m_xlines.end() || it_y == m_ylines.end() ||
      it_z == m_zlines.end() || pos[0] < m_xlines.at(0) ||
      pos[1] < m_ylines.at(0) ||
      pos[2] < m_zlines.at(0)) {
    if (m_debug) {
      std::cerr << m_className << "::ElectricFieldBinary:" << std::endl;
      std::cerr << "    Could not find the given coordinate!" << std::endl;
      std::cerr << "    You ask for the following position: " << xin << ", "
                << yin << ", " << zin << std::endl;
      std::cerr << "    The mapped position is: " << pos[0] << ", "
                << pos[1] << ", " << pos[2]
                << std::endl;
    }
    return false;
  }
  /* Lower bound returns the next mesh line behind the position in question.
   * If the position in question is on a mesh line this mesh line is returned.
   * Since we are interested in the mesh line before the position in question we
   * need to move the
   * iterator to the left except for the very first mesh line!
   */
  if (it_x == m_xlines.begin())
    i = 0;
  else
    i = std::distance(m_xlines.begin(), it_x - 1);
  if (it_y == m_ylines.begin())
    j = 0;
  else
    j = std::distance(m_ylines.begin(), it_y - 1);
  if (it_z == m_zlines.begin())
    k = 0;
  else
    k = std::distance(m_zlines.begin(), it_z - 1);
  return true;
}

int ComponentCST::Index2Element(const unsigned int i, const unsigned int j,
                                const unsigned int k) const {
  if (i > m_nx - 2 || j > m_ny - 2 || k > m_nz - 2) {
    throw "ComponentCST::Index2Element: Error. Element indices out of bounds.";
  }
  return i + j * (m_nx - 1) + k * (m_nx - 1) * (m_ny - 1);
}

void ComponentCST::GetAspectRatio(const size_t element, double& dmin,
                                  double& dmax) const {
  if (element >= m_nElements) {
    dmin = dmax = 0.;
    return;
  }
  unsigned int i, j, k;
  Element2Index(element, i, j, k);
  const double dx = fabs(m_xlines.at(i + 1) - m_xlines.at(i));
  const double dy = fabs(m_ylines.at(j + 1) - m_ylines.at(j));
  const double dz = fabs(m_zlines.at(k + 1) - m_zlines.at(k));
  dmin = std::min({dx, dy, dz});
  dmax = std::max({dx, dy, dz});
}

double ComponentCST::GetElementVolume(const size_t element) const {
  if (element >= m_nElements) return 0.;
  unsigned int i, j, k;
  Element2Index(element, i, j, k);
  const double dx = fabs(m_xlines.at(i + 1) - m_xlines.at(i));
  const double dy = fabs(m_ylines.at(j + 1) - m_ylines.at(j));
  const double dz = fabs(m_zlines.at(k + 1) - m_zlines.at(k));
  return dx * dy * dz;
}

void ComponentCST::ElectricFieldBinary(const double xin, const double yin,
                                       const double zin, double& ex, double& ey,
                                       double& ez, double& volt, Medium*& m,
                                       int& status, bool calculatePotential) const {
  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  ex = ey = ez = 0;

  bool mirrored[3];
  unsigned int i, j, k;
  double pos[3] = {0., 0., 0.};
  if (!Coordinate2Index(x, y, z, i, j, k, pos, mirrored)) {
    return;
  }
  double rx = (pos[0] - m_xlines.at(i)) / (m_xlines.at(i + 1) - m_xlines.at(i));
  double ry = (pos[1] - m_ylines.at(j)) / (m_ylines.at(j + 1) - m_ylines.at(j));
  double rz = (pos[2] - m_zlines.at(k)) / (m_zlines.at(k + 1) - m_zlines.at(k));

  float fex = GetFieldComponent(i, j, k, rx, ry, rz, 'x', m_potential);
  float fey = GetFieldComponent(i, j, k, rx, ry, rz, 'y', m_potential);
  float fez = GetFieldComponent(i, j, k, rx, ry, rz, 'z', m_potential);

  if (m_elementMaterial.size() > 0 && doShaping) {
    ShapeField(fex, fey, fez, rx, ry, rz, i, j, k, m_potential);
  }
  if (mirrored[0]) fex *= -1.f;
  if (mirrored[1]) fey *= -1.f;
  if (mirrored[2]) fez *= -1.f;
  if (m_debug) {
    std::cout << m_className << "::ElectricFieldBinary:" << std::endl;
    std::cout << "    Found position (" << x << ", " << y << ", " << z
              << "): " << std::endl;
    std::cout << "    Indices are: x: " << i << "/" << m_xlines.size()
              << "\t y: " << j << "/" << m_ylines.size() << "\t z: " << k << "/"
              << m_zlines.size() << std::endl;
    if (i != 0 && j != 0 && k != 0) {
      std::cout << "    index: " << i << "\t x before: " << m_xlines.at(i - 1)
                << "\t x behind: " << m_xlines.at(i) << "\t r = " << rx
                << "\n    index: " << j << "\t y before: " << m_ylines.at(j - 1)
                << "\t y behind: " << m_ylines.at(j) << "\t r = " << ry
                << "\n    index: " << k << "\t z before: " << m_zlines.at(k - 1)
                << "\t z behind: " << m_zlines.at(k) << "\t r = " << rz
                << std::endl;
    }
    std::cout << "    Electric field is: " << fex << ", " << fey << ", " << fez
              << "): " << std::endl;
  }
  // Get the material index of the element and return the medium taken from the
  // materials (since the material id is equal to the material vector position)
  const auto imat = m_elementMaterial.at(Index2Element(i, j, k));
  m = m_materials.at(imat).medium;
  //  m = materials[elements[imap].matmap].medium;
  status = -5;
  if (m_materials.at(imat).driftmedium) {
    if (m) {
      if (m->IsDriftable()) status = 0;
    }
  }
  if (!disableFieldComponent[0]) ex = fex;
  if (!disableFieldComponent[1]) ey = fey;
  if (!disableFieldComponent[2]) ez = fez;
  if (calculatePotential)
    volt = GetPotential(i, j, k, rx, ry, rz, m_potential);
}

float ComponentCST::GetFieldComponent(
    const unsigned int i, const unsigned int j, const unsigned int k, 
    const double rx, const double ry, const double rz,
    const char component, const std::vector<float>& potentials) const {
  float e = 0.;
  if (component == 'x') {
    const float dv1 = potentials.at(Index2Node(i + 1, j, k)) -
                      potentials.at(Index2Node(i, j, k));
    const float dv2 = potentials.at(Index2Node(i + 1, j + 1, k)) -
                      potentials.at(Index2Node(i, j + 1, k));
    const float dv3 = potentials.at(Index2Node(i + 1, j + 1, k + 1)) -
                      potentials.at(Index2Node(i, j + 1, k + 1));
    const float dv4 = potentials.at(Index2Node(i + 1, j, k + 1)) -
                      potentials.at(Index2Node(i, j, k + 1));

    const float dv11 = dv1 + (dv4 - dv1) * rz;
    const float dv21 = dv2 + (dv3 - dv2) * rz;
    const float dv = dv11 + (dv21 - dv11) * ry;
    e = -1 * dv / (m_xlines.at(i + 1) - m_xlines.at(i));
  }
  if (component == 'y') {
    const float dv1 = potentials.at(Index2Node(i, j + 1, k)) -
                      potentials.at(Index2Node(i, j, k));
    const float dv2 = potentials.at(Index2Node(i, j + 1, k + 1)) -
                      potentials.at(Index2Node(i, j, k + 1));
    const float dv3 = potentials.at(Index2Node(i + 1, j + 1, k + 1)) -
                      potentials.at(Index2Node(i + 1, j, k + 1));
    const float dv4 = potentials.at(Index2Node(i + 1, j + 1, k)) -
                      potentials.at(Index2Node(i + 1, j, k));

    const float dv11 = dv1 + (dv4 - dv1) * rx;
    const float dv21 = dv2 + (dv3 - dv2) * rx;
    const float dv = dv11 + (dv21 - dv11) * rz;
    e = -1 * dv / (m_ylines.at(j + 1) - m_ylines.at(j));
  }
  if (component == 'z') {
    const float dv1 = potentials.at(Index2Node(i, j, k + 1)) -
                      potentials.at(Index2Node(i, j, k));
    const float dv2 = potentials.at(Index2Node(i + 1, j, k + 1)) -
                      potentials.at(Index2Node(i + 1, j, k));
    const float dv3 = potentials.at(Index2Node(i + 1, j + 1, k + 1)) -
                      potentials.at(Index2Node(i + 1, j + 1, k));
    const float dv4 = potentials.at(Index2Node(i, j + 1, k + 1)) -
                      potentials.at(Index2Node(i, j + 1, k));

    const float dv11 = dv1 + (dv4 - dv1) * ry;
    const float dv21 = dv2 + (dv3 - dv2) * ry;
    const float dv = dv11 + (dv21 - dv11) * rx;
    e = -1 * dv / (m_zlines.at(k + 1) - m_zlines.at(k));
  }
  return e;
}

float ComponentCST::GetPotential(
    const unsigned int i, const unsigned int j, const unsigned int k, 
    const double rx, const double ry, const double rz,
    const std::vector<float>& potentials) const {
  double t1 = rx * 2. - 1;
  double t2 = ry * 2. - 1;
  double t3 = rz * 2. - 1;
  return (potentials.at(Index2Node(i + 1, j, k)) * (1 - t1) * (1 - t2) *
              (1 - t3) +
          potentials.at(Index2Node(i + 1, j + 1, k)) * (1 + t1) * (1 - t2) *
              (1 - t3) +
          potentials.at(Index2Node(i, j + 1, k)) * (1 + t1) * (1 + t2) *
              (1 - t3) +
          potentials.at(Index2Node(i, j, k)) * (1 - t1) * (1 + t2) * (1 - t3) +
          potentials.at(Index2Node(i + 1, j, k + 1)) * (1 - t1) * (1 - t2) *
              (1 + t3) +
          potentials.at(Index2Node(i + 1, j + 1, k + 1)) * (1 + t1) *
              (1 - t2) * (1 + t3) +
          potentials.at(Index2Node(i, j + 1, k + 1)) * (1 + t1) * (1 + t2) *
              (1 + t3) +
          potentials.at(Index2Node(i, j, k + 1)) * (1 - t1) * (1 + t2) *
              (1 + t3)) /
         8.;
}

void ComponentCST::ShapeField(float& ex, float& ey, float& ez, 
    const double rx, const double ry, const double rz,
    const unsigned int i, const unsigned int j, const unsigned int k,
    const std::vector<float>& potentials) const {

  const auto m1 = m_elementMaterial.at(Index2Element(i, j, k));
  const auto imax = m_xlines.size() - 2;
  if ((i == 0 && rx >= 0.5) || (i == imax && rx < 0.5) || (i > 0 && i < imax)) {
    if (rx >= 0.5) {
      const auto m2 = m_elementMaterial.at(Index2Element(i + 1, j, k));
      if (m1 == m2) {
        float ex_next =
            GetFieldComponent(i + 1, j, k, 0.5, ry, rz, 'x', potentials);
        ex = ex +
             (rx - 0.5) * (ex_next - ex) *
                 (m_xlines.at(i + 1) - m_xlines.at(i)) /
                 (m_xlines.at(i + 2) - m_xlines.at(i + 1));
      }
    } else {
      const auto m2 = m_elementMaterial.at(Index2Element(i - 1, j, k));
      if (m1 == m2) {
        float ex_before =
            GetFieldComponent(i - 1, j, k, 0.5, ry, rz, 'x', potentials);
        ex = ex_before +
             (rx + 0.5) * (ex - ex_before) *
                 (m_xlines.at(i) - m_xlines.at(i - 1)) /
                 (m_xlines.at(i + 1) - m_xlines.at(i));
      }
    }
  }

  const auto jmax = m_ylines.size() - 2;
  if ((j == 0 && ry >= 0.5) || (j == jmax && ry < 0.5) || (j > 0 && j < jmax)) {
    if (ry >= 0.5) {
      const auto m2 = m_elementMaterial.at(Index2Element(i, j + 1, k));
      if (m1 == m2) {
        float ey_next =
            GetFieldComponent(i, j + 1, k, rx, 0.5, rz, 'y', potentials);
        ey = ey +
             (ry - 0.5) * (ey_next - ey) *
                 (m_ylines.at(j + 1) - m_ylines.at(j)) /
                 (m_ylines.at(j + 2) - m_ylines.at(j + 1));
      }
    } else {
      const auto m2 = m_elementMaterial.at(Index2Element(i, j - 1, k));
      if (m1 == m2) {
        float ey_next =
            GetFieldComponent(i, j - 1, k, rx, 0.5, rz, 'y', potentials);
        ey = ey_next +
             (ry + 0.5) * (ey - ey_next) *
                 (m_ylines.at(j) - m_ylines.at(j - 1)) /
                 (m_ylines.at(j + 1) - m_ylines.at(j));
      }
    }
  }
  const auto kmax = m_zlines.size() - 2;
  if ((k == 0 && rz >= 0.5) || (k == kmax && rz < 0.5) || (k > 0 && k < kmax)) {
    if (rz >= 0.5) {
      const auto m2 = m_elementMaterial.at(Index2Element(i, j, k + 1));
      if (m1 == m2) {
        float ez_next =
            GetFieldComponent(i, j, k + 1, rx, ry, 0.5, 'z', potentials);
        ez = ez +
             (rz - 0.5) * (ez_next - ez) *
                 (m_zlines.at(k + 1) - m_zlines.at(k)) /
                 (m_zlines.at(k + 2) - m_zlines.at(k + 1));
      }
    } else {
      const auto m2 = m_elementMaterial.at(Index2Element(i, j, k - 1));
      if (m1 == m2) {
        float ez_next =
            GetFieldComponent(i, j, k - 1, rx, ry, 0.5, 'z', potentials);
        ez = ez_next +
             (rz + 0.5) * (ez - ez_next) *
                 (m_zlines.at(k) - m_zlines.at(k - 1)) /
                 (m_zlines.at(k + 1) - m_zlines.at(k));
      }
    }
  }
}

void ComponentCST::Element2Index(const size_t element, 
    unsigned int& i, unsigned int& j, unsigned int& k) const {
  const auto nx = m_xlines.size() - 1;
  const auto ny = m_ylines.size() - 1;
  const auto nxy = nx * ny;
  k = element / nxy;
  const auto tmp = element - k * nxy;
  j = tmp / nx;
  i = tmp - j * nx;
}

int ComponentCST::Index2Node(const unsigned int i, const unsigned int j,
                             const unsigned int k) const {
  if (i > m_nx - 1 || j > m_ny - 1 || k > m_nz - 1) {
    throw "ComponentCST::Index2Node: Error. Node indices out of bounds.";
  }
  return i + j * m_nx + k * m_nx * m_ny;
}

void ComponentCST::Node2Index(const size_t node, 
    unsigned int& i, unsigned int& j, unsigned int& k) const {

  const auto nx = m_xlines.size();
  const auto ny = m_ylines.size(); 
  const auto nxy = nx * ny;
  k = node / nxy;
  const auto tmp = node - k * nxy;
  j = tmp / nx;
  i = tmp - j * nx;
} 

}
