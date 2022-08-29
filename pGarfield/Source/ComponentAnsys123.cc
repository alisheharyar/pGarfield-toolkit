#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "Garfield/ComponentAnsys123.hh"

namespace Garfield {

ComponentAnsys123::ComponentAnsys123() : ComponentFieldMap("Ansys123") {}

bool ComponentAnsys123::Initialise(std::string elist, std::string nlist,
                                   std::string mplist, std::string prnsol,
                                   std::string unit) {
  Reset();
  // Keep track of the success.
  bool ok = true;

  // Buffer for reading
  constexpr int size = 100;
  char line[size];

  // Open the material list.
  std::ifstream fmplist(mplist);
  if (!fmplist) {
    PrintCouldNotOpen("Initialise", mplist);
    return false;
  }

  // Read the material list.
  long il = 0;
  unsigned int icurrmat = 0;
  bool readerror = false;
  while (fmplist.getline(line, size, '\n')) {
    il++;
    // Skip page feed
    if (strcmp(line, "1") == 0) {
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      fmplist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        strcmp(token, "TEMPERATURE") == 0 || strcmp(token, "PROPERTY=") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13)
      continue;
    // Read number of materials and initialise the list.
    if (strcmp(token, "LIST") == 0) {
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      const unsigned int nMaterials = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
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
        std::cout << m_className << "::Initialise: " << nMaterials 
                  << " materials.\n";
      }
    } else if (strcmp(token, "MATERIAL") == 0) {
      // Version 12 format: read material number
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      const int imat = ReadInteger(token, -1, readerror);
      if (readerror || imat < 0) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      }
      icurrmat = imat;
    } else if (strcmp(token, "TEMP") == 0) {
      // Version 12 format: read property tag and value
      token = strtok(NULL, " ");
      int itype = 0;
      if (strncmp(token, "PERX", 4) == 0) {
        itype = 1;
      } else if (strncmp(token, "RSVX", 4) == 0) {
        itype = 2;
      } else {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Unknown material property flag " << token << "\n"
                  << "    in material properties file " << mplist << " (line "
                  << il << ").\n";
        ok = false;
      }
      fmplist.getline(line, size, '\n');
      il++;
      token = NULL;
      token = strtok(line, " ");
      if (icurrmat < 1 || icurrmat > m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found out-of-range current material index "
                  << icurrmat << "\n";
        std::cerr << "    in material properties file " << mplist << ".\n";
        ok = false;
        readerror = false;
      } else if (itype == 1) {
        m_materials[icurrmat - 1].eps = ReadDouble(token, -1, readerror);
      } else if (itype == 2) {
        m_materials[icurrmat - 1].ohm = ReadDouble(token, -1, readerror);
      }
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " line " << il
                  << ").\n";
        fmplist.close();
        return false;
      }
    } else if (strcmp(token, "PROPERTY") == 0) {
      // Version 11 format
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      int itype = 0;
      if (strcmp(token, "PERX") == 0) {
        itype = 1;
      } else if (strcmp(token, "RSVX") == 0) {
        itype = 2;
      } else {
        std::cerr << m_className << "::Initialise:\n"
                  << "    Unknown material property flag " << token << "\n"
                  << "    in material properties file " << mplist << " (line "
                  << il << ").\n";
        ok = false;
      }
      token = strtok(NULL, " ");
      token = strtok(NULL, " ");
      int imat = ReadInteger(token, -1, readerror);
      if (readerror) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Error reading file " << mplist << " (line " << il
                  << ").\n";
        fmplist.close();
        return false;
      } else if (imat < 1 || imat > (int)m_materials.size()) {
        std::cerr << m_className << "::Initialise:\n";
        std::cerr << "    Found out-of-range material index " << imat << "\n";
        std::cerr << "    in material properties file " << mplist << ".\n";
        ok = false;
      } else {
        fmplist.getline(line, size, '\n');
        il++;
        fmplist.getline(line, size, '\n');
        il++;
        token = NULL;
        token = strtok(line, " ");
        token = strtok(NULL, " ");
        if (itype == 1) {
          m_materials[imat - 1].eps = ReadDouble(token, -1, readerror);
        } else if (itype == 2) {
          m_materials[imat - 1].ohm = ReadDouble(token, -1, readerror);
        }
        if (readerror) {
          std::cerr << m_className << "::Initialise:\n";
          std::cerr << "    Error reading file " << mplist << " (line " << il
                    << ").\n";
          fmplist.close();
          return false;
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
            << "    Read properties of " << m_materials.size()
            << " materials from file " << mplist << ".\n";
  if (m_debug) PrintMaterials();

  // Open the element list
  std::ifstream felist(elist);
  if (!felist) {
    PrintCouldNotOpen("Initialise", elist);
    return false;
  }

  // Read the element list
  m_elements.clear();
  int nbackground = 0;
  il = 0;
  int highestnode = 0;
  while (felist.getline(line, size, '\n')) {
    il++;
    // Skip page feed
    if (strcmp(line, "1") == 0) {
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Skip page feed if Ansys > v15.x
    if (strstr(line,"***") != NULL) {
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      felist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens
    char* token = NULL;
    // Split into tokens
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "ELEM") == 0) {
      continue;
    }
    // Read the element
    int ielem = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int imat = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    int in0 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in1 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in2 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in3 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in4 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in5 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in6 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in7 = ReadInteger(token, -1, readerror);
    if (!felist.getline(line, size, '\n')) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading element " << ielem << ".\n";
      ok = false;
      break;
    }
    ++il;
    token = NULL;
    token = strtok(line, " ");
    int in8 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in9 = ReadInteger(token, -1, readerror);

    // Check synchronisation
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << elist << " (line " << il
                << ").\n";
      felist.close();
      return false;
    } else if (ielem - 1 != (int)m_elements.size() + nbackground) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Synchronisation lost on file " << elist << " (line "
                << il << ").\n";
      std::cerr << "    Element: " << ielem << " (expected " 
                << m_elements.size() << "), material: " << imat << ",\n"
                << "    nodes: (" << in0 << ", " << in1 << ", " << in2 << ", "
                << in3 << ", " << in4 << ", " << in5 << ", " << in6 << ", "
                << in7 << ", " << in8 << ", " << in9 << ")\n";
      ok = false;
    }

    // Check the material number and ensure that epsilon is non-negative
    if (imat < 1 || imat > (int)m_materials.size()) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Out-of-range material number on file " << elist
                << " (line " << il << ").\n";
      std::cerr << "    Element: " << ielem << ", material: " << imat << ",\n";
      std::cerr << "    nodes: (" << in0 << ", " << in1 << ", " << in2 << ", "
                << in3 << ", " << in4 << ", " << in5 << ", " << in6 << ", "
                << in7 << ", " << in8 << ", " << in9 << ")\n";
      ok = false;
    }
    if (m_materials[imat - 1].eps < 0) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Element " << ielem << " in element list " << elist
                << "\n";
      std::cerr << "    uses material " << imat
                << " which has not been assigned\n";
      std::cerr << "    a positive permittivity in material list " << mplist
                << ".\n";
      ok = false;
    }

    // Check the node numbers
    if (in0 < 1 || in1 < 1 || in2 < 1 || in3 < 1 || in4 < 1 || in5 < 1 ||
        in6 < 1 || in7 < 1 || in8 < 1 || in9 < 1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Found a node number < 1 on file " << elist << " (line "
                << il << ").\n";
      std::cerr << "    Element: " << ielem << ", material: " << imat << ",\n";
      std::cerr << "    nodes: (" << in0 << ", " << in1 << ", " << in2 << ", "
                << in3 << ", " << in4 << ", " << in5 << ", " << in6 << ", "
                << in7 << ", " << in8 << ", " << in9 << ")\n";
      ok = false;
    }
    if (in0 > highestnode) highestnode = in0;
    if (in1 > highestnode) highestnode = in1;
    if (in2 > highestnode) highestnode = in2;
    if (in3 > highestnode) highestnode = in3;
    if (in4 > highestnode) highestnode = in4;
    if (in5 > highestnode) highestnode = in5;
    if (in6 > highestnode) highestnode = in6;
    if (in7 > highestnode) highestnode = in7;
    if (in8 > highestnode) highestnode = in8;
    if (in9 > highestnode) highestnode = in9;

    // Skip elements in conductors.
    if (m_deleteBackground && m_materials[imat - 1].ohm == 0) {
      nbackground++;
      continue;
    }

    // These elements must not be degenerate.
    if (in0 == in1 || in0 == in2 || in0 == in3 || in0 == in4 || in0 == in5 ||
        in0 == in6 || in0 == in7 || in0 == in8 || in0 == in9 || in1 == in2 ||
        in1 == in3 || in1 == in4 || in1 == in5 || in1 == in6 || in1 == in7 ||
        in1 == in8 || in1 == in9 || in2 == in3 || in2 == in4 || in2 == in5 ||
        in2 == in6 || in2 == in7 || in2 == in8 || in2 == in9 || in3 == in4 ||
        in3 == in5 || in3 == in6 || in3 == in7 || in3 == in8 || in3 == in9 ||
        in4 == in5 || in4 == in6 || in4 == in7 || in4 == in8 || in4 == in9 ||
        in5 == in6 || in5 == in7 || in5 == in8 || in5 == in9 || in6 == in7 ||
        in6 == in8 || in6 == in9 || in7 == in8 || in7 == in9 || in8 == in9) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Element " << ielem << " of file " << elist
                << " is degenerate,\n";
      std::cerr << "    no such elements allowed in this type of map.\n";
      ok = false;
    }
    Element newElement;
    newElement.degenerate = false;

    // Store the material reference
    newElement.matmap = imat - 1;

    // Node references
    newElement.emap[0] = in0 - 1;
    newElement.emap[1] = in1 - 1;
    newElement.emap[2] = in2 - 1;
    newElement.emap[3] = in3 - 1;
    newElement.emap[4] = in4 - 1;
    newElement.emap[7] = in5 - 1;
    newElement.emap[5] = in6 - 1;
    newElement.emap[6] = in7 - 1;
    newElement.emap[8] = in8 - 1;
    newElement.emap[9] = in9 - 1;
    m_elements.push_back(std::move(newElement));
  }
  // Close the file
  felist.close();
 
  if (m_elements.empty()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Found no valid elements in file " << elist << ".\n";
    return false;
  }

  // Tell how many lines read.
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << m_elements.size() << " elements from file "
            << elist << ",\n";
  std::cout << "    highest node number: " << highestnode << ",\n";
  std::cout << "    background elements skipped: " << nbackground << "\n";
  // Check the value of the unit
  double funit = ScalingFactor(unit);
  if (funit <= 0.) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Unknown length unit " << unit << ".\n";
    ok = false;
    funit = 1.0;
  }
  if (m_debug) {
    std::cout << m_className << ":Initialise: Unit scaling factor = " 
              << funit << ".\n";
  }

  // Open the node list
  std::ifstream fnlist(nlist);
  if (!fnlist) {
    PrintCouldNotOpen("Initialise", nlist);
    return false;
  }

  // Read the node list
  m_nodes.clear();
  il = 0;
  while (fnlist.getline(line, size, '\n')) {
    il++;
    // Skip page feed
    if (strcmp(line, "1") == 0) {
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Skip page feed if ansys > v15.x
    if (strstr(line,"***") != NULL) {
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      fnlist.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "LIST") == 0 || strcmp(token, "NODE") == 0 ||
        strcmp(token, "SORT") == 0)
      continue;
    // Read the node.
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    double xnode = ReadDouble(token, -1, readerror);
    token = strtok(NULL, " ");
    double ynode = ReadDouble(token, -1, readerror);
    token = strtok(NULL, " ");
    double znode = ReadDouble(token, -1, readerror);
    // Check syntax
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << nlist << " (line " << il
                << ").\n";
      fnlist.close();
      return false;
    }
    // Check synchronisation
    if (inode - 1 != (int)m_nodes.size()) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Synchronisation lost on file " << nlist << " (line "
                << il << ").\n";
      std::cerr << "    Node: " << inode << " (expected " << m_nodes.size()
                << "), (x,y,z) = (" << xnode << ", " << ynode << ", " << znode
                << ")\n";
      ok = false;
    }
    // Store the point coordinates
    Node newNode;
    newNode.w.clear();
    newNode.x = xnode * funit;
    newNode.y = ynode * funit;
    newNode.z = znode * funit;
    m_nodes.push_back(std::move(newNode));
  }
  // Close the file
  fnlist.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n"
            << "    Read " << m_nodes.size() << " nodes from file " 
            << nlist << ".\n";
  // Check number of nodes
  if ((int)m_nodes.size() != highestnode) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Number of nodes read (" << m_nodes.size() << ") on " << nlist
              << "\n";
    std::cerr << "    does not match element list (" << highestnode << ").\n";
    ok = false;
  }

  // Open the voltage list
  std::ifstream fprnsol(prnsol);
  if (!fprnsol) {
    PrintCouldNotOpen("Initialise", prnsol);
    return false;
  }

  // Read the voltage list
  il = 0;
  unsigned int nread = 0;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Skip page feed
    if (strcmp(line, "1") == 0) {
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      continue;
    }
    // Skip page feed if ansys > v15.x
    if (strstr(line,"***") != NULL) {
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      continue;
    }
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
    // Check syntax
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << prnsol << " (line << " << il
                << ").\n";
      fprnsol.close();
      return false;
    }
    // Check node number and store if OK
    if (inode < 1 || inode > highestnode) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Node number " << inode << " out of range\n";
      std::cerr << "    on potential file " << prnsol << " (line " << il
                << ").\n";
      ok = false;
    } else {
      m_nodes[inode - 1].v = volt;
      nread++;
    }
  }
  // Close the file
  fprnsol.close();
  // Tell how many lines read
  std::cout << m_className << "::Initialise:\n    Read "
            << nread << " potentials from file " << prnsol << ".\n";
  // Check number of nodes
  if (nread != m_nodes.size()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Number of nodes read (" << nread << ") on potential file "
              << prnsol << "\n";
    std::cerr << "    does not match the node list (" << m_nodes.size() << ").\n";
    ok = false;
  }

  // Set the ready flag
  if (ok) {
    m_ready = true;
  } else {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr
        << "    Field map could not be read and can not be interpolated.\n";
    return false;
  }
  Prepare();
  return true;
}

bool ComponentAnsys123::SetWeightingField(std::string prnsol,
                                          std::string label) {
  if (!m_ready) {
    PrintNotReady("SetWeightingField");
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream fprnsol(prnsol);
  if (!fprnsol) {
    PrintCouldNotOpen("SetWeightingField", prnsol);
    return false;
  }

  // Check if a weighting field with the same label already exists.
  const size_t iw = GetOrCreateWeightingFieldIndex(label);
  if (iw + 1 != m_wfields.size()) {
    std::cout << m_className << "::SetWeightingField:\n"
              << "    Replacing existing weighting field " << label << ".\n";
  }
  m_wfieldsOk[iw] = false;

  // Buffer for reading
  constexpr int size = 100;
  char line[size];

  bool ok = true;
  // Read the voltage list.
  int il = 0;
  unsigned int nread = 0;
  bool readerror = false;
  while (fprnsol.getline(line, size, '\n')) {
    il++;
    // Skip page feed
    if (strcmp(line, "1") == 0) {
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      continue;
    }
    // Skip page feed (Ansys > v15.x).
    if (strstr(line, "***") != NULL) {
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      fprnsol.getline(line, size, '\n');
      il++;
      continue;
    }
    // Split the line in tokens.
    char* token = NULL;
    token = strtok(line, " ");
    // Skip blank lines and headers.
    if (!token || strcmp(token, " ") == 0 || strcmp(token, "\n") == 0 ||
        int(token[0]) == 10 || int(token[0]) == 13 ||
        strcmp(token, "PRINT") == 0 || strcmp(token, "*****") == 0 ||
        strcmp(token, "LOAD") == 0 || strcmp(token, "TIME=") == 0 ||
        strcmp(token, "MAXIMUM") == 0 || strcmp(token, "VALUE") == 0 ||
        strcmp(token, "NODE") == 0)
      continue;
    // Read the element.
    int inode = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    double volt = ReadDouble(token, -1, readerror);
    // Check the syntax.
    if (readerror) {
      std::cerr << m_className << "::SetWeightingField:\n";
      std::cerr << "    Error reading file " << prnsol << " (line "
                << il << ").\n";
      fprnsol.close();
      return false;
    }
    // Check node number and store if OK.
    if (inode < 1 || inode > (int)m_nodes.size()) {
      std::cerr << m_className << "::SetWeightingField:\n";
      std::cerr << "    Node number " << inode << " out of range\n";
      std::cerr << "    on potential file " << prnsol << " (line " << il
                << ").\n";
      ok = false;
    } else {
      m_nodes[inode - 1].w[iw] = volt;
      nread++;
    }
  }
  // Close the file.
  fprnsol.close();

  std::cout << m_className << "::SetWeightingField:\n"
            << "    Read " << nread << " potentials from file "
            << prnsol << ".\n";
  // Check the number of nodes.
  if (nread != m_nodes.size()) {
    std::cerr << m_className << "::SetWeightingField:\n"
              << "    Number of nodes read from potential file " << prnsol 
              << " (" << nread << ")\n    does not match the node list ("
              << m_nodes.size() << ").\n";
    ok = false;
  }

  // Set the ready flag.
  m_wfieldsOk[iw] = ok;
  if (!ok) {
    std::cerr << m_className << "::SetWeightingField:\n";
    std::cerr << "    Field map could not be read "
              << "and cannot be interpolated.\n";
    return false;
  }
  return true;
}

}
