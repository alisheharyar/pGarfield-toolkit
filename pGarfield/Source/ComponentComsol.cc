#include "Garfield/ComponentComsol.hh"

#include <math.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "Garfield/KDTree.hh"

namespace {

bool ends_with(std::string s, std::string t) {
  if (!s.empty() && s.back() == '\r') s.pop_back();
  return s.size() >= t.size() && s.substr(s.size() - t.size(), t.size()) == t;
}

bool isComment(const std::string &line) {
  if (line.empty()) return false;
  if (line[0] == '%') return true;
  return false;
}

int readInt(std::string s) {
  std::istringstream iss(s);
  int ret;
  iss >> ret;
  return ret;
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

}  // namespace

namespace Garfield {

ComponentComsol::ComponentComsol() : ComponentFieldMap("Comsol") {}

ComponentComsol::ComponentComsol(const std::string &mesh,
                                 const std::string &mplist,
                                 const std::string &field,
                                 const std::string &unit)
    : ComponentComsol() {
  Initialise(mesh, mplist, field, unit);
}

bool ComponentComsol::Initialise(const std::string &mesh,
                                 const std::string &mplist,
                                 const std::string &field,
                                 const std::string &unit) {
  Reset();

  std::vector<int> nodeIndices;

  // Get the conversion factor to be applied to the coordinates.
  m_unit = ScalingFactor(unit);
  if (m_unit <= 0.) {
    std::cerr << m_className << "::Initialise:\n    Unknown length unit "
              << unit << ". Will use default (m).\n";
    m_unit = 100.;
  }
  // Open the materials file.
  m_materials.clear();
  std::ifstream fmplist(mplist);
  if (!fmplist) {
    PrintCouldNotOpen("Initialise", mplist);
    return false;
  }
  unsigned int nMaterials;
  fmplist >> nMaterials;
  for (unsigned int i = 0; i < nMaterials; ++i) {
    Material material;
    material.driftmedium = false;
    material.medium = nullptr;
    material.ohm = -1;
    fmplist >> material.eps;
    m_materials.push_back(std::move(material));
  }
  if (m_materials.empty()) {
    // Add default material
    Material material;
    material.driftmedium = false;
    material.medium = nullptr;
    material.eps = material.ohm = -1;
    m_materials.push_back(std::move(material));
    nMaterials = 1;
  }
  std::map<int, int> domain2material;
  int d2msize;
  fmplist >> d2msize;
  for (int i = 0; i < d2msize; ++i) {
    int domain;
    fmplist >> domain;
    fmplist >> domain2material[domain];
  }
  fmplist.close();

  // Find lowest epsilon, check for eps = 0, set default drift medium.
  if (!SetDefaultDriftMedium()) return false;

  m_nodes.clear();
  std::ifstream fmesh(mesh);
  if (!fmesh) {
    PrintCouldNotOpen("Initialise", mesh);
    return false;
  }

  std::string line;
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Could not read number of nodes from " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (!ends_with(line, "# number of mesh points") &&
           !ends_with(line, "# number of mesh vertices"));
  const int nNodes = readInt(line);
  int nInRange = 0;
  std::cout << m_className << "::Initialise: " << nNodes << " nodes.\n";
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("# Mesh point coordinates") == std::string::npos &&
           line.find("# Mesh vertex coordinates") == std::string::npos);
  for (int i = 0; i < nNodes; ++i) {
    Node newNode;
    fmesh >> newNode.x >> newNode.y >> newNode.z;
    newNode.x *= m_unit;
    newNode.y *= m_unit;
    newNode.z *= m_unit;

    if (m_range.set) {
      m_nodesHolder.push_back(std::move(newNode));
      if (CheckInRange(newNode.x, newNode.y, newNode.z)) nInRange++;
    } else {
      m_nodes.push_back(std::move(newNode));
    }
  }

  std::vector<Element> elementsHolder;

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("4 tet2 # type name") == std::string::npos);
  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (!ends_with(line, "# number of elements"));
  const int nElements = readInt(line);
  m_elements.clear();
  std::cout << m_className << "::Initialise: " << nElements << " elements.\n";
  std::getline(fmesh, line);
  // Elements 6 & 7 are swapped due to differences in COMSOL and ANSYS
  // representation
  int perm[10] = {0, 1, 2, 3, 4, 5, 7, 6, 8, 9};
  for (int i = 0; i < nElements; ++i) {
    Element newElement;
    newElement.degenerate = false;
    for (int j = 0; j < 10; ++j) {
      fmesh >> newElement.emap[perm[j]];
    }
    elementsHolder.push_back(std::move(newElement));
  }

  do {
    if (!std::getline(fmesh, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << mesh << ".\n";
      fmesh.close();
      return false;
    }
  } while (line.find("# Geometric entity indices") == std::string::npos);
  for (auto &element : elementsHolder) {
    int domain;
    fmesh >> domain;
    element.matmap = domain2material.count(domain) ? domain2material[domain]
                                                   : nMaterials - 1;
  }
  fmesh.close();

  for (auto &takeElement : elementsHolder) {
    if (ElementInRange(takeElement)) {
      for (int j = 0; j < 10; j++) {
        nodeIndices.push_back(takeElement.emap[j]);
      }
      m_elements.push_back(std::move(takeElement));
    }
  }

  if (m_range.set) {
    std::vector<int> m_nodeMap(nNodes, -1);
    // Rearange nodeIndices and delete duplicates
    sort(nodeIndices.begin(), nodeIndices.end());
    nodeIndices.erase(std::unique(nodeIndices.begin(), nodeIndices.end()),
                      nodeIndices.end());
    // Go over nodeIndices and add the corresponding m_nodesHolder node to
    // m_nodes
    for (int &i : nodeIndices) {
      m_nodes.push_back(m_nodesHolder[i]);
      // Update m_nodeMap to get correct node idex.
      m_nodeMap[i] = m_nodes.size() - 1;
    }
    // Go over m_elements and update the node idex using the map you just
    // created
    for (Element &takeElement : m_elements) {
      for (int j = 0; j < 10; ++j) {
        takeElement.emap[j] = m_nodeMap[takeElement.emap[j]];
        if (takeElement.emap[j] == -1) {
          return false;
        }
      }
    }
  }

  std::ifstream ffield(field);
  if (!ffield) {
    PrintCouldNotOpen("Initialise", field);
    return false;
  }

  const std::string hdr1 =
      "% x                       y                        z                    "
      "    V (V)";

  const std::string hdr2 =
      "% x             y              z              V (V)";
  do {
    if (!std::getline(ffield, line)) {
      std::cerr << m_className << "::Initialise:\n"
                << "    Error parsing " << field << ".\n";
      ffield.close();
      return false;
    }
  } while (line.find(hdr1) == std::string::npos &&
           line.find(hdr2) == std::string::npos);
  std::istringstream sline(line);
  std::string token;
  sline >> token;  // %
  sline >> token;  // x
  sline >> token;  // y
  sline >> token;  // z
  sline >> token;  // V
  sline >> token;  // (V)
  m_wfields.clear();
  m_wfieldsOk.clear();
  while (sline >> token) {
    std::cout << m_className << "::Initialise:\n"
              << "    Reading data for weighting field " << token << ".\n";
    m_wfields.push_back(token);
    m_wfieldsOk.push_back(true);
    sline >> token;  // (V)
  }
  const size_t nWeightingFields = m_wfields.size();

  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nNodes)) - 1, 1.)));
  std::cout << m_className << "::Initialise: Reading potentials.\n";
  PrintProgress(0.);
  // Build a k-d tree from the node coordinates.
  std::vector<std::vector<double>> points;
  for (const auto &node : m_nodes) {
    std::vector<double> point = {node.x, node.y, node.z};
    points.push_back(std::move(point));
  }
  KDTree kdtree(points);

  const int usedSize = m_range.set ? m_nodes.size() : nNodes;
  std::vector<bool> used(usedSize, false);
  for (int i = 0; i < nNodes; ++i) {
    double x, y, z, v;
    ffield >> x >> y >> z >> v;
    x *= m_unit;
    y *= m_unit;
    z *= m_unit;
    std::vector<double> w;
    for (size_t j = 0; j < nWeightingFields; ++j) {
      double p;
      ffield >> p;
      w.push_back(p);
    }
    const std::vector<double> pt = {x, y, z};
    std::vector<KDTreeResult> res;
    kdtree.n_nearest(pt, 1, res);
    if (res.empty()) {
      std::cerr << std::endl
                << m_className << "::Initialise:\n"
                << "    Could not find a matching mesh node for point (" << x
                << ", " << y << ", " << z << ")\n.";
      ffield.close();
      return false;
    }
    if (!CheckInRange(x, y, z) && res[0].dis > maxNodeDistance) continue;
    const size_t k = res[0].idx;
    used[k] = true;
    m_nodes[k].v = v;
    m_nodes[k].w = w;
    if ((i + 1) % nPrint == 0) PrintProgress(double(i + 1) / nNodes);
  }
  PrintProgress(1.);
  ffield.close();
  auto nMissing = std::count(used.begin(), used.end(), false);
  if (m_range.set) nMissing = nMissing - m_nodes.size() + nInRange;
  if (nMissing > 0) {
    std::cerr << std::endl
              << m_className << "::Initialise:\n"
              << "    Missing potentials for " << nMissing << " nodes.\n";
    // return false;
  }

  m_ready = true;
  Prepare();
  std::cout << std::endl << m_className << "::Initialise: Done.\n";
  return true;
}

bool ComponentComsol::SetWeightingPotential(const std::string &field,
                                            const std::string &label) {
  if (!m_ready) {
    std::cerr << m_className << "::SetWeightingPotential:\n"
              << "    No valid field map is present.\n"
              << "    Weighting fields cannot be added.\n";
    return false;
  }

  double x, y, z;

  // Open the voltage list.
  std::ifstream ffield(field);
  if (!ffield) {
    PrintCouldNotOpen("SetWeightingPotential", field);
    return false;
  }

  // Check if a weighting field with the same label already exists.
  const size_t iw = GetOrCreateWeightingFieldIndex(label);

  if (iw + 1 != m_wfields.size()) {
    std::cout << m_className << "::SetWeightingPotential:\n"
              << "    Replacing existing weighting field " << label << ".\n";
  }

  m_wfieldsOk[iw] = false;

  // Build a k-d tree from the node coordinates.

  std::vector<std::vector<double>> points;
  for (const auto &node : m_nodes) {
    std::vector<double> point = {node.x, node.y, node.z};
    points.push_back(std::move(point));
  }

  KDTree kdtree(points);

  std::string line;
  int nLines = 1;
  const int nNodes = m_nodes.size();
  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nNodes)) - 1, 1.)));
  std::cout << m_className << "::SetWeightingPotential:\n"
            << "    Reading weighting potentials for " << label << ".\n";
  PrintProgress(0.);

  while (std::getline(ffield, line)) {
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip lines that are not comments.
    if (isComment(line)) continue;

    std::vector<double> pvect;

    std::istringstream data(line);
    data >> x >> y >> z;
    x *= m_unit;
    y *= m_unit;
    z *= m_unit;
    if (!CheckInRange(x, y, z)) continue;
    const std::vector<double> pt = {x, y, z};
    std::vector<KDTreeResult> res;
    kdtree.n_nearest(pt, 1, res);
    if (!CheckInRange(x, y, z) && res[0].dis > maxNodeDistance) {
      continue;
    }
    if (res.empty()) {
      std::cerr << m_className << "::SetWeightingPotential:\n"
                << "    Could not find a matching mesh node for point (" << x
                << ", " << y << ", " << z << ")\n.";
      ffield.close();
      return false;
    }

    double p = 0.;
    data >> p;
    const size_t k = res[0].idx;
    m_nodes[k].w[iw] = p;
    m_wfieldsOk[iw] = true;

    if ((nLines + 1) % nPrint == 0) {
      PrintProgress(double(nLines + 1) / nNodes);
    }
    nLines++;
  }

  PrintProgress(1.);
  std::cout << std::endl << m_className << "::SetWeightingPotential: Done.\n";
  ffield.close();
  return true;
}

bool ComponentComsol::SetDynamicWeightingPotential(const std::string &field,
                                                   const std::string &label) {
  if (!m_ready) {
    std::cerr << m_className << "::SetDelayedWeightingPotential:\n"
              << "    No valid field map is present.\n"
              << "    Weighting fields cannot be added.\n";
    return false;
  }

  if (!m_timeset && !GetTimeInterval(field)) return false;

  if (!m_timeset) {
    std::cerr << m_className << "::SetDelayedWeightingPotential:\n"
              << "    No valid times slices of potential set.\n"
              << "    Please add the time slices.\n";
    return false;
  }

  const int T = m_wdtimes.size();

  double x, y, z;

  // Open the voltage list.
  std::ifstream ffield(field);
  if (!ffield) {
    PrintCouldNotOpen("SetDelayedWeightingPotential", field);
    return false;
  }

  // Check if a weighting field with the same label already exists.
  const size_t iw = GetOrCreateWeightingFieldIndex(label);

  if (iw + 1 != m_wfields.size()) {
    std::cout << m_className << "::SetDelayedWeightingPotential:\n"
              << "    Replacing existing weighting field " << label << ".\n";
  }

  m_wfieldsOk[iw] = false;

  // Build a k-d tree from the node coordinates.

  std::vector<std::vector<double>> points;
  for (const auto &node : m_nodes) {
    std::vector<double> point = {node.x, node.y, node.z};
    points.push_back(std::move(point));
  }

  KDTree kdtree(points);

  std::string line;
  int nLines = 1;
  const int nNodes = m_nodes.size();
  const unsigned int nPrint =
      std::pow(10, static_cast<unsigned int>(
                       std::max(std::floor(std::log10(nNodes)) - 1, 1.)));
  std::cout << m_className << "::SetDelayedWeightingPotential:\n"
            << "    Reading weighting potentials for " << label << ".\n";
  PrintProgress(0.);

  while (std::getline(ffield, line)) {
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip lines that are not comments.
    if (isComment(line)) continue;

    std::vector<double> pvect;

    std::istringstream data(line);
    data >> x >> y >> z;
    x *= m_unit;
    y *= m_unit;
    z *= m_unit;
    if (!CheckInRange(x, y, z)) continue;
    const std::vector<double> pt = {x, y, z};
    std::vector<KDTreeResult> res;
    kdtree.n_nearest(pt, 1, res);
    if (!CheckInRange(x, y, z) && res[0].dis > maxNodeDistance) {
      continue;
    }
    if (res.empty()) {
      std::cerr << m_className << "::SetDelayedWeightingPotential:\n"
                << "    Could not find a matching mesh node for point (" << x
                << ", " << y << ", " << z << ")\n.";
      ffield.close();
      return false;
    }

    double p = 0.;
    double p0 = 0.;
    for (int i = 0; i < T; i++) {
      data >> p;
      if (i == 0) p0 = p;
      pvect.push_back(p - p0);
    }
    const size_t k = res[0].idx;
    m_nodes[k].dw[iw] = pvect;
    m_nodes[k].w[iw] = p0;
    m_wfieldsOk[iw] = true;

    if ((nLines + 1) % nPrint == 0) {
      PrintProgress(double(nLines + 1) / nNodes);
    }
    nLines++;
  }

  PrintProgress(1.);
  std::cout << std::endl
            << m_className << "::SetDelayedWeightingPotential: Done.\n";
  ffield.close();
  return true;
}

void ComponentComsol::SetTimeInterval(const double mint, const double maxt,
                                      const double stept) {
  std::cout << std::endl
            << m_className
            << "::SetTimeInterval: Overwriting time interval of weighting "
               "potential.\n";

  if (m_wdtimes.empty()) {
    double t = mint;
    while (t <= maxt) {
      m_wdtimes.push_back(t);
      t += stept;
    }
  }
  m_timeset = true;

  std::cout << std::endl
            << m_className
            << "::SetTimeInterval: Time of weighting potential set for t in ["
            << mint << "," << maxt << "].\n";
}

bool ComponentComsol::GetTimeInterval(const std::string &field) {
  if (!m_wdtimes.empty())
    std::cout << std::endl
              << m_className
              << "::GetTimeInterval: Overwriting time interval of weighting "
                 "potential.\n";

  std::ifstream ffield(field);
  if (!ffield) {
    PrintCouldNotOpen("GetTimeInterval", field);
    return false;
  }

  std::string strtime = "t=";

  std::string line;
  // Find first occurrence of "geeks"
  size_t found = 0;
  bool searching = true;
  while (std::getline(ffield, line)) {
    // Skip empty lines.
    if (line.empty()) continue;
    // Skip lines that are not comments.
    if (line[0] == '%' && line[2] != 'x') continue;

    while (searching) {
      found = line.find(strtime, found + 1);
      searching = false;
      if (found != std::string::npos) {
        searching = true;

        int i = 2;

        std::string holder = "";

        while (true) {
          holder += line[found + i];
          i++;

          if (found + i == line.size()) break;
          if (line[found + i] == ' ') break;
        }
        m_wdtimes.push_back(stod(holder));
      }
    }
    break;
  }

  m_timeset = true;

  std::cout << std::endl
            << m_className
            << "::GetTimeInterval: Time of weighting potential set for t in ["
            << m_wdtimes.front() << "," << m_wdtimes.back() << "].\n";

  ffield.close();

  return true;
}

}  // namespace Garfield
