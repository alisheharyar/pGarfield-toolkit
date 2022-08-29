#ifndef G_COMPONENT_FIELD_MAP_H
#define G_COMPONENT_FIELD_MAP_H

#include <TMatrix.h>
#include <TVector.h>

#include <array>
#include <iostream>
#include <memory>
#include <vector>

#include "Component.hh"
#include "TMatrixD.h"
#include "TetrahedralTree.hh"

namespace Garfield {

/// Base class for components based on finite-element field maps.

class ComponentFieldMap : public Component {
 public:
  /// Default constructor.
  ComponentFieldMap() = delete;
  /// Constructor
  ComponentFieldMap(const std::string& name);
  /// Destructor
  virtual ~ComponentFieldMap();

  /// Check element aspect ratio.
  bool Check();

  /// Show x, y, z, V and angular ranges
  void PrintRange();

  /// List all currently defined materials
  void PrintMaterials();
  /// Flag a field map material as a drift medium.
  void DriftMedium(const size_t imat);
  /// Flag a field map materials as a non-drift medium.
  void NotDriftMedium(const size_t imat);
  /// Return the number of materials in the field map.
  size_t GetNumberOfMaterials() const { return m_materials.size(); }
  /// Return the relative permittivity of a field map material.
  double GetPermittivity(const size_t imat) const;
  /// Return the conductivity of a field map material.
  double GetConductivity(const size_t imat) const;
  /// Associate a field map material with a Medium object.
  void SetMedium(const size_t imat, Medium* medium);
  /// Return the Medium associated to a field map material.
  Medium* GetMedium(const size_t imat) const;
  using Component::GetMedium;
  /// Associate all field map materials with a relative permittivity
  /// of unity to a given Medium class.
  void SetGas(Medium* medium);

  /// Return the number of mesh elements.
  virtual size_t GetNumberOfElements() const { return m_elements.size(); }
  /// Return the volume and aspect ratio of a mesh element.
  bool GetElement(const size_t i, double& vol, double& dmin,
                  double& dmax) const;
  /// Return the material and node indices of a mesh element.
  virtual bool GetElement(const size_t i, size_t& mat, bool& drift,
                          std::vector<size_t>& nodes) const;
  virtual size_t GetNumberOfNodes() const { return m_nodes.size(); }
  virtual bool GetNode(const size_t i, double& x, double& y, double& z) const;
  double GetPotential(const size_t i) const;

  // Options
  void EnableCheckMapIndices(const bool on = true) {
    m_checkMultipleElement = on;
  }
  /// Option to eliminate mesh elements in conductors (default: on).
  void EnableDeleteBackgroundElements(const bool on = true) {
    m_deleteBackground = on;
  }

  /// Enable or disable the usage of the tetrahedral tree
  /// for searching the element in the mesh.
  void EnableTetrahedralTreeForElementSearch(const bool on = true) {
    m_useTetrahedralTree = on;
  }

  /// Enable or disable warnings that the calculation of the local
  /// coordinates did not achieve the requested precision.
  void EnableConvergenceWarnings(const bool on = true) {
    m_printConvergenceWarnings = on;
  }

  Medium* GetMedium(const double x, const double y, const double z) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  double DelayedWeightingPotential(double x, double y, double z, const double t,
                                   const std::string& label) override;

  bool IsInBoundingBox(const double x, const double y, const double z) const {
    return x >= m_minBoundingBox[0] && x <= m_maxBoundingBox[0] &&
           y >= m_minBoundingBox[1] && y <= m_maxBoundingBox[1] &&
           z >= m_minBoundingBox[2] && y <= m_maxBoundingBox[2];
  }

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;
  bool GetElementaryCell(double& xmin, double& ymin, double& zmin, double& xmax,
                         double& ymax, double& zmax) override;

  bool GetVoltageRange(double& vmin, double& vmax) override {
    vmin = m_mapvmin;
    vmax = m_mapvmax;
    return true;
  }

  /** Makes a weighting potential copy of a imported map which can be translated
   * and rotated. \param label name of new electrode \param labelSource name of
   * the source electrode that will be copied \param x translation in the
   * x-direction. \param y translation in the y-direction. \param z translation
   * in the z-direction. \param alpha rotation around the x-axis. \param beta
   * rotation around the y-axis. \param gamma rotation around the z-axis.
   */
  void CopyWeightingPotential(const std::string& label,
                              const std::string& labelSource, const double x,
                              const double y, const double z,
                              const double alpha, const double beta,
                              const double gamma);

  friend class ViewFEMesh;

 protected:
  bool m_is3d = true;

  enum class ElementType {
    Unknown = 0,
    Serendipity = 5,
    CurvedTetrahedron = 13
  };
  ElementType m_elementType = ElementType::CurvedTetrahedron;

  // Elements
  struct Element {
    // Nodes
    int emap[10];
    // Material
    unsigned int matmap;
    bool degenerate;
    // Bounding box of the element
    std::array<float, 3> bbMin;
    std::array<float, 3> bbMax;
  };
  std::vector<Element> m_elements;

  // Nodes
  struct Node {
    // Coordinates
    double x, y, z;
    // Potential
    double v;
    // Weighting potentials
    std::vector<double> w;
    // Delayed weighting potentials
    std::vector<std::vector<double>> dw;
  };
  std::vector<Node> m_nodes;

  // Materials
  struct Material {
    // Permittivity
    double eps;
    // Resistivity
    double ohm;
    bool driftmedium;
    // Associated medium
    Medium* medium;
  };
  std::vector<Material> m_materials;

  std::vector<std::string> m_wfields;
  std::vector<bool> m_wfieldsOk;
  std::vector<bool> m_dwfieldsOk;

  // Weighting potential copy
  struct WeightingFieldCopy {
    // Name
    std::string name;
    // Source
    size_t iSource;
    TMatrix rotMatrix;
    TVector transVector;
  };

  // Weighting potential copies
  std::vector<WeightingFieldCopy> m_wfieldCopies;

  std::vector<double> m_wdtimes;

  // Bounding box
  bool m_hasBoundingBox = false;
  std::array<double, 3> m_minBoundingBox = {{0., 0., 0.}};
  std::array<double, 3> m_maxBoundingBox = {{0., 0., 0.}};

  // Ranges and periodicities
  std::array<double, 3> m_mapmin = {{0., 0., 0.}};
  std::array<double, 3> m_mapmax = {{0., 0., 0.}};
  std::array<double, 3> m_mapamin = {{0., 0., 0.}};
  std::array<double, 3> m_mapamax = {{0., 0., 0.}};
  std::array<double, 3> m_mapna = {{0., 0., 0.}};
  std::array<double, 3> m_cells = {{0., 0., 0.}};

  double m_mapvmin = 0.;
  double m_mapvmax = 0.;

  std::array<bool, 3> m_setang;

  // Option to delete meshing in conductors
  bool m_deleteBackground = true;

  // Warnings flag
  bool m_warning = false;
  unsigned int m_nWarnings = 0;

  // Print warnings about failed convergence when refining
  // isoparametric coordinates.
  bool m_printConvergenceWarnings = true;

  // Get the scaling factor for a given length unit.
  static double ScalingFactor(std::string unit);

  // Reset the component.
  void Reset() override;

  void Prepare();

  // Calculate x, y, z, V and angular ranges.
  virtual void SetRange();

  // Periodicities
  void UpdatePeriodicity() override {
    if (!m_is3d) UpdatePeriodicity2d();
    UpdatePeriodicityCommon();
  }
  void UpdatePeriodicity2d();
  void UpdatePeriodicityCommon();

  /// Find lowest epsilon, check for eps = 0, set default drift media flags.
  bool SetDefaultDriftMedium();

  /// Interpolate the potential in a triangle.
  static double Potential3(const std::array<double, 6>& v,
                           const std::array<double, 3>& t);
  /// Interpolate the field in a triangle.
  static void Field3(const std::array<double, 6>& v,
                     const std::array<double, 3>& t, double jac[4][4],
                     const double det, double& ex, double& ey);
  /// Interpolate the potential in a curved quadrilateral.
  static double Potential5(const std::array<double, 8>& v,
                           const std::array<double, 2>& t);
  /// Interpolate the field in a curved quadrilateral.
  static void Field5(const std::array<double, 8>& v,
                     const std::array<double, 2>& t, double jac[4][4],
                     const double det, double& ex, double& ey);
  /// Interpolate the potential in a curved quadratic tetrahedron.
  static double Potential13(const std::array<double, 10>& v,
                            const std::array<double, 4>& t);
  /// Interpolate the field in a curved quadratic tetrahedron.
  static void Field13(const std::array<double, 10>& v,
                      const std::array<double, 4>& t, double jac[4][4],
                      const double det, double& ex, double& ey, double& ez);
  /// Find the element for a point in curved quadratic quadrilaterals.
  int FindElement5(const double x, const double y, const double z, double& t1,
                   double& t2, double& t3, double& t4, double jac[4][4],
                   double& det) const;
  /// Find the element for a point in curved quadratic tetrahedra.
  int FindElement13(const double x, const double y, const double z, double& t1,
                    double& t2, double& t3, double& t4, double jac[4][4],
                    double& det) const;
  /// Find the element for a point in a cube.
  int FindElementCube(const double x, const double y, const double z,
                      double& t1, double& t2, double& t3, TMatrixD*& jac,
                      std::vector<TMatrixD*>& dN) const;

  /// Move (xpos, ypos, zpos) to field map coordinates.
  void MapCoordinates(double& xpos, double& ypos, double& zpos, bool& xmirrored,
                      bool& ymirrored, bool& zmirrored, double& rcoordinate,
                      double& rotation) const;
  /// Move (ex, ey, ez) to global coordinates.
  void UnmapFields(double& ex, double& ey, double& ez, double& xpos,
                   double& ypos, double& zpos, bool& xmirrored, bool& ymirrored,
                   bool& zmirrored, double& rcoordinate,
                   double& rotation) const;

  static int ReadInteger(char* token, int def, bool& error);
  static double ReadDouble(char* token, double def, bool& error);

  virtual double GetElementVolume(const size_t i) const;
  virtual void GetAspectRatio(const size_t i, double& dmin, double& dmax) const;

  size_t GetWeightingFieldIndex(const std::string& label) const;
  size_t GetOrCreateWeightingFieldIndex(const std::string& label);

  void PrintWarning(const std::string& header);
  void PrintNotReady(const std::string& header) const;
  void PrintCouldNotOpen(const std::string& header,
                         const std::string& filename) const;
  void PrintElement(const std::string& header, const double x, const double y,
                    const double z, const double t1, const double t2,
                    const double t3, const double t4, const Element& element,
                    const unsigned int n, const int iw = -1) const;
  /// Interpolation of potential between two time slices.
  void TimeInterpolation(const double t, double& f0, double& f1, int& i0,
                         int& i1);

 private:
  /// Scan for multiple elements that contain a point
  bool m_checkMultipleElement = false;

  // Tetrahedral tree
  bool m_useTetrahedralTree = true;
  std::unique_ptr<TetrahedralTree> m_octree;

  /// Flag to check if bounding boxes of elements are cached
  bool m_cacheElemBoundingBoxes = false;

  /// Calculate local coordinates for curved quadratic triangles.
  int Coordinates3(double x, double y, double z, double& t1, double& t2,
                   double& t3, double& t4, double jac[4][4], double& det,
                   const Element& element) const;
  /// Calculate local coordinates for linear quadrilaterals.
  int Coordinates4(const double x, const double y, const double z, double& t1,
                   double& t2, double& t3, double& t4, double& det,
                   const Element& element) const;
  /// Calculate local coordinates for curved quadratic quadrilaterals.
  int Coordinates5(const double x, const double y, const double z, double& t1,
                   double& t2, double& t3, double& t4, double jac[4][4],
                   double& det, const Element& element) const;
  /// Calculate local coordinates in linear tetrahedra.
  void Coordinates12(const double x, const double y, const double z, double& t1,
                     double& t2, double& t3, double& t4,
                     const Element& element) const;
  /// Calculate local coordinates for curved quadratic tetrahedra.
  int Coordinates13(const double x, const double y, const double z, double& t1,
                    double& t2, double& t3, double& t4, double jac[4][4],
                    double& det, const Element& element) const;
  /// Calculate local coordinates for a cube.
  int CoordinatesCube(const double x, const double y, const double z,
                      double& t1, double& t2, double& t3, TMatrixD*& jac,
                      std::vector<TMatrixD*>& dN, const Element& element) const;

  /// Calculate Jacobian for curved quadratic triangles.
  void Jacobian3(const Element& element, const double u, const double v,
                 const double w, double& det, double jac[4][4]) const;
  /// Calculate Jacobian for curved quadratic quadrilaterals.
  void Jacobian5(const Element& element, const double u, const double v,
                 double& det, double jac[4][4]) const;
  /// Calculate Jacobian for curved quadratic tetrahedra.
  void Jacobian13(const Element& element, const double t, const double u,
                  const double v, const double w, double& det,
                  double jac[4][4]) const;
  /// Calculate Jacobian for a cube.
  void JacobianCube(const Element& element, const double t1, const double t2,
                    const double t3, TMatrixD*& jac,
                    std::vector<TMatrixD*>& dN) const;

  /// Calculate the bounding boxes of all elements after initialization.
  void CalculateElementBoundingBoxes();

  /// Initialize the tetrahedral tree.
  bool InitializeTetrahedralTree();

  /// Coordinate transformation matrix calculator.
  void CoordinateTransformMatrix(TMatrix& rotMatrix, TVector& transVector,
                                 const double x = 0, const double y = 0,
                                 const double z = 0, const double alpha = 0,
                                 const double beta = 0, const double gamma = 0);

  /// Get index of copy of weighting potential structure.
  size_t GetCopyWeightingPotential(const std::string& label);

  /// Transfor to coordinate system to find weighting potential value of the
  /// source field.
  void FromCopyToSourceWeightingPotential(const size_t& iwc, double& x,
                                          double& y, double& z);
};
}  // namespace Garfield

#endif
