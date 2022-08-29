#ifndef G_COMPONENT_TCAD_3D_H
#define G_COMPONENT_TCAD_3D_H

#include <memory>

#include "ComponentTcadBase.hh"
#include "TetrahedralTree.hh"

namespace Garfield {

/// Interpolation in a three-dimensional field map created by Sentaurus Device.

class ComponentTcad3d : public ComponentTcadBase<3> {
 public:
  /// Constructor
  ComponentTcad3d();
  /// Destructor
  ~ComponentTcad3d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, 
                      double& xmax, double& ymax, double& zmax) override;
  bool GetElementaryCell(double& xmin, double& ymin, double& zmin, 
                         double& xmax, double& ymax, double& zmax) override;

  /** Retrieve the properties of an element.
    * \param i index of the element
    * \param vol volume
    * \param dmin smallest length in the element
    * \param dmax largest length in the element
    * \param type element type
    * \param nodes indices of the constituent vertices
    * \param reg region
    */
  bool GetElement(const size_t i, double& vol, double& dmin, double& dmax,
                  int& type, std::vector<size_t>& nodes, int& reg) const;
  /// Get the coordinates of a mesh node and the potential and 
  /// electric field at this node.
  bool GetNode(const size_t i, double& x, double& y, double& z, double& v,
               double& ex, double& ey, double& ez) const;

 private:

  // Tetrahedral tree.
  std::unique_ptr<TetrahedralTree> m_tree;

  void Reset() override {
    Cleanup();
    m_ready = false;
  }

  size_t FindElement(const double x, const double y, const double z,
                     std::array<double, nMaxVertices>& w) const;
  bool InElement(const double x, const double y, const double z,
                 const Element& element, 
                 std::array<double, nMaxVertices>& w) const {
    if (x < element.bbMin[0] || x > element.bbMax[0] || 
        y < element.bbMin[1] || y > element.bbMax[1] || 
        z < element.bbMin[2] || z > element.bbMax[2]) {
      return false;
    }
    bool inside = false;
    switch (element.type) {
      case 2:
        if (InTriangle(x, y, z, element, w)) inside = true;
        break;
      case 5:
        if (InTetrahedron(x, y, z, element, w)) inside = true;
        break;
      default:
        std::cerr << m_className << "::InElement:\n"
                  << "    Invalid element type (" << element.type << ").\n";
        break;
    }
    return inside;
  }
  bool InTetrahedron(const double x, const double y, const double z,
                     const Element& element, 
                     std::array<double, nMaxVertices>& w) const;
  bool InTriangle(const double x, const double y, const double z,
                  const Element& element, 
                  std::array<double, nMaxVertices>& w) const;
  
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<double>& field, double& f) override;
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<std::array<double, 3> >& field, 
                   double& fx, double& fy, double& fz) override;
  void FillTree() override;

};
}
#endif
