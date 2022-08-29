#ifndef G_COMPONENT_TCAD_2D_H
#define G_COMPONENT_TCAD_2D_H

#include <memory>

#include "ComponentTcadBase.hh"
#include "QuadTree.hh"

namespace Garfield {

/// Interpolation in a two-dimensional field map created by Sentaurus Device.

class ComponentTcad2d : public ComponentTcadBase<2> {
 public:
  /// Constructor
  ComponentTcad2d();
  /// Destructor
  ~ComponentTcad2d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override {
    double v = 0.;
    ElectricField(x, y, z, ex, ey, ez, v, m, status);
  }

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, 
                      double& xmax, double& ymax, double& zmax) override;
  bool GetElementaryCell(double& xmin, double& ymin, double& zmin, 
                         double& xmax, double& ymax, double& zmax) override;
  /// Set the z-extent of the bounding box (default: unlimited).
  void SetRangeZ(const double zmin, const double zmax);

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
  /// Get the coordinates of a mesh node and the potential 
  /// and electric field at this node.
  bool GetNode(const size_t i, double& x, double& y, double& v,
               double& ex, double& ey) const;

 private:
  // Bounding box
  bool m_hasRangeZ = false;

  // Tetrahedral tree.
  std::unique_ptr<QuadTree> m_tree;

  void Reset() override {
    Cleanup();
    m_hasRangeZ = false;
    m_ready = false;
  }

  size_t FindElement(const double x, const double y,
                     std::array<double, nMaxVertices>& w) const;
  // Check whether a point is inside a given element and calculate the
  // shape functions if it is.
  bool InElement(const double x, const double y, const Element& element,
                 std::array<double, nMaxVertices>& w) const;
  bool InRectangle(const double x, const double y, const Element& element,
                   std::array<double, nMaxVertices>& w) const;
  bool InTriangle(const double x, const double y, const Element& element,
                  std::array<double, nMaxVertices>& w) const;
  bool OnLine(const double x, const double y, const Element& element,
              std::array<double, nMaxVertices>& w) const;
  bool AtPoint(const double x, const double y, const Element& element,
               std::array<double, nMaxVertices>& w) const;

  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<double>& field, double& f) override;
  bool Interpolate(const double x, const double y, const double z,
                   const std::vector<std::array<double, 2> >& field, 
                   double& fx, double& fy, double& fz) override;
  void FillTree() override;
};
}
#endif
