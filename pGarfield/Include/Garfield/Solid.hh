#ifndef G_SOLID_H
#define G_SOLID_H

#include <vector>
#include <array>

namespace Garfield {

/// Surface panel.

struct Panel {
  /// Perpendicular vector
  double a, b, c;
  /// X-coordinates of vertices
  std::vector<double> xv;
  /// Y-coordinates of vertices
  std::vector<double> yv;
  /// Z-coordinates of vertices
  std::vector<double> zv;
  /// Colour index
  int colour;
  /// Reference to solid to which the panel belongs
  int volume;
};

/// Abstract base class for solids.

class Solid {
 public:
  /// Default constructor.
  Solid() = delete;
  /// Constructor.
  Solid(const double cx, const double cy, const double cz,
        const std::string& name)
      : m_cX(cx), m_cY(cy), m_cZ(cz), m_className(name) {
    m_id = s_id++;
  }

  /// Destructor
  virtual ~Solid() {}

  /// Check whether a given point is inside the solid. If requested, 
  /// use the tesselated approximation of the solid (if applicable).
  virtual bool IsInside(const double x, const double y, const double z,
                        const bool tesselated = false) const = 0;
  /// Return the bounding box of the solid.
  virtual bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                              double& xmax, double& ymax,
                              double& zmax) const = 0;
  /// Return true if the solid is a box.
  virtual bool IsBox() const { return false; }
  /// Return true if the solid is a tube.
  virtual bool IsTube() const { return false; }
  /// Return true if the solid is a sphere.
  virtual bool IsSphere() const { return false; }
  /// Return true if the solid is a hole.
  virtual bool IsHole() const { return false; }
  /// Return true if the solid is a ridge.
  virtual bool IsRidge() const { return false; }
  /// Return true if the solid is an extrusion.
  virtual bool IsExtrusion() const { return false; }
  /// Return true if the solid is a wire.
  virtual bool IsWire() const { return false; }
 
  /// Set a label.
  void SetLabel(const std::string& label) { m_label = label; }
  /// Return the label.
  std::string GetLabel() const { return m_label; }
 
  /// Retrieve the centre point of the solid.
  bool GetCentre(double& x, double& y, double& z) const {
    x = m_cX;
    y = m_cY;
    z = m_cZ;
    return true;
  }
  /// Retrieve the direction vector.
  bool GetDirection(double& dx, double& dy, double& dz) const {
    dx = m_dX;
    dy = m_dY;
    dz = m_dZ;
    return true;
  }
  /// Retrieve the orientation (azimuthal and polar angles) of the solid.
  bool GetOrientation(double& ctheta, double& stheta, double& cphi,
                      double& sphi) const {
    ctheta = m_cTheta;
    stheta = m_sTheta;
    cphi = m_cPhi;
    sphi = m_sPhi;
    return true;
  }

  /// Return the half-length along x.
  virtual double GetHalfLengthX() const {
    return NotImplemented("GetHalfLengthX");
  }
  /// Return the half-length along y.
  virtual double GetHalfLengthY() const {
    return NotImplemented("GetHalfLengthY");
  }
  /// Return the half-length along z.
  virtual double GetHalfLengthZ() const {
    return NotImplemented("GetHalfLengthZ");
  }
  /// Return the inner radius.
  virtual double GetInnerRadius() const {
    return NotImplemented("GetInnerRadius");
  }
  /// Return the outer radius.
  virtual double GetOuterRadius() const {
    return NotImplemented("GetOuterRadius");
  }
  /// Return the radius.
  virtual double GetRadius() const { return NotImplemented("GetRadius"); }
  /// Return the lower radius (of a hole).
  virtual double GetLowerRadius() const { 
    return NotImplemented("GetLowerRadius");
  } 
  /// Return the upper radius (of a hole).
  virtual double GetUpperRadius() const { 
    return NotImplemented("GetUpperRadius");
  } 
  /// Return the x-offset of a ridge.
  virtual double GetRidgeOffset() const {
    return NotImplemented("GetRidgeOffset");
  }
  /// Return the height of a ridge.
  virtual double GetRidgeHeight() const {
    return NotImplemented("GetRidgeHeight");
  }
  /// Get the vertices defining an extrusion.
  virtual bool GetProfile(std::vector<double>& xv, 
                          std::vector<double>& yv) const;

  /// Get the ID of the solid.
  unsigned int GetId() const { return m_id; }

  /// Retrieve the surface panels of the solid.
  virtual bool SolidPanels(std::vector<Panel>& panels) = 0;

  /// Set the discretisation level (for all panels).
  virtual void SetDiscretisationLevel(const double dis) = 0;
  /// Retrieve the discretisation level of a panel.
  virtual double GetDiscretisationLevel(const Panel& panel) = 0;

  virtual void Cut(const double x0, const double y0, const double z0,
                   const double xn, const double yn, const double zn,
                   std::vector<Panel>& panels) = 0;

  enum BoundaryCondition {
    Unknown = 0,
    Voltage,
    Charge,
    Float,
    Dielectric,
    DielectricCharge,
    ParallelField,
    PerpendicularField
  };

  /// Apply Dirichlet boundary conditions (fixed voltage).
  void SetBoundaryPotential(const double v) {
    m_volt = v;
    m_bctype = Voltage;
  }
  /// Apply fixed-charge boundary conditions.
  void SetBoundaryChargeDensity(const double q) {
    m_charge = q;
    m_bctype = Charge;
  }
  /// Make the potential at the surface of the solid floating.
  void SetBoundaryFloat() { m_bctype = Float; }
  /// Make the surfaces of the solid dielectric-dielectric interfaces.
  void SetBoundaryDielectric() { m_bctype = Dielectric; }
  void SetBoundaryParallelField() { m_bctype = ParallelField; }
  void SetBoundaryPerpendicularField() { m_bctype = PerpendicularField; }

  /// Retrieve the type of boundary condition.
  BoundaryCondition GetBoundaryConditionType() const { return m_bctype; }
  /// Retrieve the potential.
  double GetBoundaryPotential() const { return m_volt; }
  /// Retrieve the surface charge density.
  double GetBoundaryChargeDensity() const { return m_charge; }

  /// Switch debugging messages on/off.
  void EnableDebugging(const bool on = true) { m_debug = on; }

  /// Set the colour of the solid.
  void SetColour(const int col) { m_colour = col; }
  /// Get the colour of the solid.
  int GetColour() const { return m_colour; }
 
  static bool Intersect(const double x1, const double y1, const double z1,
                        const double x2, const double y2, const double z2,
                        const double x0, const double y0, const double z0,
                        const double a, const double b, const double c,
                        double& xc, double& yc, double& zc);

 protected:
  /// Centre of the solid.
  double m_cX = 0., m_cY = 0., m_cZ = 0.;

  /// Direction vector.
  double m_dX = 0., m_dY = 0., m_dZ = 1.;
  /// Azimuthal angle.
  double m_cPhi = 1., m_sPhi = 0.;
  /// Polar angle.
  double m_cTheta = 1., m_sTheta = 0.;

  /// Class name.
  std::string m_className = "Solid";
 
  /// Label.
  std::string m_label = "";

  /// Debug flag.
  bool m_debug = false;

  /// Type of boundary condition.
  BoundaryCondition m_bctype = Unknown;
  /// Potential at the surface.
  double m_volt = 0.;
  /// Surface charge density.
  double m_charge = 0.;
  /// Dielectric constant.
  double m_eps = 0.;

  /// Colour.
  int m_colour = -1;

  /// Transform a point from global coordinates (x, y, z) 
  /// to local coordinates (u, v, w).
  void ToLocal(const double x, const double y, const double z, double& u,
               double& v, double& w) const {
    const double dx = x - m_cX;
    const double dy = y - m_cY;
    const double dz = z - m_cZ;

    u = m_cPhi * m_cTheta * dx + m_sPhi * m_cTheta * dy - m_sTheta * dz;
    v = -m_sPhi * dx + m_cPhi * dy;
    w = m_cPhi * m_sTheta * dx + m_sPhi * m_sTheta * dy + m_cTheta * dz;
  }
  /// Transform a point from local coordinates (u, v, w) 
  /// to global coordinates (x, y, z).
  void ToGlobal(const double u, const double v, const double w, double& x,
                double& y, double& z) const {
    x = m_cX + m_cPhi * m_cTheta * u - m_sPhi * v + m_cPhi * m_sTheta * w;
    y = m_cY + m_sPhi * m_cTheta * u + m_cPhi * v + m_sPhi * m_sTheta * w;
    z = m_cZ - m_sTheta * u + m_cTheta * w;
  }
  /// Transform a vector from global to local coordinates.
  void VectorToLocal(const double x, const double y, const double z,
                     double& u, double& v, double& w) {
    u = m_cPhi * m_cTheta * x + m_sPhi * m_cTheta * y - m_sTheta * z;
    v = -m_sPhi * x + m_cPhi * y;
    w = m_cPhi * m_sTheta * x + m_sPhi * m_sTheta * y + m_cTheta * z;
  }

  void SetDirection(const double dx, const double dy, const double dz);

 private:
  double NotImplemented(const std::string& fcn) const;
  /// ID counter.
  static unsigned int s_id;
  /// ID of the solid.
  unsigned int m_id;

};
}

#endif
