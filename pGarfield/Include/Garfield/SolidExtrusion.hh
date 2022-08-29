#ifndef G_SOLID_EXTRUSION_H
#define G_SOLID_EXTRUSION_H

#include "Solid.hh"

namespace Garfield {

/// Extrusion.

class SolidExtrusion : public Solid {
 public:
  /// Constructor from half-length and profile.
  SolidExtrusion(const double lz,
                 const std::vector<double>& xp, const std::vector<double>& yp);
  /// Constructor from half-length, profile, offset and orientation.
  SolidExtrusion(const double lz,
                 const std::vector<double>& xp, const std::vector<double>& yp,
                 const double cx, const double cy, const double cz, 
                 const double dx, const double dy, const double dz);
  /// Destructor
  ~SolidExtrusion() {}

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsExtrusion() const override { return true; }

  double GetHalfLengthZ() const override { return m_lZ; }
  bool GetProfile(std::vector<double>& xp, 
                  std::vector<double>& yp) const override {
    xp = m_xp;
    yp = m_yp;
    return !m_xp.empty();
  }
  /// Set the half-length of the extrusion.
  void SetHalfLengthZ(const double lz);
  /// Set the coordinates of the extrusion profile.
  void SetProfile(const std::vector<double>& xp,
                  const std::vector<double>& yp);
  /// Request the extrusion to be closed with a lid at +z.
  void SetTopLid(const bool closed) { m_toplid = closed; }
  /// Request the extrusion to be closed with a lid at -z.
  void SetBottomLid(const bool closed) { m_botlid = closed; }

  bool SolidPanels(std::vector<Panel>& panels) override;
  void SetDiscretisationLevel(const double dis) override {
    m_dis.fill(dis);
  }
  double GetDiscretisationLevel(const Panel& panel) override;

  void Cut(const double x0, const double y0, const double z0,
           const double xn, const double yn, const double zn,
           std::vector<Panel>& panels) override;

 private:
  /// Half length.
  double m_lZ = 0.;
  /// X coordinates of the profile.
  std::vector<double> m_xp;
  /// Y coordinates of the profile.
  std::vector<double> m_yp;

  /// Have a top lid?
  bool m_toplid = true;
  /// Have a bottom lid?
  bool m_botlid = true;

  /// Orientation of the polygon.
  bool m_clockwise = true;

  /// Discretisation levels.
  std::array<double, 3> m_dis{{-1, -1, -1}};
};
}

#endif
