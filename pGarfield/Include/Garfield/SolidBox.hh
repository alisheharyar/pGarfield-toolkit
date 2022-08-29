#ifndef G_SOLID_BOX_H
#define G_SOLID_BOX_H

#include "Solid.hh"

namespace Garfield {

/// Box.

class SolidBox : public Solid {
 public:
  /// Constructor from centre and half-widths.
  SolidBox(const double cx, const double cy, const double cz, const double lx,
           const double ly, const double lz);
  /// Constructor from centre, half-widths, and orientation.
  SolidBox(const double cx, const double cy, const double cz, const double lx,
           const double ly, const double lz, const double dx, const double dy,
           const double dz);
  /// Destructor
  ~SolidBox() {}

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsBox() const override { return true; }

  double GetHalfLengthX() const override { return m_lX; }
  double GetHalfLengthY() const override { return m_lY; }
  double GetHalfLengthZ() const override { return m_lZ; }

  void SetHalfLengthX(const double lx);
  void SetHalfLengthY(const double ly);
  void SetHalfLengthZ(const double lz);

  bool SolidPanels(std::vector<Panel>& panels) override;
  void SetDiscretisationLevel(const double dis) override {
    m_dis.fill(dis);
  } 
  double GetDiscretisationLevel(const Panel& panel) override;

  void Cut(const double x0, const double y0, const double z0,
           const double xn, const double yn, const double zn,
           std::vector<Panel>& panels) override;

 private:
  /// Half lengths.
  double m_lX = 0., m_lY = 0., m_lZ = 0.;
  /// Discretisation levels.
  std::array<double, 6> m_dis{{-1, -1, -1, -1, -1, -1}};
};
}

#endif
