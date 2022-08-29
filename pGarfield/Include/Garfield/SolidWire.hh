#ifndef G_SOLID_WIRE_H
#define G_SOLID_WIRE_H

#include "Solid.hh"

namespace Garfield {

/// Wire.

class SolidWire : public Solid {
 public:
  /// Constructor from centre, radius, and half-length.
  SolidWire(const double cx, const double cy, const double cz, const double r,
            const double lz);
  /// Constructor from centre, radius, half-length and orientation.
  SolidWire(const double cx, const double cy, const double cz, const double r,
            const double lz, const double dx, const double dy, const double dz);
  /// Destructor
  ~SolidWire() {}

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsWire() const override { return true; }

  void SetHalfLength(const double lz);
  void SetRadius(const double r);

  double GetHalfLengthZ() const override { return m_lZ; }
  double GetRadius() const override { return m_r; }

  bool SolidPanels(std::vector<Panel>& panels) override;
  void SetDiscretisationLevel(const double dis) override { m_dis = dis; }
  double GetDiscretisationLevel(const Panel& panel) override;

  void Cut(const double x0, const double y0, const double z0,
           const double xn, const double yn, const double zn,
           std::vector<Panel>& panels) override;

 private:
  /// Radius.
  double m_r;
  /// Half-length
  double m_lZ;

  /// Discretisation level.
  double m_dis = -1.;
};
}

#endif
