#ifndef G_SOLID_RIDGE_H
#define G_SOLID_RIDGE_H

#include "Solid.hh"

namespace Garfield {

/// Triangular prism (Toblerone bar).

class SolidRidge : public Solid {
 public:
  /// Constructor from centre, half-lengths, height and x-offset. 
  SolidRidge(const double cx, const double cy, const double cz,
             const double lx, const double ly, const double hz, 
             const double offsetx);
  /// Constructor from centre, half-lengths, height, x-offset and orientation. 
  SolidRidge(const double cx, const double cy, const double cz,
             const double lx, const double ly, const double hz, 
             const double offsetx,
             const double dx, const double dy, const double dz);
  /// Destructor
  ~SolidRidge() {}

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsRidge() const override { return true; }

  /// Set the half-length along x.
  void SetHalfLengthX(const double lx);
  /// Set the half-length along y.
  void SetHalfLengthY(const double ly);
  /// Set the height of the ridge.
  void SetRidgeHeight(const double hz);
  /// Set the x-offset of the ridge.
  void SetRidgeOffset(const double dx) { m_hx = dx; }

  double GetHalfLengthX() const override { return m_lX; }
  double GetHalfLengthY() const override { return m_lY; }
  double GetRidgeHeight() const override { return m_hz; }
  double GetRidgeOffset() const override { return m_hx; }

  bool SolidPanels(std::vector<Panel>& panels) override;
  void SetDiscretisationLevel(const double dis) override {
    m_dis.fill(dis);
  }
  double GetDiscretisationLevel(const Panel& panel) override;

  void Cut(const double x0, const double y0, const double z0,
           const double xn, const double yn, const double zn,
           std::vector<Panel>& panels) override;

 private:
  /// Half-length in x.
  double m_lX;
  /// Half-length in y.
  double m_lY;
  /// Height of the ridge.
  double m_hz;
  /// Offset of the ridge in x.
  double m_hx;
  
  /// Discretisation levels.
  std::array<double, 5> m_dis{{-1., -1., -1., -1., -1.}}; 
};
}

#endif
