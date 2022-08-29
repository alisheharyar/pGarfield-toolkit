#ifndef G_SOLID_TUBE_H
#define G_SOLID_TUBE_H

#include <mutex>

#include "Solid.hh"

namespace Garfield {

/// Cylindrical tube.

class SolidTube : public Solid {
 public:
  /// Constructor from centre, outer radius, and half-length.
  SolidTube(const double cx, const double cy, const double cz, const double r,
            const double lz);
  /// Constructor from centre, outer radius, half-length and orientation.
  SolidTube(const double cx, const double cy, const double cz, const double r,
            const double lz, const double dx, const double dy, const double dz);
  /// Destructor
  ~SolidTube() {}

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsTube() const override { return true; }

  void SetHalfLength(const double lz);
  void SetRadius(const double r);

  double GetHalfLengthZ() const override { return m_lZ; }
  double GetRadius() const override { return m_rMax; }

  /// When calculating the surface panels, the cylinder is
  /// approximated as a polygon with a finite number of panels.
  /// The number of corners of the polygon equals \f$4(n - 1)\f$.
  /// Thus, \f$n = 2\f$ will produce a square, \f$n = 3\f$ an octagon etc.
  void SetSectors(const unsigned int n);
  /// Specify a rotation angle (radian) of the cylinder.
  /// Such a rotation is meaningful only if the number of sectors
  /// (when approximating the circle with a polygon) has been chosen small.
  void SetRotation(const double angle) { m_rot = angle; }
  /// By default, the polygon used for approximating the cylinder when
  /// calculating surface panels is inscribed in a circle
  /// of the specified radius. If the "average-radius" flag is activated,
  /// then the radius will be interpreted as the mean radius of the polygon
  /// that approximates the cylinder.
  void SetAverageRadius(const bool average) { m_average = average; }
  /// Request the cylinder to be closed with a (polygonal) lid at +z.
  void SetTopLid(const bool closed) { m_toplid = closed; }
  /// Request the cylinder to be closed with a (polygonal) lid at -z.
  void SetBottomLid(const bool closed) { m_botlid = closed; }

  /// Return the number of sectors.
  unsigned int GetSectors() const { return m_n; }
  /// Return the current rotation angle.
  double GetRotation() const { return m_rot; }
  /// Return the status of the "average-radius" flag.
  bool GetAverage() const { return m_average; }

  bool SolidPanels(std::vector<Panel>& panels) override;
  void SetDiscretisationLevel(const double dis) override { m_dis.fill(dis); }
  double GetDiscretisationLevel(const Panel& panel) override;

  void Cut(const double x0, const double y0, const double z0,
           const double xn, const double yn, const double zn,
           std::vector<Panel>& panels) override;

 private:
  /// Mutex.
  std::mutex m_mutex;

  /// Outer radius.
  double m_rMax;
  /// Half-length
  double m_lZ;

  /// Rotation angle
  double m_rot = 0.;
  /// Number of sectors
  unsigned int m_n = 2;
  /// Average chord over the sectors.
  bool m_average = false;
  /// Radius of the approximating polygon.
  double m_rp;
  /// Inradius of the approximating polygon.
  double m_ri;
  /// X-coordinates of the approximating polygon.
  std::vector<double> m_xp;
  /// Y-coordinates of the approximating polygon.
  std::vector<double> m_yp;

  /// Have a top lid?
  bool m_toplid = true;
  /// Have a bottom lid?
  bool m_botlid = true;

  /// Discretisation levels.
  std::array<double, 3> m_dis{{-1.,-1.,-1.}};

  void UpdatePolygon();
};
}

#endif
