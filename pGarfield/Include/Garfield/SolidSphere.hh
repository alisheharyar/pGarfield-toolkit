#ifndef G_SOLID_SPHERE_H
#define G_SOLID_SPHERE_H

#include <mutex>

#include "Solid.hh"

namespace Garfield {

/// Sphere.

class SolidSphere : public Solid {
 public:
  /// Constructor from centre and outer radius.
  SolidSphere(const double cx, const double cy, const double cz,
              const double r);
  /// Constructor from centre and inner/outer radii.
  SolidSphere(const double cx, const double cy, const double cz,
              const double rmin, const double rmax);
  /// Destructor
  ~SolidSphere() {}

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) const override;
  bool IsSphere() const override { return true; }

  /// Set the radius of the sphere.
  void SetRadius(const double r);
  /// Set the inner and outer radius of the sphere.
  void SetRadii(const double rmin, const double rmax);

  double GetRadius() const override { return m_rMax; }
  double GetInnerRadius() const override { return m_rMin; }
  double GetOuterRadius() const override { return m_rMax; }

  /// When calculating surface panels, the sphere is approximated by a set of
  /// parallelograms, much the same way maps are drawn ("UV sphere"). 
  /// N specifies the number of meridians and also the number of parallels.
  void SetMeridians(const unsigned int n);

  bool SolidPanels(std::vector<Panel>& panels) override;
  void SetDiscretisationLevel(const double dis) override { m_dis = dis; }
  double GetDiscretisationLevel(const Panel& panel) override;

  void Cut(const double x0, const double y0, const double z0,
           const double xn, const double yn, const double zn,
           std::vector<Panel>& panels) override;

 private:
  /// Mutex. 
  std::mutex m_mutex;

  /// Inner and outer radii.
  double m_rMin = 0., m_rMax = 1.;

  /// Number of meridians.
  unsigned int m_n = 10;

  /// Discretisation level.
  double m_dis = -1.;

  /// Surface panels.
  std::vector<Panel> m_panelsO;
  std::vector<Panel> m_panelsI;

  void UpdatePanels();
  void MakePanels(const int vol, const double r, const bool out,
                  std::vector<Panel>& panels) const; 
};
}

#endif
