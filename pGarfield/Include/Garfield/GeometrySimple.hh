#ifndef G_GEOMETRY_SIMPLE_H
#define G_GEOMETRY_SIMPLE_H

#include <array>
#include <vector>

#include "Geometry.hh"

namespace Garfield {

/// "Native" geometry, using simple shapes.

class GeometrySimple : public Geometry {
 public:
  /// Constructor
  GeometrySimple();
  /// Destructor
  virtual ~GeometrySimple() {}

  Medium* GetMedium(const double x, const double y, const double z,
                    const bool tesselated = false) const override;

  size_t GetNumberOfSolids() const override { return m_solids.size(); }
  Solid* GetSolid(const size_t i) const override;
  Solid* GetSolid(const size_t i, Medium*& medium) const override;

  /// Add a solid to the geometry, together with the medium inside.
  void AddSolid(Solid* s, Medium* m);
  /// Get the solid at a given location (x, y, z).
  Solid* GetSolid(const double x, const double y, const double z,
                  const bool tesselated = false) const;

  /// Set a background medium.
  void SetMedium(Medium* medium) { m_medium = medium; }

  /// Reset the geometry.
  void Clear();
  /// Print a summary of the solids present in the geometry.
  void PrintSolids();

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated = false) const override;
  /// Determine whether a point is inside the envelope of the geometry.
  bool IsInBoundingBox(const double x, const double y, const double z) const;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override {
    xmin = m_bbMin[0];
    ymin = m_bbMin[1];
    zmin = m_bbMin[2];
    xmax = m_bbMax[0];
    ymax = m_bbMax[1];
    zmax = m_bbMax[2];
    return true;
  }

  /// Switch on/off debugging and warning messages.
  void EnableDebugging(const bool on = true) { m_debug = on; }

 protected:
  /// List of solids and associated media.
  std::vector<std::pair<Solid*, Medium*> > m_solids;
  /// Background medium.
  Medium* m_medium = nullptr;

  // Bounding box ranges
  bool m_hasBoundingBox = false;
  std::array<double, 3> m_bbMin = {{0., 0., 0.}};
  std::array<double, 3> m_bbMax = {{0., 0., 0.}};

  /// Switch on/off debugging messages.
  bool m_debug = false;
};
}

#endif
