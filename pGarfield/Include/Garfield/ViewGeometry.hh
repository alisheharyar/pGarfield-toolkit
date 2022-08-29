#ifndef G_VIEW_GEOMETRY
#define G_VIEW_GEOMETRY

#include <memory>
#include <string>
#include <vector>

#include <TGeoManager.h>

#include "ViewBase.hh"

namespace Garfield {

class GeometrySimple;

/// Visualize a geometry defined using the "native" shapes.

class ViewGeometry : public ViewBase {
 public:
  /// Constructor.
  ViewGeometry();
  /// Destructor.
  ~ViewGeometry();

  /// Set the geometry to be drawn.
  void SetGeometry(GeometrySimple* geo);
  /// Draw the geometry.
  void Plot(const bool twod = false);
  /// Draw a cut through the geometry at the current viewing plane. 
  void Plot2d();
  /// Draw a three-dimensional view of the geometry.
  void Plot3d();

 private:
  GeometrySimple* m_geometry = nullptr;

  std::vector<TGeoVolume*> m_volumes;
  std::vector<TGeoMedium*> m_media;

  std::unique_ptr<TGeoManager> m_geoManager;

  void Reset();
};
}
#endif
