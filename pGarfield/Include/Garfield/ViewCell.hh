#ifndef G_VIEW_CELL
#define G_VIEW_CELL

#include <memory>
#include <string>

#include <TGeoManager.h>

#include "ViewBase.hh"

namespace Garfield {

class ComponentAnalyticField;
class ComponentNeBem2d;

/// Visualize the "cell" defined in an analytic-field component.

class ViewCell : public ViewBase {
 public:
  /// Constructor
  ViewCell();
  /// Destructor
  ~ViewCell() = default;

  /// Set the component for which to draw the cell geometry.
  void SetComponent(ComponentAnalyticField* comp);
  void SetComponent(ComponentNeBem2d* comp);

  /// Make a two-dimensional drawing of the cell layout.
  void Plot2d();
  /// Make a three-dimensional drawing of the cell layout (using TGeo).
  void Plot3d();

  /// Visualize wirers using markers (default setting) 
  /// or as a circle with the actual wire radius.
  void EnableWireMarkers(const bool on = true) { m_useWireMarker = on; }
  void DisableWireMarkers() { EnableWireMarkers(false); }

 private:
  bool m_useWireMarker = true;

  ComponentAnalyticField* m_component = nullptr;
  ComponentNeBem2d* m_nebem = nullptr;

  // 3D geometry.
  std::unique_ptr<TGeoManager> m_geo;

  bool Plot(const bool twod);
  // Draw a wire in 2D.
  void PlotWire(const double x, const double y, const double d, const int type);
  // Draw a wire in 3D.
  void PlotWire(const double x, const double y, const double d, const int type,
                const double lz);
  // Draw a tube in 2D.
  void PlotTube(const double x0, const double y0, const double r, const int n);
  // Draw a tube in 3D.
  void PlotTube(const double x0, const double y0,
                const double r1, const double r2, const int n,
                const double lz); 
  // Draw a plane in 2D.
  void PlotPlane(const double x0, const double y0, 
                 const double x1, const double y1);
  // Draw a plane in 3D.
  void PlotPlane(const double dx, const double dy, const double dz,
                 const double x0, const double y0);
  // Draw a neBEM 2D layout.
  bool PlotNeBem(const bool twod);
  // Setup the TGeoManager. 
  void SetupGeo(const double dx, const double dy, const double dz);

};
}
#endif
