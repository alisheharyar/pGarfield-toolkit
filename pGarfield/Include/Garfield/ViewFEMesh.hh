#ifndef G_VIEW_FE_MESH
#define G_VIEW_FE_MESH

#include <memory>
#include <string>
#include <map>

#include <TArrayD.h>
#include <TGaxis.h>
#include <TGeoManager.h>
#include <TMatrixD.h>

#include "ComponentCST.hh"
#include "ComponentFieldMap.hh"
#include "ViewBase.hh"
#include "ViewDrift.hh"

namespace Garfield {

/// Draw the mesh of a field-map component.

class ViewFEMesh : public ViewBase {
 public:
  /// Constructor.
  ViewFEMesh();
  /// Destructor.
  ~ViewFEMesh();

  /// Set the component from which to retrieve the mesh and field.
  void SetComponent(ComponentFieldMap* cmp);

  void SetPlane(const double fx, const double fy, const double fz, 
                const double x0, const double y0, const double z0) override;
  void SetPlane(const double fx, const double fy, const double fz, 
                const double x0, const double y0, const double z0,
                const double hx, const double hy, const double hz) override;

  // Axes
  void SetXaxis(TGaxis* ax);
  void SetYaxis(TGaxis* ay);
  void SetXaxisTitle(const std::string& xtitle) { m_xaxisTitle = xtitle; }
  void SetYaxisTitle(const std::string& ytitle) { m_yaxisTitle = ytitle; }
  void EnableAxes() { m_drawAxes = true; }
  void DisableAxes() { m_drawAxes = false; }

  /// Plot method to be called by user
  bool Plot(const bool twod = true);

  /// Element fill switch; 2D only, set false for wireframe mesh
  void SetFillMesh(const bool f) { m_fillMesh = f; }

  /// Display intersection of projection plane with viewing area
  void SetDrawViewRegion(bool do_draw) { m_drawViewRegion = do_draw; }
  bool GetDrawViewRegion(void) const { return m_drawViewRegion; }

  /// Associate a color with each element material map ID;
  /// Uses ROOT color numberings
  void SetColor(int matID, int colorID) { m_colorMap[matID] = colorID; }
  void SetFillColor(int matID, int colorID) {
    m_colorMap_fill[matID] = colorID;
  }

  /// Set the optional associated ViewDrift
  void SetViewDrift(ViewDrift* vd) { m_viewDrift = vd; }

  /// Show filled mesh elements
  void SetFillMeshWithBorders() {
    m_plotMeshBorders = true;
    m_fillMesh = true;
  }

  /// Create a default set of custom-made axes.
  void CreateDefaultAxes();

  /// Disable a material so that its mesh cells are not drawn
  void DisableMaterial(int materialID) {
    m_disabledMaterial[materialID] = true;
  }

 private:
  // Options
  bool m_fillMesh = false;

  // Intersection of viewing plane with plotted area in planar coordinates
  bool m_drawViewRegion = false;
  std::vector<double> m_viewRegionX;
  std::vector<double> m_viewRegionY;

  // The field map object
  ComponentFieldMap* m_cmp = nullptr;

  // Optional associated ViewDrift object
  ViewDrift* m_viewDrift = nullptr;
  bool m_plotMeshBorders = false;

  // Axes
  TGaxis* m_xaxis = nullptr;
  TGaxis* m_yaxis = nullptr;
  std::string m_xaxisTitle = "";
  std::string m_yaxisTitle = "";
  bool m_drawAxes = false;

  // The color map
  std::map<int, int> m_colorMap;
  std::map<int, int> m_colorMap_fill;

  // Disabled materials -> not shown in the mesh view
  std::map<int, bool> m_disabledMaterial;

  std::vector<TGeoVolume*> m_volumes;
  std::vector<TGeoMedium*> m_media;
  std::unique_ptr<TGeoManager> m_geoManager;

  // Element plotting methods
  void DrawElements2d();
  void DrawElements3d();
  void DrawCST(ComponentCST* componentCST);

  void DrawDriftLines2d();
  void DrawDriftLines3d();

  bool GetPlotLimits();

  /// Return true if the specified point is in the view region.
  bool InView(const double x, const double y) const;

  bool LinesCrossed(double x1, double y1, double x2, double y2, double u1,
                    double v1, double u2, double v2, double& xc,
                    double& yc) const;
  bool IntersectPlaneArea(double& xmin, double& ymin,
                          double& xmax, double& ymax);
  bool OnLine(double x1, double y1, double x2, double y2, double u,
              double v) const;
  void RemoveCrossings(std::vector<double>& x, std::vector<double>& y);
  bool PlaneCut(double x1, double y1, double z1, double x2, double y2,
                double z2, TMatrixD& xMat);
  void ClipToView(std::vector<double>& px, std::vector<double>& py,
                  std::vector<double>& cx, std::vector<double>& cy);
  bool IsInPolygon(double x, double y, const std::vector<double>& px,
                   const std::vector<double>& py, bool& edge) const;

  void Reset();
};
}  // namespace Garfield
#endif
