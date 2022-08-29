#ifndef G_VIEW_BASE
#define G_VIEW_BASE

#include <memory>
#include <array>
#include <string>

#include <TPad.h>
#include <TCanvas.h>

namespace Garfield {

class Sensor;
class Component;

/// Base class for visualization classes.

class ViewBase {
 public:
  /// Default constructor.
  ViewBase() = delete;
  /// Constructor.
  ViewBase(const std::string& name);
  /// Destructor.
  virtual ~ViewBase() = default;

  /// Set the canvas to be painted on.
  void SetCanvas(TPad* pad) { m_pad = pad; }
  /// Unset an external canvas.
  void SetCanvas() { m_pad = nullptr; }
  /// Retrieve the canvas.
  TPad* GetCanvas();

  /// Set the x- and y-axis limits
  /// (in local coordinates of the current viewing plane, if applicable).
  void SetArea(const double xmin, const double ymin, const double xmax,
               const double ymax);
  /// Set a bounding box (if applicable).
  virtual void SetArea(const double xmin, const double ymin, const double zmin,
                       const double xmax, const double ymax, const double zmax);
  /// Use default x- and y-axis limits
  /// (based on the bounding box of the sensor/component, if applicable).
  void SetArea() {
    m_userBox = false; 
    m_userPlotLimits = false; 
  }

  /** Set the projection (viewing plane), if applicable.
    * \param fx,fy,fz normal vector
    * \param x0,y0,z0 in-plane point
    */
  virtual void SetPlane(const double fx, const double fy, const double fz,
                        const double x0, const double y0, const double z0);
  /// Set the projection plane specifying a hint for the in-plane x axis. 
  virtual void SetPlane(const double fx, const double fy, const double fz,
                        const double x0, const double y0, const double z0,
                        const double hx, const double hy, const double hz); 
  /// Rotate the viewing plane (angle in radian).
  void Rotate(const double angle);
  /// Set the viewing plane to x-y. 
  void SetPlaneXY();
  /// Set the viewing plane to x-z. 
  void SetPlaneXZ();
  /// Set the viewing plane to y-z. 
  void SetPlaneYZ();
  /// Set the viewing plane to z-x. 
  void SetPlaneZX();
  /// Set the viewing plane to z-y. 
  void SetPlaneZY();

  /// Switch on/off debugging output.
  void EnableDebugging(const bool on = true) { m_debug = on; }

  /// Find an unused function name.
  static std::string FindUnusedFunctionName(const std::string& s);
  /// Find an unused histogram name.
  static std::string FindUnusedHistogramName(const std::string& s);
  /// Find an unused canvas name.
  static std::string FindUnusedCanvasName(const std::string& s);

 protected:
  std::string m_className = "ViewBase";

  // Options
  bool m_debug = false;

  // Plot axis limits. 
  bool m_userPlotLimits = false;
  double m_xMinPlot = -1., m_xMaxPlot = 1.;
  double m_yMinPlot = -1., m_yMaxPlot = 1.;

  // Bounding box.
  bool m_userBox = false;
  double m_xMinBox = -1., m_xMaxBox = 1.;
  double m_yMinBox = -1., m_yMaxBox = 1.; 
  double m_zMinBox = -1., m_zMaxBox = 1.;

  // Viewing plane (FPROJ).
  // Default projection: x-y at z = 0.
  std::array<std::array<double, 3>, 3> m_proj{{
    {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 0}}
  }};
  std::array<double, 4> m_plane{{0, 0, 1, 0}};
  // Matrix used for projections (FPRMAT).
  std::array<std::array<double, 3>, 3> m_prmat{{
    {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}
  }};

  // Update and invert the projection matrix.
  void UpdateProjectionMatrix();
  // Determine plane coordinates.
  template<typename T> 
  void ToPlane(const T x, const T y, const T z, T& xp, T& yp) const {
    xp = m_prmat[0][0] * x + m_prmat[0][1] * y + m_prmat[0][2] * z;
    yp = m_prmat[1][0] * x + m_prmat[1][1] * y + m_prmat[1][2] * z;
  }
  // Determine whether a point is inside the bounding box.
  template<typename T> 
  bool InBox(const std::array<T, 3>& x) const {
    if (!m_userBox) return true;
    if (x[0] < m_xMinBox || x[0] > m_xMaxBox || 
        x[1] < m_yMinBox || x[1] > m_yMaxBox ||
        x[2] < m_yMinBox || x[2] > m_zMaxBox) return false;
    return true; 
  }
  // Clip the line x0 - x1 to the extent of the bounding box.
  void Clip(const std::array<float, 3>& x0,
            const std::array<float, 3>& x1, std::array<float, 3>& xc) const;
  // Draw the projection of a line onto the current viewing plane.
  void DrawLine(const std::vector<std::array<float, 3> >& xl,
                const short col, const short lw);

  // X-axis label for the current viewing plane.
  std::string LabelX();
  // Y-axis label for the current viewing plane.
  std::string LabelY();
  // Description of the current viewing plane.
  std::string PlaneDescription();

  bool PlotLimits(Sensor* sensor, double& xmin, double& ymin,
                  double& xmax, double& ymax) const;
  bool PlotLimits(Component* cmp, double& xmin, double& ymin,
                  double& xmax, double& ymax) const;
  bool PlotLimitsFromUserBox(double& xmin, double& ymin, 
                             double& xmax, double& ymax) const;
  bool PlotLimits(std::array<double, 3>& bbmin,
                  std::array<double, 3>& bbmax,
                  double& xmin, double& ymin,
                  double& xmax, double& ymax) const;

  static bool RangeSet(TVirtualPad*);
  static void SetRange(TVirtualPad* pad, const double x0, const double y0,
                       const double x1, const double y1);
 private:
  // Current pad.
  TPad* m_pad = nullptr;
  std::unique_ptr<TCanvas> m_canvas;
};
}
#endif
