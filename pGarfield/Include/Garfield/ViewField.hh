#ifndef G_VIEW_FIELD
#define G_VIEW_FIELD

#include <Rtypes.h>

#include "ViewBase.hh"

namespace Garfield {

class Sensor;
class Component;

/// Visualize the potential or electric field of a component or sensor.

class ViewField : public ViewBase {
 public:
  /// Constructor.
  ViewField();
  /// Destructor.
  ~ViewField() = default;

  /// Set the sensor for which to plot the field.
  void SetSensor(Sensor* s);
  /// Set the component for which to plot the field.
  void SetComponent(Component* c);

  /// Set the plot limits for the potential.
  void SetVoltageRange(const double vmin, const double vmax);
  /// Set the plot limits for the electric field.
  void SetElectricFieldRange(const double emin, const double emax);
  /// Set the plot limits for the weighting field.
  void SetWeightingFieldRange(const double wmin, const double wmax);
  /// Set the plot limits for the magnetic field.
  void SetMagneticFieldRange(const double bmin, const double bmax);

  /// Set the number of contour levels.
  void SetNumberOfContours(const unsigned int n);
  /// Set the number of points used for drawing 1D functions.
  void SetNumberOfSamples1d(const unsigned int n);
  /// Set the number of points used for drawing 2D functions.
  void SetNumberOfSamples2d(const unsigned int nx, const unsigned int ny);

  /** Make a contour plot of the electric potential, electric field,
    * or magnetic field.
    * \param option quantity to be plotted
    * - potential: "v", "voltage", "p", "potential"
    * - magnitude of the electric field: "emag", "field"
    * - x-component of the electric field: "ex"
    * - y-component of the electric field: "ey"
    * - z-component of the electric field: "ez"
    * - magnitude of the magnetic field: "bmag"
    * - x-component of the magnetic field: "bx"
    * - y-component of the magnetic field: "by"
    * - z-component of the magnetic field: "bz"
    **/
  void PlotContour(const std::string& option = "v");
  /** Make a 2D plot of the electric potential, electric field or 
    * magnetic field.
    * \param option quantity to be plotted (see PlotContour)
    * \param drawopt option string passed to TF2::Draw
    **/
  void Plot(const std::string& option = "v",
            const std::string& drawopt = "arr");
  /** Make a 1D plot of the potential or field along a line.
    * \param x0,y0,z0 starting point
    * \param x1,y1,z1 end point
    * \param option quantity to be plotted (see PlotContour)
    * \param normalised flag whether to use normalised x-axis coordinates
    **/
  void PlotProfile(const double x0, const double y0, const double z0,
                   const double x1, const double y1, const double z1,
                   const std::string& option = "v",
                   const bool normalised = true);

  /** Make a contour plot of the weighting potential or field.
    * \param label identifier of the electrode
    * \param option quantity to be plotted (see PlotContour)
    **/
  void PlotContourWeightingField(const std::string& label,
                                 const std::string& option);
    /** Make a 2D plot of the weighting potential or field.
      * \param label identifier of the electrode
      * \param option quantity to be plotted (see PlotContour)
      * \param drawopt option string passed to TF2::Draw
      *\param t time slice of dynamic weighting potential [ns].
      **/
    void PlotWeightingField(const std::string& label, const std::string& option,
                            const std::string& drawopt, const double t = 0.);

  /** Make a 1D plot of the weighting potential or field along a line.
    * \param label identifier of the electrode
    * \param x0,y0,z0 starting point
    * \param x1,y1,z1 end point
    * \param option quantity to be plotted (see PlotContour)
    * \param normalised flag whether to use normalised x-axis coordinates
    **/
  void PlotProfileWeightingField(const std::string& label, const double x0,
                                 const double y0, const double z0,
                                 const double x1, const double y1,
                                 const double z1,
                                 const std::string& option = "v",
                                 const bool normalised = true);
  /// Determine the range of the potential/field automatically (true)
  /// or set it explicitly (false). 
  void EnableAutoRange(const bool on = true, 
                       const bool samplePotential = true) { 
    m_useAutoRange = on;
    m_samplePotential = samplePotential; 
  }

  /** Make use (or not) of the status flag returned by the sensor/component.
    * \param on Take status flag into account (true) or ignore it (false).
    * \param v0 Value to be used for regions with status != 0.
    */
  void AcknowledgeStatus(const bool on, const double v0 = 0.) {
    m_useStatus = on;
    m_vBkg = v0;
  }

  /// Draw electric field lines from a set of starting points. 
  void PlotFieldLines(const std::vector<double>& x0,
                      const std::vector<double>& y0,
                      const std::vector<double>& z0, 
                      const bool electron = true, const bool axis = true,
                      const short col = kOrange - 3);
  /// Generates point along a line, spaced by equal flux intervals. 
  bool EqualFluxIntervals(const double x0, const double y0, const double z0,
                          const double x1, const double y1, const double z1,
                          std::vector<double>& xf, std::vector<double>& yf,
                          std::vector<double>& zf, 
                          const unsigned int nPoints = 20) const;
  /// Generate points along a line, spaced by a given flux interval.
  bool FixedFluxIntervals(const double x0, const double y0, const double z0,
                          const double x1, const double y1, const double z1,
                          std::vector<double>& xf, std::vector<double>& yf,
                          std::vector<double>& zf, 
                          const double interval = 10.) const;

 private:
  enum class Parameter { 
    Potential = 0, 
    Emag,
    Ex, 
    Ey, 
    Ez, 
    Bmag,
    Bx, 
    By, 
    Bz, 
    Unknown };

  bool m_useAutoRange = true;
  bool m_samplePotential = true;
  bool m_useStatus = false;
  double m_vBkg = 0.;

  // Sensor
  Sensor* m_sensor = nullptr;
  Component* m_component = nullptr;

  // Function range
  double m_vmin = 0., m_vmax = 100.;
  double m_emin = 0., m_emax = 10000.;
  double m_wmin = 0., m_wmax = 100.;
  double m_bmin = 0., m_bmax = 10.;

  // Number of contours
  unsigned int m_nContours = 20;
  // Number of points used to draw the functions
  unsigned int m_nSamples1d = 1000;
  unsigned int m_nSamples2dX = 200;
  unsigned int m_nSamples2dY = 200;

  bool SetPlotLimits();
  void Draw2d(const std::string& option, const bool contour,
              const bool wfield, const std::string& electrode,
              const std::string& drawopt, const double t = 0.);
  void DrawProfile(const double x0, const double y0, const double z0,
                   const double x1, const double y1, const double z1,
                   const std::string& option, const bool wfield,
                   const std::string& electrode, const bool normalised);
  Parameter GetPar(const std::string& option, std::string& title,
                   bool& bfield) const;
  double Efield(const double x, const double y, const double z,
                const Parameter par) const;
  double Wfield(const double x, const double y, const double z,
                const Parameter par, const std::string& electrode, const double t = 0.) const;
  double Bfield(const double x, const double y, const double z,
                const Parameter par) const;

};
}
#endif
