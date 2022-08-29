#ifndef G_VIEW_MEDIUM
#define G_VIEW_MEDIUM

#include <string>
#include <vector>
#include <array>

#include "ViewBase.hh"
#include "FundamentalConstants.hh"

namespace Garfield {

class Medium;

/// Plot transport coefficients as function of electric and magnetic field.

class ViewMedium : public ViewBase {
 public:
  /// Constructor.
  ViewMedium();
  /// Destructor.
  ~ViewMedium() = default;

  /// Set the medium from which to retrieve the transport coefficients.
  void SetMedium(Medium* m);

  /// Try to choose the x-axis range based on the field grid. 
  void EnableAutoRangeX(const bool on = true) { m_autoRangeX = on; } 
  /// Set the limits of the electric field.
  void SetRangeE(const double emin, const double emax, const bool logscale);
  /// Set the limits of the magnetic field.
  void SetRangeB(const double bmin, const double bmax, const bool logscale);
  /// Set the limits of the angle between electric and magnetic field.
  void SetRangeA(const double amin, const double amax, const bool logscale);
  /// Choose the y-axis range based on the function's minima/maxima.
  void EnableAutoRangeY(const bool on = true) { m_autoRangeY = on; }
  /// Set the range of the function (velocity etc.) to be plotted.
  void SetRangeY(const double ymin, const double ymax, 
                 const bool logscale = false);

  /// Set the electric field to use when plotting as function of B or angle.
  void SetElectricField(const double efield) { m_efield = efield; }
  /// Set the magnetic field to use when plotting as function of E or angle.
  void SetMagneticField(const double bfield) { m_bfield = bfield; }
  /// Set the angle to use when plotting as function of E or B.
  void SetAngle(const double angle) { m_angle = angle; }

  void EnableExport(const std::string& txtfile) { m_outfile = txtfile; }
  void DisableExport() { m_outfile = ""; }

  /** Plot the drift velocity components of electrons in the medium.
    * \param xaxis abscissa.
    *   - 'e': electric field, 
    *   - 'b': magnetic field, 
    *   - 'a': angle between E and B.
    * \param same flag to keep existing plots (true) or not.
    */
  void PlotElectronVelocity(const char xaxis = 'e', const bool same = false) {
    PlotVelocity(GetAxis(xaxis), Charge::Electron, same);
  }
  /// Plot the drift velocity components of holes in the medium.
  void PlotHoleVelocity(const char xaxis = 'e', const bool same = false) {
    PlotVelocity(GetAxis(xaxis), Charge::Hole, same);
  }
  /// Plot the ion drift velocity in the gas.
  void PlotIonVelocity(const char xaxis = 'e', const bool same = false) {
    PlotVelocity(GetAxis(xaxis), Charge::Ion, same);
  }
  /// Plot the diffusion coefficients in the medium.
  void PlotElectronDiffusion(const char xaxis = 'e', const bool same = false) {
    PlotDiffusion(GetAxis(xaxis), Charge::Electron, same);
  }
  /// Plot the diffusion coefficients of holes in the medium.
  void PlotHoleDiffusion(const char xaxis = 'e', const bool same = false) {
    PlotDiffusion(GetAxis(xaxis), Charge::Hole, same);
  }
  /// Plot the diffusion coefficients of ions in the gas.
  void PlotIonDiffusion(const char xaxis = 'e', const bool same = false) {
    PlotDiffusion(GetAxis(xaxis), Charge::Ion, same);
  }
  /// Plot the Townsend coefficient for electrons.
  void PlotElectronTownsend(const char xaxis = 'e', const bool same = false) {
    Plot(GetAxis(xaxis), Charge::Electron, Parameter::Townsend, same);
  }
  /// Plot the Townsend coefficient for holes.
  void PlotHoleTownsend(const char xaxis = 'e', const bool same = false) {
    Plot(GetAxis(xaxis), Charge::Hole, Parameter::Townsend, same);
  }
  /// Plot the attachment coefficient for electrons.
  void PlotElectronAttachment(const char xaxis = 'e', const bool same = false) {
    Plot(GetAxis(xaxis), Charge::Electron, Parameter::Attachment, same);
  }
  /// Plot the attachment coefficient for holes.
  void PlotHoleAttachment(const char xaxis = 'e', const bool same = false) {
    Plot(GetAxis(xaxis), Charge::Hole, Parameter::Attachment, same);
  }
  /// Plot the angle between drift velocity and field.
  void PlotElectronLorentzAngle(const char xaxis = 'e', const bool same = false) {
    PlotLorentzAngle(GetAxis(xaxis), Charge::Electron, same);
  }

  /// Set the (ROOT) colours to be used in the plots.
  void SetColours(const std::vector<short>& cols) { m_colours = cols; }
  /// Set user-defined plot labels.
  void SetLabels(const std::vector<std::string>& labels) { m_labels = labels; }

 private:

  enum class Parameter {
    VelocityE,
    VelocityB,
    VelocityExB,
    TransverseDiffusion,
    LongitudinalDiffusion,
    Townsend,
    Attachment,
    LorentzAngle
  };
 
  enum class Charge {
    Electron,
    Hole,
    Ion
  };

  enum class Axis {
    E,
    B,
    Angle,
    None
  };

  Medium* m_medium = nullptr;

  // X axis
  double m_eMin = 100., m_eMax = 100000.;
  double m_bMin = 0., m_bMax = 2.;
  double m_aMin = 0., m_aMax = Pi;
  bool m_logE = true;
  bool m_logB = false;
  bool m_logA = false;
  bool m_logX = true;
  bool m_autoRangeX = true;
  Axis m_xaxis = Axis::None;

  // Y axis
  double m_yMin = 0., m_yMax = 1.;
  bool m_logY = false;
  bool m_autoRangeY = true;

  // E-field to use when plotting as function of B-field or angle.
  double m_efield = 1000.;
  // B-field to use when plotting as function of E-field or angle.
  double m_bfield = 0.;
  // Angle to use when plotting as function of E-field or B-field.
  double m_angle = HalfPi;

  std::vector<double> m_xPlot;
  std::vector<std::vector<double> > m_yPlot;
  std::vector<Parameter> m_par;
  std::vector<Charge> m_q;

  std::vector<std::vector<double> > m_xGraph;
  std::vector<std::vector<double> > m_yGraph;

  std::vector<short> m_colours;
  std::vector<std::string> m_labels;

  std::string m_outfile;

  void PlotVelocity(const Axis xaxis, const Charge particle,
                    const bool same);
  void PlotDiffusion(const Axis xaxis, const Charge particle,
                     const bool same);
  void Plot(const Axis xaxis, const Charge particle,
            const Parameter par, const bool same);
  void PlotLorentzAngle(const Axis xaxis, const Charge particle,
                        const bool same);

  void ResetY();
  void ResetX(const Axis xaxis);
  void Draw();
  void Export();

  Axis GetAxis(const char xaxis) const;
  bool GetGrid(std::array<std::vector<double>, 3>& grid,
               int& ie, int& ib, int& ia, const Axis xaxis) const;
};
}
#endif
