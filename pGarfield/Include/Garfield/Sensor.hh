#ifndef G_SENSOR_H
#define G_SENSOR_H

#include <fstream>
#include <functional>
#include <mutex>
#include <tuple>
#include <utility>
#include <vector>

#include "Component.hh"
#include "Shaper.hh"

namespace Garfield {

/// %Sensor

class Sensor {
 public:
  /// Constructor
  Sensor() = default;
  /// Destructor
  ~Sensor() {}

  /// Add a component.
  void AddComponent(Component* comp);
  /// Get the number of components attached to the sensor.
  size_t GetNumberOfComponents() const { return m_components.size(); }
  /// Retrieve the pointer to a given component.
  Component* GetComponent(const unsigned int i);
  /// Activate/deactivate a given component.
  void EnableComponent(const unsigned int i, const bool on);
  /// Activate/deactivate use of the magnetic field of a given component.
  void EnableMagneticField(const unsigned int i, const bool on);
  /// Does the sensor have a non-zero magnetic field?
  bool HasMagneticField() const;

  /// Add an electrode.
  void AddElectrode(Component* comp, const std::string& label);
  /// Get the number of electrodes attached to the sensor.
  size_t GetNumberOfElectrodes() const { return m_electrodes.size(); }
  /// Remove all components, electrodes and reset the sensor.
  void Clear();

  /// Get the drift field and potential at (x, y, z).
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& medium,
                     int& status);
  /// Get the drift field at (x, y, z).
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& medium, int& status);

  /// Get the magnetic field at (x, y, z).
  void MagneticField(const double x, const double y, const double z, double& bx,
                     double& by, double& bz, int& status);

  /// Get the weighting field at (x, y, z).
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label);
  /// Get the weighting potential at (x, y, z).
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label);

  /// Get the delayed weighting potential at (x, y, z).
  double DelayedWeightingPotential(const double x, const double y,
                                   const double z, const double t,
                                   const std::string& label);

  /// Get the medium at (x, y, z).
  Medium* GetMedium(const double x, const double y, const double z);

  /// Switch debugging messages on/off.
  void EnableDebugging(const bool on = true) { m_debug = on; }

  /// Set the user area to the default.
  bool SetArea();
  /// Set the user area explicitly.
  bool SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  /// Return the current user area.
  bool GetArea(double& xmin, double& ymin, double& zmin, double& xmax,
               double& ymax, double& zmax);
  /// Check if a point is inside the user area.
  bool IsInArea(const double x, const double y, const double z);

  /// Check if a point is inside an active medium and inside the user area.
  bool IsInside(const double x, const double y, const double z);

  /// Return the voltage range.
  bool GetVoltageRange(double& vmin, double& vmax);

  /// Start a new event, when computing the average signal over multiple events.
  void NewSignal() { ++m_nEvents; }
  /// Reset signals and induced charges of all electrodes.
  void ClearSignal();

  /** Set the time window and binning for the signal calculation.
   * \param tstart start time [ns]
   * \param tstep bin width [ns]
   * \param nsteps number of bins
   */
  void SetTimeWindow(const double tstart, const double tstep,
                     const unsigned int nsteps);
  /// Retrieve the time window and binning.
  void GetTimeWindow(double& tstart, double& tstep,
                     unsigned int& nsteps) const {
    tstart = m_tStart;
    tstep = m_tStep;
    nsteps = m_nTimeBins;
  }

  /// Compute the component of the signal due to the delayed weighting field.
  void EnableDelayedSignal(const bool on = true) { m_delayedSignal = on; }
  /// Set the points in time at which to evaluate the delayed weighting field.
  void SetDelayedSignalTimes(const std::vector<double>& ts);
  /// Set the number of points to be used when averaging the delayed
  /// signal vector over a time bin (default: 0).
  /// The averaging is done with a \f$2\times navg + 1\f$ point
  /// Newton-Raphson integration.
  void SetDelayedSignalAveragingOrder(const unsigned int navg) {
    m_nAvgDelayedSignal = navg;
  }

  /// Set/override the signal in a given time bin explicitly.
  void SetSignal(const std::string& label, const unsigned int bin,
                 const double signal);
  /// Retrieve the total signal for a given electrode and time bin.
  double GetSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the electron signal for a given electrode and time bin.
  double GetSignal(const std::string& label, const unsigned int bin,
                   const int comp);
  /// Retrieve the prompt signal for a given electrode and time bin.
  double GetPromptSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the delayed signal for a given electrode and time bin.
  double GetDelayedSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the electron signal for a given electrode and time bin.
  double GetElectronSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the ion or hole signal for a given electrode and time bin.
  double GetIonSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the delayed electron signal for a given electrode and time bin.
  double GetDelayedElectronSignal(const std::string& label,
                                  const unsigned int bin);
  /// Retrieve the delayed ion/hole signal for a given electrode and time bin.
  double GetDelayedIonSignal(const std::string& label, const unsigned int bin);
  /// Calculated using the weighting potentials at the start and end points.
  double GetInducedCharge(const std::string& label);

  /// Set the function to be used for evaluating the transfer function.
  void SetTransferFunction(std::function<double(double)>);
  /// Set the points to be used for interpolating the transfer function.
  void SetTransferFunction(const std::vector<double>& times,
                           const std::vector<double>& values);
  /// Set the transfer function using a Shaper object.
  void SetTransferFunction(Shaper& shaper);
  /// Evaluate the transfer function at a given time.
  double GetTransferFunction(const double t);
  /// Print some information about the presently set transfer function.
  void PrintTransferFunction();
  /// Cache integral and FFT of the transfer function
  /// instead of recomputing it at every call (default: on).
  void EnableTransferFunctionCache(const bool on = true) {
    m_cacheTransferFunction = on;
  }
  /// Convolute the induced current on a given electrode
  /// with the transfer function.
  bool ConvoluteSignal(const std::string& label, const bool fft = false);
  /// Convolute all induced currents with the transfer function.
  bool ConvoluteSignals(const bool fft = false);
  /// Replace the signal on a given electrode by its integral.
  bool IntegrateSignal(const std::string& label);
  /// Replace the signals on all electrodes curve by their integrals.
  bool IntegrateSignals();
  /// Return whether the signal has been integrated/convoluted.
  bool IsIntegrated(const std::string& label) const;

  /// Delay the signal and subtract an attenuated copy
  /// (modelling a constant fraction discriminator).
  /// \f[
  ///   v_{out} = v_{in}\left(t - t_d\right) - f v_{in}.
  /// \f]
  bool DelayAndSubtractFraction(const double td, const double f);
  /// Set the function to be used for evaluating the noise component.
  void SetNoiseFunction(double (*f)(double t));
  /// Add noise to the induced signal.
  void AddNoise(const bool total = true, const bool electron = false,
                const bool ion = false);
  /** Add white noise to the induced signal, given a desired output ENC.
   * \param label name of the electrode
   * \param enc Equivalent Noise Charge, in electrons.
   * \param poisson flag to sample the number of noise pulses from a
   *                Poisson distribution. Otherwise the noise charge in each
   *                bin is sampled from a Gaussian distribution.
   * \param q0      amplitude of the noise delta pulses (in electrons).
   */
  void AddWhiteNoise(const std::string& label, const double enc,
                     const bool poisson = true, const double q0 = 1.);
  /// Add white noise to the induced signals on all electrodes.
  void AddWhiteNoise(const double enc, const bool poisson = true,
                     const double q0 = 1.);

  /** Determine the threshold crossings of the current signal curve.
   * \param thr threshold value
   * \param label electrode for which to compute the threshold crossings
   * \param n number of threshold crossings
   */
  bool ComputeThresholdCrossings(const double thr, const std::string& label,
                                 int& n);
  /// Get the number of threshold crossings
  /// (after having called ComputeThresholdCrossings).
  size_t GetNumberOfThresholdCrossings() const {
    return m_thresholdCrossings.size();
  }
  /** Retrieve the time and type of a given threshold crossing (after having
   * called ComputeThresholdCrossings.
   * \param i index
   * \param time threshold crossing time [ns]
   * \param level threshold (should correspond to the value requested).
   * \param rise flag whether the crossing is on a rising or falling slope.
   */
  bool GetThresholdCrossing(const unsigned int i, double& time, double& level,
                            bool& rise) const;

  /// Add the signal from a charge-carrier step.
  void AddSignal(const double q, const double t0, const double t1,
                 const double x0, const double y0, const double z0,
                 const double x1, const double y1, const double z1,
                 const bool integrateWeightingField,
                 const bool useWeightingPotential = false);

  /// Add the signal from a drift line.
  void AddSignal(const double q, const std::vector<double>& ts,
                 const std::vector<std::array<double, 3> >& xs,
                 const std::vector<std::array<double, 3> >& vs,
                 const std::vector<double>& ns, const int navg,
                 const bool useWeightingPotential = false);

  /// Exporting induced signal to a csv file.
  void ExportSignal(const std::string& label, const std::string& filename,
                    const bool chargeCariers = false);

  /// Add the induced charge from a charge carrier drift between two points.
  void AddInducedCharge(const double q, const double x0, const double y0,
                        const double z0, const double x1, const double y1,
                        const double z1);

  /// Determine whether a line between two points crosses a wire,
  /// calls Component::CrossedWire.
  bool CrossedWire(const double x0, const double y0, const double z0,
                   const double x1, const double y1, const double z1,
                   double& xc, double& yc, double& zc, const bool centre,
                   double& rc);
  /// Determine whether a point is in the trap radius of a wire.
  bool InTrapRadius(const double q0, const double x0, const double y0,
                    const double z0, double& xw, double& yw, double& rw);
  /// Determine whether a line between two points crosses a plane,
  /// calls Component::CrossedPlane.
  bool CrossedPlane(const double x0, const double y0, const double z0,
                    const double x1, const double y1, const double z1,
                    double& xc, double& yc, double& zc);

  /// Integrate the electric field flux through a line from
  /// (x0,y0,z0) to (x1,y1,z1) along a direction (xp,yp,zp).
  double IntegrateFluxLine(const double x0, const double y0, const double z0,
                           const double x1, const double y1, const double z1,
                           const double xp, const double yp, const double zp,
                           const unsigned int nI, const int isign = 0);
  // TODO!
  double GetTotalInducedCharge(const std::string& label);

 private:
  std::string m_className = "Sensor";
  /// Mutex.
  std::mutex m_mutex;

  /// Components
  std::vector<std::tuple<Component*, bool, bool> > m_components;

  struct Electrode {
    Component* comp;
    std::string label;
    std::vector<double> signal;
    std::vector<double> delayedSignal;
    std::vector<double> electronSignal;
    std::vector<double> ionSignal;
    std::vector<double> delayedElectronSignal;
    std::vector<double> delayedIonSignal;
    double charge;
    bool integrated;
  };
  /// Electrodes
  std::vector<Electrode> m_electrodes;

  // Time window for signals
  double m_tStart = 0.;
  double m_tStep = 10.;
  unsigned int m_nTimeBins = 200;
  unsigned int m_nEvents = 0;
  bool m_delayedSignal = false;
  std::vector<double> m_delayedSignalTimes;
  unsigned int m_nAvgDelayedSignal = 0;

  // Transfer function
  std::function<double(double)> m_fTransfer;
  Shaper* m_shaper = nullptr;
  std::vector<std::pair<double, double> > m_fTransferTab;
  bool m_cacheTransferFunction = true;
  // Integral of the transfer function squared.
  double m_fTransferSq = -1.;
  // FFT of the transfer function.
  std::vector<double> m_fTransferFFT;

  // Noise
  double (*m_fNoise)(double t) = nullptr;

  std::vector<std::pair<double, bool> > m_thresholdCrossings;
  double m_thresholdLevel = 0.;

  // User bounding box
  double m_xMinUser = 0., m_yMinUser = 0., m_zMinUser = 0.;
  double m_xMaxUser = 0., m_yMaxUser = 0., m_zMaxUser = 0.;
  bool m_hasUserArea = false;

  // Switch on/off debugging messages
  bool m_debug = false;

  // Return the current sensor size
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax);

  void FillSignal(Electrode& electrode, const double q,
                  const std::vector<double>& ts, const std::vector<double>& is,
                  const int navg, const bool delayed = false);
  void FillBin(Electrode& electrode, const unsigned int bin,
               const double signal, const bool electron, const bool delayed) {
    std::lock_guard<std::mutex> guard(m_mutex);
    electrode.signal[bin] += signal;
    if (delayed) electrode.delayedSignal[bin] += signal;
    if (electron) {
      electrode.electronSignal[bin] += signal;
      if (delayed) electrode.delayedElectronSignal[bin] += signal;
    } else {
      electrode.ionSignal[bin] += signal;
      if (delayed) electrode.delayedIonSignal[bin] += signal;
    }
  }

  void IntegrateSignal(Electrode& electrode);
  void ConvoluteSignal(Electrode& electrode, const std::vector<double>& tab);
  bool ConvoluteSignalFFT();
  bool ConvoluteSignalFFT(const std::string& label);
  void ConvoluteSignalFFT(Electrode& electrode, const std::vector<double>& tab,
                          const unsigned int nn);
  // Evaluate the integral over the transfer function squared.
  double TransferFunctionSq();
  double InterpolateTransferFunctionTable(const double t) const;
  void MakeTransferFunctionTable(std::vector<double>& tab);
  void FFT(std::vector<double>& data, const bool inverse, const int nn);
};
}  // namespace Garfield

#endif
