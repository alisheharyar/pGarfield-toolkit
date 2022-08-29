#ifndef G_AVALANCHE_MC_H
#define G_AVALANCHE_MC_H

#include <array>
#include <string>
#include <vector>

#include "FundamentalConstants.hh"
#include "GarfieldConstants.hh"
#include "Sensor.hh"
#include "ViewDrift.hh"

namespace Garfield {

/// Calculate drift lines and avalanches based on macroscopic transport
/// coefficients, using Monte Carlo integration.

class AvalancheMC {
 public:
  /// Constructor
  AvalancheMC();
  /// Destructor
  ~AvalancheMC() {}

  /// Set the sensor.
  void SetSensor(Sensor* s);

  /// Switch on drift line plotting.
  void EnablePlotting(ViewDrift* view);
  /// Switch off drift line plotting.
  void DisablePlotting() { m_viewer = nullptr; }

  /// Switch on calculation of induced currents (default: disabled).
  void EnableSignalCalculation(const bool on = true) { m_doSignal = on; }
  /// Set the number of points to be used when averaging the
  /// signal vector over a time bin in the Sensor class.
  /// The averaging is done with a \f$2\times navg + 1\f$ point
  /// Newton-Raphson integration. Default: 1.
  void SetSignalAveragingOrder(const unsigned int navg) { m_navg = navg; }
  /// Use the weighting potential (as opposed to the weighting field)
  /// for calculating the induced signal.
  void UseWeightingPotential(const bool on = true) {
    m_useWeightingPotential = on;
  }

  /// Switch on calculation of induced charge (default: disabled).
  void EnableInducedChargeCalculation(const bool on = true) {
    m_doInducedCharge = on;
  }

  /** Switch on Runge-Kutta-Fehlberg stepping (as opposed to simple
   * straight-line steps. */
  void EnableRKFSteps(const bool on = true) { m_doRKF = on; }

  /** Switch on equilibration of multiplication and attachment
   * over the drift line (default: enabled). */
  void EnableProjectedPathIntegration(const bool on = true) {
    m_doEquilibration = on;
  }

  /// Switch on diffusion (default: enabled).
  void EnableDiffusion() { m_useDiffusion = true; }
  /// Switch off diffusion.
  void DisableDiffusion() { m_useDiffusion = false; }

  /** Switch on attachment (and multiplication) for drift line calculation
   * (default: enabled). For avalanches the flag is ignored. */
  void EnableAttachment() { m_useAttachment = true; }
  /// Switch off attachment and multiplication.
  void DisableAttachment() { m_useAttachment = false; }

  /// Retrieve the Townsend coefficient from the component.
  void EnableTownsendMap(const bool on = true) { m_useTownsendMap = on; }
  /// Retrieve the attachment coefficient from the component.
  void EnableAttachmentMap(const bool on = true) { m_useAttachmentMap = on; }
  /// Retrieve the drift velocity from the component.
  void EnableVelocityMap(const bool on = true) { m_useVelocityMap = on; }

  /** Set a maximum avalanche size (ignore further multiplication
      once this size has been reached). */
  void EnableAvalancheSizeLimit(const unsigned int size) { m_sizeCut = size; }
  /// Do not limit the maximum avalanche size.
  void DisableAvalancheSizeLimit() { m_sizeCut = 0; }
  /// Return the currently set avalanche size limit.
  int GetAvalancheSizeLimit() const { return m_sizeCut; }

  /// Use fixed-time steps (default 20 ps).
  void SetTimeSteps(const double d = 0.02);
  /// Use fixed distance steps (default 10 um).
  void SetDistanceSteps(const double d = 0.001);
  /** Use exponentially distributed time steps with mean equal
   * to the specified multiple of the collision time (default model).*/
  void SetCollisionSteps(const unsigned int n = 100);
  /// Retrieve the step distance from a user-supplied function.
  void SetStepDistanceFunction(double (*f)(double x, double y, double z));

  /// Define a time interval (only carriers inside the interval are drifted).
  void SetTimeWindow(const double t0, const double t1);
  /// Do not limit the time interval within which carriers are drifted.
  void UnsetTimeWindow() { m_hasTimeWindow = false; }

  /// Set multiplication factor for the signal induced by electrons.
  void SetElectronSignalScalingFactor(const double scale) { m_scaleE = scale; }
  /// Set multiplication factor for the signal induced by holes.
  void SetHoleSignalScalingFactor(const double scale) { m_scaleH = scale; }
  /// Set multiplication factor for the signal induced by ions.
  void SetIonSignalScalingFactor(const double scale) { m_scaleI = scale; }

  /// Return the number of electrons and ions/holes in the avalanche.
  void GetAvalancheSize(unsigned int& ne, unsigned int& ni) const {
    ne = m_nElectrons;
    ni = std::max(m_nIons, m_nHoles);
  }

  /// Return the number of points along the last simulated drift line.
  size_t GetNumberOfDriftLinePoints() const { return m_drift.size(); }
  /// Return the coordinates and time of a point along the last drift line.
  void GetDriftLinePoint(const size_t i, double& x, double& y, double& z,
                         double& t) const;

  /** Return the number of electron trajectories in the last
   * simulated avalanche (including captured electrons). */
  size_t GetNumberOfElectronEndpoints() const {
    return m_endpointsElectrons.size();
  }
  /** Return the number of hole trajectories in the last
   * simulated avalanche (including captured holes). */
  size_t GetNumberOfHoleEndpoints() const { return m_endpointsHoles.size(); }
  /// Return the number of ion trajectories.
  size_t GetNumberOfIonEndpoints() const { return m_endpointsIons.size(); }

  /** Return the coordinates and time of start and end point of a given
   * electron drift line.
   * \param i index of the drift line
   * \param x0,y0,z0,t0 coordinates and time of the starting point
   * \param x1,y1,z1,t1 coordinates and time of the end point
   * \param status status code (see GarfieldConstants.hh)
   */
  void GetElectronEndpoint(const size_t i, double& x0, double& y0,
                           double& z0, double& t0, double& x1, double& y1,
                           double& z1, double& t1, int& status) const;
  void GetHoleEndpoint(const size_t i, double& x0, double& y0, double& z0,
                       double& t0, double& x1, double& y1, double& z1,
                       double& t1, int& status) const;
  void GetIonEndpoint(const size_t i, double& x0, double& y0, double& z0,
                      double& t0, double& x1, double& y1, double& z1,
                      double& t1, int& status) const;

  /// Simulate the drift line of an electron from a given starting point.
  bool DriftElectron(const double x, const double y, const double z,
                     const double t);
  /// Simulate the drift line of a hole from a given starting point.
  bool DriftHole(const double x, const double y, const double z,
                 const double t);
  /// Simulate the drift line of an ion from a given starting point.
  bool DriftIon(const double x, const double y, const double z,
                const double t);
  /// Simulate the drift line of a negative ion from a given starting point.
  bool DriftNegativeIon(const double x, const double y, const double z,
                        const double t);
  /** Simulate an avalanche initiated by an electron at a given starting point.
   * \param x,y,z,t coordinates and time of the initial electron
   * \param hole simulate the hole component of the avalanche or not
   */
  bool AvalancheElectron(const double x, const double y, const double z,
                         const double t, const bool hole = false);
  /// Simulate an avalanche initiated by a hole at a given starting point.
  bool AvalancheHole(const double x, const double y, const double z,
                     const double t, const bool electron = false);
  /// Simulate an avalanche initiated by an electron-hole pair.
  bool AvalancheElectronHole(const double x, const double y, const double z,
                             const double t);

  /// Add an electron to the list of particles to be transported.
  void AddElectron(const double x, const double y, const double z,
                   const double t);
  /// Add a hole to the list of particles to be transported.
  void AddHole(const double x, const double y, const double z, const double t);
  /// Add an ion to the list of particles to be transported.
  void AddIon(const double x, const double y, const double z, const double t);
  /// Resume the simulation from the current set of charge carriers.
  bool ResumeAvalanche(const bool electron = true, const bool hole = true);

  /// Switch debugging messages on/off (default: off).
  void EnableDebugging(const bool on = true) { m_debug = on; }

 private:
  std::string m_className = "AvalancheMC";

  Sensor* m_sensor = nullptr;

  struct DriftPoint {
    std::array<double, 3> x;  ///< Position.
    double t;                 ///< Time.
    Particle particle;        ///< Charge carrier type.
    unsigned int n;           ///< Number of charge carriers.
  };
  /// Current drift line
  std::vector<DriftPoint> m_drift;

  enum class StepModel {
    FixedTime,
    FixedDistance,
    CollisionTime,
    UserDistance
  };
  /// Step size model.
  StepModel m_stepModel = StepModel::CollisionTime;

  /// Fixed time step
  double m_tMc = 0.02;
  /// Fixed distance step
  double m_dMc = 0.001;
  /// Sample step size according to collision time
  int m_nMc = 100;
  /// User function returning the step size
  double (*m_fStep)(double x, double y, double z) = nullptr;

  /// Flag whether a time window should be used.
  bool m_hasTimeWindow = false;
  /// Lower limit of the time window.
  double m_tMin = 0.;
  /// Upper limit of the time window.
  double m_tMax = 0.;

  /// Max. avalanche size.
  unsigned int m_sizeCut = 0;

  /// Number of electrons produced
  unsigned int m_nElectrons = 0;
  /// Number of holes produced
  unsigned int m_nHoles = 0;
  /// Number of ions produced
  unsigned int m_nIons = 0;

  struct EndPoint {
    std::array<double, 3> x0;  ///< Starting point.
    std::array<double, 3> x1;  ///< End point.
    double t0, t1;             ///< Start and end time.
    int status;                ///< Status flag at the end point.
  };
  /// Endpoints of all electrons in the avalanche (including captured ones)
  std::vector<EndPoint> m_endpointsElectrons;
  /// Endpoints of all holes in the avalanche (including captured ones)
  std::vector<EndPoint> m_endpointsHoles;
  /// Endpoints of all ions in the avalanche
  std::vector<EndPoint> m_endpointsIons;

  ViewDrift* m_viewer = nullptr;

  bool m_doSignal = false;
  unsigned int m_navg = 1;
  bool m_useWeightingPotential = true;
  bool m_doInducedCharge = false;
  bool m_doEquilibration = true;
  bool m_doRKF = false;
  bool m_useDiffusion = true;
  bool m_useAttachment = true;
  /// Scaling factor for electron signals.
  double m_scaleE = 1.;
  /// Scaling factor for hole signals.
  double m_scaleH = 1.;
  /// Scaling factor for ion signals.
  double m_scaleI = 1.;

  /// Take Townsend coefficients from the component.
  bool m_useTownsendMap = false;
  /// Take attachment coefficients from the component.
  bool m_useAttachmentMap = false;
  /// Take the drift velocities from the component.
  bool m_useVelocityMap = false;

  bool m_debug = false;

  /// Compute a drift line with starting point x0.
  bool DriftLine(const std::array<double, 3>& x0, const double t0,
                 const Particle particle,
                 std::vector<DriftPoint>& secondaries,
                 const bool aval = false);
  /// Compute an avalanche.
  bool Avalanche(std::vector<DriftPoint>& aval,
                 const bool withElectrons, const bool withHoles);

  void AddPoint(const std::array<double, 3>& x, const double t,
                const Particle particle, const unsigned int n, 
                std::vector<DriftPoint>& points) {
    DriftPoint point;
    point.x = x;
    point.t = t;
    point.particle = particle;
    point.n = n;
    points.push_back(std::move(point));
  }
  void AddEndPoint(const std::array<double, 3>& x0, const double t0,
                   const std::array<double, 3>& x1, const double t1,
                   const int status, const Particle particle) { 
    EndPoint endPoint;
    endPoint.x0 = x0;
    endPoint.t0 = t0;
    endPoint.x1 = x1;
    endPoint.t1 = t1;
    endPoint.status = status;
    if (particle == Particle::Electron) {
      m_endpointsElectrons.push_back(std::move(endPoint));
    } else if (particle == Particle::Hole) {
      m_endpointsHoles.push_back(std::move(endPoint));
    } else if (particle == Particle::Ion) {
      m_endpointsIons.push_back(std::move(endPoint));
    }
  }

  /// Compute electric and magnetic field at a given position.
  int GetField(const std::array<double, 3>& x, std::array<double, 3>& e,
               std::array<double, 3>& b, Medium*& medium) const;
  /// Retrieve the low-field mobility.
  double GetMobility(const Particle particle, Medium* medium) const;
  /// Compute the drift velocity.
  bool GetVelocity(const Particle particle, Medium* medium,
                   const std::array<double, 3>& x,
                   const std::array<double, 3>& e,
                   const std::array<double, 3>& b,
                   std::array<double, 3>& v) const;
  /// Compute the diffusion coefficients.
  bool GetDiffusion(const Particle particle, Medium* medium,
                    const std::array<double, 3>& e,
                    const std::array<double, 3>& b, double& dl,
                    double& dt) const;
  /// Compute the attachment coefficient.
  double GetAttachment(const Particle particle, Medium* medium,
                       const std::array<double, 3>& x,
                       const std::array<double, 3>& e,
                       const std::array<double, 3>& b) const;
  /// Compute the Townsend coefficient.
  double GetTownsend(const Particle particle, Medium* medium,
                     const std::array<double, 3>& x,
                     const std::array<double, 3>& e,
                     const std::array<double, 3>& b) const;
  /// Compute end point and effective velocity for a step.
  void StepRKF(const Particle particle, const std::array<double, 3>& x0,
               const std::array<double, 3>& v0, const double dt,
               std::array<double, 3>& xf, std::array<double, 3>& vf,
               int& status) const;
  /// Add a diffusion step.
  void AddDiffusion(const double step, const double dl, const double dt,
                    std::array<double, 3>& x,
                    const std::array<double, 3>& v) const;
  /// Terminate a drift line close to the boundary.
  void Terminate(const std::array<double, 3>& x0, const double t0,
                 std::array<double, 3>& x, double& t) const;
  /// Compute multiplication and losses along the current drift line.
  bool ComputeGainLoss(const Particle particle,
                       std::vector<DriftPoint>& driftLine, int& status,
                       std::vector<DriftPoint>& secondaries,
                       const bool semiconductor = false);
  /// Compute Townsend and attachment coefficients along the current drift line.
  bool ComputeAlphaEta(const Particle particle,
                       std::vector<DriftPoint>& driftLine,
                       std::vector<double>& alphas,
                       std::vector<double>& etas) const;
  bool Equilibrate(std::vector<double>& alphas) const;
  /// Compute the induced signal for the current drift line.
  void ComputeSignal(const Particle particle, const double q,
                     const std::vector<DriftPoint>& driftLine) const;
  /// Compute the induced charge for the current drift line.
  void ComputeInducedCharge(const double q,
                            const std::vector<DriftPoint>& driftLine) const;
  void PrintError(const std::string& fcn, const std::string& par,
                  const Particle particle,
                  const std::array<double, 3>& x) const;
};
}  // namespace Garfield

#endif
