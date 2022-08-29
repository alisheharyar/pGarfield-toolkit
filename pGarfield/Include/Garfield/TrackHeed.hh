#ifndef G_TRACK_HEED_H
#define G_TRACK_HEED_H

#include <list>
#include <memory>
#include <vector>

#include "Track.hh"

namespace Heed {
class gparticle;
class particle_def;
class HeedParticle;
class HeedCondElectron;
class HeedMatterDef;
class GasDef;
class MatterDef;
class EnergyMesh;
class EnTransfCS;
class ElElasticScat;
class ElElasticScatLowSigma;
class PairProd;
class HeedDeltaElectronCS;
class HeedFieldMap;
}

namespace Garfield {

class HeedChamber;
class Medium;

/// Generate tracks using Heed++.

class TrackHeed : public Track {
 public:
  /// Constructor
  TrackHeed();
  /// Destructor
  virtual ~TrackHeed();

  bool NewTrack(const double x0, const double y0, const double z0,
                const double t0, const double dx0, const double dy0,
                const double dz0) override;
  bool GetCluster(double& xcls, double& ycls, double& zcls, double& tcls,
                  int& n, double& e, double& extra) override;
  bool GetCluster(double& xcls, double& ycls, double& zcls, double& tcls,
                  int& ne, int& ni, double& e, double& extra);
  /** Get the next "cluster" (ionising collision of the charged particle).
    * \param xcls,ycls,zcls coordinates of the collision
    * \param tcls time of the collision
    * \param ne number of electrons
    * \param ni number of ions
    * \param np number of fluorescence photons
    * \param e deposited energy
    * \param extra additional information (not always implemented)
    */
  bool GetCluster(double& xcls, double& ycls, double& zcls, double& tcls,
                  int& ne, int& ni, int& np, double& e, double& extra);
  /** Retrieve the properties of a conduction or delta electron
    * in the current cluster.
    * \param i index of the electron
    * \param x,y,z coordinates of the electron
    * \param t time
    * \param e kinetic energy (only meaningful for delta-electrons)
    * \param dx,dy,dz direction vector (only meaningful for delta-electrons)
    **/
  bool GetElectron(const unsigned int i, double& x, double& y, double& z,
                   double& t, double& e, double& dx, double& dy, double& dz);
  /** Retrieve the properties of an ion in the current cluster.
    * \param i index of the ion
    * \param x,y,z coordinates of the ion
    * \param t time
    **/
  bool GetIon(const unsigned int i, double& x, double& y, double& z,
              double& t) const;

  /** Retrieve the properties of an unabsorbed photon.
   * \param i index of the photon
   * \param x,y,z coordinates of the photon
   * \param t time
   * \param e photon energy
   * \param dx,dy,dz direction vector
   **/
  bool GetPhoton(const unsigned int i, double& x, double& y, double& z,
                 double& t, double& e, double& dx, double& dy,
                 double& dz) const;

  double GetClusterDensity() override;
  double GetStoppingPower() override;
  /// Return the W value of the medium (of the last simulated track).
  double GetW() const;
  /// Return the Fano factor of the medium (of the last simulated track).
  double GetFanoFactor() const;
  /// Return the photoabsorption cross-section at a given energy.
  double GetPhotoAbsorptionCrossSection(const double e) const;

  /// Compute the differential cross-section for a given medium.
  bool Initialise(Medium* medium, const bool verbose = false);

  /** Simulate a delta electron.
    * \param x0,y0,z0 initial position of the delta electron
    * \param t0 initial time
    * \param e0 initial kinetic energy of the delta electron
    * \param dx0,dy0,dz0 initial direction of the delta electron
    * \param ne,ni number of electrons/ions produced by the delta electron
    **/
  void TransportDeltaElectron(const double x0, const double y0, const double z0,
                              const double t0, const double e0,
                              const double dx0, const double dy0,
                              const double dz0, int& ne, int& ni);
  /** Simulate a delta electron.
    * \param x0,y0,z0 initial position of the delta electron
    * \param t0 initial time
    * \param e0 initial kinetic energy of the delta electron
    * \param dx0,dy0,dz0 initial direction of the delta electron
    * \param ne number of electrons produced by the delta electron
    **/
  void TransportDeltaElectron(const double x0, const double y0, const double z0,
                              const double t0, const double e0,
                              const double dx0, const double dy0,
                              const double dz0, int& ne);

  /** Simulate a photon.
    * \param x0,y0,z0 initial position of the photon
    * \param t0 initial time
    * \param e0 initial energy of the photon
    * \param dx0,dy0,dz0 initial direction of the photon
    * \param ne number of electrons produced by the photon
    * \param ni number of ions produced by the photon
    * \param np number of fluorescence photons
    **/
  void TransportPhoton(const double x0, const double y0, const double z0,
                       const double t0, const double e0, const double dx0,
                       const double dy0, const double dz0, 
                       int& ne, int& ni, int& np);

  /** Simulate a photon.
    * \param x0,y0,z0 initial position of the photon
    * \param t0 initial time
    * \param e0 initial energy of the photon
    * \param dx0,dy0,dz0 initial direction of the photon
    * \param ne number of electrons produced by the photon
    * \param ni number of ions produced by the photon
    **/
  void TransportPhoton(const double x0, const double y0, const double z0,
                       const double t0, const double e0, const double dx0,
                       const double dy0, const double dz0, int& ne, int& ni);
  /** Simulate a photon.
    * \param x0,y0,z0 initial position of the photon
    * \param t0 initial time
    * \param e0 initial energy of the photon
    * \param dx0,dy0,dz0 initial direction of the photon
    * \param ne number of electrons produced by the photon
    **/
  void TransportPhoton(const double x0, const double y0, const double z0,
                       const double t0, const double e0, const double dx0,
                       const double dy0, const double dz0, int& ne);

  /// Take the electric field into account in the stepping algorithm.
  void EnableElectricField();
  /// Do not take the electric field into account in the stepping algorithm.
  void DisableElectricField();
  /// Take the magnetic field into account in the stepping algorithm.
  void EnableMagneticField();
  /// Do not take the magnetic field into account in the stepping algorithm.
  void DisableMagneticField();

  /** Set parameters for calculating the particle trajectory.
    * \param maxStep
    *        maximum step length
    * \param radStraight
    *        radius beyond which to approximate circles by polylines.
    * \param stepAngleStraight
    *        max. angular step (in radian) when using polyline steps.
    * \param stepAngleCurved
    *        max. angular step (in radian) when using circular steps.
    **/
  void SetSteppingLimits(const double maxStep, const double radStraight,
                         const double stepAngleStraight,
                         const double stepAngleCurved) {
    m_maxStep = maxStep;
    m_radStraight = radStraight;
    m_stepAngleStraight = stepAngleStraight;
    m_stepAngleCurved = stepAngleCurved;
  }
  void GetSteppingLimits(double& maxStep, double& radStraight,
                         double& stepAngleStraight, double& stepAngleCurved) {
    maxStep = m_maxStep;
    radStraight = m_radStraight;
    stepAngleStraight = m_stepAngleStraight;
    stepAngleCurved = m_stepAngleCurved;
  }

  void EnableCoulombScattering(const bool on = true) { 
    m_coulombScattering = on;
  } 
  /// Switch simulation of delta electrons on.
  void EnableDeltaElectronTransport() { m_doDeltaTransport = true; }
  /// Switch simulation of delta electrons off.
  void DisableDeltaElectronTransport() { m_doDeltaTransport = false; }

  /// Simulate (or not) the photons produced in the atomic relaxation cascade.
  void EnablePhotonReabsorption(const bool on = true) {
    m_doPhotonReabsorption = on;
  }

  /// Write the photoabsorption cross-sections used to a text file.
  void EnablePhotoAbsorptionCrossSectionOutput(const bool on) {
    m_usePacsOutput = on;
  }
  /** Specify the energy mesh to be used.
    * \param e0,e1 lower/higher limit of the energy range [eV]
    * \param nsteps number of intervals
    **/
  void SetEnergyMesh(const double e0, const double e1, const int nsteps);

  /// Define particle mass and charge (for exotic particles).
  /// For standard particles Track::SetParticle should be used.
  void SetParticleUser(const double m, const double z);

  void EnableOneStepFly(const bool on) { m_oneStepFly = on; }
 private:
  // Prevent usage of copy constructor and assignment operator
  TrackHeed(const TrackHeed& heed);
  TrackHeed& operator=(const TrackHeed& heed);

  bool m_oneStepFly = false;

  bool m_ready = false;
  bool m_hasActiveTrack = false;

  double m_mediumDensity = -1.;
  std::string m_mediumName = "";

  bool m_usePacsOutput = false;

  bool m_doPhotonReabsorption = true;
  struct SimplifiedParticle {
    double x, y, z, t;
    double e;
    double dx, dy, dz;
  };
  typedef SimplifiedParticle Photon;
  std::vector<Photon> m_photons;

  bool m_coulombScattering = false;

  bool m_doDeltaTransport = true;
  typedef SimplifiedParticle DeltaElectron;
  std::vector<DeltaElectron> m_deltaElectrons;
  std::vector<Heed::HeedCondElectron> m_conductionElectrons;
  std::vector<Heed::HeedCondElectron> m_conductionIons;

  // Particle properties
  std::unique_ptr<Heed::particle_def> m_particle_def; 
  // Material properties
  std::unique_ptr<Heed::HeedMatterDef> m_matter;
  std::unique_ptr<Heed::GasDef> m_gas;
  std::unique_ptr<Heed::MatterDef> m_material;

  // Energy mesh
  double m_emin = 2.e-6;
  double m_emax = 2.e-1;
  unsigned int m_nEnergyIntervals = 200;
  std::unique_ptr<Heed::EnergyMesh> m_energyMesh;

  // Cross-sections
  std::unique_ptr<Heed::EnTransfCS> m_transferCs;
  std::unique_ptr<Heed::ElElasticScat> m_elScat;
  std::unique_ptr<Heed::ElElasticScatLowSigma> m_lowSigma;
  std::unique_ptr<Heed::PairProd> m_pairProd;
  std::unique_ptr<Heed::HeedDeltaElectronCS> m_deltaCs;

  // Interface classes
  std::unique_ptr<HeedChamber> m_chamber;
  std::unique_ptr<Heed::HeedFieldMap> m_fieldMap;

  // Bounding box
  double m_lX = 0., m_lY = 0., m_lZ = 0.;
  double m_cX = 0., m_cY = 0., m_cZ = 0.;

  // Stepping parameters.
  /// Max. step length.
  double m_maxStep = 100.;
  /// Bending radius beyond which to use straight-line approximation.
  double m_radStraight = 1000.;
  /// Angular step for curved trajectories approximated by straight-line steps.
  double m_stepAngleStraight = 0.1;
  /// Angular step for curved lines.
  double m_stepAngleCurved = 0.2;

  std::vector<Heed::gparticle*> m_particleBank;
  std::vector<Heed::gparticle*>::iterator m_bankIterator;

  bool SetupGas(Medium* medium);
  bool SetupMaterial(Medium* medium);
  bool SetupDelta(const std::string& databasePath);
  void ClearParticleBank();
  bool IsInside(const double x, const double y, const double z);
  bool UpdateBoundingBox(bool& update);
};
}

#endif
