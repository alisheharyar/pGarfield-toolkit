#ifndef G_TRACK_PAI
#define G_TRACK_PAI

#include <string>
#include <vector>

#include "Track.hh"

namespace Garfield {

/// Energy loss calculation using the Photoabsorption-Ionisation Model.

class TrackPAI : public Track {
 public:
  // Constructor
  TrackPAI();
  // Destructor
  virtual ~TrackPAI() {}

  virtual bool NewTrack(const double x0, const double y0, const double z0,
                        const double t0, const double dx0, const double dy0,
                        const double dz0);

  virtual bool GetCluster(double& xcls, double& ycls, double& zcls,
                          double& tcls, int& ncls, double& ecls, double& extra);

  virtual double GetClusterDensity();
  virtual double GetStoppingPower();

 private:
  bool m_ready = false;

  // Particle coordinates and direction
  double m_x = 0., m_y = 0., m_z = 0., m_t = 0.;
  double m_dx = 0., m_dy = 0., m_dz = 0.;
  // Particle energy and speed
  double m_e = 0.;
  double m_speed = 0.;
  // Max. energy transfer in a collision
  double m_emax = 0.;

  // Total inelastic mean free path
  double m_imfp = 0.;
  // Stopping power
  double m_dedx = 0.;

  // Dielectric function
  int m_nSteps = 1000;
  struct opticalData {
    double eps1, eps2;
    double integral;
  };
  std::vector<opticalData> m_opticalDataTable;

  // Tables for interpolation of cumulative distribution functions
  std::vector<double> m_energies;
  std::vector<double> m_cdf;
  std::vector<double> m_rutherford;

  struct electron {
    // Direction
    double dx, dy, dz;
    // Energy
    double energy;
    // Type (electron, hole)
    int type;
  };
  std::vector<electron> m_electrons;
  std::vector<electron> m_holes;

  // Medium properties
  std::string m_mediumName = "";
  double m_mediumDensity = 0.;
  double m_electronDensity = 0.;

  bool SetupMedium(Medium* medium);
  bool SetupCrossSectionTable();

  double ComputeMaxTransfer() const;

  double ComputeCsTail(const double emin, const double emax);
  double ComputeDeDxTail(const double emin, const double emax);

  double SampleEnergyDeposit(const double u, double& f) const;
  double SampleAsymptoticCs(double u) const;
  double SampleAsymptoticCsSpinZero(const double emin, double u) const;
  double SampleAsymptoticCsSpinHalf(const double emin, double u) const;
  double SampleAsymptoticCsSpinOne(const double emin, double u) const;
  double SampleAsymptoticCsElectron(const double emin, double u) const;
  double SampleAsymptoticCsPositron(const double emin, double u) const;

  double LossFunction(const double eps1, const double eps2) const {
    const double eps = eps1 * eps1 + eps2 * eps2;
    return eps > 0. ? eps2 / eps : 0.;
  }
};
}

#endif
