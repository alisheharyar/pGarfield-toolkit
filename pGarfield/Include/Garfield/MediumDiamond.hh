#ifndef G_MEDIUM_DIAMOND_H
#define G_MEDIUM_DIAMOND_H

#include <mutex>

#include "Medium.hh"

namespace Garfield {

/// Diamond.

class MediumDiamond : public Medium {
 public:
  /// Constructor
  MediumDiamond();
  /// Destructor
  virtual ~MediumDiamond() {}

  bool IsSemiconductor() const override { return true; }

  void GetComponent(const unsigned int i, std::string& label, 
                    double& f) override;

  // Electron transport parameters
  bool ElectronVelocity(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& vx, double& vy, double& vz) override;
  bool ElectronTownsend(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& alpha) override;
  bool ElectronAttachment(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          double& eta) override;
  double ElectronMobility() override { return m_eMobility; }
  // Hole transport parameters
  bool HoleVelocity(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    double& vx, double& vy, double& vz) override;
  bool HoleTownsend(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    double& alpha) override;
  bool HoleAttachment(const double ex, const double ey, const double ez,
                      const double bx, const double by, const double bz,
                      double& eta) override;
  double HoleMobility() override { return m_hMobility; }

  void SetLowFieldMobility(const double mue, const double muh);
  void UnsetLowFieldMobility();

  void SetSaturationVelocity(const double vsate, const double vsath);
  void UnsetSaturationVelocity();
 private:
  std::mutex m_mutex;
 
  // Low-field mobility
  double m_eMobility = 4.551e-6;
  double m_hMobility = 2.750e-6;
  // Hall factor
  double m_eHallFactor = 1.;
  double m_hHallFactor = 1.;
  // Saturation velocity
  double m_eSatVel = 2.6e-2;
  double m_hSatVel = 1.6e-2;

  bool m_userMobility = false;

  void UpdateTransportParameters();
};
}

#endif
