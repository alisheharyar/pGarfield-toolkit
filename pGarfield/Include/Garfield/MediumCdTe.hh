#ifndef G_MEDIUM_CDTE_H
#define G_MEDIUM_CDTE_H

#include "Medium.hh"

namespace Garfield {

/// Cadmium-Telluride.

class MediumCdTe : public Medium {
 public:
  /// Constructor
  MediumCdTe();
  /// Destructor
  virtual ~MediumCdTe() {}

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

 private:
  // Band gap energy [eV]
  // m_bandGap = 1.44;
  // Low-field mobility
  double m_eMobility = 1.05e-6;
  double m_hMobility = 0.1e-6;
  // Hall factor
  double m_eHallFactor = 1.15;
  double m_hHallFactor = 0.7;

  bool m_userMobility = false;

  void UpdateTransportParameters();
};
}

#endif
