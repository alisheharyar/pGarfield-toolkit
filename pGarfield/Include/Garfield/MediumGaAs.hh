#ifndef G_MEDIUM_GAAS_H
#define G_MEDIUM_GAAS_H

#include "Medium.hh"

namespace Garfield {

/// Gallium-Arsenide.

class MediumGaAs : public Medium {
 public:
  /// Constructor
  MediumGaAs();
  /// Destructor
  virtual ~MediumGaAs() {}

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
  // Band-gap energy [eV]
  // double m_bandGap = 1.42;
  // Low-field mobility
  double m_eMobility = 8.0e-6;
  double m_hMobility = 0.4e-6;
  // Saturation velocity
  double m_eSatVel = 7.7e-3;
  double m_hSatVel = 7.7e-3;
  // Hall factor
  double m_eHallFactor = 1.05;
  double m_hHallFactor = 1.25;
  // Impact ionization parameters
  double m_eImpactA = 1.889e5;
  double m_hImpactA = 2.215e5;
  double m_eImpactB = 5.75e5;
  double m_hImpactB = 6.57e5;

  bool m_userMobility = false;
  void UpdateTransportParameters();
};
}

#endif
