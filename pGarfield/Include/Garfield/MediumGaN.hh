#ifndef G_MEDIUM_GAN_H
#define G_MEDIUM_GAN_H

#include "Medium.hh"

namespace Garfield {

/// Gallium-Nitride.

class MediumGaN : public Medium {
 public:
  /// Constructor
  MediumGaN();
  /// Destructor
  virtual ~MediumGaN() {}

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

  /// Set the electron concentration [cm-3].
  void SetElectronConcentration(const double c);

  /// Set the low-field mobility values [cm2 / (V ns)] explicitly.
  void SetLowFieldMobility(const double mue, const double muh);
  /// Use the default mobility models. 
  void UnsetLowFieldMobility();

 private:
  // Low-field mobility.
  double m_eMobility = 1.405e-6;
  double m_hMobility = 0.170e-6;

  // Electron concentration.
  double m_eDensity = 7.78e16;

  // Hall factors.
  double m_eHallFactor = 1.;
  double m_hHallFactor = 1.;

  // Impact ionization parameters.
  // J. Baliga, Semicond. Sci. Technol. 28 (2013) 074011,
  // https://doi-org.ezproxy.cern.ch/10.1088/0268-1242/28/7/074011 
  double m_eImpactA = 1.5e5;
  double m_hImpactA = 6.4e5;
  double m_eImpactB = 1.41e7;
  double m_hImpactB = 1.46e7;

  bool m_userMobility = false;
  void UpdateTransportParameters();
};
}

#endif
