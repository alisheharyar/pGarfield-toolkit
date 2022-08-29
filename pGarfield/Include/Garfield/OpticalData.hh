#ifndef G_OPTICAL_DATA_H
#define G_OPTICAL_DATA_H

#include <string>
#include <vector>

namespace Garfield {

/// Photoabsorption cross-sections for some gases.

class OpticalData {
 public:
  /// Constructor
  OpticalData() = default;
  /// Destructor
  ~OpticalData() = default;

  /// Check whether optical data have been implemented for a given gas.
  static bool IsAvailable(const std::string& material);
  /// Photo-absorption cross-section and ionisation yield at a given energy.
  static bool PhotoabsorptionCrossSection(const std::string& material,
                                          const double energy, double& cs, 
                                          double& eta);
  /// Photo-absorption cross-section at a given energy.
  static double PhotoabsorptionCrossSection(const std::string& material,
                                            const double energy);
  /// Photo-ionisation yield at a given energy.
  static double PhotoionisationYield(const std::string& material,
                                     const double energy);

 private:
  static constexpr double OscToPacs = 8.067283e-18;
  static constexpr double Mbarn = 1.e-18;

  static bool PhotoAbsorptionCsNeon(const double e, double& cs, double& eta);
  static bool PhotoAbsorptionCsArgon(const double e, double& cs, double& eta);

  static bool PhotoAbsorptionCsCO2(const double e, double& cs, double& eta);

  static bool PhotoAbsorptionCsMethane(const double e, double& cs, double& eta);
  static bool PhotoAbsorptionCsEthane(const double e, double& cs, double& eta);
  static bool PhotoAbsorptionCsButane(const double e, double& cs, double& eta);
  static bool PhotoAbsorptionCsAcetylene(const double e, double& cs, double& eta);
  static bool PhotoAbsorptionCsCF4(const double e, double& cs, double& eta);

  static bool PhotoAbsorptionCsNitrogen(const double e, double& cs, double& eta);
};
}

#endif
