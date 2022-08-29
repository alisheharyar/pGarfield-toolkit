#ifndef G_COMPONENT_USER_H
#define G_COMPONENT_USER_H

#include <functional>

#include "Component.hh"

namespace Garfield {

/// Simple component with electric field given by a user function.

class ComponentUser : public Component {
 public:
  /// Constructor
  ComponentUser();
  /// Destructor
  ~ComponentUser() {}

  Medium* GetMedium(const double x, const double y, const double z) override {
    return !m_hasArea ? Component::GetMedium(x, y, z) : 
                        InArea(x, y, z) ? m_medium : nullptr;
  }
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;
  void MagneticField(const double x, const double y, const double z, double& bx,
                     double& by, double& bz, int& status) override;
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;
  void DelayedWeightingField(const double x, const double y, const double z,
                             const double t, double& wx, double& wy, double& wz,
                             const std::string& label) override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                      double& xmax, double& ymax, double& zmax) override;

  bool HasMagneticField() const override;

  /// Set the function to be called for calculating the electric field.
  void SetElectricField(
    std::function<void(const double, const double, const double,
                       double&, double&, double&)>);
  /// Set the function to be called for calculating the potential.
  void SetPotential(
    std::function<void(const double, const double, const double, double&)>);
  /// Set the function to be called for calculating the weighting field.
  void SetWeightingField(
    std::function<void(const double, const double, const double,
                       double&, double&, double&, const std::string&)>);
  /// Set the function to be called for calculating the weighting potential.
  void SetWeightingPotential(
    std::function<void(const double, const double, const double, 
                       double&, const std::string&)>);
  /// Set the function to be called for calculating the delayed weighting field.
  void SetDelayedWeightingField(
    std::function<void(const double, const double, const double, const double,
                       double&, double&, double&, const std::string&)>);
  /// Set the function to be called for calculating the magnetic field.
  void SetMagneticField(
    std::function<void(const double, const double, const double,
                       double&, double&, double&)>);

  /// Set the limits of the active area explicitly 
  /// (instead of using a Geometry object).
  void SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  /// Remove the explicit limits of the active area. 
  void UnsetArea();
  /// Set the medium in the active area.
  void SetMedium(Medium* medium) { m_medium = medium; }
 private:
  /// Electric field function
  std::function<void(const double, const double, const double,
                     double&, double&, double&)> m_efield;
  /// Potential function
  std::function<void(const double, const double, const double,
                     double&)> m_potential;

  /// Weighting field function
  std::function<void(const double, const double, const double, 
                     double&, double&, double&, const std::string&)> m_wfield;

  /// Weighting potential function
  std::function<void(const double, const double, const double, 
                     double&, const std::string&)> m_wpot;

  /// Delayed weighting field function
  std::function<void(const double, const double, const double, const double,
                     double&, double&, double&, const std::string&)> m_dwfield;

  /// Magnetic field function
  std::function<void(const double, const double, const double, 
                     double&, double&, double&)> m_bfield;

  // Active area.
  std::array<double, 3> m_xmin = {{0., 0., 0.}};
  std::array<double, 3> m_xmax = {{0., 0., 0.}};
  // Did we specify the active area explicitly?
  bool m_hasArea = false;
  // Medium in the active area.
  Medium* m_medium = nullptr;

  /// Reset the component
  void Reset() override;
  // Verify periodicities
  void UpdatePeriodicity() override;

  bool InArea(const double x, const double y, const double z) {
    return (x >= m_xmin[0] && x <= m_xmax[0] &&
            y >= m_xmin[1] && y <= m_xmax[1] && 
            z >= m_xmin[2] && z <= m_xmax[2]);
  }

};
}
#endif
