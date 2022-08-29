#ifndef G_COMPONENT_CONSTANT_H
#define G_COMPONENT_CONSTANT_H

#include <array>
#include <string>

#include "Component.hh"
#include "Medium.hh"

namespace Garfield {

/// Component with constant electric field.

class ComponentConstant : public Component {
 public:
  /// Constructor
  ComponentConstant();
  /// Destructor
  ~ComponentConstant() {}

  Medium* GetMedium(const double x, const double y, const double z) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                      double& xmax, double& ymax, double& zmax) override;

  /// Set the components of the electric field [V / cm].
  void SetElectricField(const double ex, const double ey, const double ez);
  /// Specify the potential at a given point.
  void SetPotential(const double x, const double y, const double z,
                    const double v = 0.);

  /// Set the components of the weighting field [1 / cm].
  void SetWeightingField(const double wx, const double wy, const double wz,
                         const std::string label);
  /// Specify the weighting potential at a given point.
  void SetWeightingPotential(const double x, const double y, const double z,
                             const double v = 0.);

  /// Set the limits of the active area explicitly 
  /// (instead of using a Geometry object).
  void SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  /// Remove the explicit limits of the active area. 
  void UnsetArea();
  /// Set the medium in the active area.
  void SetMedium(Medium* medium) { m_medium = medium; }

 private:
  // Electric field.
  std::array<double, 3> m_efield = {{0., 0., 0.}};

  // Is the potential defined?
  bool m_hasPotential = false;
  // Point where potential was specified.
  double m_x0 = 0., m_y0 = 0., m_z0 = 0.;
  // Potential at this point.
  double m_v0 = 0.;

  // Is the weighting field defined?
  bool m_hasWeightingField = false;
  // Identifier of the weighting field.
  std::string m_label = "";
  // Weighting field.
  std::array<double, 3> m_wfield = {{0., 0., 0.}};
  // Is the weighting potential defined?
  bool m_hasWeightingPotential = false;
  // Point where the weighting potential was specified.
  double m_wx0 = 0., m_wy0 = 0., m_wz0 = 0.;
  // Weighting potential at this point.
  double m_w0 = 0.;

  // Active area.
  std::array<double, 3> m_xmin = {{0., 0., 0.}};
  std::array<double, 3> m_xmax = {{0., 0., 0.}};
  // Did we specify the active area explicitly?
  bool m_hasArea = false;
  // Medium in the active area.
  Medium* m_medium = nullptr;

  void Reset() override;
  void UpdatePeriodicity() override;

  bool InArea(const double x, const double y, const double z) {
    if (x < m_xmin[0] || x > m_xmax[0] ||
        y < m_xmin[1] || y > m_xmax[1] ||
        z < m_xmin[2] || z > m_xmax[2]) {
      return false;
    }
    return true; 
  }
};
}
#endif
