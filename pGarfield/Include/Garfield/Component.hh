#ifndef G_COMPONENT_H
#define G_COMPONENT_H

#include <array>
#include <string>

#include "Geometry.hh"

namespace Garfield {

/// Abstract base class for components.

class Component {
 public:
  /// Default constructor.
  Component() = delete;
  /// Constructor
  Component(const std::string& name);
  /// Destructor
  virtual ~Component() {}

  /// Define the geometry.
  virtual void SetGeometry(Geometry* geo);
  /// Reset.
  virtual void Clear();

  /// Get the medium at a given location (x, y, z).
  virtual Medium* GetMedium(const double x, const double y, const double z);

  /** Calculate the drift field at given point.
   *
   * \param x,y,z coordinates [cm].
   * \param ex,ey,ez components of the electric field [V/cm].
   * \param m pointer to the medium at this location.
   * \param status status flag
   *
   * Status flags:
   *
   *             0: Inside an active medium
   *           > 0: Inside a wire of type X
   *     -4 ... -1: On the side of a plane where no wires are
   *            -5: Inside the mesh but not in an active medium
   *            -6: Outside the mesh
   *           -10: Unknown potential type (should not occur)
   *         other: Other cases (should not occur)
   */
  virtual void ElectricField(const double x, const double y, const double z,
                             double& ex, double& ey, double& ez, Medium*& m,
                             int& status) = 0;
  /// Calculate the drift field [V/cm] and potential [V] at (x, y, z).
  virtual void ElectricField(const double x, const double y, const double z,
                             double& ex, double& ey, double& ez, double& v,
                             Medium*& m, int& status) = 0;
  /// Calculate the voltage range [V].
  virtual bool GetVoltageRange(double& vmin, double& vmax) = 0;

  /** Calculate the weighting field at a given point and for a given electrode.
   * \param x,y,z coordinates [cm].
   * \param wx,wy,wz components of the weighting field [1/cm].
   * \param label name of the electrode
   */
  virtual void WeightingField(const double x, const double y, const double z,
                              double& wx, double& wy, double& wz,
                              const std::string& label);
  /** Calculate the weighting potential at a given point.
   * \param x,y,z coordinates [cm].
   * \param label name of the electrode.
   * \return weighting potential [dimensionless].
   */
  virtual double WeightingPotential(const double x, const double y,
                                    const double z, const std::string& label);
  /** Calculate the delayed weighting field at a given point and time
   * and for a given electrode.
   * \param x,y,z coordinates [cm].
   * \param t time [ns].
   * \param wx,wy,wz components of the weighting field [1/cm].
   * \param label name of the electrode
   */
  virtual void DelayedWeightingField(const double x, const double y,
                                     const double z, const double t, double& wx,
                                     double& wy, double& wz,
                                     const std::string& label);

  /** Calculate the delayed weighting potential at a given point and time
   * and for a given electrode.
   * \param x,y,z coordinates [cm].
   * \param t time [ns].
   * \param label name of the electrode
   */
  virtual double DelayedWeightingPotential(const double x, const double y,
                                           const double z, const double t,
                                           const std::string& label);

  /** Calculate the magnetic field at a given point.
   *
   * \param x,y,z coordinates [cm].
   * \param bx,by,bz components of the magnetic field [Tesla].
   * \param status status flag.
   */
  virtual void MagneticField(const double x, const double y, const double z,
                             double& bx, double& by, double& bz, int& status);
  /// Set a constant magnetic field.
  void SetMagneticField(const double bx, const double by, const double bz);

  /// Ready for use?
  virtual bool IsReady() { return m_ready; }

  /// Get the bounding box coordinates.
  virtual bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                              double& xmax, double& ymax, double& zmax);

  /// Get the coordinates of the elementary cell.
  virtual bool GetElementaryCell(double& xmin, double& ymin, double& zmin,
                                 double& xmax, double& ymax, double& zmax);

  /** Integrate the normal component of the electric field over a circle.
   * \param xc,yc centre of the circle [cm]
   * \param r radius [cm]
   * \param nI number of intervals for the integration
   *
   * \return charge enclosed in the circle [fC / cm]
   */
  double IntegrateFluxCircle(const double xc, const double yc, const double r,
                             const unsigned int nI = 50);
  /** Integrate the normal component of the electric field over a sphere.
   * \param xc,yc,zc centre of the sphere [cm]
   * \param r radius of the sphere [cm]
   * \param nI number of integration intervals in phi and theta
   *
   * \return charge enclosed in the sphere [fC]
   */
  double IntegrateFluxSphere(const double xc, const double yc, const double zc,
                             const double r, const unsigned int nI = 20);

  /** Integrate the normal component of the electric field over a parallelogram.
   * \param x0,y0,z0 coordinates of one of the corners [cm]
   * \param dx1,dy1,dz1 vector to one of the adjacent corners [cm]
   * \param dx2,dy2,dz2 vector to the other adjacent corner [cm]
   * \param nU,nV number of integration points in the two directions
   *
   * \return flux [V cm]
   */
  double IntegrateFluxParallelogram(
      const double x0, const double y0, const double z0, const double dx1,
      const double dy1, const double dz1, const double dx2, const double dy2,
      const double dz2, const unsigned int nU = 20, const unsigned int nV = 20);

  /// Integrate the normal component of the weighting field 
  /// over a parallelogram.
  double IntegrateWeightingFluxParallelogram(const std::string& label,
      const double x0, const double y0, const double z0, const double dx1,
      const double dy1, const double dz1, const double dx2, const double dy2,
      const double dz2, const unsigned int nU = 20, const unsigned int nV = 20);

  /** Integrate the electric field flux through a line from
    * (x0,y0,z0) to (x1,y1,z1) along a direction (xp,yp,zp).
    * \param x0,y0,z0 coordinates of the starting point
    * \param x1,y1,z1 coordinates of the end point
    * \param xp,yp,zp normal vector
    * \param nI number of intervals for the integration
    * \param isign include both negative and positive contributions (0)
                   or only contributions with a given polarity (+1,-1)
    */
  double IntegrateFluxLine(const double x0, const double y0, const double z0,
                           const double x1, const double y1, const double z1,
                           const double xp, const double yp, const double zp,
                           const unsigned int nI, const int isign = 0);

  /** Determine whether the line between two points crosses a wire.
    * \param x0,y0,z0 first point [cm].
    * \param x1,y1,z1 second point [cm]
    * \param xc,yc,zc point [cm] where the line crosses the wire or the
             coordinates of the wire centre.
    * \param centre flag whether to return the coordinates of the line-wire
    *        crossing point or of the wire centre.
    * \param rc radius [cm] of the wire.
    */
  virtual bool CrossedWire(const double x0, const double y0, const double z0,
                           const double x1, const double y1, const double z1,
                           double& xc, double& yc, double& zc,
                           const bool centre, double& rc);
  /** Determine whether a particle is inside the trap radius of a wire.
   * \param q0 charge of the particle [in elementary charges].
   * \param x0,y0,z0 position [cm] of the particle.
   * \param xw,yw coordinates of the wire (if applicable).
   * \param rw radius of the wire (if applicable).
   */
  virtual bool InTrapRadius(const double q0, const double x0, const double y0,
                            const double z0, double& xw, double& yw,
                            double& rw);
  /** Determine whether the line between two points crosses a plane.
    */
  virtual bool CrossedPlane(const double x0, const double y0, const double z0,
                            const double x1, const double y1, const double z1,
                            double& xc, double& yc, double& zc);

  /// Enable simple periodicity in the \f$x\f$ direction.
  void EnablePeriodicityX(const bool on = true) {
    m_periodic[0] = on;
    UpdatePeriodicity();
  }
  /// Enable simple periodicity in the \f$y\f$ direction.
  void EnablePeriodicityY(const bool on = true) {
    m_periodic[1] = on;
    UpdatePeriodicity();
  }
  /// Enable simple periodicity in the \f$z\f$ direction.
  void EnablePeriodicityZ(const bool on = true) {
    m_periodic[2] = on;
    UpdatePeriodicity();
  }
  /// Return periodicity flags.
  void IsPeriodic(bool& perx, bool& pery, bool& perz) {
    perx = m_periodic[0];
    pery = m_periodic[1];
    perz = m_periodic[2];
  }

  /// Enable mirror periodicity in the \f$x\f$ direction.
  void EnableMirrorPeriodicityX(const bool on = true) {
    m_mirrorPeriodic[0] = on;
    UpdatePeriodicity();
  }
  /// Enable mirror periodicity in the \f$y\f$ direction.
  void EnableMirrorPeriodicityY(const bool on = true) {
    m_mirrorPeriodic[1] = on;
    UpdatePeriodicity();
  }
  /// Enable mirror periodicity in the \f$y\f$ direction.
  void EnableMirrorPeriodicityZ(const bool on = true) {
    m_mirrorPeriodic[2] = on;
    UpdatePeriodicity();
  }
  /// Return mirror periodicity flags.
  void IsMirrorPeriodic(bool& perx, bool& pery, bool& perz) {
    perx = m_mirrorPeriodic[0];
    pery = m_mirrorPeriodic[1];
    perz = m_mirrorPeriodic[2];
  }

  /// Enable axial periodicity in the \f$x\f$ direction.
  void EnableAxialPeriodicityX(const bool on = true) {
    m_axiallyPeriodic[0] = on;
    UpdatePeriodicity();
  }
  /// Enable axial periodicity in the \f$y\f$ direction.
  void EnableAxialPeriodicityY(const bool on = true) {
    m_axiallyPeriodic[1] = on;
    UpdatePeriodicity();
  }
  /// Enable axial periodicity in the \f$z\f$ direction.
  void EnableAxialPeriodicityZ(const bool on = true) {
    m_axiallyPeriodic[2] = on;
    UpdatePeriodicity();
  }
  /// Return axial periodicity flags.
  void IsAxiallyPeriodic(bool& perx, bool& pery, bool& perz) {
    perx = m_axiallyPeriodic[0];
    pery = m_axiallyPeriodic[1];
    perz = m_axiallyPeriodic[2];
  }

  /// Enable rotation symmetry around the \f$x\f$ axis.
  void EnableRotationSymmetryX(const bool on = true) {
    m_rotationSymmetric[0] = on;
    UpdatePeriodicity();
  }
  /// Enable rotation symmetry around the \f$y\f$ axis.
  void EnableRotationSymmetryY(const bool on = true) {
    m_rotationSymmetric[1] = on;
    UpdatePeriodicity();
  }
  /// Enable rotation symmetry around the \f$z\f$ axis.
  void EnableRotationSymmetryZ(const bool on = true) {
    m_rotationSymmetric[2] = on;
    UpdatePeriodicity();
  }
  /// Return rotation symmetry flags.
  void IsRotationSymmetric(bool& rotx, bool& roty, bool& rotz) {
    rotx = m_rotationSymmetric[0];
    roty = m_rotationSymmetric[1];
    rotz = m_rotationSymmetric[2];
  }

  /// Switch on debugging messages.
  void EnableDebugging() { m_debug = true; }
  /// Switch off debugging messages.
  void DisableDebugging() { m_debug = false; }

  /// Does the component have a non-zero magnetic field?
  virtual bool HasMagneticField() const;

  /// Does the component have maps of the Townsend coefficient?
  virtual bool HasTownsendMap() const { return false; }
  /// Does the component have attachment maps?
  virtual bool HasAttachmentMap() const { return false; }
  /// Does the component have velocity maps?
  virtual bool HasVelocityMap() const { return false; }

  /// Get the electron attachment coefficient.
  virtual bool ElectronAttachment(const double /*x*/, const double /*y*/,
                                  const double /*z*/, double& eta) {
    eta = 0;
    return false;
  }
  /// Get the hole attachment coefficient.
  virtual bool HoleAttachment(const double /*x*/, const double /*y*/,
                              const double /*z*/, double& eta) {
    eta = 0;
    return false;
  }
  /// Get the electron Townsend coefficient.
  virtual bool ElectronTownsend(const double /*x*/, const double /*y*/,
                                const double /*z*/, double& alpha) {
    alpha = 0;
    return false;
  }
  /// Get the hole Townsend coefficient.
  virtual bool HoleTownsend(const double /*x*/, const double /*y*/,
                            const double /*z*/, double& alpha) {
    alpha = 0;
    return false;
  }
  /// Get the electron drift velocity.
  virtual bool ElectronVelocity(const double /*x*/, const double /*y*/,
                                const double /*z*/, double& vx, double& vy,
                                double& vz) {
    vx = vy = vz = 0;
    return false;
  }
  /// Get the hole drift velocity.
  virtual bool HoleVelocity(const double /*x*/, const double /*y*/,
                            const double /*z*/, double& vx, double& vy,
                            double& vz) {
    vx = vy = vz = 0;
    return false;
  }
  virtual bool GetElectronLifetime(const double /*x*/, const double /*y*/,
                                   const double /*z*/, double& etau) {
    etau = -1;
    return false;
  }
  virtual bool GetHoleLifetime(const double /*x*/, const double /*y*/,
                               const double /*z*/, double& htau) {
    htau = -1;
    return false;
  }

 protected:
  /// Class name.
  std::string m_className = "Component";

  /// Pointer to the geometry.
  Geometry* m_geometry = nullptr;

  /// Constant magnetic field.
  std::array<double, 3> m_b0 = {{0., 0., 0.}};

  /// Ready for use?
  bool m_ready = false;

  /// Switch on/off debugging messages
  bool m_debug = false;

  /// Simple periodicity in x, y, z.
  std::array<bool, 3> m_periodic = {{false, false, false}};
  /// Mirror periodicity in x, y, z.
  std::array<bool, 3> m_mirrorPeriodic = {{false, false, false}};
  /// Axial periodicity in x, y, z.
  std::array<bool, 3> m_axiallyPeriodic = {{false, false, false}};
  /// Rotation symmetry around x-axis, y-axis, z-axis.
  std::array<bool, 3> m_rotationSymmetric = {{false, false, false}};

  /// Reset the component.
  virtual void Reset() = 0;
  /// Verify periodicities.
  virtual void UpdatePeriodicity() = 0;
 private:

  double IntegrateFluxParallelogram(
      const double x0, const double y0, const double z0, const double dx1,
      const double dy1, const double dz1, const double dx2, const double dy2,
      const double dz2, const unsigned int nU, const unsigned int nV,
      const bool wfield, const std::string& label);

};
}  // namespace Garfield

#endif
