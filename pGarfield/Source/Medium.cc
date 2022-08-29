#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Medium.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/Random.hh"

namespace {

void PrintNotImplemented(const std::string& cls, const std::string& fcn) {
  std::cerr << cls << "::" << fcn << ": Function is not implemented.\n";
}

void PrintOutOfRange(const std::string& cls, const std::string& fcn,
                     const size_t i, const size_t j, const size_t k) {
  std::cerr << cls << "::" << fcn << ": Index (" << i << ", " << j << ", " << k
            << ") out of range.\n";
}

void PrintDataNotAvailable(const std::string& cls, const std::string& fcn) {
  std::cerr << cls << "::" << fcn << ": Data not available.\n";
}

bool CheckFields(const std::vector<double>& fields, const std::string& hdr,
                 const std::string& lbl) {
  if (fields.empty()) {
    std::cerr << hdr << ": Number of " << lbl << " must be > 0.\n";
    return false;
  }

  // Make sure the values are not negative.
  if (fields.front() < 0.) {
    std::cerr << hdr << ": " << lbl << " must be >= 0.\n";
    return false;
  }

  const size_t nEntries = fields.size();
  // Make sure the values are in strictly monotonic, ascending order.
  if (nEntries > 1) {
    for (size_t i = 1; i < nEntries; ++i) {
      if (fields[i] <= fields[i - 1]) {
        std::cerr << hdr << ": " << lbl << " are not in ascending order.\n";
        return false;
      }
    }
  }
  return true;
}
}

namespace Garfield {

int Medium::m_idCounter = -1;

Medium::Medium() : m_id(++m_idCounter) {
  // Initialise the tables.
  m_bFields.assign(1, 0.);
  m_bAngles.assign(1, HalfPi);

  // Set the default grid.
  SetFieldGrid(100., 100000., 20, true, 0., 0., 1, HalfPi, HalfPi, 1);
}

Medium::~Medium() {}

void Medium::SetTemperature(const double t) {
  if (t <= 0.) {
    std::cerr << m_className << "::SetTemperature:\n"
              << "    Temperature [K] must be greater than zero.\n";
    return;
  }
  m_temperature = t;
  m_isChanged = true;
}

void Medium::SetPressure(const double p) {
  if (p <= 0.) {
    std::cerr << m_className << "::SetPressure:\n"
              << "    Pressure [Torr] must be greater than zero.\n";
    return;
  }
  m_pressure = p;
  m_isChanged = true;
}

void Medium::SetDielectricConstant(const double eps) {
  if (eps < 1.) {
    std::cerr << m_className << "::SetDielectricConstant:\n"
              << "    Dielectric constant must be >= 1.\n";
    return;
  }
  m_epsilon = eps;
  m_isChanged = true;
}

double Medium::GetMassDensity() const {
  return m_density * AtomicMassUnit * m_a;
}

void Medium::GetComponent(const unsigned int i, std::string& label, double& f) {
  if (i >= m_nComponents) {
    std::cerr << m_className << "::GetComponent: Index out of range.\n";
  }

  label = m_name;
  f = 1.;
}

void Medium::SetAtomicNumber(const double z) {
  if (z < 1.) {
    std::cerr << m_className << "::SetAtomicNumber:\n"
              << "    Atomic number must be >= 1.\n";
    return;
  }
  m_z = z;
  m_isChanged = true;
}

void Medium::SetAtomicWeight(const double a) {
  if (a <= 0.) {
    std::cerr << m_className << "::SetAtomicWeight:\n"
              << "    Atomic weight must be greater than zero.\n";
    return;
  }
  m_a = a;
  m_isChanged = true;
}

void Medium::SetNumberDensity(const double n) {
  if (n <= 0.) {
    std::cerr << m_className << "::SetNumberDensity:\n"
              << "    Density [cm-3] must be greater than zero.\n";
    return;
  }
  m_density = n;
  m_isChanged = true;
}

void Medium::SetMassDensity(const double rho) {
  if (rho <= 0.) {
    std::cerr << m_className << "::SetMassDensity:\n"
              << "    Density [g/cm3] must be greater than zero.\n";
    return;
  }

  if (m_a <= 0.) {
    std::cerr << m_className << "::SetMassDensity:\n"
              << "    Atomic weight is not defined.\n";
    return;
  }
  m_density = rho / (AtomicMassUnit * m_a);
  m_isChanged = true;
}

bool Medium::Velocity(const double ex, const double ey, const double ez,
    const double bx, const double by, const double bz,
    const std::vector<std::vector<std::vector<double> > >& velE,
    const std::vector<std::vector<std::vector<double> > >& velB,
    const std::vector<std::vector<std::vector<double> > >& velX,
    const double q, double& vx, double& vy, double& vz) const {

  vx = vy = vz = 0.;
  // Make sure there is at least a table of velocities along E.
  if (velE.empty()) return false;

  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return false;

  // Compute the magnitude of the magnetic field.
  const double b = sqrt(bx * bx + by * by + bz * bz);
  // Compute the angle between B field and E field.
  const double ebang = GetAngle(ex, ey, ez, bx, by, bz, e, b);

  // Calculate the velocity along E.
  double ve = 0.;
  if (!Interpolate(e0, b, ebang, velE, ve, m_intpVel, m_extrVel)) {
    std::cerr << m_className << "::Velocity: Interpolation along E failed.\n";
    return false;
  }
  if (b < Small) {
    // No magnetic field.
    const double mu = q * ve / e;
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
    return true;
  } else if (velX.empty() || velB.empty() ||
             (m_bFields.size() == 1 && fabs(m_bFields[0]) < Small)) {
    // Magnetic field, velocities along ExB, Bt not available.
    Langevin(ex, ey, ez, bx, by, bz, q * ve / e, vx, vy, vz);
    return true;
  }

  // Magnetic field, velocities along ExB and Bt available.
  // Compute unit vectors along E, E x B and Bt.
  double ue[3] = {ex / e, ey / e, ez / e};
  double uexb[3] = {ey * bz - ez * by, ez * bx - ex * bz, ex * by - ey * bx};
  const double exb =
      sqrt(uexb[0] * uexb[0] + uexb[1] * uexb[1] + uexb[2] * uexb[2]);
  if (exb > 0.) {
    uexb[0] /= exb;
    uexb[1] /= exb;
    uexb[2] /= exb;
  } else {
    uexb[0] = ue[0];
    uexb[1] = ue[1];
    uexb[2] = ue[2];
  }

  double ubt[3] = {uexb[1] * ez - uexb[2] * ey, 
                   uexb[2] * ex - uexb[0] * ez,
                   uexb[0] * ey - uexb[1] * ex};
  const double bt = sqrt(ubt[0] * ubt[0] + ubt[1] * ubt[1] + ubt[2] * ubt[2]);
  if (bt > 0.) {
    ubt[0] /= bt;
    ubt[1] /= bt;
    ubt[2] /= bt;
  } else {
    ubt[0] = ue[0];
    ubt[1] = ue[1];
    ubt[2] = ue[2];
  }

  if (m_debug) {
    std::cout << m_className << "::Velocity:\n";
    std::printf("    unit vector along E:     (%15.5f, %15.5f, %15.5f)\n",
                ue[0], ue[1], ue[2]);
    std::printf("    unit vector along E x B: (%15.5f, %15.5f, %15.5f)\n",
                uexb[0], uexb[1], uexb[2]);
    std::printf("    unit vector along Bt:    (%15.5f, %15.5f, %15.5f)\n",
                ubt[0], ubt[1], ubt[2]);
  }

  // Calculate the velocities in all directions.
  double vexb = 0.;
  if (!Interpolate(e0, b, ebang, velX, vexb, m_intpVel, m_extrVel)) {
    std::cerr << m_className << "::Velocity: Interpolation along ExB failed.\n";
    return false;
  }
  vexb *= q;
  double vbt = 0.;
  if (!Interpolate(e0, b, ebang, velB, vbt, m_intpVel, m_extrVel)) {
    std::cerr << m_className << "::Velocity: Interpolation along Bt failed.\n";
    return false;
  }
  if (ex * bx + ey * by + ez * bz > 0.) {
    vbt = fabs(vbt);
  } else {
    vbt = -fabs(vbt);
  }
  vbt *= q * q;
  vx = q * (ve * ue[0] + vbt * ubt[0] + vexb * uexb[0]);
  vy = q * (ve * ue[1] + vbt * ubt[1] + vexb * uexb[1]);
  vz = q * (ve * ue[2] + vbt * ubt[2] + vexb * uexb[2]);
  return true;
}

void Medium::Langevin(const double ex, const double ey, const double ez,
                      double bx, double by, double bz, const double mu,
                      double& vx, double& vy, double& vz) {

  bx *= Tesla2Internal;
  by *= Tesla2Internal;
  bz *= Tesla2Internal;
  const double b2 = bx * bx + by * by + bz * bz;
  const double mu2 = mu * mu;
  const double eb = bx * ex + by * ey + bz * ez;
  const double f = mu / (1. + mu2 * b2);
  vx = f * (ex + mu * (ey * bz - ez * by) + mu2 * bx * eb);
  vy = f * (ey + mu * (ez * bx - ex * bz) + mu2 * by * eb);
  vz = f * (ez + mu * (ex * by - ey * bx) + mu2 * bz * eb);
}

void Medium::Langevin(const double ex, const double ey, const double ez,
                      double bx, double by, double bz, 
                      const double mu, const double muH,
                      double& vx, double& vy, double& vz) {

  bx *= Tesla2Internal;
  by *= Tesla2Internal;
  bz *= Tesla2Internal;
  const double b2 = bx * bx + by * by + bz * bz;
  const double mu2 = muH * muH;
  const double f = mu / (1. + mu2 * b2);
  const double eb = bx * ex + by * ey + bz * ez;
  vx = f * (ex + muH * (ey * bz - ez * by) + mu2 * bx * eb);
  vy = f * (ey + muH * (ez * bx - ex * bz) + mu2 * by * eb);
  vz = f * (ez + muH * (ex * by - ey * bx) + mu2 * bz * eb);
}

bool Medium::Diffusion(const double ex, const double ey, const double ez,
                       const double bx, const double by, const double bz,
                       const std::vector<std::vector<std::vector<double> > >& difL,
                       const std::vector<std::vector<std::vector<double> > >& difT,
                       double& dl, double& dt) const { 

  dl = dt = 0.;
  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return true;

  // Compute the magnitude of the magnetic field.
  const double b = m_tab2d ? sqrt(bx * bx + by * by + bz * bz) : 0.;
  // Compute the angle between B field and E field.
  const double ebang = m_tab2d ? GetAngle(ex, ey, ez, bx, by, bz, e, b) : 0.;

  // Interpolate.
  if (!difL.empty()) {
    if (!Interpolate(e0, b, ebang, difL, dl, m_intpDif, m_extrDif)) dl = 0.;
  }
  if (!difT.empty()) {
    if (!Interpolate(e0, b, ebang, difT, dt, m_intpDif, m_extrDif)) dt = 0.;
  }

  // If no data available, calculate
  // the diffusion coefficients using the Einstein relation
  if (difL.empty() || difT.empty()) {
    const double d = sqrt(2. * BoltzmannConstant * m_temperature / e);
    if (difL.empty()) dl = d;
    if (difT.empty()) dt = d;
  }
  // Verify values and apply scaling.
  dl = ScaleDiffusion(std::max(dl, 0.));
  dt = ScaleDiffusion(std::max(dt, 0.));
  return true;
}

bool Medium::Diffusion(const double ex, const double ey, const double ez,
  const double bx, const double by, const double bz,
  const std::vector<std::vector<std::vector<std::vector<double> > > >& diff,
  double cov[3][3]) const {

  // Initialise the tensor.
  cov[0][0] = cov[0][1] = cov[0][2] = 0.;
  cov[1][0] = cov[1][1] = cov[1][2] = 0.;
  cov[2][0] = cov[2][1] = cov[2][2] = 0.;

  if (diff.empty()) return false;

  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return true;

  // Compute the magnitude of the magnetic field.
  const double b = m_tab2d ? sqrt(bx * bx + by * by + bz * bz) : 0.;
  // Compute the angle between B field and E field.
  const double ebang = m_tab2d ? GetAngle(ex, ey, ez, bx, by, bz, e, b) : 0.;

  for (int j = 0; j < 6; ++j) {
    // Interpolate.
    double y = 0.;
    if (!Interpolate(e0, b, ebang, diff[j], y, m_intpDif, m_extrDif)) y = 0.;
    // Apply scaling.
    y = ScaleDiffusionTensor(y);
    if (j < 3) {
      cov[j][j] = y;
    } else if (j == 3) {
      cov[0][1] = cov[1][0] = y;
    } else if (j == 4) {
      cov[0][2] = cov[2][0] = y;
    } else if (j == 5) {
      cov[1][2] = cov[2][1] = y;
    }
  }
  return true;
}

bool Medium::Alpha(const double ex, const double ey, const double ez,
                   const double bx, const double by, const double bz,
                   const std::vector<std::vector<std::vector<double> > >& tab,
                   unsigned int intp, const unsigned int thr, 
                   const std::pair<unsigned int, unsigned int>& extr, 
                   double& alpha) const {
                  
  alpha = 0.;
  if (tab.empty()) return false;

  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return true;

  // Compute the magnitude of the magnetic field.
  const double b = m_tab2d ? sqrt(bx * bx + by * by + bz * bz) : 0.;
  // Compute the angle between B field and E field.
  const double ebang = m_tab2d ? GetAngle(ex, ey, ez, bx, by, bz, e, b) : 0.;

  // Interpolate.
  if (e0 < m_eFields[thr]) intp = 1;
  if (!Interpolate(e0, b, ebang, tab, alpha, intp, extr)) alpha = -30.;
  if (alpha < -20.) {
    alpha = 0.;
  } else {
    alpha = exp(alpha);
  }
  return true;
}

bool Medium::ElectronVelocity(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& vx, double& vy, double& vz) {

  return Velocity(ex, ey, ez, bx, by, bz, m_eVelE, m_eVelB, m_eVelX, -1., 
                  vx, vy, vz);
}

bool Medium::ElectronDiffusion(const double ex, const double ey,
                               const double ez, const double bx,
                               const double by, const double bz, double& dl,
                               double& dt) {

  return Diffusion(ex, ey, ez, bx, by, bz, m_eDifL, m_eDifT, dl, dt);
}

bool Medium::ElectronDiffusion(const double ex, const double ey,
                               const double ez, const double bx,
                               const double by, const double bz,
                               double cov[3][3]) {

  return Diffusion(ex, ey, ez, bx, by, bz, m_eDifM, cov);
}

bool Medium::ElectronTownsend(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& alpha) {

  if (!Alpha(ex, ey, ez, bx, by, bz, m_eAlp, m_intpAlp, m_eThrAlp, m_extrAlp, 
             alpha)) {
    return false;
  } 
  // Apply scaling.
  alpha = ScaleTownsend(alpha);
  return true;
}

bool Medium::ElectronAttachment(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& eta) {

  if (!Alpha(ex, ey, ez, bx, by, bz, m_eAtt, m_intpAtt, m_eThrAtt, m_extrAtt, 
             eta)) {
    return false;
  } 
  // Apply scaling.
  eta = ScaleAttachment(eta);
  return true;
}

bool Medium::ElectronLorentzAngle(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz,
                                  double& lor) {
  lor = 0.;
  if (m_eLor.empty()) return false;

  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return true;

  // Compute the magnitude of the magnetic field.
  const double b = m_tab2d ? sqrt(bx * bx + by * by + bz * bz) : 0.;
  // Compute the angle between B field and E field.
  const double ebang = m_tab2d ? GetAngle(ex, ey, ez, bx, by, bz, e, b) : 0.;

  // Interpolate.
  if (!Interpolate(e0, b, ebang, m_eLor, lor, m_intpLor, m_extrLor)) lor = 0.;
  // Apply scaling.
  lor = ScaleLorentzAngle(lor);
  return true;
}

double Medium::ElectronMobility() {
  if (m_eVelE.empty()) return -1.;
  return m_eVelE[0][0][0] / UnScaleElectricField(m_eFields[0]);
}

double Medium::GetElectronEnergy(const double px, const double py,
                                 const double pz, double& vx, double& vy,
                                 double& vz, const int band) {
  if (band > 0) {
    std::cerr << m_className << "::GetElectronEnergy:\n";
    std::cerr << "    Unknown band index.\n";
  }

  vx = SpeedOfLight * px / ElectronMass;
  vy = SpeedOfLight * py / ElectronMass;
  vz = SpeedOfLight * pz / ElectronMass;

  return 0.5 * (px * px + py * py + pz * pz) / ElectronMass;
}

void Medium::GetElectronMomentum(const double e, double& px, double& py,
                                 double& pz, int& band) {
  const double p = sqrt(2. * ElectronMass * e) / SpeedOfLight;
  RndmDirection(px, py, pz, p);
  band = -1;
}

double Medium::GetElectronNullCollisionRate(const int /*band*/) {
  if (m_debug) PrintNotImplemented(m_className, "GetElectronNullCollisionRate");
  return 0.;
}

double Medium::GetElectronCollisionRate(const double /*e*/,
                                        const int /*band*/) {
  if (m_debug) PrintNotImplemented(m_className, "GetElectronCollisionRate");
  return 0.;
}

bool Medium::ElectronCollision(const double e, int& type, int& level, 
    double& e1, double& dx, double& dy, double& dz, 
    std::vector<std::pair<Particle, double> >& /*secondaries*/,
    int& ndxc, int& band) {
  type = level = -1;
  e1 = e;
  ndxc = band = 0;
  RndmDirection(dx, dy, dz);

  if (m_debug) PrintNotImplemented(m_className, "GetElectronCollision");
  return false;
}

bool Medium::GetDeexcitationProduct(const unsigned int /*i*/, double& t,
                                    double& s, int& type,
                                    double& energy) const {
  if (m_debug) PrintNotImplemented(m_className, "GetDeexcitationProduct");
  t = s = energy = 0.;
  type = 0;
  return false;
}

bool Medium::HoleVelocity(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          double& vx, double& vy, double& vz) {

  return Velocity(ex, ey, ez, bx, by, bz, m_hVelE, m_hVelB, m_hVelX, +1., 
                  vx, vy, vz);
}

bool Medium::HoleDiffusion(const double ex, const double ey, const double ez,
                           const double bx, const double by, const double bz,
                           double& dl, double& dt) {
  return Diffusion(ex, ey, ez, bx, by, bz, m_hDifL, m_hDifT, dl, dt);
}

bool Medium::HoleDiffusion(const double ex, const double ey, const double ez,
                           const double bx, const double by, const double bz,
                           double cov[3][3]) {

  return Diffusion(ex, ey, ez, bx, by, bz, m_hDifM, cov);
}

bool Medium::HoleTownsend(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          double& alpha) {

  if (!Alpha(ex, ey, ez, bx, by, bz, m_hAlp, m_intpAlp, m_hThrAlp, m_extrAlp, 
             alpha)) {
    return false;
  } 
  // Apply scaling.
  alpha = ScaleTownsend(alpha);
  return true;
}

bool Medium::HoleAttachment(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& eta) {

  if (!Alpha(ex, ey, ez, bx, by, bz, m_hAtt, m_intpAtt, m_hThrAtt, m_extrAtt, 
             eta)) {
    return false;
  } 
  // Apply scaling.
  eta = ScaleAttachment(eta);
  return true;
}

double Medium::HoleMobility() {
  if (m_hVelE.empty()) return -1.;
  return m_hVelE[0][0][0] / UnScaleElectricField(m_eFields[0]);
}

bool Medium::IonVelocity(const double ex, const double ey, const double ez,
                         const double bx, const double by, const double bz,
                         double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_iMob.empty()) return false;
  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return true;
  // Compute the magnitude of the electric field.
  const double b = sqrt(bx * bx + by * by + bz * bz);

  // Compute the angle between B field and E field.
  const double ebang = m_tab2d ? GetAngle(ex, ey, ez, bx, by, bz, e, b) : 0.;
  double mu = 0.;
  if (!Interpolate(e0, b, ebang, m_iMob, mu, m_intpMob, m_extrMob)) mu = 0.;

  constexpr double q = 1.;
  mu *= q;
  if (b < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, vx, vy, vz);
  }

  return true;
}

bool Medium::IonDiffusion(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          double& dl, double& dt) {

  return Diffusion(ex, ey, ez, bx, by, bz, m_iDifL, m_iDifT, dl, dt);
}

bool Medium::IonDissociation(const double ex, const double ey, const double ez,
                             const double bx, const double by, const double bz,
                             double& diss) {

  if (!Alpha(ex, ey, ez, bx, by, bz, m_iDis, m_intpDis, m_iThrDis, m_extrDis, 
             diss)) {
    return false;
  } 
  // Apply scaling.
  diss = ScaleDissociation(diss);
  return true;
}

double Medium::IonMobility() {
  return m_iMob.empty() ? -1. : m_iMob[0][0][0];
}

bool Medium::NegativeIonVelocity(
    const double ex, const double ey, const double ez,
    const double bx, const double by, const double bz,
    double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_iMob.empty() && m_nMob.empty()) return false;

  // Compute the magnitude of the electric field.
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  const double e0 = ScaleElectricField(e);
  if (e < Small || e0 < Small) return true;
  // Compute the magnitude of the electric field.
  const double b = sqrt(bx * bx + by * by + bz * bz);

  // Compute the angle between B field and E field.
  const double ebang = m_tab2d ? GetAngle(ex, ey, ez, bx, by, bz, e, b) : 0.;
  double mu = 0.;
  if (m_nMob.empty()) {
    if (!Interpolate(e0, b, ebang, m_iMob, mu, m_intpMob, m_extrMob)) mu = 0.;
  } else {
    if (!Interpolate(e0, b, ebang, m_nMob, mu, m_intpMob, m_extrMob)) mu = 0.;
  }
  constexpr double q = -1.;
  mu *= q;
  if (b < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, vx, vy, vz);
  }
  return true;
}

double Medium::NegativeIonMobility() {
  return m_nMob.empty() ? IonMobility() : m_nMob[0][0][0];
}

bool Medium::GetOpticalDataRange(double& emin, double& emax,
                                 const unsigned int i) {
  if (i >= m_nComponents) {
    std::cerr << m_className << "::GetOpticalDataRange: Index out of range.\n";
    return false;
  }

  if (m_debug) PrintNotImplemented(m_className, "GetOpticalDataRange");
  emin = emax = 0.;
  return false;
}

bool Medium::GetDielectricFunction(const double e, double& eps1, double& eps2,
                                   const unsigned int i) {
  if (i >= m_nComponents) {
    std::cerr << m_className << "::GetDielectricFunction: Index out of range.\n";
    return false;
  }

  if (e < 0.) {
    std::cerr << m_className << "::GetDielectricFunction: Energy must be > 0.\n";
    return false;
  }

  if (m_debug) PrintNotImplemented(m_className, "GetDielectricFunction");
  eps1 = 1.;
  eps2 = 0.;
  return false;
}

bool Medium::GetPhotoAbsorptionCrossSection(const double e, double& sigma,
                                            const unsigned int i) {
  if (i >= m_nComponents) {
    std::cerr << m_className << "::GetPhotoAbsorptionCrossSection:\n";
    std::cerr << "    Component " << i << " does not exist.\n";
    return false;
  }

  if (e < 0.) {
    std::cerr << m_className << "::GetPhotoAbsorptionCrossSection:\n";
    std::cerr << "    Energy must be > 0.\n";
    return false;
  }

  if (m_debug) {
    PrintNotImplemented(m_className, "GetPhotoAbsorptionCrossSection");
  }
  sigma = 0.;
  return false;
}

double Medium::GetPhotonCollisionRate(const double e) {
  double sigma = 0.;
  if (!GetPhotoAbsorptionCrossSection(e, sigma)) return 0.;

  return sigma * m_density * SpeedOfLight;
}

bool Medium::GetPhotonCollision(const double e, int& type, int& level,
                                double& e1, double& ctheta, int& nsec,
                                double& esec) {
  type = level = -1;
  e1 = e;
  ctheta = 1.;
  nsec = 0;
  esec = 0.;
  return false;
}

void Medium::SetFieldGrid(double emin, double emax, const size_t ne, bool logE,
                          double bmin, double bmax, const size_t nb,
                          double amin, double amax, const size_t na) {
  // Check if the requested E-field range makes sense.
  if (ne <= 0) {
    std::cerr << m_className << "::SetFieldGrid:\n"
              << "    Number of E-fields must be > 0.\n";
    return;
  }

  if (emin < 0. || emax < 0.) {
    std::cerr << m_className << "::SetFieldGrid:\n"
              << "    Electric fields must be positive.\n";
    return;
  }

  if (emax < emin) {
    std::cerr << m_className << "::SetFieldGrid: Swapping min./max. E-field.\n";
    std::swap(emin, emax);
  }

  double estep = 0.;
  if (logE) {
    // Logarithmic scale
    if (emin < Small) {
      std::cerr << m_className << "::SetFieldGrid:\n"
                << "    Min. E-field must be non-zero for log. scale.\n";
      return;
    }
    if (ne == 1) {
      std::cerr << m_className << "::SetFieldGrid:\n"
                << "    Number of E-fields must be > 1 for log. scale.\n";
      return;
    }
    estep = pow(emax / emin, 1. / (ne - 1.));
  } else {
    // Linear scale
    if (ne > 1) estep = (emax - emin) / (ne - 1.);
  }

  // Check if the requested B-field range makes sense.
  if (nb <= 0) {
    std::cerr << m_className << "::SetFieldGrid:\n"
              << "    Number of B-fields must be > 0.\n";
    return;
  }
  if (bmax < 0. || bmin < 0.) {
    std::cerr << m_className << "::SetFieldGrid:\n"
              << "    Magnetic fields must be positive.\n";
    return;
  }
  if (bmax < bmin) {
    std::cerr << m_className << "::SetFieldGrid: Swapping min./max. B-field.\n";
    std::swap(bmin, bmax);
  }

  const double bstep = nb > 1 ? (bmax - bmin) / (nb - 1.) : 0.;

  // Check if the requested angular range makes sense.
  if (na <= 0) {
    std::cerr << m_className << "::SetFieldGrid:\n"
              << "    Number of angles must be > 0.\n";
    return;
  }
  if (amax < 0. || amin < 0.) {
    std::cerr << m_className << "::SetFieldGrid:"
              << "    Angles must be positive.\n";
    return;
  }
  if (amax < amin) {
    std::cerr << m_className << "::SetFieldGrid: Swapping min./max. angle.\n";
    std::swap(amin, amax);
  }
  const double astep = na > 1 ? (amax - amin) / (na - 1.) : 0;

  // Setup the field grids.
  std::vector<double> eFields(ne);
  std::vector<double> bFields(nb);
  std::vector<double> bAngles(na);
  for (size_t i = 0; i < ne; ++i) {
    eFields[i] = logE ? emin * pow(estep, i) : emin + i * estep;
  }
  for (size_t i = 0; i < nb; ++i) {
    bFields[i] = bmin + i * bstep;
  }
  for (size_t i = 0; i < na; ++i) {
    bAngles[i] = amin + i * astep;
  }
  SetFieldGrid(eFields, bFields, bAngles);
}

void Medium::SetFieldGrid(const std::vector<double>& efields,
                          const std::vector<double>& bfields,
                          const std::vector<double>& angles) {
  const std::string hdr = m_className + "::SetFieldGrid";
  if (!CheckFields(efields, hdr, "E-fields")) return;
  if (!CheckFields(bfields, hdr, "B-fields")) return;
  if (!CheckFields(angles, hdr, "angles")) return;

  if (m_debug) {
    std::cout << m_className << "::SetFieldGrid:\n    E-fields:\n";
    for (const auto efield : efields) std::cout << "      " << efield << "\n";
    std::cout << "    B-fields:\n";
    for (const auto bfield : bfields) std::cout << "      " << bfield << "\n";
    std::cout << "    Angles:\n";
    for (const auto angle : angles) std::cout << "      " << angle << "\n";
  }

  // Clone the existing tables.
  // Electrons
  Clone(m_eVelE, efields, bfields, angles, m_intpVel, m_extrVel, 0.,
        "electron velocity along E");
  Clone(m_eVelB, efields, bfields, angles, m_intpVel, m_extrVel, 0.,
        "electron velocity along Bt");
  Clone(m_eVelX, efields, bfields, angles, m_intpVel, m_extrVel, 0.,
        "electron velocity along ExB");
  Clone(m_eDifL, efields, bfields, angles, m_intpDif, m_extrDif, 0.,
        "electron longitudinal diffusion");
  Clone(m_eDifT, efields, bfields, angles, m_intpDif, m_extrDif, 0.,
        "electron transverse diffusion");
  Clone(m_eAlp, efields, bfields, angles, m_intpAlp, m_extrAlp, -30., 
        "electron Townsend coefficient");
  Clone(m_eAtt, efields, bfields, angles, m_intpAtt, m_extrAtt, -30., 
        "electron attachment coefficient");
  Clone(m_eLor, efields, bfields, angles, m_intpLor, m_extrLor, 0., 
        "electron Lorentz angle");
  if (!m_eDifM.empty()) {
    Clone(m_eDifM, 6, efields, bfields, angles, m_intpDif, m_extrDif, 0., 
          "electron diffusion tensor");
  }

  // Holes
  Clone(m_hVelE, efields, bfields, angles, m_intpVel, m_extrVel, 0.,
        "hole velocity along E");
  Clone(m_hVelB, efields, bfields, angles, m_intpVel, m_extrVel, 0.,
        "hole velocity along Bt");
  Clone(m_hVelX, efields, bfields, angles, m_intpVel, m_extrVel, 0.,
        "hole velocity along ExB");
  Clone(m_hDifL, efields, bfields, angles, m_intpDif, m_extrDif, 0.,
        "hole longitudinal diffusion");
  Clone(m_hDifT, efields, bfields, angles, m_intpDif, m_extrDif, 0.,
        "hole transverse diffusion");
  Clone(m_hAlp, efields, bfields, angles, m_intpAlp, m_extrAlp, -30., 
        "hole Townsend coefficient");
  Clone(m_hAtt, efields, bfields, angles, m_intpAtt, m_extrAtt, -30., 
        "hole attachment coefficient");
  if (!m_hDifM.empty()) {
    Clone(m_hDifM, 6, efields, bfields, angles, m_intpDif, m_extrDif, 0., 
          "hole diffusion tensor");
  }

  // Ions
  Clone(m_iMob, efields, bfields, angles, m_intpMob, m_extrMob, 0.,
        "ion mobility");
  Clone(m_iDifL, efields, bfields, angles, m_intpDif, m_extrDif, 0., 
        "ion longitudinal diffusion");
  Clone(m_iDifT, efields, bfields, angles, m_intpDif, m_extrDif, 0., 
        "ion transverse diffusion");
  Clone(m_iDis, efields, bfields, angles, m_intpDis, m_extrDis, -30., 
        "ion dissociation");

  if (bfields.size() > 1 || angles.size() > 1) m_tab2d = true;
  m_eFields = efields;
  m_bFields = bfields;
  m_bAngles = angles;
}

void Medium::GetFieldGrid(std::vector<double>& efields,
                          std::vector<double>& bfields,
                          std::vector<double>& angles) {
  efields = m_eFields;
  bfields = m_bFields;
  angles = m_bAngles;
}

bool Medium::SetEntry(const size_t i, const size_t j, const size_t k,
                      const std::string& fcn,
                      std::vector<std::vector<std::vector<double> > >& tab,
                      const double val) {

  if (i >= m_eFields.size() || j >= m_bFields.size() || k >= m_bAngles.size()) {
    PrintOutOfRange(m_className, "Set" + fcn, i, j, k);
    return false;
  }
  if (tab.empty()) {
    Init(m_eFields.size(), m_bFields.size(), m_bAngles.size(), tab, val);
  }
  tab[k][j][i] = val;
  return true;
}

bool Medium::GetEntry(const size_t i, const size_t j, const size_t k, 
                      const std::string& fcn, 
                      const std::vector<std::vector<std::vector<double> > >& tab,
                      double& val) const {
  val = 0.;
  if (i >= m_eFields.size() || j >= m_bFields.size() || k >= m_bAngles.size()) {
    PrintOutOfRange(m_className, "Get" + fcn, i, j, k);
    return false;
  }
  if (tab.empty()) {
    if (m_debug) PrintDataNotAvailable(m_className, "Get" + fcn);
    return false;
  }
  val = tab[k][j][i];
  return true;
}

void Medium::ResetTables() {
  ResetElectronVelocity();
  ResetElectronDiffusion();
  ResetElectronTownsend();
  ResetElectronAttachment();
  ResetElectronLorentzAngle();
  
  ResetHoleVelocity();
  ResetHoleDiffusion();
  ResetHoleTownsend();
  ResetHoleAttachment();
 
  ResetIonMobility();
  ResetIonDiffusion();
  ResetIonDissociation();

  ResetNegativeIonMobility();
}

void Medium::Clone(std::vector<std::vector<std::vector<double> > >& tab,
                   const std::vector<double>& efields,
                   const std::vector<double>& bfields,
                   const std::vector<double>& angles,
                   const unsigned int intp,
                   const std::pair<unsigned int, unsigned int>& extr,
                   const double init, const std::string& lbl) {
  if (m_debug) {
    std::cout << m_className << "::Clone: Copying " << lbl << " to new grid.\n";
  }

  if (tab.empty()) {
    if (m_debug) std::cout << m_className << "::Clone: Table is empty.\n";
    return;
  }
  // Get the dimensions of the new grid.
  const auto nE = efields.size();
  const auto nB = bfields.size();
  const auto nA = angles.size();

  // Create a temporary table to store the values at the new grid points.
  std::vector<std::vector<std::vector<double> > > tabClone;
  Init(nE, nB, nA, tabClone, init);

  // Fill the temporary table.
  for (size_t i = 0; i < nE; ++i) {
    const double e = efields[i];
    for (size_t j = 0; j < nB; ++j) {
      const double b = bfields[j];
      for (size_t k = 0; k < nA; ++k) {
        const double a = angles[k];
        double val = 0.;
        if (!Interpolate(e, b, a, tab, val, intp, extr)) {
          std::cerr << m_className << "::Clone:\n"
                    << "    Interpolation of " << lbl << " failed.\n"
                    << "    Cannot copy value to new grid at E = " << e
                    << ", B = " << b << ", angle: " << a << "\n";
          continue;
        }
        tabClone[k][j][i] = val;
      }
    }
  }
  // Copy the values to the original table.
  tab.swap(tabClone);
}

void Medium::Clone(
    std::vector<std::vector<std::vector<std::vector<double> > > >& tab,
    const size_t n, const std::vector<double>& efields,
    const std::vector<double>& bfields, const std::vector<double>& angles,
    const unsigned int intp, const std::pair<unsigned int, unsigned int>& extr,
    const double init, const std::string& lbl) {
  // If the table does not exist, do nothing.
  if (tab.empty()) return;

  // Get the dimensions of the new grid.
  const auto nE = efields.size();
  const auto nB = bfields.size();
  const auto nA = angles.size();

  // Create a temporary table to store the values at the new grid points.
  std::vector<std::vector<std::vector<std::vector<double> > > > tabClone;
  Init(nE, nB, nA, n, tabClone, init);

  // Fill the temporary table.
  for (size_t l = 0; l < n; ++l) {
    for (size_t i = 0; i < nE; ++i) {
      const double e = efields[i];
      for (size_t j = 0; j < nB; ++j) {
        const double b = bfields[j];
        for (size_t k = 0; k < nA; ++k) {
          const double a = angles[k];
          double val = 0.;
          if (!Interpolate(e, b, a, tab[l], val, intp, extr)) {
            std::cerr << m_className << "::Clone:\n"
                      << "    Interpolation of " << lbl << " failed.\n"
                      << "    Cannot copy value to new grid at index " << l
                      << ", E = " << e << ", B = " << b << ", angle: " << a
                      << "\n";
            continue;
          }
          tabClone[l][k][j][i] = val;
        }
      }
    }
  }
  // Copy the values to the original table.
  tab.swap(tabClone);
}

bool Medium::SetIonMobility(const size_t ie, const size_t ib,
                            const size_t ia, const double mu) {
  // Check the index.
  if (ie >= m_eFields.size() || ib >= m_bFields.size() ||
      ia >= m_bAngles.size()) {
    PrintOutOfRange(m_className, "SetIonMobility", ie, ib, ia);
    return false;
  }

  if (m_iMob.empty()) {
    std::cerr << m_className << "::SetIonMobility:\n"
              << "    Ion mobility table not initialised.\n";
    return false;
  }

  if (mu == 0.) {
    std::cerr << m_className << "::SetIonMobility: Zero value not allowed.\n";
    return false;
  }

  m_iMob[ia][ib][ie] = mu;
  if (m_debug) {
    std::cout << m_className << "::SetIonMobility:\n    Ion mobility at E = "
              << m_eFields[ie] << " V/cm, B = " 
              << m_bFields[ib] << " T, angle " 
              << m_bAngles[ia] << " set to " << mu << " cm2/(V ns).\n";
  }
  return true;
}

bool Medium::SetIonMobility(const std::vector<double>& efields,
                            const std::vector<double>& mobs,
                            const bool negativeIons) {
  if (efields.size() != mobs.size()) {
    std::cerr << m_className << "::SetIonMobility:\n"
              << "    E-field and mobility arrays have different sizes.\n";
    return false;
  }

  if (negativeIons) {
    ResetNegativeIonMobility();
  } else {
    ResetIonMobility();
  }
  const auto nE = m_eFields.size();
  const auto nB = m_bFields.size();
  const auto nA = m_bAngles.size();
  if (negativeIons) {
    Init(nE, nB, nA, m_nMob, 0.);
  } else {
    Init(nE, nB, nA, m_iMob, 0.);
  }
  for (size_t i = 0; i < nE; ++i) {
    const double e = m_eFields[i];
    const double mu = Interpolate1D(e, mobs, efields, m_intpMob, m_extrMob);
    if (negativeIons) {
      m_nMob[0][0][i] = mu;
    } else {
      m_iMob[0][0][i] = mu;
    }
  }
  if (!m_tab2d) return true;
  for (size_t i = 0; i < nA; ++i) {
    for (size_t j = 0; j < nB; ++j) {
      for (size_t k = 0; k < nE; ++k) {
        if (negativeIons) {
          m_nMob[i][j][k] = m_nMob[0][0][k];
        } else {
          m_iMob[i][j][k] = m_iMob[0][0][k];
        }
      }
    }
  }
  return true;
}

void Medium::SetExtrapolationMethodVelocity(const std::string& low,
                                            const std::string& high) {
  SetExtrapolationMethod(low, high, m_extrVel, "Velocity");
}

void Medium::SetExtrapolationMethodDiffusion(const std::string& low,
                                             const std::string& high) {
  SetExtrapolationMethod(low, high, m_extrDif, "Diffusion");
}

void Medium::SetExtrapolationMethodTownsend(const std::string& low,
                                            const std::string& high) {
  SetExtrapolationMethod(low, high, m_extrAlp, "Townsend");
}

void Medium::SetExtrapolationMethodAttachment(const std::string& low,
                                              const std::string& high) {
  SetExtrapolationMethod(low, high, m_extrAtt, "Attachment");
}

void Medium::SetExtrapolationMethodIonMobility(const std::string& low,
                                               const std::string& high) {
  SetExtrapolationMethod(low, high, m_extrMob, "IonMobility");
}

void Medium::SetExtrapolationMethodIonDissociation(const std::string& low,
                                                   const std::string& high) {
  SetExtrapolationMethod(low, high, m_extrDis, "IonDissociation");
}

void Medium::SetExtrapolationMethod(const std::string& low,
                                    const std::string& high,
                                    std::pair<unsigned int, unsigned int>& extr,
                                    const std::string& fcn) {
  unsigned int i = 0;
  if (GetExtrapolationIndex(low, i)) {
    extr.first = i;
  } else {
    std::cerr << m_className << "::SetExtrapolationMethod" << fcn << ":\n"
              << "    Unknown extrapolation method (" << low << ")\n";
  }
  unsigned int j = 0;
  if (GetExtrapolationIndex(high, j)) {
    extr.second = j;
  } else {
    std::cerr << m_className << "::SetExtrapolationMethod" << fcn << ":\n"
              << "    Unknown extrapolation method (" << high << ")\n";
  }
}

bool Medium::GetExtrapolationIndex(std::string str, unsigned int& nb) const {
  // Convert to upper-case.
  std::transform(str.begin(), str.end(), str.begin(), toupper);

  if (str == "CONST" || str == "CONSTANT") {
    nb = 0;
  } else if (str == "LIN" || str == "LINEAR") {
    nb = 1;
  } else if (str == "EXP" || str == "EXPONENTIAL") {
    nb = 2;
  } else {
    return false;
  }

  return true;
}

size_t Medium::SetThreshold(
    const std::vector<std::vector<std::vector<double> > >& tab) const {

  if (tab.empty()) return 0;
  const auto nE = m_eFields.size();
  const auto nB = m_bFields.size();
  const auto nA = m_bAngles.size();
  for (size_t i = 0; i < nE; ++i) {
    bool below = false;
    for (size_t k = 0; k < nA; ++k) {
      for (size_t j = 0; j < nB; ++j) {
        if (tab[k][j][i] < -20.) {
          below = true;
          break;
        } 
      }
      if (below) break; 
    }
    if (below) continue;
    return i;
  } 
  return nE - 1; 
}

void Medium::SetInterpolationMethodVelocity(const unsigned int intrp) {
  if (intrp > 0) m_intpVel = intrp;
}

void Medium::SetInterpolationMethodDiffusion(const unsigned int intrp) {
  if (intrp > 0) m_intpDif = intrp;
}

void Medium::SetInterpolationMethodTownsend(const unsigned int intrp) {
  if (intrp > 0) m_intpAlp = intrp;
}

void Medium::SetInterpolationMethodAttachment(const unsigned int intrp) {
  if (intrp > 0) m_intpAtt = intrp;
}

void Medium::SetInterpolationMethodIonMobility(const unsigned int intrp) {
  if (intrp > 0) m_intpMob = intrp;
}

void Medium::SetInterpolationMethodIonDissociation(const unsigned int intrp) {
  if (intrp > 0) m_intpDis = intrp;
}

double Medium::GetAngle(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        const double emag, const double bmag) const {
  const double eb = emag * bmag; 
  if (eb <= 0.) return m_bAngles[0];
  const double einb = fabs(ex * bx + ey * by + ez * bz);
  if (einb > 0.2 * eb) {
    double exb[3] = {ex * by - ey * bx, ex * bz - ez * bx, ez * by - ey * bz};
    return asin(
        std::min(1., sqrt(exb[0] * exb[0] + exb[1] * exb[1] + exb[2] * exb[2]) / eb));
  }
  return acos(std::min(1., einb / eb));
}

bool Medium::Interpolate(
    const double e, const double b, const double a,
    const std::vector<std::vector<std::vector<double> > >& table, double& y,
    const unsigned int intp,
    const std::pair<unsigned int, unsigned int>& extr) const {
  if (table.empty()) {
    y = 0.;
    return false;  // TODO: true!
  }

  if (m_tab2d) {
    return Numerics::Boxin3(table, m_bAngles, m_bFields, m_eFields,
                            m_bAngles.size(), m_bFields.size(), m_eFields.size(),
                            a, b, e, y, intp);
  } else {
    y = Interpolate1D(e, table[0][0], m_eFields, intp, extr);
  }
  return true;
}

double Medium::Interpolate1D(
    const double e, const std::vector<double>& table,
    const std::vector<double>& fields, const unsigned int intpMeth,
    const std::pair<unsigned int, unsigned int>& extr) const {
  // This function is a generalized version of the Fortran functions
  // GASVEL, GASVT1, GASVT2, GASLOR, GASMOB, GASDFT, and GASDFL
  // for the case of a 1D table. All variables are generic.

  const auto nSizeTable = fields.size();

  if (e < 0. || nSizeTable < 1) return 0.;

  double result = 0.;

  if (nSizeTable == 1) {
    // Only one point
    result = table[0];
  } else if (e < fields[0]) {
    // Extrapolation towards small fields
    if (fields[0] >= fields[1]) {
      if (m_debug) {
        std::cerr << m_className << "::Interpolate1D:\n";
        std::cerr << "    First two field values coincide.\n";
        std::cerr << "    No extrapolation to lower fields.\n";
      }
      result = table[0];
    } else if (extr.first == 1) {
      // Linear extrapolation
      const double extr4 = (table[1] - table[0]) / (fields[1] - fields[0]);
      const double extr3 = table[0] - extr4 * fields[0];
      result = extr3 + extr4 * e;
    } else if (extr.first == 2) {
      // Logarithmic extrapolation
      const double extr4 = log(table[1] / table[0]) / (fields[1] - fields[0]);
      const double extr3 = log(table[0] - extr4 * fields[0]);
      result = std::exp(std::min(50., extr3 + extr4 * e));
    } else {
      result = table[0];
    }
  } else if (e > fields[nSizeTable - 1]) {
    // Extrapolation towards large fields
    if (fields[nSizeTable - 1] <= fields[nSizeTable - 2]) {
      if (m_debug) {
        std::cerr << m_className << "::Interpolate1D:\n";
        std::cerr << "    Last two field values coincide.\n";
        std::cerr << "    No extrapolation to higher fields.\n";
      }
      result = table[nSizeTable - 1];
    } else if (extr.second == 1) {
      // Linear extrapolation
      const double extr2 = (table[nSizeTable - 1] - table[nSizeTable - 2]) /
                           (fields[nSizeTable - 1] - fields[nSizeTable - 2]);
      const double extr1 =
          table[nSizeTable - 1] - extr2 * fields[nSizeTable - 1];
      result = extr1 + extr2 * e;
    } else if (extr.second == 2) {
      // Logarithmic extrapolation
      const double extr2 = log(table[nSizeTable - 1] / table[nSizeTable - 2]) /
                           (fields[nSizeTable - 1] - fields[nSizeTable - 2]);
      const double extr1 =
          log(table[nSizeTable - 1]) - extr2 * fields[nSizeTable - 1];
      result = exp(std::min(50., extr1 + extr2 * e));
    } else {
      result = table[nSizeTable - 1];
    }
  } else {
    // Intermediate points, spline interpolation (not implemented).
    // Intermediate points, Newtonian interpolation
    result = Numerics::Divdif(table, fields, nSizeTable, e, intpMeth);
  }

  return result;
}

void Medium::Init(const size_t nE, const size_t nB, const size_t nA,
                  std::vector<std::vector<std::vector<double> > >& tab,
                  const double val) {
  if (nE == 0 || nB == 0 || nA == 0) {
    std::cerr << m_className << "::Init: Invalid grid.\n";
    return;
  }
  tab.assign(
      nA, std::vector<std::vector<double> >(nB, std::vector<double>(nE, val)));
}

void Medium::Init(
    const size_t nE, const size_t nB, const size_t nA, const size_t nT,
    std::vector<std::vector<std::vector<std::vector<double> > > >& tab,
    const double val) {
  if (nE == 0 || nB == 0 || nA == 0 || nT == 0) {
    std::cerr << m_className << "::Init: Invalid grid.\n";
    return;
  }

  tab.assign(nT, std::vector<std::vector<std::vector<double> > >(
                     nA, std::vector<std::vector<double> >(
                             nB, std::vector<double>(nE, val))));
}
}
