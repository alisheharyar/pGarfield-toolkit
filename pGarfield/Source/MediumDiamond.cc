#include <cmath>
#include <iostream>

#include "Garfield/MediumDiamond.hh"
#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

MediumDiamond::MediumDiamond() : Medium() {
  m_className = "MediumDiamond";
  m_name = "Diamond";

  SetTemperature(300.);
  SetDielectricConstant(5.7);
  SetAtomicNumber(6);
  SetAtomicWeight(12.01);
  SetMassDensity(3.51);

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = false;

  m_w = 13.6;
  // TODO:
  m_fano = 0.1;
}

void MediumDiamond::GetComponent(const unsigned int i, std::string& label,
                                 double& f) {
  if (i == 0) {
    label = "C";
    f = 1.;
  } else {
    std::cerr << m_className << "::GetComponent: Index out of range.\n";
  }
}

bool MediumDiamond::ElectronVelocity(
    const double ex, const double ey, const double ez, 
    const double bx, const double by, const double bz, 
    double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_isChanged) {
    UpdateTransportParameters();
    m_isChanged = false;
  }
  if (!m_eVelE.empty()) {
    // Interpolation in user table.
    return Medium::ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }
  // Calculate the mobility
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  const double mu = -m_eMobility / (1. + m_eMobility * emag / m_eSatVel);
  const double b2 = bx * bx + by * by + bz * bz;
  if (b2 < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, m_eHallFactor * mu, vx, vy, vz); 
  }
  return true;
}

bool MediumDiamond::ElectronTownsend(
    const double ex, const double ey, const double ez, 
    const double bx, const double by, const double bz, double& alpha) {
  alpha = 0.;
  if (!m_eAlp.empty()) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool MediumDiamond::ElectronAttachment(
    const double ex, const double ey, const double ez, 
    const double bx, const double by, const double bz, double& eta) {
  eta = 0.;
  if (!m_eAtt.empty()) {
    // Interpolation in user table.
    return Medium::ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

bool MediumDiamond::HoleVelocity(
    const double ex, const double ey, const double ez,
    const double bx, const double by, const double bz,
    double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_isChanged) {
    UpdateTransportParameters();
    m_isChanged = false;
  }
  if (!m_hVelE.empty()) {
    // Interpolation in user table.
    return Medium::HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }
  // Calculate the mobility
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  const double mu = m_hMobility / (1. + m_hMobility * emag / m_hSatVel);
  const double b2 = bx * bx + by * by + bz * bz;
  if (b2 < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, m_hHallFactor * mu, vx, vy, vz);
  }
  return true;
}

bool MediumDiamond::HoleTownsend(
    const double ex, const double ey, const double ez,
    const double bx, const double by, const double bz, double& alpha) {
  alpha = 0.;
  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool MediumDiamond::HoleAttachment(
    const double ex, const double ey, const double ez, 
    const double bx, const double by, const double bz, double& eta) {
  eta = 0.;
  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

void MediumDiamond::SetLowFieldMobility(const double mue, const double muh) {
  if (mue <= 0. || muh <= 0.) {
    std::cerr << m_className << "::SetLowFieldMobility:\n"
              << "    Mobility must be greater than zero.\n";
    return;
  }
  m_eMobility = mue;
  m_hMobility = muh;
  m_userMobility = true;
  m_isChanged = true;
}

void MediumDiamond::UnsetLowFieldMobility() {
  m_userMobility = false;
  m_isChanged = true;
}

void MediumDiamond::SetSaturationVelocity(const double vsate, 
                                          const double vsath) {
  std::lock_guard<std::mutex> guard(m_mutex);
  if (vsate <= 0. || vsath <= 0.) {
    std::cerr << m_className << "::SetSaturationVelocity:\n"
              << "    Velocity must be greater than zero.\n";
    return;
  }
  m_eSatVel = vsate;
  m_hSatVel = vsath;
}

void MediumDiamond::UnsetSaturationVelocity() {
  std::lock_guard<std::mutex> guard(m_mutex);
  m_eSatVel = 2.6e-2;
  m_hSatVel = 1.6e-2;
}
 
void MediumDiamond::UpdateTransportParameters() {
  std::lock_guard<std::mutex> guard(m_mutex);

  // M. Pomorski, Electronic Properties of Single Crystal CVD Diamond and  
  // Its Suitability for Particle Detection in Hadron Physics Experiments. 
  // Wolfgang von Goethe University, Frankfurt am Main
  // E. Bossini, N. Miafra, Front. Phys. 8 (2020), 
  // https://doi.org/10.3389/fphy.2020.00248
  if (!m_userMobility) {
    const double t = m_temperature / 300.;
    m_eMobility = 4.551e-6 * pow(t, -1.5);
    m_hMobility = 2.750e-6 * pow(t, -1.5);
  }
}

}
