#include <cmath>
#include <iostream>

#include "Garfield/MediumCdTe.hh"
#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

MediumCdTe::MediumCdTe() : Medium() {
  m_className = "MediumCdTe";
  m_name = "CdTe";

  SetTemperature(300.);
  SetDielectricConstant(10.9);
  SetAtomicNumber(48.52);
  SetAtomicWeight(240.01);
  SetMassDensity(5.85);

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = false;

  m_w = 4.43;
  m_fano = 0.1;
}

void MediumCdTe::GetComponent(const unsigned int i, std::string& label,
                              double& f) {
  if (i == 0) {
    label = "Cd";
    f = 0.5;
  } else if (i == 1) {
    label = "Te";
    f = 0.5;
  } else {
    std::cerr << m_className << "::GetComponent: Index out of range.\n";
  }
}

bool MediumCdTe::ElectronVelocity(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz, double& vx,
                                  double& vy, double& vz) {
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
  const double mu = -m_eMobility;
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

bool MediumCdTe::ElectronTownsend(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz,
                                  double& alpha) {
  alpha = 0.;
  if (!m_eAlp.empty()) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool MediumCdTe::ElectronAttachment(const double ex, const double ey,
                                    const double ez, const double bx,
                                    const double by, const double bz,
                                    double& eta) {
  eta = 0.;
  if (!m_eAtt.empty()) {
    // Interpolation in user table.
    return Medium::ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

bool MediumCdTe::HoleVelocity(const double ex, const double ey, const double ez,
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
  const double mu = m_hMobility;
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

bool MediumCdTe::HoleTownsend(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& alpha) {
  alpha = 0.;
  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return false;
}

bool MediumCdTe::HoleAttachment(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& eta) {
  eta = 0.;
  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

void MediumCdTe::SetLowFieldMobility(const double mue, const double muh) {
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

void MediumCdTe::UnsetLowFieldMobility() {
  m_userMobility = false;
  m_isChanged = true;
}

void MediumCdTe::UpdateTransportParameters() {

  if (!m_userMobility) {
    const double t = m_temperature / 300.;
    m_eMobility = 1.05e-6 * pow(t, -1.7);
    m_hMobility = 0.1e-6 * pow(t, 1.67);
  }
}

}
