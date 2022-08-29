#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumGaAs.hh"
#include "Garfield/Random.hh"

namespace Garfield {

MediumGaAs::MediumGaAs() : Medium() {
  m_className = "MediumGaAs";
  m_name = "GaAs";

  SetTemperature(300.);
  SetDielectricConstant(12.9);
  SetAtomicNumber(32);
  SetAtomicWeight(72.32);
  SetMassDensity(5.317);

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = false;

  m_w = 4.35;
  m_fano = 0.1;
}

void MediumGaAs::GetComponent(const unsigned int i, std::string& label,
                              double& f) {
  if (i == 0) {
    label = "Ga";
    f = 0.5;
  } else if (i == 1) {
    label = "As";
    f = 0.5;
  }
}

bool MediumGaAs::ElectronVelocity(const double ex, const double ey,
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
  // Calculate the mobility.
  // - J. J. Barnes, R. J. Lomax, G. I. Haddad, 
  //   IEEE Trans. Electron Devices ED-23 (1976), 1042.
  const double e2 = ex * ex + ey * ey + ez * ez;
  // Inverse of the critical field.
  constexpr double r = 1. / 4000.;
  constexpr double r4 = r * r * r * r;
  const double er4 = e2 * e2 * r4;
  const double mu = -(m_eMobility + er4 * m_eSatVel / sqrt(e2)) / (1. + er4);
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

bool MediumGaAs::ElectronTownsend(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz,
                                  double& alpha) {
  alpha = 0.;
  if (m_isChanged) {
    UpdateTransportParameters();
    m_isChanged = false;
  }
  if (!m_eAlp.empty()) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  if (emag > Small) {
    alpha = m_eImpactA * exp(-pow(m_eImpactB / emag, 1.82));
  } 
  return true;
}

bool MediumGaAs::ElectronAttachment(const double ex, const double ey,
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

bool MediumGaAs::HoleVelocity(const double ex, const double ey, const double ez,
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
  // Calculate the mobility.
  // - J. J. Barnes, R. J. Lomax, G. I. Haddad, 
  //   IEEE Trans. Electron Devices ED-23 (1976), 1042–1048.
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  // Inverse of the critical field.
  constexpr double r = 1. / 4000.;
  const double mu = (m_hMobility + m_hSatVel * r) / (1. + emag * r);
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

bool MediumGaAs::HoleTownsend(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& alpha) {
  alpha = 0.;
  if (m_isChanged) {
    UpdateTransportParameters();
    m_isChanged = false;
  }
  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  if (emag > Small) {
    // alpha = m_hImpactA * exp(-m_hImpactB / emag);
    alpha = m_hImpactA * exp(-pow(m_hImpactB / emag, 1.75));
  } 
  return true;
}

bool MediumGaAs::HoleAttachment(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& eta) {
  eta = 0.;
  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

void MediumGaAs::SetLowFieldMobility(const double mue, const double muh) {

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

void MediumGaAs::UnsetLowFieldMobility() {
  m_userMobility = false;
  m_isChanged = true;
}

void MediumGaAs::UpdateTransportParameters() {

  const double t = m_temperature / 300.;
  // Update the low field lattice mobility.
  if (!m_userMobility) {
    // Temperature dependence as in Sentaurus Device and Silvaco Atlas.
    constexpr double eMu0 = 8.0e-6;
    constexpr double hMu0 = 0.4e-6;
    m_eMobility = eMu0 / t;
    m_hMobility = hMu0 * pow(t, -2.1);
  }
  //  - J. S. Blakemore, Journal of Applied Physics 53, R123 (1982)
  //    https://doi.org/10.1063/1.331665
  // m_eMobility = 8.0e-6 * pow(t, -2.3);

  // Update the saturation velocity.
  //  - M. J. Littlejohn, J. R. Hauser, T. H. Glisson, 
  //    J. Appl. Phys. 48 (1977), 4587
  m_eSatVel = m_hSatVel = std::max(1.13e-2 - 3.6e-3 * t, 5.e-4);

  // Update the impact ionization parameters.
  // Selberherr model parameters from Silvaco Atlas.
  m_eImpactA = 1.889e5 * (1. + 0.588 * (t  - 1));
  m_hImpactA = 2.215e5 * (1. + 0.588 * (t  - 1));
  m_eImpactB = 5.75e5 * (1. + 0.248 * (t  - 1));
  m_hImpactB = 6.57e5 * (1. + 0.248 * (t  - 1));
}
}
