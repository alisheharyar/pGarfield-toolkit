#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumGaN.hh"
#include "Garfield/Random.hh"

namespace Garfield {

MediumGaN::MediumGaN() : Medium() {
  m_className = "MediumGaN";
  m_name = "GaN";

  // J. Wang et al., 
  // Review of using gallium nitride for ionizing radiation detection,
  // Appl. Phys. Rev. 2 (2015), 
  // http://dx.doi.org/10.1063/1.4929913

  SetTemperature(300.);
  SetDielectricConstant(8.9);
  SetAtomicNumber(19);
  SetAtomicWeight(41.865);
  SetMassDensity(6.15);

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = false;

  m_w = 8.9;
  m_fano = 0.1;
}

void MediumGaN::GetComponent(const unsigned int i, std::string& label,
                             double& f) {
  if (i == 0) {
    label = "Ga";
    f = 0.5;
  } else if (i == 1) {
    label = "N";
    f = 0.5;
  }
}

bool MediumGaN::ElectronVelocity(const double ex, const double ey,
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
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  constexpr double vsat = 1.27e-2;
  constexpr double ec = 172.e3;
  const double e0 = emag / ec;
  const double e1 = pow(e0, 4.19);
  const double den = 1. + e1 + 3.24 * pow(e0, 0.885);
  const double mu = -(m_eMobility + vsat * e1 / emag) / den; 
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

bool MediumGaN::ElectronTownsend(const double ex, const double ey,
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
    alpha = m_eImpactA * exp(-m_eImpactB / emag);
  } 
  return true;
}

bool MediumGaN::ElectronAttachment(const double ex, const double ey,
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

bool MediumGaN::HoleVelocity(const double ex, const double ey, const double ez,
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
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  // Values for saturation velocity and exponent from Sentaurus Synopsys. 
  constexpr double vsat = 7.e-3;
  constexpr double beta = 0.725;
  constexpr double invbeta = 1. / beta;
  const double r = m_hMobility * emag / vsat;
  const double mu = m_hMobility / pow(1. + pow(r, beta), invbeta);
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

bool MediumGaN::HoleTownsend(const double ex, const double ey, const double ez,
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
    alpha = m_hImpactA * exp(-m_hImpactB / emag);
  } 
  return true;
}

bool MediumGaN::HoleAttachment(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& eta) {
  eta = 0.;
  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return true;
}

void MediumGaN::SetElectronConcentration(const double c) {

  if (c < 0.) {
    std::cerr << m_className << "::SetElectronConcentration:\n"
              << "    Concentration cannot be negative.\n";
    return;
  }
  m_eDensity = c;
}

void MediumGaN::SetLowFieldMobility(const double mue, const double muh) {

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

void MediumGaN::UnsetLowFieldMobility() {
  m_userMobility = false;
  m_isChanged = true;
}

void MediumGaN::UpdateTransportParameters() {

  if (m_userMobility) return;
  const double t = m_temperature / 300.;

  // Electron low field mobility.
  // - F. Schwierz, Solid-State Electronics 49 (2005), 889
  //   https://doi.org/10.1016/j.sse.2005.03.006
  const double eMuMin = 0.080e-6 * pow(t, -0.2);
  const double eMuMax = 1.405e-6 * pow(t, -2.85);
  const double cRef = 7.78e16 * pow(t, 1.3);
  const double alpha = 0.71 * pow(t, 0.31);
  const double den = 1. + pow(m_eDensity / cRef, alpha);
  m_eMobility = eMuMin + (eMuMax - eMuMin) / den; 
  // Low-field mobility for holes (using only the lattice mobility).
  // - T. TMnatsakanov et al., Solid-State Electronics 47 (2003), 111
  //   https://doi.org/10.1016/S0038-1101(02)00256-3 
  m_hMobility = 0.170e-6 * pow(t, -5.);

}
}
