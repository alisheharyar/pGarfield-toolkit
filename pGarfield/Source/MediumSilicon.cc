#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Random.hh"
#include "Garfield/Utilities.hh"

namespace Garfield {

MediumSilicon::MediumSilicon()
    : Medium(),
      m_eStepXL(m_eFinalXL / nEnergyStepsXL),
      m_eStepG(m_eFinalG / nEnergyStepsG),
      m_eStepV(m_eFinalV / nEnergyStepsV) {
  m_className = "MediumSilicon";
  m_name = "Si";

  SetTemperature(300.);
  SetDielectricConstant(11.9);
  SetAtomicNumber(14.);
  SetAtomicWeight(28.0855);
  SetMassDensity(2.329);

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = true;

  m_w = 3.6;
  m_fano = 0.11;

  m_ieMinL = int(m_eMinL / m_eStepXL) + 1;
  m_ieMinG = int(m_eMinG / m_eStepG) + 1;

  // Load the density of states table.
  InitialiseDensityOfStates();
}

void MediumSilicon::SetDoping(const char type, const double c) {
  if (toupper(type) == 'N') {
    m_dopingType = 'n';
    if (c > Small) {
      m_dopingConcentration = c;
    } else {
      std::cerr << m_className << "::SetDoping:\n"
                << "    Doping concentration must be greater than zero.\n"
                << "    Using default value for n-type silicon "
                << "(10^12 cm-3) instead.\n";
      m_dopingConcentration = 1.e12;
    }
  } else if (toupper(type) == 'P') {
    m_dopingType = 'p';
    if (c > Small) {
      m_dopingConcentration = c;
    } else {
      std::cerr << m_className << "::SetDoping:\n"
                << "    Doping concentration must be greater than zero.\n"
                << "    Using default value for p-type silicon "
                << "(10^18 cm-3) instead.\n";
      m_dopingConcentration = 1.e18;
    }
  } else if (toupper(type) == 'I') {
    m_dopingType = 'i';
    m_dopingConcentration = 0.;
  } else {
    std::cerr << m_className << "::SetDoping:\n"
              << "    Unknown dopant type (" << type << ").\n"
              << "    Available types are n, p and i (intrinsic).\n";
    return;
  }

  m_isChanged = true;
}

void MediumSilicon::GetDoping(char& type, double& c) const {
  type = m_dopingType;
  c = m_dopingConcentration;
}

void MediumSilicon::SetTrapCrossSection(const double ecs, const double hcs) {
  if (ecs < 0.) {
    std::cerr << m_className << "::SetTrapCrossSection:\n"
              << "    Capture cross-section [cm2] must non-negative.\n";
  } else {
    m_eTrapCs = ecs;
  }

  if (hcs < 0.) {
    std::cerr << m_className << "::SetTrapCrossSection:\n"
              << "    Capture cross-section [cm2] must be non-negative.n";
  } else {
    m_hTrapCs = hcs;
  }

  m_trappingModel = 0;
  m_isChanged = true;
}

void MediumSilicon::SetTrapDensity(const double n) {
  if (n < 0.) {
    std::cerr << m_className << "::SetTrapDensity:\n"
              << "    Trap density [cm-3] must be non-negative.\n";
  } else {
    m_eTrapDensity = n;
    m_hTrapDensity = n;
  }

  m_trappingModel = 0;
  m_isChanged = true;
}

void MediumSilicon::SetTrappingTime(const double etau, const double htau) {
  if (etau <= 0.) {
    std::cerr << m_className << "::SetTrappingTime:\n"
              << "    Trapping time [ns-1] must be positive.\n";
  } else {
    m_eTrapTime = etau;
  }

  if (htau <= 0.) {
    std::cerr << m_className << "::SetTrappingTime:\n"
              << "    Trapping time [ns-1] must be positive.\n";
  } else {
    m_hTrapTime = htau;
  }

  m_trappingModel = 1;
  m_isChanged = true;
}

bool MediumSilicon::ElectronVelocity(const double ex, const double ey,
                                     const double ez, const double bx,
                                     const double by, const double bz,
                                     double& vx, double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronVelocity:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_eVelE.empty()) {
    // Interpolation in user table.
    return Medium::ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }

  // Calculate the mobility.
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  const double mu = -ElectronMobility(emag);

  if (fabs(bx) < Small && fabs(by) < Small && fabs(bz) < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, m_eHallFactor * mu, vx, vy, vz);
  }
  return true;
}

bool MediumSilicon::ElectronTownsend(const double ex, const double ey,
                                     const double ez, const double bx,
                                     const double by, const double bz,
                                     double& alpha) {
  alpha = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronTownsend:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_eAlp.empty()) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }

  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  alpha = ElectronAlpha(emag);
  return true;
}

bool MediumSilicon::ElectronAttachment(const double ex, const double ey,
                                       const double ez, const double bx,
                                       const double by, const double bz,
                                       double& eta) {
  eta = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronAttachment:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_eAtt.empty()) {
    // Interpolation in user table.
    return Medium::ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  }

  switch (m_trappingModel) {
    case 0:
      eta = m_eTrapCs * m_eTrapDensity;
      break;
    case 1:
      double vx, vy, vz;
      ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      eta = m_eTrapTime * sqrt(vx * vx + vy * vy + vz * vz);
      if (eta > 0.) eta = -1. / eta;
      break;
    default:
      std::cerr << m_className << "::ElectronAttachment: Unknown model. Bug!\n";
      return false;
  }

  return true;
}

bool MediumSilicon::HoleVelocity(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz, double& vx,
                                 double& vy, double& vz) {
  vx = vy = vz = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::HoleVelocity:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_hVelE.empty()) {
    // Interpolation in user table.
    return Medium::HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }

  // Calculate the mobility.
  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  const double mu = HoleMobility(emag);

  if (fabs(bx) < Small && fabs(by) < Small && fabs(bz) < Small) {
    vx = mu * ex;
    vy = mu * ey;
    vz = mu * ez;
  } else {
    Langevin(ex, ey, ez, bx, by, bz, mu, m_hHallFactor * mu, vx, vy, vz);
  }
  return true;
}

bool MediumSilicon::HoleTownsend(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz,
                                 double& alpha) {
  alpha = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::HoleTownsend:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_hAlp.empty()) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }

  const double emag = sqrt(ex * ex + ey * ey + ez * ez);
  alpha = HoleAlpha(emag);
  return true;
}

bool MediumSilicon::HoleAttachment(const double ex, const double ey,
                                   const double ez, const double bx,
                                   const double by, const double bz,
                                   double& eta) {
  eta = 0.;
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::HoleAttachment:\n"
                << "    Error calculating the transport parameters.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (!m_hAtt.empty()) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }

  switch (m_trappingModel) {
    case 0:
      eta = m_hTrapCs * m_hTrapDensity;
      break;
    case 1:
      double vx, vy, vz;
      HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      eta = m_hTrapTime * sqrt(vx * vx + vy * vy + vz * vz);
      if (eta > 0.) eta = -1. / eta;
      break;
    default:
      std::cerr << m_className << "::HoleAttachment: Unknown model. Bug!\n";
      return false;
  }

  return true;
}

void MediumSilicon::SetLowFieldMobility(const double mue, const double muh) {
  if (mue <= 0. || muh <= 0.) {
    std::cerr << m_className << "::SetLowFieldMobility:\n"
              << "    Mobility must be greater than zero.\n";
    return;
  }

  m_eMobility = mue;
  m_hMobility = muh;
  m_hasUserMobility = true;
  m_isChanged = true;
}

void MediumSilicon::SetLatticeMobilityModelMinimos() {
  m_latticeMobilityModel = LatticeMobility::Minimos;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetLatticeMobilityModelSentaurus() {
  m_latticeMobilityModel = LatticeMobility::Sentaurus;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetLatticeMobilityModelReggiani() {
  m_latticeMobilityModel = LatticeMobility::Reggiani;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetDopingMobilityModelMinimos() {
  m_dopingMobilityModel = DopingMobility::Minimos;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetDopingMobilityModelMasetti() {
  m_dopingMobilityModel = DopingMobility::Masetti;
  m_hasUserMobility = false;
  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocity(const double vsate,
                                          const double vsath) {
  if (vsate <= 0. || vsath <= 0.) {
    std::cout << m_className << "::SetSaturationVelocity:\n"
              << "    Restoring default values.\n";
    m_hasUserSaturationVelocity = false;
  } else {
    m_eSatVel = vsate;
    m_hSatVel = vsath;
    m_hasUserSaturationVelocity = true;
  }

  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocityModelMinimos() {
  m_saturationVelocityModel = SaturationVelocity::Minimos;
  m_hasUserSaturationVelocity = false;
  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocityModelCanali() {
  m_saturationVelocityModel = SaturationVelocity::Canali;
  m_hasUserSaturationVelocity = false;
  m_isChanged = true;
}

void MediumSilicon::SetSaturationVelocityModelReggiani() {
  m_saturationVelocityModel = SaturationVelocity::Reggiani;
  m_hasUserSaturationVelocity = false;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelMinimos() {
  m_highFieldMobilityModel = HighFieldMobility::Minimos;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelCanali() {
  m_highFieldMobilityModel = HighFieldMobility::Canali;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelReggiani() {
  m_highFieldMobilityModel = HighFieldMobility::Reggiani;
  m_isChanged = true;
}

void MediumSilicon::SetHighFieldMobilityModelConstant() {
  m_highFieldMobilityModel = HighFieldMobility::Constant;
}

void MediumSilicon::SetImpactIonisationModelVanOverstraetenDeMan() {
  m_impactIonisationModel = ImpactIonisation::VanOverstraeten;
  m_isChanged = true;
}

void MediumSilicon::SetImpactIonisationModelGrant() {
  m_impactIonisationModel = ImpactIonisation::Grant;
  m_isChanged = true;
}

void MediumSilicon::SetImpactIonisationModelMassey() {
  m_impactIonisationModel = ImpactIonisation::Massey;
  m_isChanged = true;
}

void MediumSilicon::SetImpactIonisationModelOkutoCrowell() {
  m_impactIonisationModel = ImpactIonisation::Okuto;
  m_isChanged = true;
}

bool MediumSilicon::SetMaxElectronEnergy(const double e) {
  if (e <= m_eMinG + Small) {
    std::cerr << m_className << "::SetMaxElectronEnergy:\n"
              << "    Requested upper electron energy limit (" << e
              << " eV) is too small.\n";
    return false;
  }

  m_eFinalG = e;
  // Determine the energy interval size.
  m_eStepG = m_eFinalG / nEnergyStepsG;

  m_isChanged = true;

  return true;
}

double MediumSilicon::GetElectronEnergy(const double px, const double py,
                                        const double pz, double& vx, double& vy,
                                        double& vz, const int band) {
  // Effective masses
  double mx = ElectronMass, my = ElectronMass, mz = ElectronMass;
  // Energy offset
  double e0 = 0.;
  if (band >= 0 && band < m_nValleysX) {
    // X valley
    if (m_anisotropic) {
      switch (band) {
        case 0:
        case 1:
          // X 100, -100
          mx *= m_mLongX;
          my *= m_mTransX;
          mz *= m_mTransX;
          break;
        case 2:
        case 3:
          // X 010, 0-10
          mx *= m_mTransX;
          my *= m_mLongX;
          mz *= m_mTransX;
          break;
        case 4:
        case 5:
          // X 001, 00-1
          mx *= m_mTransX;
          my *= m_mTransX;
          mz *= m_mLongX;
          break;
        default:
          std::cerr << m_className << "::GetElectronEnergy:\n"
                    << "    Unexpected band index " << band << "!\n";
          break;
      }
    } else {
      // Conduction effective mass
      const double mc = 3. / (1. / m_mLongX + 2. / m_mTransX);
      mx *= mc;
      my *= mc;
      mz *= mc;
    }
  } else if (band < m_nValleysX + m_nValleysL) {
    // L valley, isotropic approximation
    e0 = m_eMinL;
    // Effective mass
    const double mc = 3. / (1. / m_mLongL + 2. / m_mTransL);
    mx *= mc;
    my *= mc;
    mz *= mc;
  } else if (band == m_nValleysX + m_nValleysL) {
    // Higher band(s)
  }

  if (m_nonParabolic) {
    // Non-parabolicity parameter
    double alpha = 0.;
    if (band < m_nValleysX) {
      // X valley
      alpha = m_alphaX;
    } else if (band < m_nValleysX + m_nValleysL) {
      // L valley
      alpha = m_alphaL;
    }

    const double p2 = 0.5 * (px * px / mx + py * py / my + pz * pz / mz);
    if (alpha > 0.) {
      const double e = 0.5 * (sqrt(1. + 4 * alpha * p2) - 1.) / alpha;
      const double a = SpeedOfLight / (1. + 2 * alpha * e);
      vx = a * px / mx;
      vy = a * py / my;
      vz = a * pz / mz;
      return e0 + e;
    }
  }

  const double e = 0.5 * (px * px / mx + py * py / my + pz * pz / mz);
  vx = SpeedOfLight * px / mx;
  vy = SpeedOfLight * py / my;
  vz = SpeedOfLight * pz / mz;
  return e0 + e;
}

void MediumSilicon::GetElectronMomentum(const double e, double& px, double& py,
                                        double& pz, int& band) {
  int nBands = m_nValleysX;
  if (e > m_eMinL) nBands += m_nValleysL;
  if (e > m_eMinG) ++nBands;

  // If the band index is out of range, choose one at random.
  if (band < 0 || band > m_nValleysX + m_nValleysL ||
      (e < m_eMinL || band >= m_nValleysX) ||
      (e < m_eMinG || band == m_nValleysX + m_nValleysL)) {
    if (e < m_eMinL) {
      band = int(m_nValleysX * RndmUniform());
      if (band >= m_nValleysX) band = m_nValleysX - 1;
    } else {
      const double dosX = GetConductionBandDensityOfStates(e, 0);
      const double dosL = GetConductionBandDensityOfStates(e, m_nValleysX);
      const double dosG =
          GetConductionBandDensityOfStates(e, m_nValleysX + m_nValleysL);
      const double dosSum = m_nValleysX * dosX + m_nValleysL * dosL + dosG;
      if (dosSum < Small) {
        band = m_nValleysX + m_nValleysL;
      } else {
        const double r = RndmUniform() * dosSum;
        if (r < dosX) {
          band = int(m_nValleysX * RndmUniform());
          if (band >= m_nValleysX) band = m_nValleysX - 1;
        } else if (r < dosX + dosL) {
          band = m_nValleysX + int(m_nValleysL * RndmUniform());
          if (band >= m_nValleysX + m_nValleysL) band = m_nValleysL - 1;
        } else {
          band = m_nValleysX + m_nValleysL;
        }
      }
    }
    if (m_debug) {
      std::cout << m_className << "::GetElectronMomentum:\n"
                << "    Randomised band index: " << band << "\n";
    }
  }
  if (band < m_nValleysX) {
    // X valleys
    double pstar = sqrt(2. * ElectronMass * e);
    if (m_nonParabolic) {
      const double alpha = m_alphaX;
      pstar *= sqrt(1. + alpha * e);
    }

    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    const double phi = TwoPi * RndmUniform();

    if (m_anisotropic) {
      const double pl = pstar * sqrt(m_mLongX);
      const double pt = pstar * sqrt(m_mTransX);
      switch (band) {
        case 0:
        case 1:
          // 100
          px = pl * ctheta;
          py = pt * cos(phi) * stheta;
          pz = pt * sin(phi) * stheta;
          break;
        case 2:
        case 3:
          // 010
          px = pt * sin(phi) * stheta;
          py = pl * ctheta;
          pz = pt * cos(phi) * stheta;
          break;
        case 4:
        case 5:
          // 001
          px = pt * cos(phi) * stheta;
          py = pt * sin(phi) * stheta;
          pz = pl * ctheta;
          break;
        default:
          // Other band; should not occur.
          std::cerr << m_className << "::GetElectronMomentum:\n"
                    << "    Unexpected band index (" << band << ").\n";
          px = pstar * stheta * cos(phi);
          py = pstar * stheta * sin(phi);
          pz = pstar * ctheta;
          break;
      }
    } else {
      pstar *= sqrt(3. / (1. / m_mLongX + 2. / m_mTransX));
      px = pstar * cos(phi) * stheta;
      py = pstar * sin(phi) * stheta;
      pz = pstar * ctheta;
    }
  } else if (band < m_nValleysX + m_nValleysL) {
    // L valleys
    double pstar = sqrt(2. * ElectronMass * (e - m_eMinL));
    if (m_nonParabolic) {
      const double alpha = m_alphaL;
      pstar *= sqrt(1. + alpha * (e - m_eMinL));
    }
    pstar *= sqrt(3. / (1. / m_mLongL + 2. / m_mTransL));
    RndmDirection(px, py, pz, pstar);
  } else if (band == m_nValleysX + m_nValleysL) {
    // Higher band
    const double pstar = sqrt(2. * ElectronMass * e);
    RndmDirection(px, py, pz, pstar);
  }
}

double MediumSilicon::GetElectronNullCollisionRate(const int band) {
  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::GetElectronNullCollisionRate:\n"
                << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }

  if (band >= 0 && band < m_nValleysX) {
    return m_cfNullElectronsX;
  } else if (band >= m_nValleysX && band < m_nValleysX + m_nValleysL) {
    return m_cfNullElectronsL;
  } else if (band == m_nValleysX + m_nValleysL) {
    return m_cfNullElectronsG;
  }
  std::cerr << m_className << "::GetElectronNullCollisionRate:\n"
            << "    Band index (" << band << ") out of range.\n";
  return 0.;
}

double MediumSilicon::GetElectronCollisionRate(const double e, const int band) {
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n"
              << "    Electron energy must be positive.\n";
    return 0.;
  }

  if (e > m_eFinalG) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n"
              << "    Collision rate at " << e << " eV (band " << band
              << ") is not included in the current table.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::GetElectronCollisionRate:\n"
                << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }

  if (band >= 0 && band < m_nValleysX) {
    int iE = int(e / m_eStepXL);
    if (iE >= nEnergyStepsXL)
      iE = nEnergyStepsXL - 1;
    else if (iE < 0)
      iE = 0;
    return m_cfTotElectronsX[iE];
  } else if (band >= m_nValleysX && band < m_nValleysX + m_nValleysL) {
    int iE = int(e / m_eStepXL);
    if (iE >= nEnergyStepsXL)
      iE = nEnergyStepsXL - 1;
    else if (iE < m_ieMinL)
      iE = m_ieMinL;
    return m_cfTotElectronsL[iE];
  } else if (band == m_nValleysX + m_nValleysL) {
    int iE = int(e / m_eStepG);
    if (iE >= nEnergyStepsG)
      iE = nEnergyStepsG - 1;
    else if (iE < m_ieMinG)
      iE = m_ieMinG;
    return m_cfTotElectronsG[iE];
  }

  std::cerr << m_className << "::GetElectronCollisionRate:\n"
            << "    Band index (" << band << ") out of range.\n";
  return 0.;
}

bool MediumSilicon::ElectronCollision(const double e, int& type, 
    int& level, double& e1, double& px, double& py, double& pz, 
    std::vector<std::pair<Particle, double> >& secondaries, int& ndxc,
    int& band) {
  if (e > m_eFinalG) {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Requested electron energy (" << e << " eV) exceeds the "
              << "current energy range (" << m_eFinalG << " eV).\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  } else if (e <= 0.) {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Electron energy must be greater than zero.\n";
    return false;
  }

  if (m_isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << m_className << "::ElectronCollision:\n"
                << "    Error calculating the collision rates table.\n";
      return false;
    }
    m_isChanged = false;
  }

  // Energy loss
  double loss = 0.;
  // Sample the scattering process.
  if (band >= 0 && band < m_nValleysX) {
    // X valley
    // Get the energy interval.
    int iE = int(e / m_eStepXL);
    if (iE >= nEnergyStepsXL) iE = nEnergyStepsXL - 1;
    if (iE < 0) iE = 0;
    // Select the scattering process.
    const double r = RndmUniform();
    if (r <= m_cfElectronsX[iE][0]) {
      level = 0;
    } else if (r >= m_cfElectronsX[iE][m_nLevelsX - 1]) {
      level = m_nLevelsX - 1;
    } else {
      const auto begin = m_cfElectronsX[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nLevelsX, r) - begin;
    }

    // Get the collision type.
    type = m_scatTypeElectronsX[level];
    // Fill the collision counters.
    ++m_nCollElectronDetailed[level];
    ++m_nCollElectronBand[band];
    if (type == ElectronCollisionTypeAcousticPhonon) {
      ++m_nCollElectronAcoustic;
    } else if (type == ElectronCollisionTypeIntervalleyG) {
      // Intervalley scattering (g type)
      ++m_nCollElectronIntervalley;
      // Final valley is in opposite direction.
      switch (band) {
        case 0:
          band = 1;
          break;
        case 1:
          band = 0;
          break;
        case 2:
          band = 3;
          break;
        case 3:
          band = 2;
          break;
        case 4:
          band = 5;
          break;
        case 5:
          band = 4;
          break;
        default:
          break;
      }
    } else if (type == ElectronCollisionTypeIntervalleyF) {
      // Intervalley scattering (f type)
      ++m_nCollElectronIntervalley;
      // Final valley is perpendicular.
      switch (band) {
        case 0:
        case 1:
          band = int(RndmUniform() * 4) + 2;
          break;
        case 2:
        case 3:
          band = int(RndmUniform() * 4);
          if (band > 1) band += 2;
          break;
        case 4:
        case 5:
          band = int(RndmUniform() * 4);
          break;
      }
    } else if (type == ElectronCollisionTypeInterbandXL) {
      // XL scattering
      ++m_nCollElectronIntervalley;
      // Final valley is in L band.
      band = m_nValleysX + int(RndmUniform() * m_nValleysL);
      if (band >= m_nValleysX + m_nValleysL)
        band = m_nValleysX + m_nValleysL - 1;
    } else if (type == ElectronCollisionTypeInterbandXG) {
      ++m_nCollElectronIntervalley;
      band = m_nValleysX + m_nValleysL;
    } else if (type == ElectronCollisionTypeImpurity) {
      ++m_nCollElectronImpurity;
    } else if (type == ElectronCollisionTypeIonisation) {
      ++m_nCollElectronIonisation;
    } else {
      std::cerr << m_className << "::ElectronCollision:\n";
      std::cerr << "    Unexpected collision type (" << type << ").\n";
    }

    // Get the energy loss.
    loss = m_energyLossElectronsX[level];

  } else if (band >= m_nValleysX && band < m_nValleysX + m_nValleysL) {
    // L valley
    // Get the energy interval.
    int iE = int(e / m_eStepXL);
    if (iE >= nEnergyStepsXL) iE = nEnergyStepsXL - 1;
    if (iE < m_ieMinL) iE = m_ieMinL;
    // Select the scattering process.
    const double r = RndmUniform();
    if (r <= m_cfElectronsL[iE][0]) {
      level = 0;
    } else if (r >= m_cfElectronsL[iE][m_nLevelsL - 1]) {
      level = m_nLevelsL - 1;
    } else {
      const auto begin = m_cfElectronsL[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nLevelsL, r) - begin;
    }

    // Get the collision type.
    type = m_scatTypeElectronsL[level];
    // Fill the collision counters.
    ++m_nCollElectronDetailed[m_nLevelsX + level];
    ++m_nCollElectronBand[band];
    if (type == ElectronCollisionTypeAcousticPhonon) {
      ++m_nCollElectronAcoustic;
    } else if (type == ElectronCollisionTypeOpticalPhonon) {
      ++m_nCollElectronOptical;
    } else if (type == ElectronCollisionTypeIntervalleyG ||
               type == ElectronCollisionTypeIntervalleyF) {
      // Equivalent intervalley scattering
      ++m_nCollElectronIntervalley;
      // Randomise the final valley.
      band = m_nValleysX + int(RndmUniform() * m_nValleysL);
      if (band >= m_nValleysX + m_nValleysL) band = m_nValleysX + m_nValleysL;
    } else if (type == ElectronCollisionTypeInterbandXL) {
      // LX scattering
      ++m_nCollElectronIntervalley;
      // Randomise the final valley.
      band = int(RndmUniform() * m_nValleysX);
      if (band >= m_nValleysX) band = m_nValleysX - 1;
    } else if (type == ElectronCollisionTypeInterbandLG) {
      // LG scattering
      ++m_nCollElectronIntervalley;
      band = m_nValleysX + m_nValleysL;
    } else if (type == ElectronCollisionTypeImpurity) {
      ++m_nCollElectronImpurity;
    } else if (type == ElectronCollisionTypeIonisation) {
      ++m_nCollElectronIonisation;
    } else {
      std::cerr << m_className << "::ElectronCollision:\n"
                << "    Unexpected collision type (" << type << ").\n";
    }

    // Get the energy loss.
    loss = m_energyLossElectronsL[level];
  } else if (band == m_nValleysX + m_nValleysL) {
    // Higher bands
    // Get the energy interval.
    int iE = int(e / m_eStepG);
    if (iE >= nEnergyStepsG) iE = nEnergyStepsG - 1;
    if (iE < m_ieMinG) iE = m_ieMinG;
    // Select the scattering process.
    const double r = RndmUniform();
    if (r <= m_cfElectronsG[iE][0]) {
      level = 0;
    } else if (r >= m_cfElectronsG[iE][m_nLevelsG - 1]) {
      level = m_nLevelsG - 1;
    } else {
      const auto begin = m_cfElectronsG[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nLevelsG, r) - begin;
    }

    // Get the collision type.
    type = m_scatTypeElectronsG[level];
    // Fill the collision counters.
    ++m_nCollElectronDetailed[m_nLevelsX + m_nLevelsL + level];
    ++m_nCollElectronBand[band];
    if (type == ElectronCollisionTypeAcousticPhonon) {
      ++m_nCollElectronAcoustic;
    } else if (type == ElectronCollisionTypeOpticalPhonon) {
      ++m_nCollElectronOptical;
    } else if (type == ElectronCollisionTypeIntervalleyG ||
               type == ElectronCollisionTypeIntervalleyF) {
      // Equivalent intervalley scattering
      ++m_nCollElectronIntervalley;
    } else if (type == ElectronCollisionTypeInterbandXG) {
      // GX scattering
      ++m_nCollElectronIntervalley;
      // Randomise the final valley.
      band = int(RndmUniform() * m_nValleysX);
      if (band >= m_nValleysX) band = m_nValleysX - 1;
    } else if (type == ElectronCollisionTypeInterbandLG) {
      // GL scattering
      ++m_nCollElectronIntervalley;
      // Randomise the final valley.
      band = m_nValleysX + int(RndmUniform() * m_nValleysL);
      if (band >= m_nValleysX + m_nValleysL)
        band = m_nValleysX + m_nValleysL - 1;
    } else if (type == ElectronCollisionTypeIonisation) {
      ++m_nCollElectronIonisation;
    } else {
      std::cerr << m_className << "::ElectronCollision:\n"
                << "    Unexpected collision type (" << type << ").\n";
    }

    // Get the energy loss.
    loss = m_energyLossElectronsG[level];
  } else {
    std::cerr << m_className << "::ElectronCollision:\n"
              << "    Band index (" << band << ") out of range.\n";
    return false;
  }

  // Secondaries
  ndxc = 0;
  // Ionising collision
  if (type == ElectronCollisionTypeIonisation) {
    double ee = 0., eh = 0.;
    ComputeSecondaries(e, ee, eh);
    loss = ee + eh + m_bandGap;
    // Add the secondary electron.
    secondaries.emplace_back(std::make_pair(Particle::Electron, ee));
    // Add the hole.
    secondaries.emplace_back(std::make_pair(Particle::Hole, eh));
  }

  if (e < loss) loss = e - 0.00001;
  // Update the energy.
  e1 = e - loss;
  if (e1 < Small) e1 = Small;

  // Update the momentum.
  if (band >= 0 && band < m_nValleysX) {
    // X valleys
    double pstar = sqrt(2. * ElectronMass * e1);
    if (m_nonParabolic) {
      const double alpha = m_alphaX;
      pstar *= sqrt(1. + alpha * e1);
    }

    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    const double phi = TwoPi * RndmUniform();

    if (m_anisotropic) {
      const double pl = pstar * sqrt(m_mLongX);
      const double pt = pstar * sqrt(m_mTransX);
      switch (band) {
        case 0:
        case 1:
          // 100
          px = pl * ctheta;
          py = pt * cos(phi) * stheta;
          pz = pt * sin(phi) * stheta;
          break;
        case 2:
        case 3:
          // 010
          px = pt * sin(phi) * stheta;
          py = pl * ctheta;
          pz = pt * cos(phi) * stheta;
          break;
        case 4:
        case 5:
          // 001
          px = pt * cos(phi) * stheta;
          py = pt * sin(phi) * stheta;
          pz = pl * ctheta;
          break;
        default:
          return false;
      }
    } else {
      pstar *= sqrt(3. / (1. / m_mLongX + 2. / m_mTransX));
      px = pstar * cos(phi) * stheta;
      py = pstar * sin(phi) * stheta;
      pz = pstar * ctheta;
    }
    return true;

  } else if (band >= m_nValleysX && band < m_nValleysX + m_nValleysL) {
    // L valleys
    double pstar = sqrt(2. * ElectronMass * (e1 - m_eMinL));
    if (m_nonParabolic) {
      const double alpha = m_alphaL;
      pstar *= sqrt(1. + alpha * (e1 - m_eMinL));
    }
    pstar *= sqrt(3. / (1. / m_mLongL + 2. / m_mTransL));
    RndmDirection(px, py, pz, pstar);
    return true;
  }
  const double pstar = sqrt(2. * ElectronMass * e1);
  RndmDirection(px, py, pz, pstar);
  return true;
}

void MediumSilicon::ResetCollisionCounters() {
  m_nCollElectronAcoustic = m_nCollElectronOptical = 0;
  m_nCollElectronIntervalley = 0;
  m_nCollElectronImpurity = 0;
  m_nCollElectronIonisation = 0;
  const auto nLevels = m_nLevelsX + m_nLevelsL + m_nLevelsG;
  m_nCollElectronDetailed.assign(nLevels, 0);
  const auto nBands = m_nValleysX + m_nValleysL + 1;
  m_nCollElectronBand.assign(nBands, 0);
}

unsigned int MediumSilicon::GetNumberOfElectronCollisions() const {
  return m_nCollElectronAcoustic + m_nCollElectronOptical +
         m_nCollElectronIntervalley + m_nCollElectronImpurity +
         m_nCollElectronIonisation;
}

unsigned int MediumSilicon::GetNumberOfLevels() const {
  return m_nLevelsX + m_nLevelsL + m_nLevelsG;
}

unsigned int MediumSilicon::GetNumberOfElectronCollisions(
    const unsigned int level) const {
  const unsigned int nLevels = m_nLevelsX + m_nLevelsL + m_nLevelsG;
  if (level >= nLevels) {
    std::cerr << m_className << "::GetNumberOfElectronCollisions:\n"
              << "    Scattering rate term (" << level << ") does not exist.\n";
    return 0;
  }
  return m_nCollElectronDetailed[level];
}

unsigned int MediumSilicon::GetNumberOfElectronBands() const {
  return m_nValleysX + m_nValleysL + 1;
}

int MediumSilicon::GetElectronBandPopulation(const int band) {
  const int nBands = m_nValleysX + m_nValleysL + 1;
  if (band < 0 || band >= nBands) {
    std::cerr << m_className << "::GetElectronBandPopulation:\n";
    std::cerr << "    Band index (" << band << ") out of range.\n";
    return 0;
  }
  return m_nCollElectronBand[band];
}

bool MediumSilicon::GetOpticalDataRange(double& emin, double& emax,
                                        const unsigned int i) {
  if (i != 0) {
    std::cerr << m_className << "::GetOpticalDataRange: Index out of range.\n";
  }

  // Make sure the optical data table has been loaded.
  if (m_opticalDataEnergies.empty()) {
    if (!LoadOpticalData(m_opticalDataFile)) {
      std::cerr << m_className << "::GetOpticalDataRange:\n"
                << "    Optical data table could not be loaded.\n";
      return false;
    }
  }

  emin = m_opticalDataEnergies.front();
  emax = m_opticalDataEnergies.back();
  if (m_debug) {
    std::cout << m_className << "::GetOpticalDataRange:\n"
              << "    " << emin << " < E [eV] < " << emax << "\n";
  }
  return true;
}

bool MediumSilicon::GetDielectricFunction(const double e, double& eps1,
                                          double& eps2, const unsigned int i) {
  if (i != 0) {
    std::cerr << m_className + "::GetDielectricFunction: Index out of range.\n";
    return false;
  }

  // Make sure the optical data table has been loaded.
  if (m_opticalDataEnergies.empty()) {
    if (!LoadOpticalData(m_opticalDataFile)) {
      std::cerr << m_className << "::GetDielectricFunction:\n";
      std::cerr << "    Optical data table could not be loaded.\n";
      return false;
    }
  }

  // Make sure the requested energy is within the range of the table.
  const double emin = m_opticalDataEnergies.front();
  const double emax = m_opticalDataEnergies.back();
  if (e < emin || e > emax) {
    std::cerr << m_className << "::GetDielectricFunction:\n"
              << "    Requested energy (" << e << " eV) "
              << " is outside the range of the optical data table.\n"
              << "    " << emin << " < E [eV] < " << emax << "\n";
    eps1 = eps2 = 0.;
    return false;
  }

  // Locate the requested energy in the table.
  const auto begin = m_opticalDataEnergies.cbegin();
  const auto it1 = std::upper_bound(begin, m_opticalDataEnergies.cend(), e);
  if (it1 == begin) {
    eps1 = m_opticalDataEpsilon.front().first;
    eps2 = m_opticalDataEpsilon.front().second;
    return true;
  }
  const auto it0 = std::prev(it1);

  // Interpolate the real part of dielectric function.
  const double x0 = *it0;
  const double x1 = *it1;
  const double lnx0 = log(*it0);
  const double lnx1 = log(*it1);
  const double lnx = log(e);
  const double y0 = m_opticalDataEpsilon[it0 - begin].first;
  const double y1 = m_opticalDataEpsilon[it1 - begin].first;
  if (y0 <= 0. || y1 <= 0.) {
    // Use linear interpolation if one of the values is negative.
    eps1 = y0 + (e - x0) * (y1 - y0) / (x1 - x0);
  } else {
    // Otherwise use log-log interpolation.
    const double lny0 = log(y0);
    const double lny1 = log(y1);
    eps1 = lny0 + (lnx - lnx0) * (lny1 - lny0) / (lnx1 - lnx0);
    eps1 = exp(eps1);
  }

  // Interpolate the imaginary part of dielectric function,
  // using log-log interpolation.
  const double lnz0 = log(m_opticalDataEpsilon[it0 - begin].second);
  const double lnz1 = log(m_opticalDataEpsilon[it1 - begin].second);
  eps2 = lnz0 + (lnx - lnx0) * (lnz1 - lnz0) / (lnx1 - lnx0);
  eps2 = exp(eps2);
  return true;
}

bool MediumSilicon::Initialise() {
  if (!m_isChanged) {
    if (m_debug) {
      std::cerr << m_className << "::Initialise: Nothing changed.\n";
    }
    return true;
  }
  if (!UpdateTransportParameters()) {
    std::cerr << m_className << "::Initialise:\n    Error preparing "
              << "transport parameters/calculating collision rates.\n";
    return false;
  }
  return true;
}

bool MediumSilicon::UpdateTransportParameters() {
  std::lock_guard<std::mutex> guard(m_mutex);

  // Calculate impact ionisation coefficients.
  UpdateImpactIonisation();

  if (!m_hasUserMobility) {
    // Calculate lattice mobility.
    UpdateLatticeMobility();

    // Calculate doping mobility
    switch (m_dopingMobilityModel) {
      case DopingMobility::Minimos:
        UpdateDopingMobilityMinimos();
        break;
      case DopingMobility::Masetti:
        UpdateDopingMobilityMasetti();
        break;
      default:
        std::cerr << m_className << "::UpdateTransportParameters:\n    "
                  << "Unknown doping mobility model. Program bug!\n";
        break;
    }
  }

  // Calculate saturation velocity
  if (!m_hasUserSaturationVelocity) {
    UpdateSaturationVelocity();
  }

  // Calculate high field saturation parameters
  if (m_highFieldMobilityModel == HighFieldMobility::Canali) {
    UpdateHighFieldMobilityCanali();
  }

  if (m_debug) {
    std::cout << m_className << "::UpdateTransportParameters:\n"
              << "    Low-field mobility [cm2 V-1 ns-1]\n"
              << "      Electrons: " << m_eMobility << "\n"
              << "      Holes:     " << m_hMobility << "\n";
    if (m_highFieldMobilityModel == HighFieldMobility::Constant) {
      std::cout << "    Mobility is not field-dependent.\n";
    } else {
      std::cout << "    Saturation velocity [cm / ns]\n"
                << "      Electrons: " << m_eSatVel << "\n"
                << "      Holes:     " << m_hSatVel << "\n";
    }
  }

  if (!ElectronScatteringRates()) return false;
  if (!HoleScatteringRates()) return false;

  // Reset the collision counters.
  ResetCollisionCounters();

  return true;
}

void MediumSilicon::UpdateLatticeMobility() {

  // Temperature normalized to 300 K
  const double t = m_temperature / 300.;

  switch (m_latticeMobilityModel) {
    case LatticeMobility::Minimos:
      // - S. Selberherr, W. Haensch, M. Seavey, J. Slotboom,
      //   Solid State Electronics 33 (1990), 1425
      // - Minimos 6.1 User's Guide (1999)
      m_eLatticeMobility = 1.43e-6 * pow(t, -2.);
      m_hLatticeMobility = 0.46e-6 * pow(t, -2.18);
      break;
    case LatticeMobility::Sentaurus:
      // - C. Lombardi et al., IEEE Trans. CAD 7 (1988), 1164
      // - Sentaurus Device User Guide (2007)
      m_eLatticeMobility = 1.417e-6 * pow(t, -2.5);
      m_hLatticeMobility = 0.4705e-6 * pow(t, -2.2);
      break;
    case LatticeMobility::Reggiani:
      // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
      m_eLatticeMobility = 1.320e-6 * pow(t, -2.);
      m_hLatticeMobility = 0.460e-6 * pow(t, -2.2);
      break;
    default:
      std::cerr << m_className << "::UpdateLatticeMobility:\n"
                << "    Unknown lattice mobility model. Program bug!\n";
      break;
  }
}

void MediumSilicon::UpdateDopingMobilityMinimos() {
  // References:
  // - S. Selberherr, W. Haensch, M. Seavey, J. Slotboom,
  //   Solid State Electronics 33 (1990), 1425-1436
  // - Minimos 6.1 User's Guide (1999)

  // Mobility reduction due to ionised impurity scattering
  // Surface term not taken into account
  double eMuMin = 0.080e-6;
  double hMuMin = 0.045e-6;
  if (m_temperature > 200.) {
    const double c0 = pow(m_temperature / 300., -0.45);
    eMuMin *= c0;
    hMuMin *= c0;
  } else {
    const double c0 = pow(2. / 3., -0.45) * pow(m_temperature / 200., -0.15);
    eMuMin *= c0;
    hMuMin *= c0;
  }
  const double c1 = pow(m_temperature / 300., 3.2);
  const double eRefC = 1.12e17 * c1;
  const double hRefC = 2.23e17 * c1;
  const double alpha = 0.72 * pow(m_temperature / 300., 0.065);
  // Assume impurity concentration equal to doping concentration
  m_eMobility = eMuMin +
                (m_eLatticeMobility - eMuMin) /
                    (1. + pow(m_dopingConcentration / eRefC, alpha));
  m_hMobility = hMuMin +
                (m_hLatticeMobility - hMuMin) /
                    (1. + pow(m_dopingConcentration / hRefC, alpha));
}

void MediumSilicon::UpdateDopingMobilityMasetti() {
  // Reference:
  // - G. Masetti, M. Severi, S. Solmi,
  //   IEEE Trans. Electron Devices 30 (1983), 764-769
  // - Sentaurus Device User Guide (2007)
  // - Minimos NT User Guide (2004)

  if (m_dopingConcentration < 1.e13) {
    m_eMobility = m_eLatticeMobility;
    m_hMobility = m_hLatticeMobility;
    return;
  }

  // Parameters adopted from Minimos NT User Guide
  constexpr double eMuMin1 = 0.0522e-6;
  constexpr double eMuMin2 = 0.0522e-6;
  constexpr double eMu1 = 0.0434e-6;
  constexpr double hMuMin1 = 0.0449e-6;
  constexpr double hMuMin2 = 0.;
  constexpr double hMu1 = 0.029e-6;
  constexpr double eCr = 9.68e16;
  constexpr double eCs = 3.42e20;
  constexpr double hCr = 2.23e17;
  constexpr double hCs = 6.10e20;
  constexpr double hPc = 9.23e16;
  constexpr double eAlpha = 0.68;
  constexpr double eBeta = 2.;
  constexpr double hAlpha = 0.719;
  constexpr double hBeta = 2.;

  m_eMobility = eMuMin1 +
                (m_eLatticeMobility - eMuMin2) /
                    (1. + pow(m_dopingConcentration / eCr, eAlpha)) -
                eMu1 / (1. + pow(eCs / m_dopingConcentration, eBeta));

  m_hMobility = hMuMin1 * exp(-hPc / m_dopingConcentration) +
                (m_hLatticeMobility - hMuMin2) /
                    (1. + pow(m_dopingConcentration / hCr, hAlpha)) -
                hMu1 / (1. + pow(hCs / m_dopingConcentration, hBeta));
}

void MediumSilicon::UpdateSaturationVelocity() {

  switch (m_saturationVelocityModel) {
    case SaturationVelocity::Minimos:
      // - R. Quay, C. Moglestue, V. Palankovski, S. Selberherr,
      //   Materials Science in Semiconductor Processing 3 (2000), 149
      // - Minimos NT User Guide (2004)
      m_eSatVel = 1.e-2 / (1. + 0.74 * (m_temperature / 300. - 1.));
      m_hSatVel = 0.704e-2 / (1. + 0.37 * (m_temperature / 300. - 1.));
      break;
    case SaturationVelocity::Canali:
      // - C. Canali, G. Majni, R. Minder, G. Ottaviani,
      //   IEEE Transactions on Electron Devices 22 (1975), 1045
      // - Sentaurus Device User Guide (2007)
      m_eSatVel = 1.07e-2 * pow(300. / m_temperature, 0.87);
      m_hSatVel = 8.37e-3 * pow(300. / m_temperature, 0.52);
      break;
    case SaturationVelocity::Reggiani:
      // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
      m_eSatVel = 1.470e-2 * sqrt(tanh(150. / m_temperature));
      m_hSatVel = 0.916e-2 * sqrt(tanh(300. / m_temperature));
      break;
    default:
      std::cerr << m_className << "::UpdateSaturationVelocity:\n" 
                << "    Unknown saturation velocity model. Program bug!\n";
      break;
  }
}

void MediumSilicon::UpdateHighFieldMobilityCanali() {
  // References:
  // - C. Canali, G. Majni, R. Minder, G. Ottaviani,
  //   IEEE Transactions on Electron Devices 22 (1975), 1045-1047
  // - Sentaurus Device User Guide (2007)

  // Temperature dependent exponent in high-field mobility formula
  m_eBetaCanali = 1.109 * pow(m_temperature / 300., 0.66);
  m_hBetaCanali = 1.213 * pow(m_temperature / 300., 0.17);
  m_eBetaCanaliInv = 1. / m_eBetaCanali;
  m_hBetaCanaliInv = 1. / m_hBetaCanali;
}

void MediumSilicon::UpdateImpactIonisation() {

  if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten ||
      m_impactIonisationModel == ImpactIonisation::Grant) {

    // Temperature dependence as in Sentaurus Device
    // Optical phonon energy
    constexpr double hbarOmega = 0.063;
    // Temperature scaling coefficient
    const double gamma =
        tanh(hbarOmega / (2. * BoltzmannConstant * 300.)) /
        tanh(hbarOmega / (2. * BoltzmannConstant * m_temperature));

    if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten) {
      // - R. van Overstraeten and H. de Man,
      //   Solid State Electronics 13 (1970), 583
      // - W. Maes, K. de Meyer and R. van Overstraeten,
      //   Solid State Electronics 33 (1990), 705
      // - Sentaurus Device User Guide (2016)

      // Low field coefficients taken from Maes, de Meyer, van Overstraeten
      // eImpactA0 = gamma * 3.318e5;
      // eImpactB0 = gamma * 1.135e6;
      m_eImpactA0 = gamma * 7.03e5;
      m_eImpactB0 = gamma * 1.231e6;
      m_eImpactA1 = gamma * 7.03e5;
      m_eImpactB1 = gamma * 1.231e6;

      m_hImpactA0 = gamma * 1.582e6;
      m_hImpactB0 = gamma * 2.036e6;
      m_hImpactA1 = gamma * 6.71e5;
      m_hImpactB1 = gamma * 1.693e6;
    } else if (m_impactIonisationModel == ImpactIonisation::Grant) {
      // W. N. Grant, Solid State Electronics 16 (1973), 1189
      // Sentaurus Device User Guide (2007)
      m_eImpactA0 = 2.60e6 * gamma;
      m_eImpactB0 = 1.43e6 * gamma;
      m_eImpactA1 = 6.20e5 * gamma;
      m_eImpactB1 = 1.08e6 * gamma;
      m_eImpactA2 = 5.05e5 * gamma;
      m_eImpactB2 = 9.90e5 * gamma;

      m_hImpactA0 = 2.00e6 * gamma;
      m_hImpactB0 = 1.97e6 * gamma;
      m_hImpactA1 = 5.60e5 * gamma;
      m_hImpactB1 = 1.32e6 * gamma;
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Massey) {
    // D. J. Massey, J. P. R. David, and G. J. Rees,
    // IEEE Trans. Electron Devices 53 (2006), 2328
    m_eImpactA0 = 4.43e5;
    m_eImpactB0 = 9.66e5 + 4.99e2 * m_temperature;

    m_hImpactA0 = 1.13e6;
    m_hImpactB0 = 1.71e6 + 1.09e3 * m_temperature;
  } else if (m_impactIonisationModel == ImpactIonisation::Okuto) {
    const double dt = m_temperature - 300.;
    m_eImpactA0 = 0.426 * (1. + 3.05e-4 * dt);
    m_hImpactA0 = 0.243 * (1. + 5.35e-4 * dt);
    m_eImpactB0 = 4.81e5 * (1. + 6.86e-4 * dt);
    m_hImpactB0 = 6.53e5 * (1. + 5.67e-4 * dt);
  } else {
    std::cerr << m_className << "::UpdateImpactIonisation:\n"
              << "    Unknown impact ionisation model. Program bug!\n";
  }
}

double MediumSilicon::ElectronMobility(const double emag) const {

  if (emag < Small) return 0.;

  if (m_highFieldMobilityModel == HighFieldMobility::Minimos) {
    // Minimos User's Guide (1999)
    const double r = 2 * m_eMobility * emag / m_eSatVel;
    return 2. * m_eMobility / (1. + sqrt(1. + r * r));
  } else if (m_highFieldMobilityModel == HighFieldMobility::Canali) {
    // Sentaurus Device User Guide (2007)
    const double r = m_eMobility * emag / m_eSatVel;
    return m_eMobility / pow(1. + pow(r, m_eBetaCanali), m_eBetaCanaliInv);
  } else if (m_highFieldMobilityModel == HighFieldMobility::Reggiani) {
    // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
    const double r = m_eMobility * emag / m_eSatVel;
    constexpr double k = 1. / 1.5;
    return m_eMobility / pow(1. + pow(r, 1.5), k);
  }
  return m_eMobility;
}

double MediumSilicon::ElectronAlpha(const double emag) const {

  if (emag < Small) return 0.;
  
  if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten) {
    // - R. van Overstraeten and H. de Man,
    //   Solid State Electronics 13 (1970), 583
    // - W. Maes, K. de Meyer and R. van Overstraeten,
    //   Solid State Electronics 33 (1990), 705
    // - Sentaurus Device User Guide (2016)
    if (emag < 4e5) {
      return m_eImpactA0 * exp(-m_eImpactB0 / emag);
    } else {
      return m_eImpactA1 * exp(-m_eImpactB1 / emag);
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Grant) {
    // W. N. Grant, Solid State Electronics 16 (1973), 1189
    if (emag < 2.4e5) {
      return m_eImpactA0 * exp(-m_eImpactB0 / emag);
    } else if (emag < 5.3e5) {
      return m_eImpactA1 * exp(-m_eImpactB1 / emag);
    } else {
      return m_eImpactA2 * exp(-m_eImpactB2 / emag);
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Massey) {
    return m_eImpactA0 * exp(-m_eImpactB0 / emag);
  } else if (m_impactIonisationModel == ImpactIonisation::Okuto) {
    const double f = m_eImpactB0 / emag;
    return m_eImpactA0 * emag * exp(-f * f);
  }
  std::cerr << m_className << "::ElectronAlpha: Unknown model. Program bug!\n";
  return 0.;
}

double MediumSilicon::HoleMobility(const double emag) const {

  if (emag < Small) return 0.;

  if (m_highFieldMobilityModel == HighFieldMobility::Minimos) {
    // Minimos User's Guide (1999)
    return m_hMobility / (1. + m_hMobility * emag / m_eSatVel);
  } else if (m_highFieldMobilityModel == HighFieldMobility::Canali) {
    // Sentaurus Device User Guide (2007)
    const double r = m_hMobility * emag / m_hSatVel;
    return m_hMobility / pow(1. + pow(r, m_hBetaCanali), m_hBetaCanaliInv);
  } else if (m_highFieldMobilityModel == HighFieldMobility::Reggiani) {
    // M. A. Omar, L. Reggiani, Solid State Electronics 30 (1987), 693
    const double r = m_hMobility * emag / m_hSatVel;
    return m_hMobility / sqrt(1. + r * r);
  }
  return m_hMobility;
}

double MediumSilicon::HoleAlpha(const double emag) const {

  if (emag < Small) return 0.;

  if (m_impactIonisationModel == ImpactIonisation::VanOverstraeten) {
    // - R. van Overstraeten and H. de Man,
    //   Solid State Electronics 13 (1970), 583
    // - Sentaurus Device User Guide (2016)
    if (emag < 4e5) {
      return m_hImpactA0 * exp(-m_hImpactB0 / emag);
    } else {
      return m_hImpactA1 * exp(-m_hImpactB1 / emag);
    }
  } else if (m_impactIonisationModel == ImpactIonisation::Grant) {
    // W. N. Grant, Solid State Electronics 16 (1973), 1189
    if (emag < 5.3e5) {
      return m_hImpactA0 * exp(-m_hImpactB0 / emag);
    } else {
      return m_hImpactA1 * exp(-m_hImpactB1 / emag);
    } 
  } else if (m_impactIonisationModel == ImpactIonisation::Massey) {
    return m_hImpactA0 * exp(-m_hImpactB0 / emag);
  } else if (m_impactIonisationModel == ImpactIonisation::Okuto) {
    const double f = m_hImpactB0 / emag;
    return m_hImpactA0 * emag * exp(-f * f);
  }
  std::cerr << m_className << "::HoleAlpha: Unknown model. Program bug!\n";
  return 0.;
}

bool MediumSilicon::LoadOpticalData(const std::string& filename) {
  // Clear the optical data table.
  m_opticalDataEnergies.clear();
  m_opticalDataEpsilon.clear();

  // Get the path to the data directory.
  char* pPath = getenv("GARFIELD_HOME");
  if (pPath == 0) {
    std::cerr << m_className << "::LoadOpticalData:\n"
              << "    Environment variable GARFIELD_HOME is not set.\n";
    return false;
  }
  const std::string filepath = std::string(pPath) + "/Data/" + filename;

  // Open the file.
  std::ifstream infile(filepath);
  // Make sure the file could actually be opened.
  if (!infile) {
    std::cerr << m_className << "::LoadOpticalData:\n"
              << "    Error opening file " << filename << ".\n";
    return false;
  }

  double lastEnergy = -1.;
  double energy, eps1, eps2, loss;
  // Read the file line by line.
  std::string line;
  std::istringstream dataStream;
  int i = 0;
  while (!infile.eof()) {
    ++i;
    // Read the next line.
    std::getline(infile, line);
    // Strip white space from the beginning of the line.
    ltrim(line);
    if (line.empty()) continue;
    // Skip comments.
    if (line[0] == '#' || line[0] == '*' || (line[0] == '/' && line[1] == '/'))
      continue;
    // Extract the values.
    dataStream.str(line);
    dataStream >> energy >> eps1 >> eps2 >> loss;
    if (dataStream.eof()) break;
    // Check if the data has been read correctly.
    if (infile.fail()) {
      std::cerr << m_className << "::LoadOpticalData:\n    Error reading file "
                << filename << " (line " << i << ").\n";
      return false;
    }
    // Reset the stringstream.
    dataStream.str("");
    dataStream.clear();
    // Make sure the values make sense.
    // The table has to be in ascending order
    //  with respect to the photon energy.
    if (energy <= lastEnergy) {
      std::cerr << m_className << "::LoadOpticalData:\n    Table is not in "
                << "monotonically increasing order (line " << i << ").\n"
                << "    " << lastEnergy << "  " << energy << "  " << eps1
                << "  " << eps2 << "\n";
      return false;
    }
    // The imaginary part of the dielectric function has to be positive.
    if (eps2 < 0.) {
      std::cerr << m_className << "::LoadOpticalData:\n    Negative value "
                << "of the loss function (line " << i << ").\n";
      return false;
    }
    // Ignore negative photon energies.
    if (energy <= 0.) continue;
    // Add the values to the list.
    m_opticalDataEnergies.emplace_back(energy);
    m_opticalDataEpsilon.emplace_back(std::make_pair(eps1, eps2));
    lastEnergy = energy;
  }

  if (m_opticalDataEnergies.empty()) {
    std::cerr << m_className << "::LoadOpticalData:\n"
              << "    Import of data from file " << filepath << "failed.\n"
              << "    No valid data found.\n";
    return false;
  }

  if (m_debug) {
    std::cout << m_className << "::LoadOpticalData:\n    Read "
              << m_opticalDataEnergies.size() << " values from file "
              << filepath << "\n";
  }
  return true;
}

bool MediumSilicon::ElectronScatteringRates() {
  // Reset the scattering rates
  m_cfTotElectronsX.assign(nEnergyStepsXL, 0.);
  m_cfTotElectronsL.assign(nEnergyStepsXL, 0.);
  m_cfTotElectronsG.assign(nEnergyStepsG, 0.);
  m_cfElectronsX.assign(nEnergyStepsXL, std::vector<double>());
  m_cfElectronsL.assign(nEnergyStepsXL, std::vector<double>());
  m_cfElectronsG.assign(nEnergyStepsG, std::vector<double>());
  m_energyLossElectronsX.clear();
  m_energyLossElectronsL.clear();
  m_energyLossElectronsG.clear();
  m_scatTypeElectronsX.clear();
  m_scatTypeElectronsL.clear();
  m_scatTypeElectronsG.clear();
  m_cfNullElectronsX = 0.;
  m_cfNullElectronsL = 0.;
  m_cfNullElectronsG = 0.;

  m_nLevelsX = 0;
  m_nLevelsL = 0;
  m_nLevelsG = 0;
  // Fill the scattering rate tables.
  ElectronAcousticScatteringRates();
  ElectronOpticalScatteringRates();
  ElectronImpurityScatteringRates();
  ElectronIntervalleyScatteringRatesXX();
  ElectronIntervalleyScatteringRatesXL();
  ElectronIntervalleyScatteringRatesLL();
  ElectronIntervalleyScatteringRatesXGLG();
  ElectronIonisationRatesXL();
  ElectronIonisationRatesG();

  if (m_debug) {
    std::cout << m_className << "::ElectronScatteringRates:\n"
              << "    " << m_nLevelsX << " X-valley scattering terms\n"
              << "    " << m_nLevelsL << " L-valley scattering terms\n"
              << "    " << m_nLevelsG << " higher band scattering terms\n";
  }

  std::ofstream outfileX;
  std::ofstream outfileL;
  if (m_cfOutput) {
    outfileX.open("ratesX.txt", std::ios::out);
    outfileL.open("ratesL.txt", std::ios::out);
  }

  m_ieMinL = int(m_eMinL / m_eStepXL) + 1;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    // Sum up the scattering rates of all processes.
    for (int j = m_nLevelsX; j--;) m_cfTotElectronsX[i] += m_cfElectronsX[i][j];
    for (int j = m_nLevelsL; j--;) m_cfTotElectronsL[i] += m_cfElectronsL[i][j];

    if (m_cfOutput) {
      outfileX << i * m_eStepXL << " " << m_cfTotElectronsX[i] << " ";
      for (int j = 0; j < m_nLevelsX; ++j) {
        outfileX << m_cfElectronsX[i][j] << " ";
      }
      outfileX << "\n";
      outfileL << i * m_eStepXL << " " << m_cfTotElectronsL[i] << " ";
      for (int j = 0; j < m_nLevelsL; ++j) {
        outfileL << m_cfElectronsL[i][j] << " ";
      }
      outfileL << "\n";
    }

    if (m_cfTotElectronsX[i] > m_cfNullElectronsX) {
      m_cfNullElectronsX = m_cfTotElectronsX[i];
    }
    if (m_cfTotElectronsL[i] > m_cfNullElectronsL) {
      m_cfNullElectronsL = m_cfTotElectronsL[i];
    }

    // Make sure the total scattering rate is positive.
    if (m_cfTotElectronsX[i] <= 0.) {
      std::cerr << m_className << "::ElectronScatteringRates:\n    X-valley "
                << "scattering rate at " << i * m_eStepXL << " eV <= 0.\n";
      return false;
    }
    // Normalise the rates.
    for (int j = 0; j < m_nLevelsX; ++j) {
      m_cfElectronsX[i][j] /= m_cfTotElectronsX[i];
      if (j > 0) m_cfElectronsX[i][j] += m_cfElectronsX[i][j - 1];
    }

    // Make sure the total scattering rate is positive.
    if (m_cfTotElectronsL[i] <= 0.) {
      if (i < m_ieMinL) continue;
      std::cerr << m_className << "::ElectronScatteringRates:\n    L-valley "
                << "scattering rate at " << i * m_eStepXL << " eV <= 0.\n";
      return false;
    }
    // Normalise the rates.
    for (int j = 0; j < m_nLevelsL; ++j) {
      m_cfElectronsL[i][j] /= m_cfTotElectronsL[i];
      if (j > 0) m_cfElectronsL[i][j] += m_cfElectronsL[i][j - 1];
    }
  }

  if (m_cfOutput) {
    outfileX.close();
    outfileL.close();
  }

  std::ofstream outfileG;
  if (m_cfOutput) {
    outfileG.open("ratesG.txt", std::ios::out);
  }
  m_ieMinG = int(m_eMinG / m_eStepG) + 1;
  for (int i = 0; i < nEnergyStepsG; ++i) {
    // Sum up the scattering rates of all processes.
    for (int j = m_nLevelsG; j--;) m_cfTotElectronsG[i] += m_cfElectronsG[i][j];

    if (m_cfOutput) {
      outfileG << i * m_eStepG << " " << m_cfTotElectronsG[i] << " ";
      for (int j = 0; j < m_nLevelsG; ++j) {
        outfileG << m_cfElectronsG[i][j] << " ";
      }
      outfileG << "\n";
    }

    if (m_cfTotElectronsG[i] > m_cfNullElectronsG) {
      m_cfNullElectronsG = m_cfTotElectronsG[i];
    }

    // Make sure the total scattering rate is positive.
    if (m_cfTotElectronsG[i] <= 0.) {
      if (i < m_ieMinG) continue;
      std::cerr << m_className << "::ElectronScatteringRates:\n    Higher "
                << "band scattering rate at " << i * m_eStepG << " eV <= 0.\n";
    }
    // Normalise the rates.
    for (int j = 0; j < m_nLevelsG; ++j) {
      m_cfElectronsG[i][j] /= m_cfTotElectronsG[i];
      if (j > 0) m_cfElectronsG[i][j] += m_cfElectronsG[i][j - 1];
    }
  }

  if (m_cfOutput) {
    outfileG.close();
  }

  return true;
}

bool MediumSilicon::ElectronAcousticScatteringRates() {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Acoustic phonon intraband scattering
  // Acoustic deformation potential [eV]
  constexpr double defpot2 = 9. * 9.;
  // Longitudinal velocity of sound [cm/ns]
  constexpr double u = 9.04e-4;
  // Prefactor for acoustic deformation potential scattering
  const double cIntra = TwoPi * SpeedOfLight * SpeedOfLight * kbt * defpot2 /
                        (Hbar * u * u * rho);

  // Fill the scattering rate tables.
  double en = Small;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    const double dosX = GetConductionBandDensityOfStates(en, 0);
    const double dosL = GetConductionBandDensityOfStates(en, m_nValleysX);

    m_cfElectronsX[i].push_back(cIntra * dosX);
    m_cfElectronsL[i].push_back(cIntra * dosL);
    en += m_eStepXL;
  }

  en = Small;
  for (int i = 0; i < nEnergyStepsG; ++i) {
    const double dosG =
        GetConductionBandDensityOfStates(en, m_nValleysX + m_nValleysL);
    m_cfElectronsG[i].push_back(cIntra * dosG);
    en += m_eStepG;
  }

  // Assume that energy loss is negligible.
  m_energyLossElectronsX.push_back(0.);
  m_energyLossElectronsL.push_back(0.);
  m_energyLossElectronsG.push_back(0.);
  m_scatTypeElectronsX.push_back(ElectronCollisionTypeAcousticPhonon);
  m_scatTypeElectronsL.push_back(ElectronCollisionTypeAcousticPhonon);
  m_scatTypeElectronsG.push_back(ElectronCollisionTypeAcousticPhonon);
  ++m_nLevelsX;
  ++m_nLevelsL;
  ++m_nLevelsG;

  return true;
}

bool MediumSilicon::ElectronOpticalScatteringRates() {
  // Reference:
  //  - K. Hess (editor),
  //    Monte Carlo device simulation: full band and beyond
  //    Chapter 5
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - M. Lundstrom,
  //    Fundamentals of carrier transport

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Coupling constant [eV/cm]
  constexpr double dtk = 2.2e8;
  // Phonon energy [eV]
  constexpr double eph = 63.0e-3;
  // Phonon cccupation numbers
  const double nocc = 1. / (exp(eph / kbt) - 1);
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c = c0 * dtk * dtk / eph;

  double en = 0.;
  // L valleys
  /*
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    // Absorption
    if (en > m_eMinL) {
      double dos = GetConductionBandDensityOfStates(en + eph, m_nValleysX);
      m_cfElectronsL[i].push_back(c * nocc * dos);
    } else {
      m_cfElectronsL[i].push_back(0.);
    }
    // Emission
    if (en > m_eMinL + eph) {
      double dos = GetConductionBandDensityOfStates(en - eph, m_nValleysX);
      m_cfElectronsL[i].push_back(c * (nocc + 1) * dos);
    } else {
      m_cfElectronsL[i].push_back(0.);
    }
    en += m_eStepXL;
  }
  //*/

  en = 0.;
  // Higher band(s)
  for (int i = 0; i < nEnergyStepsG; ++i) {
    // Absorption
    if (en > m_eMinG) {
      double dos =
          GetConductionBandDensityOfStates(en + eph, m_nValleysX + m_nValleysL);
      m_cfElectronsG[i].push_back(c * nocc * dos);
    } else {
      m_cfElectronsG[i].push_back(0.);
    }
    // Emission
    if (en > m_eMinG + eph) {
      double dos =
          GetConductionBandDensityOfStates(en - eph, m_nValleysX + m_nValleysL);
      m_cfElectronsG[i].push_back(c * (nocc + 1) * dos);
    } else {
      m_cfElectronsG[i].push_back(0.);
    }
    en += m_eStepG;
  }

  // Absorption
  // m_energyLossElectronsL.push_back(-eph);
  m_energyLossElectronsG.push_back(-eph);
  // Emission
  // m_energyLossElectronsL.push_back(eph);
  m_energyLossElectronsG.push_back(eph);
  // m_scatTypeElectronsL.push_back(ElectronCollisionTypeOpticalPhonon);
  // m_scatTypeElectronsL.push_back(ElectronCollisionTypeOpticalPhonon);
  m_scatTypeElectronsG.push_back(ElectronCollisionTypeOpticalPhonon);
  m_scatTypeElectronsG.push_back(ElectronCollisionTypeOpticalPhonon);

  // m_nLevelsL += 2;
  m_nLevelsG += 2;

  return true;
}

bool MediumSilicon::ElectronIntervalleyScatteringRatesXX() {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  constexpr unsigned int nPhonons = 6;
  // f-type scattering: transition between orthogonal axes (multiplicity 4)
  // g-type scattering: transition between opposite axes (multiplicity 1)
  // Sequence of transitions in the table:
  // TA (g) - LA (g) - LO (g) - TA (f) - LA (f) - TO (f)
  // Coupling constants [eV/cm]
  constexpr double dtk[nPhonons] = {0.5e8, 0.8e8, 1.1e9, 0.3e8, 2.0e8, 2.0e8};
  // Phonon energies [eV]
  constexpr double eph[nPhonons] = {12.06e-3, 18.53e-3, 62.04e-3,
                                    18.86e-3, 47.39e-3, 59.03e-3};
  // Phonon cccupation numbers
  double nocc[nPhonons];
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c[nPhonons];

  for (unsigned int j = 0; j < nPhonons; ++j) {
    nocc[j] = 1. / (exp(eph[j] / kbt) - 1);
    c[j] = c0 * dtk[j] * dtk[j] / eph[j];
    if (j > 2) c[j] *= 4;
  }

  double en = 0.;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    for (unsigned int j = 0; j < nPhonons; ++j) {
      // Absorption
      double dos = GetConductionBandDensityOfStates(en + eph[j], 0);
      m_cfElectronsX[i].push_back(c[j] * nocc[j] * dos);
      // Emission
      if (en > eph[j]) {
        dos = GetConductionBandDensityOfStates(en - eph[j], 0);
        m_cfElectronsX[i].push_back(c[j] * (nocc[j] + 1) * dos);
      } else {
        m_cfElectronsX[i].push_back(0.);
      }
    }
    en += m_eStepXL;
  }

  for (unsigned int j = 0; j < nPhonons; ++j) {
    // Absorption
    m_energyLossElectronsX.push_back(-eph[j]);
    // Emission
    m_energyLossElectronsX.push_back(eph[j]);
    if (j <= 2) {
      m_scatTypeElectronsX.push_back(ElectronCollisionTypeIntervalleyG);
      m_scatTypeElectronsX.push_back(ElectronCollisionTypeIntervalleyG);
    } else {
      m_scatTypeElectronsX.push_back(ElectronCollisionTypeIntervalleyF);
      m_scatTypeElectronsX.push_back(ElectronCollisionTypeIntervalleyF);
    }
  }

  m_nLevelsX += 2 * nPhonons;

  return true;
}

bool MediumSilicon::ElectronIntervalleyScatteringRatesXL() {
  // Reference:
  // - M. Lundstrom, Fundamentals of carrier transport
  // - M. Martin et al.,
  //   Semicond. Sci. Technol. 8, 1291-1297

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  constexpr unsigned int nPhonons = 4;

  // Coupling constants [eV/cm]
  constexpr double dtk[nPhonons] = {2.e8, 2.e8, 2.e8, 2.e8};
  // Phonon energies [eV]
  constexpr double eph[nPhonons] = {58.e-3, 55.e-3, 41.e-3, 17.e-3};
  // Number of equivalent valleys
  constexpr unsigned int zX = 6;
  constexpr unsigned int zL = 8;

  // Phonon cccupation numbers
  double nocc[nPhonons] = {0.};
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c[nPhonons];

  for (unsigned int j = 0; j < nPhonons; ++j) {
    nocc[j] = 1. / (exp(eph[j] / kbt) - 1);
    c[j] = c0 * dtk[j] * dtk[j] / eph[j];
  }

  double en = 0.;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    for (unsigned int j = 0; j < nPhonons; ++j) {
      // XL
      // Absorption
      if (en + eph[j] > m_eMinL) {
        double dos = GetConductionBandDensityOfStates(en + eph[j], m_nValleysX);
        m_cfElectronsX[i].push_back(zL * c[j] * nocc[j] * dos);
      } else {
        m_cfElectronsX[i].push_back(0.);
      }
      // Emission
      if (en - eph[j] > m_eMinL) {
        double dos = GetConductionBandDensityOfStates(en - eph[j], m_nValleysX);
        m_cfElectronsX[i].push_back(zL * c[j] * (nocc[j] + 1) * dos);
      } else {
        m_cfElectronsX[i].push_back(0.);
      }
      // LX
      if (en > m_eMinL) {
        // Absorption
        double dos = GetConductionBandDensityOfStates(en + eph[j], 0);
        m_cfElectronsL[i].push_back(zX * c[j] * nocc[j] * dos);
        // Emission
        dos = GetConductionBandDensityOfStates(en - eph[j], 0);
        m_cfElectronsL[i].push_back(zX * c[j] * (nocc[j] + 1) * dos);
      } else {
        m_cfElectronsL[i].push_back(0.);
        m_cfElectronsL[i].push_back(0.);
      }
    }
    en += m_eStepXL;
  }

  for (unsigned int j = 0; j < nPhonons; ++j) {
    // Absorption
    m_energyLossElectronsX.push_back(-eph[j]);
    m_energyLossElectronsL.push_back(-eph[j]);
    // Emission
    m_energyLossElectronsX.push_back(eph[j]);
    m_energyLossElectronsL.push_back(eph[j]);
    m_scatTypeElectronsX.push_back(ElectronCollisionTypeInterbandXL);
    m_scatTypeElectronsX.push_back(ElectronCollisionTypeInterbandXL);
    m_scatTypeElectronsL.push_back(ElectronCollisionTypeInterbandXL);
    m_scatTypeElectronsL.push_back(ElectronCollisionTypeInterbandXL);
  }

  m_nLevelsX += 2 * nPhonons;
  m_nLevelsL += 2 * nPhonons;

  return true;
}

bool MediumSilicon::ElectronIntervalleyScatteringRatesLL() {
  // Reference:
  //  - K. Hess (editor),
  //    Monte Carlo device simulation: full band and beyond
  //    Chapter 5
  //  - M. J. Martin et al.,
  //    Semicond. Sci. Technol. 8, 1291-1297

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  const int nPhonons = 1;
  // Coupling constant [eV/cm]
  const double dtk[nPhonons] = {2.63e8};
  // Phonon energy [eV]
  const double eph[nPhonons] = {38.87e-3};
  // Phonon cccupation numbers
  double nocc[nPhonons];
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c[nPhonons];

  for (int j = 0; j < nPhonons; ++j) {
    nocc[j] = 1. / (exp(eph[j] / kbt) - 1);
    c[j] = c0 * dtk[j] * dtk[j] / eph[j];
    c[j] *= 7;
  }

  double en = 0.;
  double dos = 0.;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    for (int j = 0; j < nPhonons; ++j) {
      // Absorption
      dos = GetConductionBandDensityOfStates(en + eph[j], m_nValleysX);
      m_cfElectronsL[i].push_back(c[j] * nocc[j] * dos);
      // Emission
      if (en > m_eMinL + eph[j]) {
        dos = GetConductionBandDensityOfStates(en - eph[j], m_nValleysX);
        m_cfElectronsL[i].push_back(c[j] * (nocc[j] + 1) * dos);
      } else {
        m_cfElectronsL[i].push_back(0.);
      }
    }
    en += m_eStepXL;
  }

  for (int j = 0; j < nPhonons; ++j) {
    // Absorption
    m_energyLossElectronsL.push_back(-eph[j]);
    // Emission
    m_energyLossElectronsL.push_back(eph[j]);
    m_scatTypeElectronsL.push_back(ElectronCollisionTypeIntervalleyF);
    m_scatTypeElectronsL.push_back(ElectronCollisionTypeIntervalleyF);
  }

  m_nLevelsL += 2 * nPhonons;

  return true;
}

bool MediumSilicon::ElectronIntervalleyScatteringRatesXGLG() {
  // Reference:
  //  - K. Hess (editor),
  //    Monte Carlo device simulation: full band and beyond
  //    Chapter 5

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  const int nPhonons = 1;

  // Coupling constants [eV/cm]
  // Average of XG and LG
  const double dtk[nPhonons] = {2.43e8};
  // Phonon energies [eV]
  const double eph[nPhonons] = {37.65e-3};
  // Number of equivalent valleys
  const int zX = 6;
  const int zL = 8;
  const int zG = 1;

  // Phonon cccupation numbers
  double nocc[nPhonons] = {0.};
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c[nPhonons];

  for (int j = 0; j < nPhonons; ++j) {
    nocc[j] = 1. / (exp(eph[j] / kbt) - 1);
    c[j] = c0 * dtk[j] * dtk[j] / eph[j];
  }

  double en = 0.;
  double dos = 0.;
  // XG, LG
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    for (int j = 0; j < nPhonons; ++j) {
      // Absorption
      if (en + eph[j] > m_eMinG) {
        dos = GetConductionBandDensityOfStates(en + eph[j],
                                               m_nValleysX + m_nValleysL);
        m_cfElectronsX[i].push_back(zG * c[j] * nocc[j] * dos);
        m_cfElectronsL[i].push_back(zG * c[j] * nocc[j] * dos);
      } else {
        m_cfElectronsX[i].push_back(0.);
        m_cfElectronsL[i].push_back(0.);
      }
      // Emission
      if (en - eph[j] > m_eMinG) {
        dos = GetConductionBandDensityOfStates(en - eph[j],
                                               m_nValleysX + m_nValleysL);
        m_cfElectronsX[i].push_back(zG * c[j] * (nocc[j] + 1) * dos);
        m_cfElectronsL[i].push_back(zG * c[j] * (nocc[j] + 1) * dos);
      } else {
        m_cfElectronsX[i].push_back(0.);
        m_cfElectronsL[i].push_back(0.);
      }
    }
    en += m_eStepXL;
  }

  // GX, GL
  en = 0.;
  double dosX = 0., dosL = 0.;
  for (int i = 0; i < nEnergyStepsG; ++i) {
    for (int j = 0; j < nPhonons; ++j) {
      // Absorption
      dosX = GetConductionBandDensityOfStates(en + eph[j], 0);
      dosL = GetConductionBandDensityOfStates(en + eph[j], m_nValleysX);
      m_cfElectronsG[i].push_back(zX * c[j] * nocc[j] * dosX);
      if (en > m_eMinL) {
        m_cfElectronsG[i].push_back(zL * c[j] * nocc[j] * dosL);
      } else {
        m_cfElectronsG[i].push_back(0.);
      }
      // Emission
      dosX = GetConductionBandDensityOfStates(en - eph[j], 0);
      dosL = GetConductionBandDensityOfStates(en - eph[j], m_nValleysX);
      if (en > eph[j]) {
        m_cfElectronsG[i].push_back(zX * c[j] * (nocc[j] + 1) * dosX);
      } else {
        m_cfElectronsG[i].push_back(0.);
      }
      if (en - eph[j] > m_eMinL) {
        m_cfElectronsG[i].push_back(zL * c[j] * (nocc[j] + 1) * dosL);
      } else {
        m_cfElectronsG[i].push_back(0.);
      }
    }
    en += m_eStepG;
  }

  for (int j = 0; j < nPhonons; ++j) {
    // Absorption (XL)
    m_energyLossElectronsX.push_back(-eph[j]);
    m_energyLossElectronsL.push_back(-eph[j]);
    // Emission (XL)
    m_energyLossElectronsX.push_back(eph[j]);
    m_energyLossElectronsL.push_back(eph[j]);
    // Absorption (G)
    m_energyLossElectronsG.push_back(-eph[j]);
    m_energyLossElectronsG.push_back(-eph[j]);
    // Emission (G)
    m_energyLossElectronsG.push_back(eph[j]);
    m_energyLossElectronsG.push_back(eph[j]);

    m_scatTypeElectronsX.push_back(ElectronCollisionTypeInterbandXG);
    m_scatTypeElectronsX.push_back(ElectronCollisionTypeInterbandXG);
    m_scatTypeElectronsL.push_back(ElectronCollisionTypeInterbandLG);
    m_scatTypeElectronsL.push_back(ElectronCollisionTypeInterbandLG);

    m_scatTypeElectronsG.push_back(ElectronCollisionTypeInterbandXG);
    m_scatTypeElectronsG.push_back(ElectronCollisionTypeInterbandLG);
    m_scatTypeElectronsG.push_back(ElectronCollisionTypeInterbandXG);
    m_scatTypeElectronsG.push_back(ElectronCollisionTypeInterbandLG);
  }

  m_nLevelsX += 2 * nPhonons;
  m_nLevelsL += 2 * nPhonons;
  m_nLevelsG += 4 * nPhonons;

  return true;
}

bool MediumSilicon::ElectronIonisationRatesXL() {
  // References:
  // - E. Cartier, M. V. Fischetti, E. A. Eklund and F. R. McFeely,
  //   Appl. Phys. Lett 62, 3339-3341
  // - DAMOCLES web page: www.research.ibm.com/DAMOCLES

  // Coefficients [ns-1]
  constexpr double p[3] = {6.25e1, 3.e3, 6.8e5};
  // Threshold energies [eV]
  constexpr double eth[3] = {1.2, 1.8, 3.45};

  double en = 0.;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    double fIon = 0.;
    if (en > eth[0]) {
      fIon += p[0] * (en - eth[0]) * (en - eth[0]);
    }
    if (en > eth[1]) {
      fIon += p[1] * (en - eth[1]) * (en - eth[1]);
    }
    if (en > eth[2]) {
      fIon += p[2] * (en - eth[2]) * (en - eth[2]);
    }
    m_cfElectronsX[i].push_back(fIon);
    m_cfElectronsL[i].push_back(fIon);
    en += m_eStepXL;
  }

  m_energyLossElectronsX.push_back(eth[0]);
  m_energyLossElectronsL.push_back(eth[0]);
  m_scatTypeElectronsX.push_back(ElectronCollisionTypeIonisation);
  m_scatTypeElectronsL.push_back(ElectronCollisionTypeIonisation);
  ++m_nLevelsX;
  ++m_nLevelsL;

  return true;
}

bool MediumSilicon::ElectronIonisationRatesG() {
  // References:
  // - E. Cartier, M. V. Fischetti, E. A. Eklund and F. R. McFeely,
  //   Appl. Phys. Lett 62, 3339-3341
  // - S. Tanuma, C. J. Powell and D. R. Penn
  //   Surf. Interface Anal. (2010)

  // Coefficients [ns-1]
  constexpr double p[3] = {6.25e1, 3.e3, 6.8e5};
  // Threshold energies [eV]
  constexpr double eth[3] = {1.2, 1.8, 3.45};

  double en = 0.;
  for (int i = 0; i < nEnergyStepsG; ++i) {
    double fIon = 0.;
    if (en > eth[0]) {
      fIon += p[0] * (en - eth[0]) * (en - eth[0]);
    }
    if (en > eth[1]) {
      fIon += p[1] * (en - eth[1]) * (en - eth[1]);
    }
    if (en > eth[2]) {
      fIon += p[2] * (en - eth[2]) * (en - eth[2]);
    }
    if (en >= m_eMinG) {
      m_cfElectronsG[i].push_back(fIon);
    } else {
      m_cfElectronsG[i].push_back(0.);
    }
    en += m_eStepG;
  }

  m_energyLossElectronsG.push_back(eth[0]);
  m_scatTypeElectronsG.push_back(ElectronCollisionTypeIonisation);
  ++m_nLevelsG;

  return true;
}

bool MediumSilicon::ElectronImpurityScatteringRates() {
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Band parameters
  // Density of states effective masses
  const double mdX =
      ElectronMass * pow(m_mLongX * m_mTransX * m_mTransX, 1. / 3.);
  const double mdL =
      ElectronMass * pow(m_mLongL * m_mTransL * m_mTransL, 1. / 3.);

  // Dielectric constant
  const double eps = GetDielectricConstant();
  // Impurity concentration
  const double impurityConcentration = m_dopingConcentration;
  if (impurityConcentration < Small) return true;

  // Screening length
  const double ls = sqrt(eps * kbt / (4 * Pi * FineStructureConstant * HbarC *
                                      impurityConcentration));
  const double ebX = 0.5 * HbarC * HbarC / (mdX * ls * ls);
  const double ebL = 0.5 * HbarC * HbarC / (mdL * ls * ls);

  // Prefactor
  // const double c = pow(2., 2.5) * Pi * impurityConcentration *
  //                 pow(FineStructureConstant * HbarC, 2) *
  //                 SpeedOfLight / (eps * eps * sqrt(md) * eb * eb);
  // Use momentum-transfer cross-section
  const double cX = impurityConcentration * Pi *
                    pow(FineStructureConstant * HbarC, 2) * SpeedOfLight /
                    (sqrt(2 * mdX) * eps * eps);
  const double cL = impurityConcentration * Pi *
                    pow(FineStructureConstant * HbarC, 2) * SpeedOfLight /
                    (sqrt(2 * mdL) * eps * eps);

  double en = 0.;
  for (int i = 0; i < nEnergyStepsXL; ++i) {
    const double gammaX = en * (1. + m_alphaX * en);
    const double gammaL = (en - m_eMinL) * (1. + m_alphaL * (en - m_eMinL));
    // m_cfElectrons[i][iLevel] = c * sqrt(gamma) * (1. + 2 * alpha * en) /
    //                         (1. + 4. * gamma / eb);
    if (gammaX <= 0.) {
      m_cfElectronsX[i].push_back(0.);
    } else {
      const double b = 4 * gammaX / ebX;
      m_cfElectronsX[i].push_back((cX / pow(gammaX, 1.5)) *
                                  (log(1. + b) - b / (1. + b)));
    }
    if (en <= m_eMinL || gammaL <= 0.) {
      m_cfElectronsL[i].push_back(0.);
    } else {
      const double b = 4 * gammaL / ebL;
      m_cfElectronsL[i].push_back((cL / pow(gammaL, 1.5)) *
                                  (log(1. + b) - b / (1. + b)));
    }
    en += m_eStepXL;
  }

  m_energyLossElectronsX.push_back(0.);
  m_energyLossElectronsL.push_back(0.);
  m_scatTypeElectronsX.push_back(ElectronCollisionTypeImpurity);
  m_scatTypeElectronsL.push_back(ElectronCollisionTypeImpurity);
  ++m_nLevelsX;
  ++m_nLevelsL;

  return true;
}

bool MediumSilicon::HoleScatteringRates() {
  // Reset the scattering rates
  m_cfTotHoles.assign(nEnergyStepsV, 0.);
  m_cfHoles.assign(nEnergyStepsV, std::vector<double>());
  m_energyLossHoles.clear();
  m_scatTypeHoles.clear();
  m_cfNullHoles = 0.;

  m_nLevelsV = 0;
  // Fill the scattering rates table
  HoleAcousticScatteringRates();
  HoleOpticalScatteringRates();
  // HoleImpurityScatteringRates();
  HoleIonisationRates();

  std::ofstream outfile;
  if (m_cfOutput) {
    outfile.open("ratesV.txt", std::ios::out);
  }

  for (int i = 0; i < nEnergyStepsV; ++i) {
    // Sum up the scattering rates of all processes.
    for (int j = m_nLevelsV; j--;) m_cfTotHoles[i] += m_cfHoles[i][j];

    if (m_cfOutput) {
      outfile << i * m_eStepV << " " << m_cfTotHoles[i] << " ";
      for (int j = 0; j < m_nLevelsV; ++j) {
        outfile << m_cfHoles[i][j] << " ";
      }
      outfile << "\n";
    }

    if (m_cfTotHoles[i] > m_cfNullHoles) {
      m_cfNullHoles = m_cfTotHoles[i];
    }

    // Make sure the total scattering rate is positive.
    if (m_cfTotHoles[i] <= 0.) {
      std::cerr << m_className << "::HoleScatteringRates:\n"
                << "    Scattering rate at " << i * m_eStepV << " eV <= 0.\n";
      return false;
    }
    // Normalise the rates.
    for (int j = 0; j < m_nLevelsV; ++j) {
      m_cfHoles[i][j] /= m_cfTotHoles[i];
      if (j > 0) m_cfHoles[i][j] += m_cfHoles[i][j - 1];
    }
  }

  if (m_cfOutput) {
    outfile.close();
  }

  return true;
}

bool MediumSilicon::HoleAcousticScatteringRates() {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - DAMOCLES web page: www.research.ibm.com/DAMOCLES
  //  - M. Lundstrom, Fundamentals of carrier transport

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Acoustic phonon intraband scattering
  // Acoustic deformation potential [eV]
  // DAMOCLES: 4.6 eV; Lundstrom: 5 eV
  constexpr double defpot2 = 4.6 * 4.6;
  // Longitudinal velocity of sound [cm/ns]
  constexpr double u = 9.04e-4;
  // Prefactor for acoustic deformation potential scattering
  const double cIntra = TwoPi * SpeedOfLight * SpeedOfLight * kbt * defpot2 /
                        (Hbar * u * u * rho);

  // Fill the scattering rate tables.
  double en = Small;
  for (int i = 0; i < nEnergyStepsV; ++i) {
    const double dos = GetValenceBandDensityOfStates(en, 0);
    m_cfHoles[i].push_back(cIntra * dos);
    en += m_eStepV;
  }

  // Assume that energy loss is negligible.
  m_energyLossHoles.push_back(0.);
  m_scatTypeHoles.push_back(ElectronCollisionTypeAcousticPhonon);
  ++m_nLevelsV;

  return true;
}

bool MediumSilicon::HoleOpticalScatteringRates() {
  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705
  //  - DAMOCLES web page: www.research.ibm.com/DAMOCLES
  //  - M. Lundstrom, Fundamentals of carrier transport

  // Mass density [(eV/c2)/cm3]
  const double rho = m_density * m_a * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * m_temperature;

  // Coupling constant [eV/cm]
  // DAMOCLES: 6.6, Lundstrom: 6.0
  constexpr double dtk = 6.6e8;
  // Phonon energy [eV]
  constexpr double eph = 63.0e-3;
  // Phonon cccupation numbers
  const double nocc = 1. / (exp(eph / kbt) - 1);
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c = c0 * dtk * dtk / eph;

  double en = 0.;
  for (int i = 0; i < nEnergyStepsV; ++i) {
    // Absorption
    double dos = GetValenceBandDensityOfStates(en + eph, 0);
    m_cfHoles[i].push_back(c * nocc * dos);
    // Emission
    if (en > eph) {
      dos = GetValenceBandDensityOfStates(en - eph, 0);
      m_cfHoles[i].push_back(c * (nocc + 1) * dos);
    } else {
      m_cfHoles[i].push_back(0.);
    }
    en += m_eStepV;
  }

  // Absorption
  m_energyLossHoles.push_back(-eph);
  // Emission
  m_energyLossHoles.push_back(eph);
  m_scatTypeHoles.push_back(ElectronCollisionTypeOpticalPhonon);
  m_scatTypeHoles.push_back(ElectronCollisionTypeOpticalPhonon);

  m_nLevelsV += 2;

  return true;
}

bool MediumSilicon::HoleIonisationRates() {
  // References:
  //  - DAMOCLES web page: www.research.ibm.com/DAMOCLES

  // Coefficients [ns-1]
  constexpr double p[2] = {2., 1.e3};
  // Threshold energies [eV]
  constexpr double eth[2] = {1.1, 1.45};
  // Exponents
  constexpr double b[2] = {6., 4.};

  double en = 0.;
  for (int i = 0; i < nEnergyStepsV; ++i) {
    double fIon = 0.;
    if (en > eth[0]) {
      fIon += p[0] * pow(en - eth[0], b[0]);
    }
    if (en > eth[1]) {
      fIon += p[1] * pow(en - eth[1], b[1]);
    }
    m_cfHoles[i].push_back(fIon);
    en += m_eStepV;
  }

  m_energyLossHoles.push_back(eth[0]);
  m_scatTypeHoles.push_back(ElectronCollisionTypeIonisation);
  ++m_nLevelsV;

  return true;
}

double MediumSilicon::GetConductionBandDensityOfStates(const double e,
                                                       const int band) {
  if (band < 0) {
    int iE = int(e / m_eStepDos);
    const int nPoints = m_fbDosConduction.size();
    if (iE >= nPoints || iE < 0) {
      return 0.;
    } else if (iE == nPoints - 1) {
      return m_fbDosConduction[nPoints - 1];
    }

    const double dos = m_fbDosConduction[iE] +
                       (m_fbDosConduction[iE + 1] - m_fbDosConduction[iE]) *
                           (e / m_eStepDos - iE);
    return dos * 1.e21;

  } else if (band < m_nValleysX) {
    // X valleys
    if (e <= 0.) return 0.;
    // Density-of-states effective mass (cube)
    const double md3 = pow(ElectronMass, 3) * m_mLongX * m_mTransX * m_mTransX;

    if (m_fullBandDos) {
      if (e < m_eMinL) {
        return GetConductionBandDensityOfStates(e, -1) / m_nValleysX;
      } else if (e < m_eMinG) {
        // Subtract the fraction of the full-band density of states
        // attributed to the L valleys.
        const double dosX =
            GetConductionBandDensityOfStates(e, -1) -
            GetConductionBandDensityOfStates(e, m_nValleysX) * m_nValleysL;
        return dosX / m_nValleysX;
      } else {
        // Subtract the fraction of the full-band density of states
        // attributed to the L valleys and the higher bands.
        const double dosX =
            GetConductionBandDensityOfStates(e, -1) -
            GetConductionBandDensityOfStates(e, m_nValleysX) * m_nValleysL -
            GetConductionBandDensityOfStates(e, m_nValleysX + m_nValleysL);
        if (dosX <= 0.) return 0.;
        return dosX / m_nValleysX;
      }
    } else if (m_nonParabolic) {
      const double alpha = m_alphaX;
      return sqrt(md3 * e * (1. + alpha * e) / 2.) * (1. + 2 * alpha * e) /
             (Pi2 * pow(HbarC, 3.));
    } else {
      return sqrt(md3 * e / 2.) / (Pi2 * pow(HbarC, 3.));
    }
  } else if (band < m_nValleysX + m_nValleysL) {
    // L valleys
    if (e <= m_eMinL) return 0.;

    // Density-of-states effective mass (cube)
    const double md3 = pow(ElectronMass, 3) * m_mLongL * m_mTransL * m_mTransL;
    // Non-parabolicity parameter
    const double alpha = m_alphaL;

    if (m_fullBandDos) {
      // Energy up to which the non-parabolic approximation is used.
      const double ej = m_eMinL + 0.5;
      if (e <= ej) {
        return sqrt(md3 * (e - m_eMinL) * (1. + alpha * (e - m_eMinL))) *
               (1. + 2 * alpha * (e - m_eMinL)) /
               (Sqrt2 * Pi2 * pow(HbarC, 3.));
      } else {
        // Fraction of full-band density of states attributed to L valleys
        double fL = sqrt(md3 * (ej - m_eMinL) * (1. + alpha * (ej - m_eMinL))) *
                    (1. + 2 * alpha * (ej - m_eMinL)) /
                    (Sqrt2 * Pi2 * pow(HbarC, 3.));
        fL = m_nValleysL * fL / GetConductionBandDensityOfStates(ej, -1);

        double dosXL = GetConductionBandDensityOfStates(e, -1);
        if (e > m_eMinG) {
          dosXL -=
              GetConductionBandDensityOfStates(e, m_nValleysX + m_nValleysL);
        }
        if (dosXL <= 0.) return 0.;
        return fL * dosXL / 8.;
      }
    } else if (m_nonParabolic) {
      return sqrt(md3 * (e - m_eMinL) * (1. + alpha * (e - m_eMinL))) *
             (1. + 2 * alpha * (e - m_eMinL)) / (Sqrt2 * Pi2 * pow(HbarC, 3.));
    } else {
      return sqrt(md3 * (e - m_eMinL) / 2.) / (Pi2 * pow(HbarC, 3.));
    }
  } else if (band == m_nValleysX + m_nValleysL) {
    // Higher bands
    const double ej = 2.7;
    if (m_eMinG >= ej) {
      std::cerr << m_className << "::GetConductionBandDensityOfStates:\n"
                << "    Cannot determine higher band density-of-states.\n"
                << "    Program bug. Check offset energy!\n";
      return 0.;
    }
    if (e < m_eMinG) {
      return 0.;
    } else if (e < ej) {
      // Coexistence of XL and higher bands.
      const double dj = GetConductionBandDensityOfStates(ej, -1);
      // Assume linear increase of density-of-states.
      return dj * (e - m_eMinG) / (ej - m_eMinG);
    } else {
      return GetConductionBandDensityOfStates(e, -1);
    }
  }

  std::cerr << m_className << "::GetConductionBandDensityOfStates:\n"
            << "    Band index (" << band << ") out of range.\n";
  return ElectronMass * sqrt(ElectronMass * e / 2.) / (Pi2 * pow(HbarC, 3.));
}

double MediumSilicon::GetValenceBandDensityOfStates(const double e,
                                                    const int band) {
  if (band <= 0) {
    // Total (full-band) density of states.
    const int nPoints = m_fbDosValence.size();
    int iE = int(e / m_eStepDos);
    if (iE >= nPoints || iE < 0) {
      return 0.;
    } else if (iE == nPoints - 1) {
      return m_fbDosValence[nPoints - 1];
    }

    const double dos =
        m_fbDosValence[iE] +
        (m_fbDosValence[iE + 1] - m_fbDosValence[iE]) * (e / m_eStepDos - iE);
    return dos * 1.e21;
  }

  std::cerr << m_className << "::GetConductionBandDensityOfStates:\n"
            << "    Band index (" << band << ") out of range.\n";
  return 0.;
}

void MediumSilicon::ComputeSecondaries(const double e0, double& ee,
                                       double& eh) {
  const int nPointsValence = m_fbDosValence.size();
  const int nPointsConduction = m_fbDosConduction.size();
  const double widthValenceBand = m_eStepDos * nPointsValence;
  const double widthConductionBand = m_eStepDos * nPointsConduction;

  bool ok = false;
  while (!ok) {
    // Sample a hole energy according to the valence band DOS.
    eh = RndmUniformPos() * std::min(widthValenceBand, e0);
    int ih = std::min(int(eh / m_eStepDos), nPointsValence - 1);
    while (RndmUniform() > m_fbDosValence[ih] / m_fbDosMaxV) {
      eh = RndmUniformPos() * std::min(widthValenceBand, e0);
      ih = std::min(int(eh / m_eStepDos), nPointsValence - 1);
    }
    // Sample an electron energy according to the conduction band DOS.
    ee = RndmUniformPos() * std::min(widthConductionBand, e0);
    int ie = std::min(int(ee / m_eStepDos), nPointsConduction - 1);
    while (RndmUniform() > m_fbDosConduction[ie] / m_fbDosMaxC) {
      ee = RndmUniformPos() * std::min(widthConductionBand, e0);
      ie = std::min(int(ee / m_eStepDos), nPointsConduction - 1);
    }
    // Calculate the energy of the primary electron.
    const double ep = e0 - m_bandGap - eh - ee;
    if (ep < Small) continue;
    if (ep > 5.) return;
    // Check if the primary electron energy is consistent with the DOS.
    int ip = std::min(int(ep / m_eStepDos), nPointsConduction - 1);
    if (RndmUniform() > m_fbDosConduction[ip] / m_fbDosMaxC) continue;
    ok = true;
  }
}

void MediumSilicon::InitialiseDensityOfStates() {
  m_eStepDos = 0.1;

  m_fbDosValence = {
      {0.,      1.28083,  2.08928, 2.70763, 3.28095, 3.89162, 4.50547, 5.15043,
       5.89314, 6.72667,  7.67768, 8.82725, 10.6468, 12.7003, 13.7457, 14.0263,
       14.2731, 14.5527,  14.8808, 15.1487, 15.4486, 15.7675, 16.0519, 16.4259,
       16.7538, 17.0589,  17.3639, 17.6664, 18.0376, 18.4174, 18.2334, 16.7552,
       15.1757, 14.2853,  13.6516, 13.2525, 12.9036, 12.7203, 12.6104, 12.6881,
       13.2862, 14.0222,  14.9366, 13.5084, 9.77808, 6.15266, 3.47839, 2.60183,
       2.76747, 3.13985,  3.22524, 3.29119, 3.40868, 3.6118,  3.8464,  4.05776,
       4.3046,  4.56219,  4.81553, 5.09909, 5.37616, 5.67297, 6.04611, 6.47252,
       6.9256,  7.51254,  8.17923, 8.92351, 10.0309, 11.726,  16.2853, 18.2457,
       12.8879, 7.86019,  6.02275, 5.21777, 4.79054, 3.976,   3.11855, 2.46854,
       1.65381, 0.830278, 0.217735}};

  m_fbDosConduction = {
      {0.,      1.5114,  2.71026,  3.67114,  4.40173, 5.05025, 5.6849,  6.28358,
       6.84628, 7.43859, 8.00204,  8.80658,  9.84885, 10.9579, 12.0302, 13.2051,
       14.6948, 16.9879, 18.4492,  18.1933,  17.6747, 16.8135, 15.736,  14.4965,
       13.1193, 12.1817, 12.6109,  15.3148,  19.4936, 23.0093, 24.4106, 22.2834,
       19.521,  18.9894, 18.8015,  17.9363,  17.0252, 15.9871, 14.8486, 14.3797,
       14.2426, 14.3571, 14.7271,  14.681,   14.3827, 14.2789, 14.144,  14.1684,
       14.1418, 13.9237, 13.7558,  13.5691,  13.4567, 13.2693, 12.844,  12.4006,
       12.045,  11.7729, 11.3607,  11.14,    11.0586, 10.5475, 9.73786, 9.34423,
       9.4694,  9.58071, 9.6967,   9.84854,  10.0204, 9.82705, 9.09102, 8.30665,
       7.67306, 7.18925, 6.79675,  6.40713,  6.21687, 6.33267, 6.5223,  6.17877,
       5.48659, 4.92208, 4.44239,  4.02941,  3.5692,  3.05953, 2.6428,  2.36979,
       2.16273, 2.00627, 1.85206,  1.71265,  1.59497, 1.46681, 1.34913, 1.23951,
       1.13439, 1.03789, 0.924155, 0.834962, 0.751017}};

  auto it = std::max_element(m_fbDosValence.begin(), m_fbDosValence.end());
  m_fbDosMaxV = *it;
  it = std::max_element(m_fbDosConduction.begin(), m_fbDosConduction.end());
  m_fbDosMaxC = *it;
}
}
