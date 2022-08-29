#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MagboltzInterface.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/OpticalData.hh"
#include "Garfield/Random.hh"
#include "Garfield/Utilities.hh"
#include "Garfield/ViewBase.hh"

namespace {

bool IsComment(const std::string& line) {
  if (line.empty()) return false;
  if (line[0] == '#') return true;
  return false;
}

std::string GetDescription(const unsigned int index,
                           char scrpt[][Garfield::Magboltz::nCharDescr]) {
  return std::string(scrpt[index],
                     scrpt[index] + Garfield::Magboltz::nCharDescr);
}

std::string GetDescription(const unsigned int i1, const unsigned int i2,
                           char scrpt[][6][Garfield::Magboltz::nCharDescr]) {
  return std::string(scrpt[i1][i2],
                     scrpt[i1][i2] + Garfield::Magboltz::nCharDescr);
}

void SetScatteringParameters(const int model, const double parIn, double& cut,
                             double& parOut) {
  cut = 1.;
  parOut = 0.5;
  if (model <= 0) return;

  if (model >= 2) {
    parOut = parIn;
    return;
  }

  // Set cuts on angular distribution and
  // renormalise forward scattering probability.
  if (parIn <= 1.) {
    parOut = parIn;
    return;
  }

  const double cns = parIn - 0.5;
  const double thetac = asin(2. * sqrt(cns - cns * cns));
  const double fac = (1. - cos(thetac)) / pow(sin(thetac), 2.);
  parOut = cns * fac + 0.5;
  cut = thetac * 2. / Garfield::Pi;
}

}

namespace Garfield {

const int MediumMagboltz::DxcTypeRad = 0;
const int MediumMagboltz::DxcTypeCollIon = 1;
const int MediumMagboltz::DxcTypeCollNonIon = -1;

MediumMagboltz::MediumMagboltz()
    : MediumGas(),
      m_eMax(40.),
      m_eStep(m_eMax / Magboltz::nEnergySteps),
      m_eHigh(400.),
      m_eHighLog(log(m_eHigh)),
      m_lnStep(1.),
      m_eFinalGamma(20.),
      m_eStepGamma(m_eFinalGamma / nEnergyStepsGamma) {
  m_className = "MediumMagboltz";

  // Set physical constants in Magboltz common blocks.
  Magboltz::cnsts_.echarg = ElementaryCharge * 1.e-15;
  Magboltz::cnsts_.emass = ElectronMassGramme;
  Magboltz::cnsts_.amu = AtomicMassUnit;
  Magboltz::cnsts_.pir2 = BohrRadius * BohrRadius * Pi;
  Magboltz::inpt_.ary = RydbergEnergy;

  // Set parameters in Magboltz common blocks.
  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = Magboltz::nEnergySteps;
  // Select the scattering model.
  Magboltz::inpt_.nAniso = 2;
  // Max. energy [eV]
  Magboltz::inpt_.efinal = m_eMax;
  // Energy step size [eV]
  Magboltz::inpt_.estep = m_eStep;
  // Temperature and pressure
  Magboltz::inpt_.akt = BoltzmannConstant * m_temperature;
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;
  // Disable Penning transfer.
  Magboltz::inpt_.ipen = 0;

  m_description.assign(Magboltz::nMaxLevels,
                       std::string(Magboltz::nCharDescr, ' '));

  m_cfTot.assign(Magboltz::nEnergySteps, 0.);
  m_cfTotLog.assign(nEnergyStepsLog, 0.);
  m_cf.assign(Magboltz::nEnergySteps,
              std::vector<double>(Magboltz::nMaxLevels, 0.));
  m_cfLog.assign(nEnergyStepsLog,
                 std::vector<double>(Magboltz::nMaxLevels, 0.));

  m_isChanged = true;

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = true;

  m_scaleExc.fill(1.);
}
MediumMagboltz::MediumMagboltz(const std::string& gas1, const double f1,
                               const std::string& gas2, const double f2,
                               const std::string& gas3, const double f3,
                               const std::string& gas4, const double f4,
                               const std::string& gas5, const double f5,
                               const std::string& gas6, const double f6) 
    : MediumMagboltz() {
  SetComposition(gas1, f1, gas2, f2, gas3, f3, 
                 gas4, f4, gas5, f5, gas6, f6);
}

bool MediumMagboltz::SetMaxElectronEnergy(const double e) {
  if (e <= Small) {
    std::cerr << m_className << "::SetMaxElectronEnergy: Invalid energy.\n";
    return false;
  }
  m_eMax = e;

  std::lock_guard<std::mutex> guard(m_mutex);
  // Determine the energy interval size.
  m_eStep = std::min(m_eMax, m_eHigh) / Magboltz::nEnergySteps;

  // Force recalculation of the scattering rates table.
  m_isChanged = true;

  return true;
}

bool MediumMagboltz::SetMaxPhotonEnergy(const double e) {
  if (e <= Small) {
    std::cerr << m_className << "::SetMaxPhotonEnergy: Invalid energy.\n";
    return false;
  }
  m_eFinalGamma = e;

  // Determine the energy interval size.
  m_eStepGamma = m_eFinalGamma / nEnergyStepsGamma;

  // Force recalculation of the scattering rates table.
  m_isChanged = true;

  return true;
}

void MediumMagboltz::SetSplittingFunctionOpalBeaty() {
  m_useOpalBeaty = true;
  m_useGreenSawada = false;
}

void MediumMagboltz::SetSplittingFunctionGreenSawada() {
  m_useOpalBeaty = false;
  m_useGreenSawada = true;
  if (m_isChanged) return;

  bool allset = true;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (!m_hasGreenSawada[i]) {
      if (allset) {
        std::cout << m_className << "::SetSplittingFunctionGreenSawada:\n";
        allset = false;
      }
      std::cout << "    Fit parameters for " << m_gas[i] << " not available.\n"
                << "    Using Opal-Beaty formula instead.\n";
    }
  }
}

void MediumMagboltz::SetSplittingFunctionFlat() {
  m_useOpalBeaty = false;
  m_useGreenSawada = false;
}

void MediumMagboltz::EnableDeexcitation() {
  if (m_usePenning) {
    std::cout << m_className << "::EnableDeexcitation:\n"
              << "    Penning transfer will be switched off.\n";
  }
  // if (m_useRadTrap) {
  //   std::cout << "    Radiation trapping is switched on.\n";
  // } else {
  //   std::cout << "    Radiation trapping is switched off.\n";
  // }
  m_usePenning = false;
  m_useDeexcitation = true;
  m_isChanged = true;
  m_dxcProducts.clear();
}

void MediumMagboltz::EnableRadiationTrapping() {
  m_useRadTrap = true;
  if (!m_useDeexcitation) {
    std::cout << m_className << "::EnableRadiationTrapping:\n    "
              << "Radiation trapping is enabled but de-excitation is not.\n";
  } else {
    m_isChanged = true;
  }
}

bool MediumMagboltz::EnablePenningTransfer() {
  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);
  if (!MediumGas::EnablePenningTransfer()) return false;

  m_usePenning = true;
  return true;
}   

bool MediumMagboltz::EnablePenningTransfer(const double r,
                                           const double lambda) {
   
  if (!MediumGas::EnablePenningTransfer(r, lambda)) return false;
 
  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);

  // Make sure that the collision rate table is up to date.
  if (!Update()) return false;
  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation) {
      ++nLevelsFound;
    }
    m_rPenning[i] = m_rPenningGlobal;
    m_lambdaPenning[i] = m_lambdaPenningGlobal;
  }

  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n    "
              << "Updated Penning transfer parameters for " << nLevelsFound
              << " excitation cross-sections.\n";
    if (nLevelsFound != m_excLevels.size() && !m_excLevels.empty()) {
      std::cerr << m_className << "::EnablePenningTransfer:\n    Warning: "
                << "mismatch between number of excitation cross-sections ("
                << nLevelsFound << ")\n    and number of excitation rates in "
                << "the gas table (" << m_excLevels.size() << ").\n    "
                << "The gas table was probably calculated using a different "
                << "version of Magboltz.\n";
    }
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n    "
              << "No excitation cross-sections in the present energy range.\n";
  }

  if (m_useDeexcitation) {
    std::cout << m_className << "::EnablePenningTransfer:\n    "
              << "Deexcitation handling will be switched off.\n";
  }
  m_usePenning = true;
  return true;
}

bool MediumMagboltz::EnablePenningTransfer(const double r, const double lambda,
                                           std::string gasname) {

  if (!MediumGas::EnablePenningTransfer(r, lambda, gasname)) return false;

  // Get (again) the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) return false;

  // Look (again) for this gas in the present mixture.
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      iGas = i;
      break;
    }
  }

  // Make sure that the collision rate table is up to date.
  if (!Update()) return false;
  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (int(m_csType[i] / nCsTypes) != iGas) continue;
    if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation) {
      ++nLevelsFound;
    }
    m_rPenning[i] = m_rPenningGas[iGas];
    m_lambdaPenning[i] = m_lambdaPenningGas[iGas];
  }

  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n";
    if (m_lambdaPenningGas[iGas] > 0.) {
      std::cout << "    Penning transfer parameters for " << nLevelsFound
                << " " << gasname << " excitation levels set to:\n"
                << "      r = " << m_rPenningGas[iGas] << ", lambda = "
                << m_lambdaPenningGas[iGas] << " cm\n";
    } else {
      std::cout << "    Penning transfer probability for " << nLevelsFound
                << " " << gasname << " excitation levels set to r = "
                << m_rPenningGas[iGas] << "\n";
    }
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n    " << gasname
              << " has no excitation levels in the present energy range.\n";
  }

  m_usePenning = true;
  return true;
}

void MediumMagboltz::DisablePenningTransfer() {

  MediumGas::DisablePenningTransfer();
  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);

  m_usePenning = false;
}

bool MediumMagboltz::DisablePenningTransfer(std::string gasname) {

  if (!MediumGas::DisablePenningTransfer(gasname)) return false;
  // Get the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) return false;

  // Look (again) for this gas in the present gas mixture.
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      iGas = i;
      break;
    }
  }

  if (iGas < 0) return false;

  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (int(m_csType[i] / nCsTypes) == iGas) {
      m_rPenning[i] = 0.;
      m_lambdaPenning[i] = 0.;
    } else {
      if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation &&
          m_rPenning[i] > Small) {
        ++nLevelsFound;
      }
    }
  }

  if (nLevelsFound == 0) {
    // There are no more excitation levels with r > 0.
    std::cout << m_className << "::DisablePenningTransfer:\n"
              << "    Penning transfer switched off for all excitations.\n";
    m_usePenning = false;
  }
  return true;
}

void MediumMagboltz::SetExcitationScaling(const double r, std::string gasname) {
  if (r <= 0.) {
    std::cerr << m_className << "::SetExcitationScaling: Incorrect value.\n";
    return;
  }

  // Get the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) {
    std::cerr << m_className << "::SetExcitationScaling: Unknown gas name.\n";
    return;
  }

  // Look for this gas in the present gas mixture.
  bool found = false;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_scaleExc[i] = r;
      found = true;
      break;
    }
  }

  if (!found) {
    std::cerr << m_className << "::SetExcitationScaling:\n"
              << "    Specified gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return;
  }

  // Force re-calculation of the collision rate table.
  m_isChanged = true;
}

bool MediumMagboltz::Initialise(const bool verbose) {
  if (!m_isChanged) {
    if (m_debug) {
      std::cerr << m_className << "::Initialise: Nothing changed.\n";
    }
    return true;
  }
  return Update(verbose);
}

void MediumMagboltz::PrintGas() {
  MediumGas::PrintGas();

  if (m_isChanged) {
    if (!Initialise()) return;
  }

  std::cout << "    Electron cross-sections:\n";
  int igas = -1; 
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    // Collision type
    int type = m_csType[i] % nCsTypes;
    if (igas != int(m_csType[i] / nCsTypes)) {
      igas = int(m_csType[i] / nCsTypes);
      std::cout << "      " << m_gas[igas] << "\n";
    }
    // Description (from Magboltz)
    // Threshold energy
    double e = m_rgas[igas] * m_energyLoss[i];
    std::cout << "        Level " << i << ": " << m_description[i] << "\n";
    std::cout << "          Type " << type;
    if (type == ElectronCollisionTypeElastic) {
      std::cout << " (elastic)\n";
    } else if (type == ElectronCollisionTypeIonisation) {
      std::cout << " (ionisation). Ionisation threshold: " << e << " eV.\n";
    } else if (type == ElectronCollisionTypeAttachment) {
      std::cout << " (attachment)\n";
    } else if (type == ElectronCollisionTypeInelastic) {
      std::cout << " (inelastic). Energy loss: " << e << " eV.\n";
    } else if (type == ElectronCollisionTypeExcitation) {
      std::cout << " (excitation). Excitation energy: " << e << " eV.\n";
    } else if (type == ElectronCollisionTypeSuperelastic) {
      std::cout << " (super-elastic). Energy gain: " << -e << " eV.\n";
    } else if (type == ElectronCollisionTypeVirtual) {
      std::cout << " (virtual)\n";
    } else {
      std::cout << " (unknown)\n";
    }
    if (type == ElectronCollisionTypeExcitation && m_usePenning &&
        e > m_minIonPot) {
      std::cout << "          Penning transfer coefficient: " 
                << m_rPenning[i] << "\n";
    } else if (type == ElectronCollisionTypeExcitation && m_useDeexcitation) {
      const int idxc = m_iDeexcitation[i];
      if (idxc < 0 || idxc >= (int)m_deexcitations.size()) {
        std::cout << "          Deexcitation cascade not implemented.\n";
        continue;
      }
      const auto& dxc = m_deexcitations[idxc];
      if (dxc.osc > 0.) {
        std::cout << "          Oscillator strength: " << dxc.osc << "\n";
      }
      std::cout << "          Decay channels:\n";
      const int nChannels = dxc.type.size();
      for (int j = 0; j < nChannels; ++j) {
        if (dxc.type[j] == DxcTypeRad) {
          std::cout << "          Radiative decay to ";
          if (dxc.final[j] < 0) {
            std::cout << "ground state: ";
          } else {
            std::cout << m_deexcitations[dxc.final[j]].label << ": ";
          }
        } else if (dxc.type[j] == DxcTypeCollIon) {
          if (dxc.final[j] < 0) {
            std::cout << "          Penning ionisation: ";
          } else {
            std::cout << "          Associative ionisation: ";
          }
        } else if (dxc.type[j] == DxcTypeCollNonIon) {
          if (dxc.final[j] >= 0) {
            std::cout << "          Collision-induced transition to "
                      << m_deexcitations[dxc.final[j]].label << ": ";
          } else {
            std::cout << "          Loss: ";
          }
        }
        const double br = j == 0 ? dxc.p[j] : dxc.p[j] - dxc.p[j - 1];
        std::cout << std::setprecision(5) << br * 100. << "%\n";
      }
    }
  }
}

double MediumMagboltz::GetElectronNullCollisionRate(const int band) {
  // If necessary, update the collision rates table.
  if (!Update()) return 0.;

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronNullCollisionRate: Band > 0.\n";
  }

  return m_cfNull;
}

double MediumMagboltz::GetElectronCollisionRate(const double e,
                                                const int band) {
  // Check if the electron energy is within the currently set range.
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate: Invalid energy.\n";
    return m_cfTot[0];
  }
  if (e > m_eMax && m_useAutoAdjust) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n    Rate at " << e
              << " eV is not included in the current table.\n    "
              << "Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  // If necessary, update the collision rates table.
  if (!Update()) return 0.;

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronCollisionRate: Band > 0.\n";
  }

  // Get the energy interval.
  if (e <= m_eHigh) {
    // Linear binning
    constexpr int iemax = Magboltz::nEnergySteps - 1;
    const int iE = std::min(std::max(int(e / m_eStep), 0), iemax);
    return m_cfTot[iE];
  }

  // Logarithmic binning
  const double eLog = log(e);
  int iE = int((eLog - m_eHighLog) / m_lnStep);
  // Calculate the collision rate by log-log interpolation.
  const double fmax = m_cfTotLog[iE];
  const double fmin = iE == 0 ? log(m_cfTot.back()) : m_cfTotLog[iE - 1];
  const double emin = m_eHighLog + iE * m_lnStep;
  const double f = fmin + (eLog - emin) * (fmax - fmin) / m_lnStep;
  return exp(f);
}

double MediumMagboltz::GetElectronCollisionRate(const double e,
                                                const unsigned int level,
                                                const int band) {
  // Check if the electron energy is within the currently set range.
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate: Invalid energy.\n";
    return 0.;
  }

  // Check if the level exists.
  if (level >= m_nTerms) {
    std::cerr << m_className << "::GetElectronCollisionRate: Invalid level.\n";
    return 0.;
  }

  // Get the total scattering rate.
  double rate = GetElectronCollisionRate(e, band);
  // Get the energy interval.
  if (e <= m_eHigh) {
    // Linear binning
    constexpr int iemax = Magboltz::nEnergySteps - 1;
    const int iE = std::min(std::max(int(e / m_eStep), 0), iemax);
    if (level == 0) {
      rate *= m_cf[iE][0];
    } else {
      rate *= m_cf[iE][level] - m_cf[iE][level - 1];
    }
  } else {
    // Logarithmic binning
    const int iE = int((log(e) - m_eHighLog) / m_lnStep);
    if (level == 0) {
      rate *= m_cfLog[iE][0];
    } else {
      rate *= m_cfLog[iE][level] - m_cfLog[iE][level - 1];
    }
  }
  return rate;
}

bool MediumMagboltz::ElectronCollision(const double e, int& type, 
    int& level, double& e1, double& dx, double& dy, double& dz, 
    std::vector<std::pair<Particle, double> >& secondaries, int& ndxc,
    int& band) {
  ndxc = 0;
  if (e <= 0.) {
    std::cerr << m_className << "::ElectronCollision: Invalid energy.\n";
    return false;
  }
  // Check if the electron energy is within the currently set range.
  if (e > m_eMax && m_useAutoAdjust) {
    std::cerr << m_className << "::ElectronCollision:\n    Provided energy ("
              << e << " eV) exceeds current energy range.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  // If necessary, update the collision rates table.
  if (!Update()) return false;

  if (m_debug && band > 0) {
    std::cerr << m_className << "::ElectronCollision: Band > 0.\n";
  }

  double angCut = 1.;
  double angPar = 0.5;

  if (e <= m_eHigh) {
    // Linear binning
    // Get the energy interval.
    constexpr int iemax = Magboltz::nEnergySteps - 1;
    const int iE = std::min(std::max(int(e / m_eStep), 0), iemax);

    // Sample the scattering process.
    const double r = RndmUniform();
    if (r <= m_cf[iE][0]) {
      level = 0;
    } else if (r >= m_cf[iE][m_nTerms - 1]) {
      level = m_nTerms - 1;
    } else {
      const auto begin = m_cf[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nTerms, r) - begin;
    }
    // Get the angular distribution parameters.
    angCut = m_scatCut[iE][level];
    angPar = m_scatPar[iE][level];
  } else {
    // Logarithmic binning
    // Get the energy interval.
    const int iE = std::min(std::max(int(log(e / m_eHigh) / m_lnStep), 0),
                            nEnergyStepsLog - 1);
    // Sample the scattering process.
    const double r = RndmUniform();
    if (r <= m_cfLog[iE][0]) {
      level = 0;
    } else if (r >= m_cfLog[iE][m_nTerms - 1]) {
      level = m_nTerms - 1;
    } else {
      const auto begin = m_cfLog[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nTerms, r) - begin;
    }
    // Get the angular distribution parameters.
    angCut = m_scatCutLog[iE][level];
    angPar = m_scatParLog[iE][level];
  }

  // Extract the collision type.
  type = m_csType[level] % nCsTypes;
  const int igas = int(m_csType[level] / nCsTypes);
  // Increase the collision counters.
  ++m_nCollisions[type];
  ++m_nCollisionsDetailed[level];

  // Get the energy loss for this process.
  double loss = m_energyLoss[level];

  if (type == ElectronCollisionTypeVirtual) return true;

  if (type == ElectronCollisionTypeIonisation) {
    // Sample the secondary electron energy according to
    // the Opal-Beaty-Peterson parameterisation.
    double esec = 0.;
    if (e < loss) loss = e - 0.0001;
    if (m_useOpalBeaty) {
      // Get the splitting parameter.
      const double w = m_wOpalBeaty[level];
      esec = w * tan(RndmUniform() * atan(0.5 * (e - loss) / w));
      // Rescaling (SST)
      // esec = w * pow(esec / w, 0.9524);
    } else if (m_useGreenSawada) {
      const double gs = m_parGreenSawada[igas][0];
      const double gb = m_parGreenSawada[igas][1];
      const double w = gs * e / (e + gb);
      const double ts = m_parGreenSawada[igas][2];
      const double ta = m_parGreenSawada[igas][3];
      const double tb = m_parGreenSawada[igas][4];
      const double esec0 = ts - ta / (e + tb);
      const double r = RndmUniform();
      esec = esec0 +
             w * tan((r - 1.) * atan(esec0 / w) +
                     r * atan((0.5 * (e - loss) - esec0) / w));
    } else {
      esec = RndmUniform() * (e - loss);
    }
    if (esec <= 0) esec = Small;
    loss += esec;
    // Add the secondary electron.
    secondaries.emplace_back(std::make_pair(Particle::Electron, esec));
    // Add the ion.
    secondaries.emplace_back(std::make_pair(Particle::Ion, 0.));
    bool fluorescence = false;
    if (m_yFluorescence[level] > Small) {
      if (RndmUniform() < m_yFluorescence[level]) fluorescence = true;
    } 
    // Add Auger and photo electrons (if any).
    if (fluorescence) {
      if (m_nAuger2[level] > 0) {
        const double eav = m_eAuger2[level] / m_nAuger2[level];
        for (unsigned int i = 0; i < m_nAuger2[level]; ++i) {
          secondaries.emplace_back(std::make_pair(Particle::Electron, eav));
        }
      }
      if (m_nFluorescence[level] > 0) {
        const double eav = m_eFluorescence[level] / m_nFluorescence[level];
        for (unsigned int i = 0; i < m_nFluorescence[level]; ++i) {
          secondaries.emplace_back(std::make_pair(Particle::Electron, eav));
        }
      }
    } else if (m_nAuger1[level] > 0) {
      const double eav = m_eAuger1[level] / m_nAuger1[level];
      for (unsigned int i = 0; i < m_nAuger1[level]; ++i) {
        secondaries.emplace_back(std::make_pair(Particle::Electron, eav));
      }
    } 
  } else if (type == ElectronCollisionTypeExcitation) {
    // Follow the de-excitation cascade (if switched on).
    if (m_useDeexcitation && m_iDeexcitation[level] >= 0) {
      int fLevel = 0;
      ComputeDeexcitationInternal(m_iDeexcitation[level], fLevel);
      ndxc = m_dxcProducts.size();
    } else if (m_usePenning) {
      m_dxcProducts.clear();
      // Simplified treatment of Penning ionisation.
      // If the energy threshold of this level exceeds the
      // ionisation potential of one of the gases,
      // create a new electron (with probability rPenning).
      if (m_debug) {
        std::cout << m_className << "::ElectronCollision:\n"
                  << "    Level: " << level << "\n"
                  << "    Ionization potential: " << m_minIonPot << "\n"
                  << "    Excitation energy: " << loss * m_rgas[igas] << "\n"
                  << "    Penning probability: " << m_rPenning[level] << "\n"; 
      }
      if (loss * m_rgas[igas] > m_minIonPot &&
          RndmUniform() < m_rPenning[level]) {
        // The energy of the secondary electron is assumed to be given by
        // the difference of excitation and ionisation threshold.
        double esec = loss * m_rgas[igas] - m_minIonPot;
        if (esec <= 0) esec = Small;
        // Add the secondary electron to the list.
        dxcProd newDxcProd;
        newDxcProd.t = 0.;
        newDxcProd.s = 0.;
        if (m_lambdaPenning[level] > Small) {
          // Uniform distribution within a sphere of radius lambda
          newDxcProd.s = m_lambdaPenning[level] * std::cbrt(RndmUniformPos());
        }
        newDxcProd.energy = esec;
        newDxcProd.type = DxcProdTypeElectron;
        m_dxcProducts.push_back(std::move(newDxcProd));
        ndxc = 1;
        ++m_nPenning;
      }
    }
  }

  if (e < loss) loss = e - 0.0001;

  // Determine the scattering angle.
  double ctheta0 = 1. - 2. * RndmUniform();
  if (m_useAnisotropic) {
    switch (m_scatModel[level]) {
      case 0:
        break;
      case 1:
        ctheta0 = 1. - RndmUniform() * angCut;
        if (RndmUniform() > angPar) ctheta0 = -ctheta0;
        break;
      case 2:
        ctheta0 = (ctheta0 + angPar) / (1. + angPar * ctheta0);
        break;
      default:
        std::cerr << m_className << "::ElectronCollision:\n"
                  << "    Unknown scattering model.\n"
                  << "    Using isotropic distribution.\n";
        break;
    }
  }

  const double s1 = m_rgas[igas];
  const double s2 = (s1 * s1) / (s1 - 1.);
  const double theta0 = acos(ctheta0);
  const double arg = std::max(1. - s1 * loss / e, Small);
  const double d = 1. - ctheta0 * sqrt(arg);

  // Update the energy.
  e1 = std::max(e * (1. - loss / (s1 * e) - 2. * d / s2), Small);
  double q = std::min(sqrt((e / e1) * arg) / s1, 1.);
  const double theta = asin(q * sin(theta0));
  double ctheta = cos(theta);
  if (ctheta0 < 0.) {
    const double u = (s1 - 1.) * (s1 - 1.) / arg;
    if (ctheta0 * ctheta0 > u) ctheta = -ctheta;
  }
  const double stheta = sin(theta);
  // Calculate the direction after the collision.
  dz = std::min(dz, 1.);
  const double argZ = sqrt(dx * dx + dy * dy);

  // Azimuth is chosen at random.
  const double phi = TwoPi * RndmUniform();
  const double cphi = cos(phi);
  const double sphi = sin(phi);
  if (argZ == 0.) {
    dz = ctheta;
    dx = cphi * stheta;
    dy = sphi * stheta;
  } else {
    const double a = stheta / argZ;
    const double dz1 = dz * ctheta + argZ * stheta * sphi;
    const double dy1 = dy * ctheta + a * (dx * cphi - dy * dz * sphi);
    const double dx1 = dx * ctheta - a * (dy * cphi + dx * dz * sphi);
    dz = dz1;
    dy = dy1;
    dx = dx1;
  }

  return true;
}

bool MediumMagboltz::GetDeexcitationProduct(const unsigned int i, double& t,
                                            double& s, int& type,
                                            double& energy) const {
  if (i >= m_dxcProducts.size() || !(m_useDeexcitation || m_usePenning)) {
    return false;
  }
  t = m_dxcProducts[i].t;
  s = m_dxcProducts[i].s;
  type = m_dxcProducts[i].type;
  energy = m_dxcProducts[i].energy;
  return true;
}

double MediumMagboltz::GetPhotonCollisionRate(const double e) {
  if (e <= 0.) {
    std::cerr << m_className << "::GetPhotonCollisionRate: Invalid  energy.\n";
    return m_cfTotGamma[0];
  }
  if (e > m_eFinalGamma && m_useAutoAdjust) {
    std::cerr << m_className << "::GetPhotonCollisionRate:\n    Rate at " << e
              << " eV is not included in the current table.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxPhotonEnergy(1.05 * e);
  }

  if (!Update()) return 0.;

  const int iE =
      std::min(std::max(int(e / m_eStepGamma), 0), nEnergyStepsGamma - 1);

  double cfSum = m_cfTotGamma[iE];
  if (m_useDeexcitation && m_useRadTrap && !m_deexcitations.empty()) {
    // Loop over the excitations.
    for (const auto& dxc : m_deexcitations) {
      if (dxc.cf > 0. && fabs(e - dxc.energy) <= dxc.width) {
        cfSum += dxc.cf *
                 TMath::Voigt(e - dxc.energy, dxc.sDoppler, 2 * dxc.gPressure);
      }
    }
  }

  return cfSum;
}

bool MediumMagboltz::GetPhotonCollision(const double e, int& type, int& level,
                                        double& e1, double& ctheta, int& nsec,
                                        double& esec) {
  if (e <= 0.) {
    std::cerr << m_className << "::GetPhotonCollision: Invalid energy.\n";
    return false;
  }
  if (e > m_eFinalGamma && m_useAutoAdjust) {
    std::cerr << m_className << "::GetPhotonCollision:\n    Provided energy ("
              << e << " eV) exceeds current energy range.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxPhotonEnergy(1.05 * e);
  }

  if (!Update()) return false;

  // Energy interval
  const int iE =
      std::min(std::max(int(e / m_eStepGamma), 0), nEnergyStepsGamma - 1);

  double r = m_cfTotGamma[iE];
  if (m_useDeexcitation && m_useRadTrap && !m_deexcitations.empty()) {
    int nLines = 0;
    std::vector<double> pLine(0);
    std::vector<int> iLine(0);
    // Loop over the excitations.
    const unsigned int nDeexcitations = m_deexcitations.size();
    for (unsigned int i = 0; i < nDeexcitations; ++i) {
      const auto& dxc = m_deexcitations[i];
      if (dxc.cf > 0. && fabs(e - dxc.energy) <= dxc.width) {
        r += dxc.cf *
             TMath::Voigt(e - dxc.energy, dxc.sDoppler, 2 * dxc.gPressure);
        pLine.push_back(r);
        iLine.push_back(i);
        ++nLines;
      }
    }
    r *= RndmUniform();
    if (nLines > 0 && r >= m_cfTotGamma[iE]) {
      // Photon is absorbed by a discrete line.
      for (int i = 0; i < nLines; ++i) {
        if (r <= pLine[i]) {
          ++m_nPhotonCollisions[PhotonCollisionTypeExcitation];
          int fLevel = 0;
          ComputeDeexcitationInternal(iLine[i], fLevel);
          type = PhotonCollisionTypeExcitation;
          nsec = m_dxcProducts.size();
          return true;
        }
      }
      std::cerr << m_className << "::GetPhotonCollision:\n";
      std::cerr << "    Random sampling of deexcitation line failed.\n";
      std::cerr << "    Program bug!\n";
      return false;
    }
  } else {
    r *= RndmUniform();
  }

  if (r <= m_cfGamma[iE][0]) {
    level = 0;
  } else if (r >= m_cfGamma[iE][m_nPhotonTerms - 1]) {
    level = m_nPhotonTerms - 1;
  } else {
    const auto begin = m_cfGamma[iE].cbegin();
    level = std::lower_bound(begin, begin + m_nPhotonTerms, r) - begin;
  }

  nsec = 0;
  esec = e1 = 0.;
  type = csTypeGamma[level];
  // Collision type
  type = type % nCsTypesGamma;
  int ngas = int(csTypeGamma[level] / nCsTypesGamma);
  ++m_nPhotonCollisions[type];
  // Ionising collision
  if (type == 1) {
    esec = std::max(e - m_ionPot[ngas], Small);
    nsec = 1;
  }

  // Determine the scattering angle
  ctheta = 2 * RndmUniform() - 1.;

  return true;
}

void MediumMagboltz::ResetCollisionCounters() {
  m_nCollisions.fill(0);
  m_nCollisionsDetailed.assign(m_nTerms, 0);
  m_nPenning = 0;
  m_nPhotonCollisions.fill(0);
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions() const {
  return std::accumulate(std::begin(m_nCollisions), std::end(m_nCollisions), 0);
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions(
    unsigned int& nElastic, unsigned int& nIonisation,
    unsigned int& nAttachment, unsigned int& nInelastic,
    unsigned int& nExcitation, unsigned int& nSuperelastic) const {
  nElastic = m_nCollisions[ElectronCollisionTypeElastic];
  nIonisation = m_nCollisions[ElectronCollisionTypeIonisation];
  nAttachment = m_nCollisions[ElectronCollisionTypeAttachment];
  nInelastic = m_nCollisions[ElectronCollisionTypeInelastic];
  nExcitation = m_nCollisions[ElectronCollisionTypeExcitation];
  nSuperelastic = m_nCollisions[ElectronCollisionTypeSuperelastic];
  return nElastic + nIonisation + nAttachment + nInelastic + nExcitation +
         nSuperelastic;
}

unsigned int MediumMagboltz::GetNumberOfLevels() {
  if (!Update()) return 0;
  return m_nTerms;
}

bool MediumMagboltz::GetLevel(const unsigned int i, int& ngas, int& type,
                              std::string& descr, double& e) {
  if (!Update()) return false;

  if (i >= m_nTerms) {
    std::cerr << m_className << "::GetLevel: Index out of range.\n";
    return false;
  }

  // Collision type
  type = m_csType[i] % nCsTypes;
  ngas = int(m_csType[i] / nCsTypes);
  // Description (from Magboltz)
  descr = m_description[i];
  // Threshold energy
  e = m_rgas[ngas] * m_energyLoss[i];
  if (m_debug) {
    std::cout << m_className << "::GetLevel:\n"
              << "    Level " << i << ": " << descr << "\n"
              << "    Type " << type << "\n"
              << "    Threshold energy: " << e << " eV\n";
    if (type == ElectronCollisionTypeExcitation && m_usePenning &&
        e > m_minIonPot) {
      std::cout << "    Penning transfer coefficient: " << m_rPenning[i]
                << "\n";
    } else if (type == ElectronCollisionTypeExcitation && m_useDeexcitation) {
      const int idxc = m_iDeexcitation[i];
      if (idxc < 0 || idxc >= (int)m_deexcitations.size()) {
        std::cout << "    Deexcitation cascade not implemented.\n";
        return true;
      }
      const auto& dxc = m_deexcitations[idxc];
      if (dxc.osc > 0.) {
        std::cout << "    Oscillator strength: " << dxc.osc << "\n";
      }
      std::cout << "    Decay channels:\n";
      const int nChannels = dxc.type.size();
      for (int j = 0; j < nChannels; ++j) {
        if (dxc.type[j] == DxcTypeRad) {
          std::cout << "      Radiative decay to ";
          if (dxc.final[j] < 0) {
            std::cout << "ground state: ";
          } else {
            std::cout << m_deexcitations[dxc.final[j]].label << ": ";
          }
        } else if (dxc.type[j] == DxcTypeCollIon) {
          if (dxc.final[j] < 0) {
            std::cout << "      Penning ionisation: ";
          } else {
            std::cout << "      Associative ionisation: ";
          }
        } else if (dxc.type[j] == DxcTypeCollNonIon) {
          if (dxc.final[j] >= 0) {
            std::cout << "      Collision-induced transition to "
                      << m_deexcitations[dxc.final[j]].label << ": ";
          } else {
            std::cout << "      Loss: ";
          }
        }
        const double br = j == 0 ? dxc.p[j] : dxc.p[j] - dxc.p[j - 1];
        std::cout << std::setprecision(5) << br * 100. << "%\n";
      }
    }
  }

  return true;
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions(
    const unsigned int level) const {
  if (level >= m_nTerms) {
    std::cerr << m_className << "::GetNumberOfElectronCollisions: "
              << "Level " << level << " does not exist.\n";
    return 0;
  }
  return m_nCollisionsDetailed[level];
}

unsigned int MediumMagboltz::GetNumberOfPhotonCollisions() const {
  return std::accumulate(std::begin(m_nPhotonCollisions),
                         std::end(m_nPhotonCollisions), 0);
}

unsigned int MediumMagboltz::GetNumberOfPhotonCollisions(
    unsigned int& nElastic, unsigned int& nIonising,
    unsigned int& nInelastic) const {
  nElastic = m_nPhotonCollisions[0];
  nIonising = m_nPhotonCollisions[1];
  nInelastic = m_nPhotonCollisions[2];
  return nElastic + nIonising + nInelastic;
}

int MediumMagboltz::GetGasNumberMagboltz(const std::string& input) {

  if (input.empty()) return 0;

  if (input == "CF4") {
    return 1;
  } else if (input == "Ar") {
    return 2;
  } else if (input == "He" || input == "He-4") {
    // Helium 4
    return 3;
  } else if (input == "He-3") {
    // Helium 3
    return 4;
  } else if (input == "Ne") {
    return 5;
  } else if (input == "Kr") {
    return 6;
  } else if (input == "Xe") {
    return 7;
  } else if (input == "CH4") {
    // Methane
    return 8;
  } else if (input == "C2H6") {
    // Ethane
    return 9;
  } else if (input == "C3H8") {
    // Propane
    return 10;
  } else if (input == "iC4H10") {
    // Isobutane
    return 11;
  } else if (input == "CO2") {
    return 12;
  } else if (input == "neoC5H12") {
    // Neopentane
    return 13;
  } else if (input == "H2O") {
    return 14;
  } else if (input == "O2") {
    return 15;
  } else if (input == "N2") {
    return 16;
  } else if (input == "NO") {
    // Nitric oxide (NO)
    return 17;
  } else if (input == "N2O") {
    // Nitrous oxide (N2O)
    return 18;
  } else if (input == "C2H4") {
    // Ethene (C2H4)
    return 19;
  } else if (input == "C2H2") {
    // Acetylene (C2H2)
    return 20;
  } else if (input == "H2") {
    // Hydrogen
    return 21;
  } else if (input == "D2") {
    // Deuterium
    return 22;
  } else if (input == "CO") {
    // Carbon monoxide (CO)
    return 23;
  } else if (input == "Methylal") {
    // Methylal (dimethoxymethane, CH3-O-CH2-O-CH3, "hot" version)
    return 24;
  } else if (input == "DME") {
    return 25;
  } else if (input == "Reid-Step") {
    return 26;
  } else if (input == "Maxwell-Model") {
    return 27;
  } else if (input == "Reid-Ramp") {
    return 28;
  } else if (input == "C2F6") {
    return 29;
  } else if (input == "SF6") {
    return 30;
  } else if (input == "NH3") {
    return 31;
  } else if (input == "C3H6") {
    // Propene
    return 32;
  } else if (input == "cC3H6") {
    // Cyclopropane
    return 33;
  } else if (input == "CH3OH") {
    // Methanol
    return 34;
  } else if (input == "C2H5OH") {
    // Ethanol
    return 35;
  } else if (input == "C3H7OH") {
    // Propanol
    return 36;
  } else if (input == "Cs") {
    return 37;
  } else if (input == "F2") {
    // Fluorine
    return 38;
  } else if (input == "CS2") {
    return 39;
  } else if (input == "COS") {
    return 40;
  } else if (input == "CD4") {
    // Deuterated methane
    return 41;
  } else if (input == "BF3") {
    return 42;
  } else if (input == "C2HF5" || input == "C2H2F4") {
    return 43;
  } else if (input == "TMA") {
    return 44;
  } else if (input == "nC3H7OH") {
    // n-propanol
    return 46;
  } else if (input == "paraH2") {
    // Para hydrogen
    return 47;
  } else if (input == "orthoD2") {
    // Ortho deuterium
    return 48;
  } else if (input == "CHF3") {
    return 50;
  } else if (input == "CF3Br") {
    return 51;
  } else if (input == "C3F8") {
    return 52;
  } else if (input == "O3") {
    // Ozone
    return 53;
  } else if (input == "Hg") {
    // Mercury
    return 54;
  } else if (input == "H2S") {
    return 55;
  } else if (input == "nC4H10") {
    // n-Butane
    return 56;
  } else if (input == "nC5H12") {
    // n-Pentane
    return 57;
  } else if (input == "N2 (Phelps)") {
    return 58;
  } else if (input == "GeH4") {
    // Germane, GeH4
    return 59;
  } else if (input == "SiH4") {
    // Silane, SiH4
    return 60;
  } else if (input == "CCl4") {
    return 61;
  }

  std::cerr << "MediumMagboltz::GetGasNumberMagboltz:\n"
            << "    Gas " << input << " is not defined.\n";
  return 0;
}

bool MediumMagboltz::Update(const bool verbose) {

  std::lock_guard<std::mutex> guard(m_mutex);
  if (!m_isChanged) return true;
  if (!Mixer(verbose)) {
    std::cerr << m_className 
              << "::Update: Error calculating the collision rates table.\n";
    return false;
  }
  m_isChanged = false;
  return true;
}


bool MediumMagboltz::Mixer(const bool verbose) {

  // Set constants and parameters in Magboltz common blocks.
  Magboltz::cnsts_.echarg = ElementaryCharge * 1.e-15;
  Magboltz::cnsts_.emass = ElectronMassGramme;
  Magboltz::cnsts_.amu = AtomicMassUnit;
  Magboltz::cnsts_.pir2 = BohrRadius * BohrRadius * Pi;
  Magboltz::inpt_.ary = RydbergEnergy;

  Magboltz::inpt_.akt = BoltzmannConstant * m_temperature;
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;

  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = Magboltz::nEnergySteps;
  Magboltz::inpt_.nAniso = m_useAnisotropic ? 2 : 0;

  for (unsigned int i = 0; i < Magboltz::nEnergySteps; ++i) {
    const double en = (i + 0.5) * m_eStep;
    Magboltz::mix2_.eg[i] = en;
    Magboltz::mix2_.eroot[i] = sqrt(en);
    Magboltz::dens_.den[i] = 0.;
  }
  constexpr int iemax = Magboltz::nEnergySteps - 1;

  // Calculate the atomic density (ideal gas law).
  const double dens = GetNumberDensity();
  // Prefactor for calculation of scattering rate from cross-section.
  const double prefactor = dens * SpeedOfLight * sqrt(2. / ElectronMass);

  m_rgas.fill(1.);

  m_ionPot.fill(-1.);
  m_minIonPot = -1.;

  m_parGreenSawada.fill({1., 0., 0., 0., 0.});
  m_hasGreenSawada.fill(false);

  m_wOpalBeaty.fill(1.);
  m_energyLoss.fill(0.);
  m_csType.fill(0);

  m_yFluorescence.fill(0.);
  m_nAuger1.fill(0);
  m_eAuger1.fill(0.);
  m_nAuger2.fill(0);
  m_eAuger2.fill(0.);
  m_nFluorescence.fill(0);
  m_eFluorescence.fill(0.);

  m_scatModel.fill(0);

  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);

  m_deexcitations.clear();
  m_iDeexcitation.fill(-1);

  // Reset the collision rates.
  m_cfTot.assign(Magboltz::nEnergySteps, 0.);
  m_cfTotLog.assign(nEnergyStepsLog, 0.);

  m_cf.assign(Magboltz::nEnergySteps,
              std::vector<double>(Magboltz::nMaxLevels, 0.));
  m_cfLog.assign(nEnergyStepsLog,
                 std::vector<double>(Magboltz::nMaxLevels, 0.));

  m_scatPar.assign(Magboltz::nEnergySteps,
                   std::vector<double>(Magboltz::nMaxLevels, 0.5));
  m_scatCut.assign(Magboltz::nEnergySteps,
                   std::vector<double>(Magboltz::nMaxLevels, 1.));

  m_scatParLog.assign(nEnergyStepsLog,
                      std::vector<double>(Magboltz::nMaxLevels, 0.5));
  m_scatCutLog.assign(nEnergyStepsLog,
                      std::vector<double>(Magboltz::nMaxLevels, 1.));

  // Cross-sections
  // 0: total, 1: elastic,
  // 2: ionisation, 3: attachment,
  // 4, 5: unused
  static double q[Magboltz::nEnergySteps][6];
  // Inelastic cross-sections
  static double qIn[Magboltz::nEnergySteps][Magboltz::nMaxInelasticTerms];
  // Ionisation cross-sections
  static double qIon[Magboltz::nEnergySteps][Magboltz::nMaxIonisationTerms];
  // Attachment cross-sections
  static double qAtt[Magboltz::nEnergySteps][Magboltz::nMaxAttachmentTerms];
  // "Null-collision" cross-sections
  static double qNull[Magboltz::nEnergySteps][Magboltz::nMaxNullTerms];
  // Parameters for scattering angular distribution
  static double pEqEl[Magboltz::nEnergySteps][6];
  // Parameters for angular distribution in inelastic collisions
  static double pEqIn[Magboltz::nEnergySteps][Magboltz::nMaxInelasticTerms];
  // Parameters for angular distribution in ionising collisions
  static double pEqIon[Magboltz::nEnergySteps][Magboltz::nMaxIonisationTerms];
  // Penning transfer parameters
  static double penFra[Magboltz::nMaxInelasticTerms][3];
  // Description of cross-section terms
  static char scrpt[Magboltz::nMaxLevelsPerComponent][Magboltz::nCharDescr];
  // Description of "null-collision" cross-section terms
  static char scrptn[Magboltz::nMaxNullTerms][Magboltz::nCharDescr];

  // Check the gas composition and establish the gas numbers.
  int gasNumber[m_nMaxGases];
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const int ng = GetGasNumberMagboltz(m_gas[i]);
    if (ng <= 0) {
      std::cerr << m_className << "::Mixer:\n    Gas " << m_gas[i]
                << " does not have a gas number in Magboltz.\n";
      return false;
    }
    gasNumber[i] = ng;
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n    " << Magboltz::nEnergySteps
              << " linear energy steps between 0 and "
              << std::min(m_eMax, m_eHigh) << " eV.\n";
    if (m_eMax > m_eHigh) {
      std::cout << "    " << nEnergyStepsLog << " logarithmic steps between "
                << m_eHigh << " and " << m_eMax << " eV\n";
    }
  }
  m_nTerms = 0;

  std::ofstream outfile;
  if (m_useCsOutput) {
    outfile.open("cs.txt", std::ios::out);
    outfile << "# energy [eV] vs. cross-section [cm2]\n";
  }

  // Loop over the gases in the mixture.
  for (unsigned int iGas = 0; iGas < m_nComponents; ++iGas) {
    Magboltz::inpt_.efinal = std::min(m_eMax, m_eHigh);
    Magboltz::inpt_.estep = m_eStep;
    Magboltz::mix2_.eg[iemax] = (iemax + 0.5) * m_eStep;
    Magboltz::mix2_.eroot[iemax] = sqrt((iemax + 0.5) * m_eStep);
    char name[Magboltz::nCharName];
    // Number of inelastic cross-sections
    static std::int64_t nIn = 0;
    // Number of ionisation cross-sections
    static std::int64_t nIon = 0;
    // Number of attachment cross-sections
    static std::int64_t nAtt = 0;
    // Number of "null-collision" cross-sections
    static std::int64_t nNull = 0;
    // Virial coefficient (not used).
    static double virial = 0.;
    // Thresholds/characteristic energies.
    static double e[6];
    // Energy losses for inelastic cross-sections.
    static double eIn[Magboltz::nMaxInelasticTerms];
    // Ionisation thresholds.
    static double eIon[Magboltz::nMaxIonisationTerms];
    // Scattering algorithms
    static std::int64_t kIn[Magboltz::nMaxInelasticTerms];
    static std::int64_t kEl[6];
    // Opal-Beaty parameter
    static double eoby[Magboltz::nMaxIonisationTerms];
    // Scaling factor for "null-collision" terms
    static double scln[Magboltz::nMaxNullTerms];
    // Parameters for simulation of Auger and fluorescence processes.
    static std::int64_t nc0[Magboltz::nMaxIonisationTerms];
    static std::int64_t ng1[Magboltz::nMaxIonisationTerms];
    static std::int64_t ng2[Magboltz::nMaxIonisationTerms];
    static double ec0[Magboltz::nMaxIonisationTerms];
    static double wklm[Magboltz::nMaxIonisationTerms];
    static double efl[Magboltz::nMaxIonisationTerms];
    static double eg1[Magboltz::nMaxIonisationTerms];
    static double eg2[Magboltz::nMaxIonisationTerms];
    // Retrieve the cross-section data for this gas from Magboltz.
    std::int64_t ngs = gasNumber[iGas];
    Magboltz::gasmix_(
        &ngs, q[0], qIn[0], &nIn, &e[0], eIn, name, &virial, eoby, pEqEl[0], 
        pEqIn[0], penFra[0], kEl, kIn, qIon[0], pEqIon[0], eIon, &nIon, 
        qAtt[0], &nAtt, qNull[0], &nNull, scln, nc0, ec0, wklm, efl,
        ng1, eg1, ng2, eg2, scrpt, scrptn,  
        Magboltz::nCharName, Magboltz::nCharDescr, Magboltz::nCharDescr);
    if (m_debug || verbose) {
      const double m = (2. / e[1]) * ElectronMass / AtomicMassUnitElectronVolt;
      std::cout << "    " << name << "\n"
                << "      mass: " << m << " amu\n";
      if (nIon > 1) {
        std::cout << "      ionisation threshold: " << eIon[0] << " eV\n";
      } else {
        std::cout << "      ionisation threshold: " << e[2] << " eV\n";
      }
      if (e[3] > 0. && e[4] > 0.) {
        std::cout << "      cross-sections at minimum ionising energy:\n"
                  << "        excitation: " << e[3] * 1.e18 << " Mbarn\n"
                  << "        ionisation: " << e[4] * 1.e18 << " Mbarn\n";
      }
    }
    int np0 = m_nTerms;
    // Make sure there is still sufficient space.
    if (np0 + nIn + nIon + nAtt + nNull >= Magboltz::nMaxLevels) {
      std::cerr << m_className << "::Mixer:\n"
                << "    Max. number of levels (" << Magboltz::nMaxLevels
                << ") exceeded.\n";
      return false;
    }
    const double van = m_fraction[iGas] * prefactor;
    int np = np0;
    if (m_useCsOutput) {
      outfile << "# cross-sections for " << name << "\n";
      outfile << "# cross-section types:\n";
      outfile << "# elastic\n";
    }
    // Elastic scattering
    ++m_nTerms;
    m_scatModel[np] = kEl[1];
    const double r = 1. + 0.5 * e[1];
    m_rgas[iGas] = r;
    m_energyLoss[np] = 0.;
    m_description[np] = GetDescription(1, scrpt);
    m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeElastic;
    bool withIon = false;
    // Ionisation
    if (nIon > 1) {
      for (int j = 0; j < nIon; ++j) {
        if (m_eMax < eIon[j]) continue;
        withIon = true;
        ++m_nTerms;
        ++np;
        m_scatModel[np] = kEl[2];
        m_energyLoss[np] = eIon[j] / r;
        m_wOpalBeaty[np] = eoby[j];
        m_yFluorescence[np] = wklm[j];
        m_nAuger1[np] = nc0[j];
        m_eAuger1[np] = ec0[j];
        m_nFluorescence[np] = ng1[j];
        m_eFluorescence[np] = eg1[j];
        m_nAuger2[np] = ng2[j];
        m_eAuger2[np] = eg2[j];
        m_description[np] = GetDescription(2 + j, scrpt);
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeIonisation;
        if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
      }
      m_parGreenSawada[iGas][0] = eoby[0];
      m_parGreenSawada[iGas][4] = 2 * eIon[0];
      m_ionPot[iGas] = eIon[0];
    } else {
      if (m_eMax >= e[2]) {
        withIon = true;
        ++m_nTerms;
        ++np;
        m_scatModel[np] = kEl[2];
        m_energyLoss[np] = e[2] / r;
        m_wOpalBeaty[np] = eoby[0];
        m_parGreenSawada[iGas][0] = eoby[0];
        m_parGreenSawada[iGas][4] = 2 * e[2];
        m_ionPot[iGas] = e[2];
        m_description[np] = GetDescription(2, scrpt);
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeIonisation;
        if (m_useCsOutput) outfile << "# ionisation (gross)\n";
      }
    }
    // Attachment
    for (int j = 0; j < nAtt; ++j) {
      ++m_nTerms;
      ++np;
      m_scatModel[np] = 0;
      m_energyLoss[np] = 0.;
      m_description[np] = GetDescription(2 + nIon + j, scrpt);
      m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeAttachment;
      if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
    }
    // Inelastic terms
    int nExc = 0, nSuperEl = 0;
    for (int j = 0; j < nIn; ++j) {
      ++np;
      m_scatModel[np] = kIn[j];
      m_energyLoss[np] = eIn[j] / r;
      m_description[np] = GetDescription(4 + nIon + nAtt + j, scrpt);
      if ((m_description[np][1] == 'E' && m_description[np][2] == 'X') ||
          (m_description[np][0] == 'E' && m_description[np][1] == 'X') ||
          (m_gas[iGas] == "N2" && eIn[j] > 6.)) {
        // Excitation
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeExcitation;
        ++nExc;
      } else if (eIn[j] < 0.) {
        // Super-elastic collision
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeSuperelastic;
        ++nSuperEl;
      } else {
        // Inelastic collision
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeInelastic;
      }
      if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
    }
    m_nTerms += nIn;
    if (nNull > 0) {
      for (int j = 0; j < nNull; ++j) {
        ++m_nTerms;
        ++np;
        m_scatModel[np] = 0;
        m_energyLoss[np] = 0.;
        m_description[np] = GetDescription(j, scrptn);
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeVirtual;
        if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
      }
    }
    // Loop over the energy table.
    for (unsigned int iE = 0; iE < Magboltz::nEnergySteps; ++iE) {
      np = np0;
      if (m_useCsOutput) {
        outfile << (iE + 0.5) * m_eStep << "  " << q[iE][1] << "  ";
      }
      // Elastic scattering
      m_cf[iE][np] = q[iE][1] * van;
      SetScatteringParameters(m_scatModel[np], pEqEl[iE][1], m_scatCut[iE][np],
                              m_scatPar[iE][np]);
      // Ionisation
      if (withIon) {
        if (nIon > 1) {
          for (int j = 0; j < nIon; ++j) {
            if (m_eMax < eIon[j]) continue;
            ++np;
            m_cf[iE][np] = qIon[iE][j] * van;
            SetScatteringParameters(m_scatModel[np], pEqIon[iE][j],
                                    m_scatCut[iE][np], m_scatPar[iE][np]);
            if (m_useCsOutput) outfile << qIon[iE][j] << "  ";
          }
        } else {
          ++np;
          m_cf[iE][np] = q[iE][2] * van;
          SetScatteringParameters(m_scatModel[np], pEqEl[iE][2],
                                  m_scatCut[iE][np], m_scatPar[iE][np]);
          if (m_useCsOutput) outfile << q[iE][2] << "  ";
        }
      }
      // Attachment
      for (int j = 0; j < nAtt; ++j) {
        ++np;
        m_cf[iE][np] = qAtt[iE][j] * van;
        // m_cf[iE][np] = q[iE][3] * van;
        m_scatPar[iE][np] = 0.5;
        if (m_useCsOutput) outfile << qAtt[iE][j] << "  ";
      }
      // Inelastic terms
      for (int j = 0; j < nIn; ++j) {
        ++np;
        if (m_useCsOutput) outfile << qIn[iE][j] << "  ";
        m_cf[iE][np] = qIn[iE][j] * van;
        // Scale the excitation cross-sections (for error estimates).
        m_cf[iE][np] *= m_scaleExc[iGas];
        if (m_cf[iE][np] < 0.) {
          std::cerr << m_className << "::Mixer:\n"
                    << "    Negative inelastic cross-section at "
                    << (iE + 0.5) * m_eStep << " eV. Set to zero.\n";
          m_cf[iE][np] = 0.;
        }
        SetScatteringParameters(m_scatModel[np], pEqIn[iE][j],
                                m_scatCut[iE][np], m_scatPar[iE][np]);
      }
      if ((m_debug || verbose) && nIn > 0 && iE == iemax) {
        std::cout << "      " << nIn << " inelastic terms (" << nExc
                  << " excitations, " << nSuperEl << " superelastic, "
                  << nIn - nExc - nSuperEl << " other)\n";
      }
      if (nNull > 0) {
        for (int j = 0; j < nNull; ++j) {
          ++np;
          m_cf[iE][np] = qNull[iE][j] * van * scln[j];
          if (m_useCsOutput) outfile << qNull[iE][j] << "  ";
        }
      }
      if (m_useCsOutput) outfile << "\n";
    }
    if (m_eMax <= m_eHigh) continue;
    // Fill the high-energy part (logarithmic binning).
    // Calculate the growth factor.
    const double rLog = pow(m_eMax / m_eHigh, 1. / nEnergyStepsLog);
    m_lnStep = log(rLog);
    // Set the upper limit of the first bin.
    double emax = m_eHigh * rLog;

    for (int iE = 0; iE < nEnergyStepsLog; ++iE) {
      Magboltz::inpt_.estep = emax / (Magboltz::nEnergySteps - 0.5);
      Magboltz::inpt_.efinal = emax + 0.5 * Magboltz::inpt_.estep;
      Magboltz::mix2_.eg[iemax] = emax;
      Magboltz::mix2_.eroot[iemax] = sqrt(emax);
      Magboltz::gasmix_(
          &ngs, q[0], qIn[0], &nIn, e, eIn, name, &virial, eoby, pEqEl[0], 
          pEqIn[0], penFra[0], kEl, kIn, qIon[0], pEqIon[0], eIon, &nIon, 
          qAtt[0], &nAtt, qNull[0], &nNull, scln, nc0, ec0, wklm, efl,
          ng1, eg1, ng2, eg2, scrpt, scrptn, 
          Magboltz::nCharName, Magboltz::nCharDescr, Magboltz::nCharDescr);
      np = np0;
      if (m_useCsOutput) outfile << emax << "  " << q[iemax][1] << "  ";
      // Elastic scattering
      m_cfLog[iE][np] = q[iemax][1] * van;
      SetScatteringParameters(m_scatModel[np], pEqEl[iemax][1],
                              m_scatCutLog[iE][np], m_scatParLog[iE][np]);
      // Ionisation
      if (withIon) {
        if (nIon > 1) {
          for (int j = 0; j < nIon; ++j) {
            if (m_eMax < eIon[j]) continue;
            ++np;
            m_cfLog[iE][np] = qIon[iemax][j] * van;
            SetScatteringParameters(m_scatModel[np], pEqIon[iemax][j],
                                    m_scatCutLog[iE][np], m_scatParLog[iE][np]);
            if (m_useCsOutput) outfile << qIon[iemax][j] << "  ";
          }
        } else {
          ++np;
          // Gross cross-section
          m_cfLog[iE][np] = q[iemax][2] * van;
          // Counting cross-section
          // m_cfLog[iE][np] = q[iemax][4] * van;
          SetScatteringParameters(m_scatModel[np], pEqEl[iemax][2],
                                  m_scatCutLog[iE][np], m_scatParLog[iE][np]);
          if (m_useCsOutput) outfile << q[iemax][2] << "  ";
        }
      }
      // Attachment
      for (int j = 0; j < nAtt; ++j) {
        ++np;
        m_cfLog[iE][np] = qAtt[iemax][j] * van;
        // m_cfLog[iE][np] = q[iemax][3] * van;
        if (m_useCsOutput) outfile << qAtt[iemax][j] << "  ";
      }
      // Inelastic terms
      for (int j = 0; j < nIn; ++j) {
        ++np;
        if (m_useCsOutput) outfile << qIn[iemax][j] << "  ";
        m_cfLog[iE][np] = qIn[iemax][j] * van;
        // Scale the excitation cross-sections (for error estimates).
        m_cfLog[iE][np] *= m_scaleExc[iGas];
        if (m_cfLog[iE][np] < 0.) {
          std::cerr << m_className << "::Mixer:\n"
                    << "    Negative inelastic cross-section at " << emax
                    << " eV. Set to zero.\n";
          m_cfLog[iE][np] = 0.;
        }
        SetScatteringParameters(m_scatModel[np], pEqIn[iemax][j],
                                m_scatCutLog[iE][np], m_scatParLog[iE][np]);
      }
      if (nNull > 0) {
        for (int j = 0; j < nNull; ++j) {
          ++np;
          m_cfLog[iE][np] = qNull[iemax][j] * van * scln[j];
          if (m_useCsOutput) outfile << qNull[iemax][j] << "  ";
        }
      }
      if (m_useCsOutput) outfile << "\n";
      // Increase the energy.
      emax *= rLog;
    }
  }
  if (m_useCsOutput) outfile.close();

  // Find the smallest ionisation threshold.
  auto it = std::min_element(std::begin(m_ionPot),
                             std::begin(m_ionPot) + m_nComponents);
  m_minIonPot = *it;
  std::string minIonPotGas = m_gas[std::distance(std::begin(m_ionPot), it)];

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n"
              << "    Lowest ionisation threshold in the mixture: "
              << m_minIonPot << " eV (" << minIonPotGas << ")\n";
  }

  for (unsigned int iE = 0; iE < Magboltz::nEnergySteps; ++iE) {
    // Calculate the total collision frequency.
    for (unsigned int k = 0; k < m_nTerms; ++k) {
      if (m_cf[iE][k] < 0.) {
        std::cerr << m_className << "::Mixer:\n"
                  << "    Negative collision rate at " << (iE + 0.5) * m_eStep
                  << " eV, cross-section " << k << ". Set to zero.\n";
        std::cout << m_description[k] << "\n";
        m_cf[iE][k] = 0.;
      }
      m_cfTot[iE] += m_cf[iE][k];
    }
    // Normalise the collision probabilities.
    if (m_cfTot[iE] > 0.) {
      for (unsigned int k = 0; k < m_nTerms; ++k) m_cf[iE][k] /= m_cfTot[iE];
    }
    for (unsigned int k = 1; k < m_nTerms; ++k) {
      m_cf[iE][k] += m_cf[iE][k - 1];
    }
    const double ekin = m_eStep * (iE + 0.5);
    m_cfTot[iE] *= sqrt(ekin);
    // Use relativistic expression at high energies.
    if (ekin > 1.e3) {
      const double re = ekin / ElectronMass;
      m_cfTot[iE] *= sqrt(1. + 0.5 * re) / (1. + re);
    }
  }

  if (m_eMax > m_eHigh) {
    const double rLog = pow(m_eMax / m_eHigh, 1. / nEnergyStepsLog);
    for (int iE = 0; iE < nEnergyStepsLog; ++iE) {
      // Calculate the total collision frequency.
      for (unsigned int k = 0; k < m_nTerms; ++k) {
        if (m_cfLog[iE][k] < 0.) m_cfLog[iE][k] = 0.;
        m_cfTotLog[iE] += m_cfLog[iE][k];
      }
      // Normalise the collision probabilities.
      if (m_cfTotLog[iE] > 0.) {
        for (unsigned int k = 0; k < m_nTerms; ++k) {
          m_cfLog[iE][k] /= m_cfTotLog[iE];
        }
      }
      for (unsigned int k = 1; k < m_nTerms; ++k) {
        m_cfLog[iE][k] += m_cfLog[iE][k - 1];
      }
      const double ekin = m_eHigh * pow(rLog, iE + 1);
      const double re = ekin / ElectronMass;
      m_cfTotLog[iE] *= sqrt(ekin) * sqrt(1. + re) / (1. + re);
      // Store the logarithm (for log-log interpolation)
      m_cfTotLog[iE] = log(m_cfTotLog[iE]);
    }
  }

  // Determine the null collision frequency.
  m_cfNull = 0.;
  for (unsigned int j = 0; j < Magboltz::nEnergySteps; ++j) {
    if (m_cfTot[j] > m_cfNull) m_cfNull = m_cfTot[j];
  }
  if (m_eMax > m_eHigh) {
    for (int j = 0; j < nEnergyStepsLog; ++j) {
      const double r = exp(m_cfTotLog[j]);
      if (r > m_cfNull) m_cfNull = r;
    }
  }

  // Reset the collision counters.
  m_nCollisionsDetailed.assign(m_nTerms, 0);
  m_nCollisions.fill(0);

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n"
              << "    Energy [eV]    Collision Rate [ns-1]\n";
    const double emax = std::min(m_eHigh, m_eMax);
    for (int i = 0; i < 8; ++i) {
      const double en = (2 * i + 1) * emax / 16;
      const double cf = m_cfTot[(i + 1) * Magboltz::nEnergySteps / 16];
      std::printf("    %10.2f    %18.2f\n", en, cf);
    }
  }

  // Set up the de-excitation channels.
  if (m_useDeexcitation) {
    ComputeDeexcitationTable(verbose);
    for (const auto& dxc : m_deexcitations) {
      if (dxc.p.size() == dxc.final.size() && dxc.p.size() == dxc.type.size())
        continue;
      std::cerr << m_className << "::Mixer:\n"
                << "    Mismatch in deexcitation channel count. Program bug!\n"
                << "    Deexcitation handling is switched off.\n";
      m_useDeexcitation = false;
      break;
    }
  }

  // Fill the photon collision rates table.
  if (!ComputePhotonCollisionTable(verbose)) {
    // std::cerr << m_className << "::Mixer:\n"
    //          << "    Photon collision rates could not be calculated.\n";
    if (m_useDeexcitation) {
      std::cerr << "    Deexcitation handling is switched off.\n";
      m_useDeexcitation = false;
    }
  }

  // Reset the Penning transfer parameters.
  if (m_debug) {
    std::cout << m_className << "::Mixer: Resetting transfer probabilities.\n"
              << "    Global: " << m_rPenningGlobal << "\n";
    for (unsigned int i = 0; i < m_nMaxGases; ++i) {
      std::cout << "    Component " << i << ": " << m_rPenningGas[i] << "\n";
    }
  }
  if (m_rPenningGlobal > Small) {
    m_rPenning.fill(m_rPenningGlobal);
    m_lambdaPenning.fill(m_lambdaPenningGlobal);
  }
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    int iGas = int(m_csType[i] / nCsTypes);
    if (m_rPenningGas[iGas] > Small) {
      m_rPenning[i] = m_rPenningGas[iGas];
      m_lambdaPenning[i] = m_lambdaPenningGas[iGas];
    }
  }

  // Set the Green-Sawada splitting function parameters.
  SetupGreenSawada();

  return true;
}

void MediumMagboltz::PlotElectronCrossSections() {

  if (!Update()) return;

  std::array<float, Magboltz::nEnergySteps> en;
  for (unsigned int k = 0; k < Magboltz::nEnergySteps; ++k) {
    en[k] = (k + 0.5) * m_eStep;
  } 
  std::array<std::array<float, Magboltz::nEnergySteps>, 5> cs;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    for (size_t j = 0; j < 5; ++j) cs[j].fill(0.);
    const double density = GetNumberDensity() * m_fraction[i];
    const double scale = sqrt(0.5 * ElectronMass) / (density * SpeedOfLight);
    for (unsigned int j = 0; j < m_nTerms; ++j) {
      if (int(m_csType[j] / nCsTypes) != int(i)) continue;
      int cstype = m_csType[j] % nCsTypes;
      if (cstype >= ElectronCollisionTypeVirtual) continue;
      // Group inelastic and superelastic collisions.
      if (cstype == 5) cstype = 3;
      for (unsigned int k = 0; k < Magboltz::nEnergySteps; ++k) {
        double cf = m_cf[k][j];
        if (j > 0) cf -= m_cf[k][j - 1]; 
        cs[cstype][k] += cf * 1.e18 * scale;
      }
    }
    const std::string name = ViewBase::FindUnusedCanvasName("cCs");
    TCanvas* canvas = new TCanvas(name.c_str(), m_gas[i].c_str(), 800, 600);
    canvas->cd();
    canvas->SetLogx();
    canvas->SetLogy();
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->DrawFrame(en[0], 0.01, en.back(), 100., 
                      ";energy [eV];#sigma [Mbarn]");
    TGraph gr(Magboltz::nEnergySteps);
    gr.SetLineWidth(3);
    const std::array<short, 5> cols = {kBlack, kCyan - 2, kRed + 2, 
                                       kGreen + 3, kMagenta + 3}; 
    for (size_t j = 0; j < 5; ++j) {
      if (*std::max_element(cs[j].begin(), cs[j].end()) < 1.e-10) continue;
      gr.SetLineColor(cols[j]);
      gr.DrawGraph(Magboltz::nEnergySteps, en.data(), cs[j].data(), "lsame");
    }
    canvas->Update();
  }

}

void MediumMagboltz::SetupGreenSawada() {
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const double ta = 1000.;
    const double tb = m_parGreenSawada[i][4];
    m_hasGreenSawada[i] = true;
    if (m_gas[i] == "He" || m_gas[i] == "He-3") {
      m_parGreenSawada[i] = {15.5, 24.5, -2.25, ta, tb};
    } else if (m_gas[i] == "Ne") {
      m_parGreenSawada[i] = {24.3, 21.6, -6.49, ta, tb};
    } else if (m_gas[i] == "Ar") {
      m_parGreenSawada[i] = {6.92, 7.85, 6.87, ta, tb};
    } else if (m_gas[i] == "Kr") {
      m_parGreenSawada[i] = {7.95, 13.5, 3.90, ta, tb};
    } else if (m_gas[i] == "Xe") {
      m_parGreenSawada[i] = {7.93, 11.5, 3.81, ta, tb};
    } else if (m_gas[i] == "H2" || m_gas[i] == "D2") {
      m_parGreenSawada[i] = {7.07, 7.7, 1.87, ta, tb};
    } else if (m_gas[i] == "N2") {
      m_parGreenSawada[i] = {13.8, 15.6, 4.71, ta, tb};
    } else if (m_gas[i] == "O2") {
      m_parGreenSawada[i] = {18.5, 12.1, 1.86, ta, tb};
    } else if (m_gas[i] == "CH4") {
      m_parGreenSawada[i] = {7.06, 12.5, 3.45, ta, tb};
    } else if (m_gas[i] == "H2O") {
      m_parGreenSawada[i] = {12.8, 12.6, 1.28, ta, tb};
    } else if (m_gas[i] == "CO") {
      m_parGreenSawada[i] = {13.3, 14.0, 2.03, ta, tb};
    } else if (m_gas[i] == "C2H2") {
      m_parGreenSawada[i] = {9.28, 5.8, 1.37, ta, tb};
    } else if (m_gas[i] == "NO") {
      m_parGreenSawada[i] = {10.4, 9.5, -4.30, ta, tb};
    } else if (m_gas[i] == "CO2") {
      m_parGreenSawada[i] = {12.3, 13.8, -2.46, ta, tb};
    } else {
      m_parGreenSawada[i][3] = 0.;
      m_hasGreenSawada[i] = false;
      if (m_useGreenSawada) {
        std::cout << m_className << "::SetupGreenSawada:\n"
                  << "    Fit parameters for " << m_gas[i]
                  << " not available.\n"
                  << "    Opal-Beaty formula is used instead.\n";
      }
    }
  }
}

void MediumMagboltz::ComputeDeexcitationTable(const bool verbose) {
  m_iDeexcitation.fill(-1);
  m_deexcitations.clear();

  // Indices of "de-excitable" gases (only Ar for the time being).
  int iAr = -1;
  
  std::map<std::string, int> lvl;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    // Skip non-excitation levels.
    if (m_csType[i] % nCsTypes != ElectronCollisionTypeExcitation) continue;
    // Extract the index of the gas.
    const int ngas = int(m_csType[i] / nCsTypes);
    std::string level;
    if (m_gas[ngas] == "Ar") {
      // Argon
      if (iAr < 0) iAr = ngas;
      // Get the level description (as specified in Magboltz).
      level = "       ";
      for (int j = 0; j < 7; ++j) level[j] = m_description[i][5 + j];
      rtrim(level);     
      level = "Ar_" + level;
    } else {
      continue;
    }
    
    lvl[level] = m_deexcitations.size();
    m_iDeexcitation[i] = m_deexcitations.size();

    Deexcitation dxc;
    dxc.gas = ngas;
    dxc.level = i;
    dxc.label = level;
    // Excitation energy
    dxc.energy = m_energyLoss[i] * m_rgas[ngas];
    // Oscillator strength
    dxc.osc = dxc.cf = 0.;
    dxc.sDoppler = dxc.gPressure = dxc.width = 0.;

    if (level == "Ar_HIGH") {
      // This (artificial) level represents the sum of higher J = 1 states.
      // The deeexcitation cascade is simulated by allocating it
      // with equal probability to one of the five nearest levels below.
      dxc.type.assign(5, DxcTypeCollNonIon);
      dxc.p = {100., 100., 100., 100., 100.};
      dxc.final = {lvl["Ar_6D5"], lvl["Ar_5S1!"], lvl["Ar_4S2"], lvl["Ar_5S4"],
                   lvl["Ar_6D2"]};
    }
    m_deexcitations.push_back(std::move(dxc));
  }

  if (m_deexcitations.empty()) return;

  std::string path = ""; 
  auto installdir = std::getenv("GARFIELD_INSTALL");
  if (!installdir) {
    std::cerr << m_className << "::ComputeDeexcitationTable:\n"
              << "    Environment variable GARFIELD_INSTALL not set.\n";
  } else {
    path = std::string(installdir) + "/share/Garfield/Data/Deexcitation/";
  } 
  
  std::string filename = path + "OscillatorStrengths_Ar.txt";
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << m_className << "::ComputeDeexcitationTable:\n"
              << "    Could not open " << filename << ".\n";
    return;
  }
  for (std::string line; std::getline(infile, line);) {
    ltrim(line);
    if (line.empty() || IsComment(line)) continue;
    auto words = tokenize(line);
    if (words.size() < 2) continue; 
    std::string level = "Ar_" + words[0];
    if (lvl.count(level) == 0) {
      std::cout << "    Unexpected level " << level << "\n";
      continue;
    }
    m_deexcitations[lvl[level]].osc = std::stod(words[1]); 
  } 
  infile.close();

  filename = path + "TransitionRates_Ar.txt";
  infile.open(filename);
  if (!infile.is_open()) {
    std::cerr << m_className << "::ComputeDeexcitationTable:\n"
              << "    Could not open " << filename << ".\n";
    return;
  }
  for (std::string line; std::getline(infile, line);) {
    ltrim(line);
    if (line.empty() || IsComment(line)) continue;
    auto words = tokenize(line);
    if (words.size() < 3) continue; 
    std::string level0 = "Ar_" + words[0];
    if (lvl.count(level0) == 0) {
      std::cout << "    Unexpected level " << level0 << "\n";
      continue;
    }
    auto& dxc = m_deexcitations[lvl[level0]];
    if (words[1] == "Ground") {
      dxc.final.push_back(-1); 
    } else {
      std::string level1 = "Ar_" + words[1];
      if (lvl.count(level1) == 0) {
        std::cout << "    Unexpected level " << level1 << "\n";
        continue;
      }
      dxc.final.push_back(lvl[level1]);
    }
    dxc.p.push_back(std::stod(words[2]));
    dxc.type.push_back(DxcTypeRad);
  } 
  infile.close();
 
  if (m_debug || verbose) {
    std::cout << m_className << "::ComputeDeexcitationTable:\n";
    std::cout << "    Found " << m_deexcitations.size() << " levels "
              << "with available radiative de-excitation data.\n";
  }

  // Collisional de-excitation channels
  if (iAr >= 0) {
    // Add the Ar dimer ground state.
    Deexcitation dimer;
    dimer.label = "Ar_Dimer";
    dimer.level = -1;
    dimer.gas = iAr;
    dimer.energy = 14.71;
    dimer.osc = dimer.cf = 0.;
    dimer.sDoppler = dimer.gPressure = dimer.width = 0.;
    lvl["Ar_Dimer"] = m_deexcitations.size();
    m_deexcitations.push_back(std::move(dimer));
    // Add an Ar excimer level.
    Deexcitation excimer;
    excimer.label = "Ar_Excimer";
    excimer.level = -1;
    excimer.gas = iAr;
    excimer.energy = 14.71;
    excimer.osc = excimer.cf = 0.;
    excimer.sDoppler = excimer.gPressure = excimer.width = 0.;
    lvl["Ar_Excimer"] = m_deexcitations.size();
    m_deexcitations.push_back(std::move(excimer));
    const double nAr = GetNumberDensity() * m_fraction[iAr];

    filename = path + "RateConstants_Ar_Ar.txt";
    infile.open(filename);
    if (!infile.is_open()) {
      std::cerr << m_className << "::ComputeDeexcitationTable:\n"
                << "    Could not open " << filename << ".\n";
      return;
    }
    for (std::string line; std::getline(infile, line);) {
      ltrim(line);
      if (line.empty() || IsComment(line)) continue;
      auto words = tokenize(line);
      if (words.size() < 3) continue; 
      std::string level0 = "Ar_" + words[0];
      if (lvl.count(level0) == 0) {
        std::cout << "    Unexpected level " << level0 << "\n";
        continue;
      }
      auto& dxc = m_deexcitations[lvl[level0]];
      std::string level1 = "Ar_" + words[1];
      if (lvl.count(level1) == 0) {
        std::cout << "    Unexpected level " << level1 << "\n";
        continue;
      }
      const double k = std::stod(words[2]);
      // Three-body collisions lead to excimer formation.
      // Two-body collisions give rise to collisional mixing.
      if (level1 == "Ar_Excimer") {
        dxc.p.push_back(k * nAr * nAr);
      } else {
        dxc.p.push_back(k * nAr);
      }
      dxc.final.push_back(lvl[level1]);
      dxc.type.push_back(DxcTypeCollNonIon);
    } 
    infile.close();

    // Transfer from 3d and 5s levels to 4p levels.
    std::vector<std::string> levels3d5s = {
        "3D6", "3D5", "3D3", "3D4!", "3D4", "3D1!!", "3D1!", "3D2", 
        "3S1!!!!", "3S1!!", "3S1!!!", "3S1!", "2S5", "2S4", "2S3", "2S2"
    };
    std::vector<int> levels4p;
    for (unsigned int j = 1; j <= 10; ++j) {
      std::string level = "Ar_2P" + std::to_string(j);
      if (lvl.count(level) == 0) {
        std::cout << "    Unexpected level " << level << ".\n";
      } else {
        levels4p.push_back(lvl[level]);
      }
    }
    for (const std::string& level0 : levels3d5s) {
      if (lvl.count("Ar_" + level0) == 0) {
        std::cout << "    Unexpected level " << level0 << ".\n";
        continue;
      }
      auto& dxc = m_deexcitations[lvl["Ar_" + level0]];
      // Parameter to be tuned (order of magnitude guess).
      constexpr double k4p = 1.e-20;
      const double p4p = 0.1 * k4p * nAr; 
      for (const auto level1 : levels4p) {
        dxc.p.push_back(p4p);
        dxc.final.push_back(level1);
        dxc.type.push_back(DxcTypeCollNonIon);
      }
    }
    std::vector<std::string> levels = {
      "4D5", "3S4", "4D2", "4S1!", "3S2", "5D5", "4S4", "5D2", 
      "6D5", "5S1!", "4S2", "5S4", "6D2"};
    for (const std::string& level0 : levels) {
      if (lvl.count("Ar_" + level0) == 0) {
        std::cout << "    Unexpected level " << level0 << ".\n";
        continue;
      }
      auto& dxc = m_deexcitations[lvl["Ar_" + level0]];
      // Transfer to 4p levels.
      constexpr double k4p = 1.e-20;
      const double p4p = 0.1 * k4p * nAr; 
      for (const auto level1 : levels4p) {
        dxc.p.push_back(p4p);
        dxc.final.push_back(level1);
        dxc.type.push_back(DxcTypeCollNonIon);
      }
      // Hornbeck-Molnar ionisation
      // P. Becker and F. Lampe, J. Chem. Phys. 42 (1965), 3857-3863
      // A. Bogaerts and R. Gijbels, Phys. Rev. A 52 (1995), 3743-3751
      // This value seems high, to be checked!
      constexpr double kHM = 2.e-18;
      dxc.p.push_back(kHM * nAr);
      dxc.final.push_back(lvl["Ar_Dimer"]);
      dxc.type.push_back(DxcTypeCollIon);
    }
  }

  // Collisional deexcitation by quenching gases.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    std::string gas = m_gas[i];
    // Collision radius
    double rQ = 0.;
    if (m_gas[i] == "CO2") {
      rQ = 165.e-10;
    } else if (m_gas[i] == "CH4") {
      rQ = 190.e-10;
    } else if (m_gas[i] == "C2H6") {
      rQ = 195.e-10;
    } else if (m_gas[i] == "iC4H10") {
      rQ = 250.e-10;
      gas = "nC4H10";
    } else if (m_gas[i] == "C2H2") {
      rQ = 165.e-10;
    } else if (m_gas[i] == "CF4") {
      rQ = 235.e-10;
    } else {
      continue;
    }

    // Partial density.
    const double nQ = GetNumberDensity() * m_fraction[i];

    filename = path + "RateConstants_Ar_" + m_gas[i] + ".txt";
    infile.open(filename);
    if (!infile.is_open()) {
      std::cerr << m_className << "::ComputeDeexcitationTable:\n"
                << "    Could not open " << filename << ".\n";
      return;
    }
    for (std::string line; std::getline(infile, line);) {
      ltrim(line);
      if (line.empty() || IsComment(line)) continue;
      auto words = tokenize(line);
      if (words.size() < 2) continue; 
      std::string level0 = "Ar_" + words[0];
      if (lvl.count(level0) == 0) {
        std::cout << "    Unexpected level " << level0 << "\n";
        continue;
      }
      auto& dxc = m_deexcitations[lvl[level0]];
      if (dxc.energy < m_ionPot[i]) {
        AddPenningDeexcitation(dxc, std::stod(words[1]) * nQ, 0.);
      } else {
        const double eta = OpticalData::PhotoionisationYield(gas, dxc.energy);
        double pIon = pow(eta, 0.4);
        if (words.size() > 2 && !IsComment(words[2])) {
          pIon = std::stod(words[2]);
        }
        AddPenningDeexcitation(dxc, std::stod(words[1]) * nQ, pIon);
      }
    }

    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      if (level.find("Ar_1S") == 0 || level.find("Ar_2P") == 0) {
        continue;
      } 
      const double eta = OpticalData::PhotoionisationYield(gas, dxc.energy);
      const double pIon = pow(eta, 0.4);
      if (dxc.osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        double pacs = OpticalData::PhotoabsorptionCrossSection(gas, dxc.energy);
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, i);
        if (dxc.energy < m_ionPot[i]) {
          AddPenningDeexcitation(dxc, kQ * nQ, 0.);
        } else {
          AddPenningDeexcitation(dxc, kQ * nQ, pIon);
        }
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        constexpr double rAr3d = 436.e-10;
        const double kQ = RateConstantHardSphere(rAr3d, rQ, iAr, i);
        if (dxc.energy < m_ionPot[i]) {
          AddPenningDeexcitation(dxc, kQ * nQ, 0.);
        } else {
          AddPenningDeexcitation(dxc, kQ * nQ, pIon);
        }
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        constexpr double rAr5s = 635.e-10;
        const double kQ = RateConstantHardSphere(rAr5s, rQ, iAr, i);
        if (dxc.energy < m_ionPot[i]) {
          AddPenningDeexcitation(dxc, kQ * nQ, 0.);
        } else {
          AddPenningDeexcitation(dxc, kQ * nQ, pIon);
        }
      }
    }
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputeDeexcitationTable:\n"
              << "      Level  Energy [eV]                    Lifetimes [ns]\n"
              << "                            Total    Radiative       "
              << "     Collisional\n"
              << "                               "
              << "                Ionisation  Transfer      Loss\n";
  }

  for (auto& dxc : m_deexcitations) {
    // Calculate the total decay rate of each level.
    dxc.rate = 0.;
    double fRad = 0.;
    double fCollIon = 0., fCollTransfer = 0., fCollLoss = 0.;
    const unsigned int nChannels = dxc.type.size();
    for (unsigned int j = 0; j < nChannels; ++j) {
      dxc.rate += dxc.p[j];
      if (dxc.type[j] == DxcTypeRad) {
        fRad += dxc.p[j];
      } else if (dxc.type[j] == DxcTypeCollIon) {
        fCollIon += dxc.p[j];
      } else if (dxc.type[j] == DxcTypeCollNonIon) {
        if (dxc.final[j] < 0) {
          fCollLoss += dxc.p[j];
        } else {
          fCollTransfer += dxc.p[j];
        }
      } else {
        std::cerr << m_className << "::ComputeDeexcitationTable:\n    "
                  << "Unknown type of deexcitation channel (level " << dxc.label
                  << "). Program bug!\n";
      }
    }
    if (dxc.rate <= 0.) continue;
    // Print the radiative and collisional decay rates.
    if (m_debug || verbose) {
      std::cout << std::setw(12) << dxc.label << "  " << std::fixed
                << std::setprecision(3) << std::setw(7) << dxc.energy << "  "
                << std::setw(10) << 1. / dxc.rate << "  ";
      if (fRad > 0.) {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                  << 1. / fRad << " ";
      } else {
        std::cout << "---------- ";
      }
      if (fCollIon > 0.) {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                  << 1. / fCollIon << " ";
      } else {
        std::cout << "---------- ";
      }
      if (fCollTransfer > 0.) {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                  << 1. / fCollTransfer << " ";
      } else {
        std::cout << "---------- ";
      }
      if (fCollLoss > 0.) {
        std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                  << 1. / fCollLoss << "\n";
      } else {
        std::cout << "---------- \n";
      }
    }
    // Normalise the branching ratios.
    for (unsigned int j = 0; j < nChannels; ++j) {
      dxc.p[j] /= dxc.rate;
      if (j > 0) dxc.p[j] += dxc.p[j - 1];
    }
  }
}

double MediumMagboltz::RateConstantWK(const double energy, const double osc,
                                      const double pacs, const int igas1,
                                      const int igas2) const {
  // Calculate rate constant from Watanabe-Katsuura formula.
  const double m1 = ElectronMassGramme / (m_rgas[igas1] - 1.);
  const double m2 = ElectronMassGramme / (m_rgas[igas2] - 1.);
  // Compute the reduced mass.
  double mR = (m1 * m2 / (m1 + m2)) / AtomicMassUnit;
  const double uA = (RydbergEnergy / energy) * osc;
  const double uQ = (2 * RydbergEnergy / energy) * pacs /
                    (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
  return 2.591e-19 * pow(uA * uQ, 0.4) * pow(m_temperature / mR, 0.3);
}

double MediumMagboltz::RateConstantHardSphere(const double r1, const double r2,
                                              const int igas1,
                                              const int igas2) const {
  // Hard sphere cross-section
  const double r = r1 + r2;
  const double sigma = r * r * Pi;
  // Reduced mass
  const double m1 = ElectronMass / (m_rgas[igas1] - 1.);
  const double m2 = ElectronMass / (m_rgas[igas2] - 1.);
  const double mR = m1 * m2 / (m1 + m2);
  // Relative velocity
  const double vel =
      SpeedOfLight * sqrt(8. * BoltzmannConstant * m_temperature / (Pi * mR));
  return sigma * vel;
}

void MediumMagboltz::ComputeDeexcitation(int iLevel, int& fLevel) {
  if (!m_useDeexcitation) {
    std::cerr << m_className << "::ComputeDeexcitation: Not enabled.\n";
    return;
  }

  // Make sure that the tables are updated.
  if (!Update()) return;

  if (iLevel < 0 || iLevel >= (int)m_nTerms) {
    std::cerr << m_className << "::ComputeDeexcitation: Index out of range.\n";
    return;
  }

  iLevel = m_iDeexcitation[iLevel];
  if (iLevel < 0 || iLevel >= (int)m_deexcitations.size()) {
    std::cerr << m_className << "::ComputeDeexcitation:\n"
              << "    Level is not deexcitable.\n";
    return;
  }

  ComputeDeexcitationInternal(iLevel, fLevel);
  if (fLevel >= 0 && fLevel < (int)m_deexcitations.size()) {
    fLevel = m_deexcitations[fLevel].level;
  }
}

void MediumMagboltz::ComputeDeexcitationInternal(int iLevel, int& fLevel) {
  m_dxcProducts.clear();

  double t = 0.;
  fLevel = iLevel;
  while (iLevel >= 0 && iLevel < (int)m_deexcitations.size()) {
    const auto& dxc = m_deexcitations[iLevel];
    const int nChannels = dxc.p.size();
    if (dxc.rate <= 0. || nChannels <= 0) {
      // This level is a dead end.
      fLevel = iLevel;
      return;
    }
    // Determine the de-excitation time.
    t += -log(RndmUniformPos()) / dxc.rate;
    // Select the transition.
    fLevel = -1;
    int type = DxcTypeRad;
    const double r = RndmUniform();
    for (int j = 0; j < nChannels; ++j) {
      if (r <= dxc.p[j]) {
        fLevel = dxc.final[j];
        type = dxc.type[j];
        break;
      }
    }
    if (type == DxcTypeRad) {
      // Radiative decay
      dxcProd photon;
      photon.s = 0.;
      photon.t = t;
      photon.type = DxcProdTypePhoton;
      photon.energy = dxc.energy;
      if (fLevel >= 0) {
        // Decay to a lower lying excited state.
        photon.energy -= m_deexcitations[fLevel].energy;
        if (photon.energy < Small) photon.energy = Small;
        m_dxcProducts.push_back(std::move(photon));
        // Proceed with the next level in the cascade.
        iLevel = fLevel;
      } else {
        // Decay to ground state.
        double delta = RndmVoigt(0., dxc.sDoppler, dxc.gPressure);
        while (photon.energy + delta < Small || fabs(delta) >= dxc.width) {
          delta = RndmVoigt(0., dxc.sDoppler, dxc.gPressure);
        }
        photon.energy += delta;
        m_dxcProducts.push_back(std::move(photon));
        // Deexcitation cascade is over.
        fLevel = iLevel;
        return;
      }
    } else if (type == DxcTypeCollIon) {
      // Ionisation electron
      dxcProd electron;
      electron.s = 0.;
      electron.t = t;
      electron.type = DxcProdTypeElectron;
      electron.energy = dxc.energy;
      if (fLevel >= 0) {
        // Associative ionisation
        electron.energy -= m_deexcitations[fLevel].energy;
        if (electron.energy < Small) electron.energy = Small;
        ++m_nPenning;
        m_dxcProducts.push_back(std::move(electron));
        // Proceed with the next level in the cascade.
        iLevel = fLevel;
      } else {
        // Penning ionisation
        electron.energy -= m_minIonPot;
        if (electron.energy < Small) electron.energy = Small;
        ++m_nPenning;
        m_dxcProducts.push_back(std::move(electron));
        // Deexcitation cascade is over.
        fLevel = iLevel;
        return;
      }
    } else if (type == DxcTypeCollNonIon) {
      // Proceed with the next level in the cascade.
      iLevel = fLevel;
    } else {
      std::cerr << m_className << "::ComputeDeexcitationInternal:\n"
                << "    Unknown deexcitation type (" << type << "). Bug!\n";
      // Abort the calculation.
      fLevel = iLevel;
      return;
    }
  }
}

bool MediumMagboltz::ComputePhotonCollisionTable(const bool verbose) {

  // Atomic density
  const double dens = GetNumberDensity();

  // Reset the collision rate arrays.
  m_cfTotGamma.assign(nEnergyStepsGamma, 0.);
  m_cfGamma.assign(nEnergyStepsGamma, std::vector<double>());
  csTypeGamma.clear();

  m_nPhotonTerms = 0;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const double prefactor = dens * SpeedOfLight * m_fraction[i];
    // Check if optical data for this gas is available.
    std::string gasname = m_gas[i];
    if (gasname == "iC4H10") {
      gasname = "nC4H10";
      if (m_debug || verbose) {
        std::cout << m_className << "::ComputePhotonCollisionTable:\n"
                  << "    Photoabsorption cross-section for "
                  << "iC4H10 not available.\n"
                  << "    Using n-butane cross-section instead.\n";
      }
    }
    if (!OpticalData::IsAvailable(gasname)) return false;
    csTypeGamma.push_back(i * nCsTypesGamma + PhotonCollisionTypeIonisation);
    csTypeGamma.push_back(i * nCsTypesGamma + PhotonCollisionTypeInelastic);
    m_nPhotonTerms += 2;
    for (int j = 0; j < nEnergyStepsGamma; ++j) {
      const double en = (j + 0.5) * m_eStepGamma;
      // Retrieve total photoabsorption cross-section and ionisation yield.
      const double cs = OpticalData::PhotoabsorptionCrossSection(gasname, en);
      const double eta = OpticalData::PhotoionisationYield(gasname, en);
      m_cfTotGamma[j] += cs * prefactor;
      // Ionisation
      m_cfGamma[j].push_back(cs * prefactor * eta);
      // Inelastic absorption
      m_cfGamma[j].push_back(cs * prefactor * (1. - eta));
    }
  }

  // If requested, write the cross-sections to file.
  if (m_useCsOutput) {
    std::ofstream csfile;
    csfile.open("csgamma.txt", std::ios::out);
    for (int j = 0; j < nEnergyStepsGamma; ++j) {
      csfile << (j + 0.5) * m_eStepGamma << "  ";
      for (unsigned int i = 0; i < m_nPhotonTerms; ++i)
        csfile << m_cfGamma[j][i] << "  ";
      csfile << "\n";
    }
    csfile.close();
  }

  // Calculate the cumulative rates.
  for (int j = 0; j < nEnergyStepsGamma; ++j) {
    for (unsigned int i = 0; i < m_nPhotonTerms; ++i) {
      if (i > 0) m_cfGamma[j][i] += m_cfGamma[j][i - 1];
    }
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputePhotonCollisionTable:\n";
    std::cout << "    Energy [eV]      Mean free path [um]\n";
    for (int i = 0; i < 10; ++i) {
      const int j = (2 * i + 1) * nEnergyStepsGamma / 20;
      const double en = (2 * i + 1) * m_eFinalGamma / 20;
      const double imfp = m_cfTotGamma[j] / SpeedOfLight;
      if (imfp > 0.) {
        printf("    %10.2f    %18.4f\n", en, 1.e4 / imfp);
      } else {
        printf("    %10.2f          ------------\n", en);
      }
    }
  }

  if (!m_useDeexcitation) return true;

  // Conversion factor from oscillator strength to cross-section
  constexpr double f2cs =
      FineStructureConstant * 2 * Pi2 * HbarC * HbarC / ElectronMass;
  // Discrete absorption lines
  int nResonanceLines = 0;
  for (auto& dxc : m_deexcitations) {
    if (dxc.osc < Small) continue;
    const double prefactor = dens * SpeedOfLight * m_fraction[dxc.gas];
    dxc.cf = prefactor * f2cs * dxc.osc;
    // Compute the line width due to Doppler broadening.
    const double mgas = ElectronMass / (m_rgas[dxc.gas] - 1.);
    const double wDoppler = sqrt(BoltzmannConstant * m_temperature / mgas);
    dxc.sDoppler = wDoppler * dxc.energy;
    // Compute the half width at half maximum due to resonance broadening.
    //   A. W. Ali and H. R. Griem, Phys. Rev. 140, 1044
    //   A. W. Ali and H. R. Griem, Phys. Rev. 144, 366
    const double kResBroad = 1.92 * Pi * sqrt(1. / 3.);
    dxc.gPressure = kResBroad * FineStructureConstant * pow(HbarC, 3) *
                    dxc.osc * dens * m_fraction[dxc.gas] /
                    (ElectronMass * dxc.energy);
    // Make an estimate for the width within which a photon can be
    // absorbed by the line
    constexpr double nWidths = 1000.;
    // Calculate the FWHM of the Voigt distribution according to the
    // approximation formula given in
    // Olivero and Longbothum, J. Quant. Spectr. Rad. Trans. 17, 233-236
    const double fwhmGauss = dxc.sDoppler * sqrt(2. * log(2.));
    const double fwhmLorentz = dxc.gPressure;
    const double fwhmVoigt =
        0.5 * (1.0692 * fwhmLorentz + sqrt(0.86639 * fwhmLorentz * fwhmLorentz +
                                           4 * fwhmGauss * fwhmGauss));
    dxc.width = nWidths * fwhmVoigt;
    ++nResonanceLines;
  }

  if (nResonanceLines <= 0) {
    std::cerr << m_className << "::ComputePhotonCollisionTable:\n"
              << "    No resonance lines found.\n";
    return true;
  }

  if (!(m_debug || verbose)) return true;
  std::cout << m_className << "::ComputePhotonCollisionTable:\n    "
            << "Discrete absorption lines:\n   Energy [eV]   "
            << "Line width (FWHM) [eV]    Mean free path [um]\n        "
            << "              Doppler    Pressure         (peak)\n";
  for (const auto& dxc : m_deexcitations) {
    if (dxc.osc < Small) continue;
    const double wp = 2 * dxc.gPressure;
    const double wd = 2 * sqrt(2 * log(2.)) * dxc.sDoppler;
    const double imfpP =
        (dxc.cf / SpeedOfLight) * TMath::Voigt(0., dxc.sDoppler, wp);
    if (imfpP > 0.) {
      printf("  %6.3f +/- %6.1e  %6.2e  %6.3e  %14.4f\n", dxc.energy,
             dxc.width, wd, wp, 1.e4 / imfpP);
    } else {
      printf("  %6.3f +/- %6.1e  %6.2e  %6.3e  -------------\n", dxc.energy,
             dxc.width, wd, wp);
    }
  }
  return true;
}

void MediumMagboltz::RunMagboltz(
    const double emag, const double bmag, const double btheta, const int ncoll,
    bool verbose, double& vx, double& vy, double& vz, double& dl, double& dt,
    double& alpha, double& eta, double& lor, double& vxerr, double& vyerr,
    double& vzerr, double& dlerr, double& dterr, double& alphaerr,
    double& etaerr, double& lorerr, double& alphatof,
    std::array<double, 6>& difftens) {
  // Initialize the values.
  vx = vy = vz = 0.;
  dl = dt = 0.;
  alpha = eta = alphatof = 0.;
  lor = 0.;
  vxerr = vyerr = vzerr = 0.;
  dlerr = dterr = 0.;
  alphaerr = etaerr = 0.;
  lorerr = 0.;

  // Set the input parameters in the Magboltz common blocks.
  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nAniso = 2;
  if (m_autoEnergyLimit) {
    Magboltz::inpt_.efinal = 0.;
  } else {
    Magboltz::inpt_.efinal = std::min(m_eMax, m_eHigh);
  }
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;
  Magboltz::inpt_.ipen = 0;
  Magboltz::setp_.nmax = ncoll;

  Magboltz::thrm_.ithrm = m_useGasMotion ? 1 : 0;

  Magboltz::setp_.efield = emag;
  // Convert from Tesla to kGauss.
  Magboltz::bfld_.bmag = bmag * 10.;
  // Convert from radians to degree.
  Magboltz::bfld_.btheta = btheta * RadToDegree;

  // Set the gas composition in Magboltz.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const int ng = GetGasNumberMagboltz(m_gas[i]);
    if (ng <= 0) {
      std::cerr << m_className << "::RunMagboltz:\n    Gas " << m_gas[i]
                << " does not have a gas number in Magboltz.\n";
      return;
    }
    Magboltz::gasn_.ngasn[i] = ng;
    Magboltz::ratio_.frac[i] = 100 * m_fraction[i];
  }

  // Run Magboltz.
  Magboltz::magboltz_();

  // Velocities. Convert to cm / ns.
  vx = Magboltz::vel_.wx * 1.e-9;
  vxerr = Magboltz::velerr_.dwx;
  vy = Magboltz::vel_.wy * 1.e-9;
  vyerr = Magboltz::velerr_.dwy;
  vz = Magboltz::vel_.wz * 1.e-9;
  vzerr = Magboltz::velerr_.dwz;

  // Calculate the Lorentz angle.
  const double vt = sqrt(vx * vx + vy * vy);
  const double v2 = (vx * vx + vy * vy + vz * vz);
  lor = atan2(vt, vz);
  if (vt > 0. && v2 > 0. && fabs(lor) > 0.) {
    const double dvx = vx * vxerr;
    const double dvy = vy * vyerr;
    const double dvz = vz * vzerr;
    const double a = vz / vt;
    lorerr = sqrt(a * a * (vx * vx * dvx * dvx + vy * vy * dvy * dvy) +
                  vt * vt * dvz * dvz) /
             v2;
    lorerr /= lor;
  }

  // Diffusion coefficients.
  dt = sqrt(0.2 * 0.5 * (Magboltz::diflab_.difxx + Magboltz::diflab_.difyy) /
            vz) *
       1.e-4;
  dterr = 0.5 * sqrt(Magboltz::diferl_.dfter * Magboltz::diferl_.dfter +
                     vzerr * vzerr);
  dl = sqrt(0.2 * Magboltz::diflab_.difzz / vz) * 1.e-4;
  dlerr = 0.5 * sqrt(Magboltz::diferl_.dfler * Magboltz::diferl_.dfler +
                     vzerr * vzerr);
  // Diffusion tensor.
  difftens[0] = 0.2e-4 * Magboltz::diflab_.difzz / vz;
  difftens[1] = 0.2e-4 * Magboltz::diflab_.difxx / vz;
  difftens[2] = 0.2e-4 * Magboltz::diflab_.difyy / vz;
  difftens[3] = 0.2e-4 * Magboltz::diflab_.difxz / vz;
  difftens[4] = 0.2e-4 * Magboltz::diflab_.difyz / vz;
  difftens[5] = 0.2e-4 * Magboltz::diflab_.difxy / vz;
  // Townsend and attachment coefficients.
  alpha = Magboltz::ctowns_.alpha;
  alphaerr = Magboltz::ctwner_.alper;
  eta = Magboltz::ctowns_.att;
  etaerr = Magboltz::ctwner_.atter;

  // Calculate effective Townsend SST coefficient from TOF results. 
  if (fabs(Magboltz::tofout_.tofdl) > 0.) {
    const double wrzn = 1.e5 * Magboltz::tofout_.tofwr;
    const double fc1 = 0.5 * wrzn / Magboltz::tofout_.tofdl;
    const double fc2 = (Magboltz::tofout_.ralpha - 
                        Magboltz::tofout_.rattof) * 1.e12 / 
                       Magboltz::tofout_.tofdl;
    alphatof = fc1 - sqrt(fc1 * fc1 - fc2);
  }
  // Print the results.
  if (!(m_debug || verbose)) return;
  std::cout << m_className << "::RunMagboltz: Results:\n";
  printf("    Drift velocity along E:   %12.8f cm/ns +/- %5.2f%%\n", vz, vzerr);
  printf("    Drift velocity along Bt:  %12.8f cm/ns +/- %5.2f%%\n", vx, vxerr);
  printf("    Drift velocity along ExB: %12.8f cm/ns +/- %5.2f%%\n", vy, vyerr);
  printf("    Lorentz angle:            %12.3f degree\n", lor * RadToDegree);
  printf("    Longitudinal diffusion:   %12.8f cm1/2 +/- %5.2f%%\n", dl, dlerr);
  printf("    Transverse diffusion:     %12.8f cm1/2 +/- %5.2f%%\n", dt, dterr);
  printf("    Townsend coefficient:     %12.4f cm-1  +/- %5.2f%%\n", alpha,
         alphaerr);
  printf("    Attachment coefficient:   %12.4f cm-1  +/- %5.2f%%\n", eta,
         etaerr);
  if (alphatof > 0.) {
    printf("    TOF effective Townsend:   %12.4f cm-1 (alpha - eta)\n",
           alphatof);
  }
}

void MediumMagboltz::GenerateGasTable(const int numColl, const bool verbose) {
  // Set the reference pressure and temperature.
  m_pressureTable = m_pressure;
  m_temperatureTable = m_temperature;

  // Initialize the parameter arrays.
  const unsigned int nEfields = m_eFields.size();
  const unsigned int nBfields = m_bFields.size();
  const unsigned int nAngles = m_bAngles.size();
  Init(nEfields, nBfields, nAngles, m_eVelE, 0.);
  Init(nEfields, nBfields, nAngles, m_eVelB, 0.);
  Init(nEfields, nBfields, nAngles, m_eVelX, 0.);
  Init(nEfields, nBfields, nAngles, m_eDifL, 0.);
  Init(nEfields, nBfields, nAngles, m_eDifT, 0.);
  Init(nEfields, nBfields, nAngles, m_eLor, 0.);
  Init(nEfields, nBfields, nAngles, m_eAlp, -30.);
  Init(nEfields, nBfields, nAngles, m_eAlp0, -30.);
  Init(nEfields, nBfields, nAngles, m_eAtt, -30.);
  Init(nEfields, nBfields, nAngles, 6, m_eDifM, 0.);

  m_excRates.clear();
  m_ionRates.clear();
  // Retrieve the excitation and ionisation cross-sections in the gas mixture.
  GetExcitationIonisationLevels();
  std::cout << m_className << "::GenerateGasTable: Found "
            << m_excLevels.size() << " excitations and "
            << m_ionLevels.size() << " ionisations.\n";
  for (const auto& exc : m_excLevels) {
    std::cout << "    " << exc.label << ", energy = " << exc.energy << " eV.\n";
  }
  for (const auto& ion : m_ionLevels) {
    std::cout << "    " << ion.label << ", energy = " << ion.energy << " eV.\n";
  }
  if (!m_excLevels.empty()) {
    Init(nEfields, nBfields, nAngles, m_excLevels.size(), m_excRates, 0.);
  }
  if (!m_ionLevels.empty()) {
    Init(nEfields, nBfields, nAngles, m_ionLevels.size(), m_ionRates, 0.);
  }
  double vx = 0., vy = 0., vz = 0.;
  double difl = 0., dift = 0.;
  double alpha = 0., eta = 0.;
  double lor = 0.;
  double vxerr = 0., vyerr = 0., vzerr = 0.;
  double diflerr = 0., difterr = 0.;
  double alphaerr = 0., etaerr = 0.;
  double alphatof = 0.;
  double lorerr = 0.;
  std::array<double, 6> difftens;

  // Run through the grid of E- and B-fields and angles.
  for (unsigned int i = 0; i < nEfields; ++i) {
    const double e = m_eFields[i];
    for (unsigned int j = 0; j < nAngles; ++j) {
      const double a = m_bAngles[j];
      for (unsigned int k = 0; k < nBfields; ++k) {
        const double b = m_bFields[k];
        std::cout << m_className << "::GenerateGasTable: E = " << e
                  << " V/cm, B = " << b << " T, angle: " << a << " rad\n";
        RunMagboltz(e, b, a, numColl, verbose, vx, vy, vz, difl, dift, alpha,
                    eta, lor, vxerr, vyerr, vzerr, diflerr, difterr, alphaerr,
                    etaerr, lorerr, alphatof, difftens);
        m_eVelE[j][k][i] = vz;
        m_eVelX[j][k][i] = vy;
        m_eVelB[j][k][i] = vx;
        m_eDifL[j][k][i] = difl;
        m_eDifT[j][k][i] = dift;
        m_eLor[j][k][i] = lor;
        m_eAlp[j][k][i] = alpha > 0. ? log(alpha) : -30.;
        m_eAlp0[j][k][i] = alpha > 0. ? log(alpha) : -30.;
        m_eAtt[j][k][i] = eta > 0. ? log(eta) : -30.;
        for (unsigned int l = 0; l < 6; ++l) {
          m_eDifM[l][j][k][i] = difftens[l];
        }
        // Retrieve the excitation and ionisation rates.
        unsigned int nNonZero = 0;
        if (m_useGasMotion) {
          // Retrieve the total collision frequency and number of collisions.
          double ftot = 0., fel = 0., fion = 0., fatt = 0., fin = 0.;
          std::int64_t ntotal = 0;
          Magboltz::colft_(&ftot, &fel, &fion, &fatt, &fin, &ntotal);
          if (ntotal == 0) continue;
          // Convert from ps-1 to ns-1.
          const double scale = 1.e3 * ftot / ntotal;
          for (unsigned int ig = 0; ig < m_nComponents; ++ig) {
            const auto nL = Magboltz::larget_.last[ig];
            for (std::int64_t il = 0; il < nL; ++il) {
              if (Magboltz::larget_.iarry[il][ig] <= 0) break;
              // Skip levels that are not ionisations or inelastic collisions.
              const int cstype = (Magboltz::larget_.iarry[il][ig] - 1) % 5;
              if (cstype != 1 && cstype != 3) continue;
              // const int igas = int((Magboltz::larget_.iarry[il][ig] - 1) / 5);
              auto descr = GetDescription(il, ig, Magboltz::script_.dscrpt);
              descr = m_gas[ig] + descr;
              if (cstype == 3) {
                const unsigned int nExc = m_excLevels.size();
                for (unsigned int ie = 0; ie < nExc; ++ie) {
                  if (descr != m_excLevels[ie].label) continue;
                  const auto ncoll = Magboltz::outptt_.icoln[il][ig];
                  m_excRates[ie][j][k][i] = scale * ncoll;
                  if (ncoll > 0) ++nNonZero;
                  break;
                }
              } else if (cstype == 1) {
                const unsigned int nIon = m_ionLevels.size();
                for (unsigned int ii = 0; ii < nIon; ++ii) {
                  if (descr != m_ionLevels[ii].label) continue;
                  const auto ncoll = Magboltz::outptt_.icoln[il][ig];
                  m_ionRates[ii][j][k][i] = scale * ncoll;
                  if (ncoll > 0) ++nNonZero;
                  break;
                }
              }
            }
          }
        } else {
          // Retrieve the total collision frequency and number of collisions.
          double ftot = 0., fel = 0., fion = 0., fatt = 0., fin = 0.;
          std::int64_t ntotal = 0;
          Magboltz::colf_(&ftot, &fel, &fion, &fatt, &fin, &ntotal);
          if (ntotal == 0) continue;
          // Convert from ps-1 to ns-1.
          const double scale = 1.e3 * ftot / ntotal;
          for (std::int64_t il = 0; il < Magboltz::nMaxLevels; ++il) {
            if (Magboltz::large_.iarry[il] <= 0) break;
            // Skip levels that are not ionisations or inelastic collisions.
            const int cstype = (Magboltz::large_.iarry[il] - 1) % 5;
            if (cstype != 1 && cstype != 3) continue;
            const int igas = int((Magboltz::large_.iarry[il] - 1) / 5);
            std::string descr = GetDescription(il, Magboltz::scrip_.dscrpt);
            descr = m_gas[igas] + descr;
            if (cstype == 3) {
              const unsigned int nExc = m_excLevels.size();
              for (unsigned int ie = 0; ie < nExc; ++ie) {
                if (descr != m_excLevels[ie].label) continue;
                m_excRates[ie][j][k][i] = scale * Magboltz::outpt_.icoln[il];
                if (Magboltz::outpt_.icoln[il] > 0) ++nNonZero;
                break;
              }
            } else if (cstype == 1) {
              const unsigned int nIon = m_ionLevels.size();
              for (unsigned int ii = 0; ii < nIon; ++ii) {
                if (descr != m_ionLevels[ii].label) continue;
                m_ionRates[ii][j][k][i] = scale * Magboltz::outpt_.icoln[il];
                if (Magboltz::outpt_.icoln[il] > 0) ++nNonZero;
                break;
              }
            }
          }
        }
        if (nNonZero > 0) {
          std::cout << "    Excitation and ionisation rates:\n";
          std::cout << "                         Level                         "
                    << "          Rate [ns-1]\n";
          const unsigned int nExc = m_excLevels.size();
          for (unsigned int ie = 0; ie < nExc; ++ie) {
            if (m_excRates[ie][j][k][i] <= 0) continue;
            std::cout << std::setw(60) << m_excLevels[ie].label;
            std::printf(" %15.8f\n", m_excRates[ie][j][k][i]);
          }
          const unsigned int nIon = m_ionLevels.size();
          for (unsigned int ii = 0; ii < nIon; ++ii) {
            if (m_ionRates[ii][j][k][i] <= 0) continue;
            std::cout << std::setw(60) << m_ionLevels[ii].label;
            std::printf(" %15.8f\n", m_ionRates[ii][j][k][i]);
          }
        }
      }
    }
  }
  // Set the threshold indices.
  SetThreshold(m_eAlp);
  SetThreshold(m_eAtt);
}

void MediumMagboltz::GetExcitationIonisationLevels() {
 
  // Reset.
  m_excLevels.clear();
  m_ionLevels.clear();
  // Cross-sections.
  static double q[Magboltz::nEnergySteps][6];
  static double qIn[Magboltz::nEnergySteps][Magboltz::nMaxInelasticTerms];
  static double qIon[Magboltz::nEnergySteps][Magboltz::nMaxIonisationTerms];
  static double qAtt[Magboltz::nEnergySteps][Magboltz::nMaxAttachmentTerms];
  static double qNull[Magboltz::nEnergySteps][Magboltz::nMaxNullTerms];
  // Parameters for angular distributions.
  static double pEqEl[Magboltz::nEnergySteps][6];
  static double pEqIn[Magboltz::nEnergySteps][Magboltz::nMaxInelasticTerms];
  static double pEqIon[Magboltz::nEnergySteps][Magboltz::nMaxIonisationTerms];
  // Penning transfer parameters
  static double penFra[Magboltz::nMaxInelasticTerms][3];
  // Description of cross-section terms
  static char scrpt[Magboltz::nMaxLevelsPerComponent][Magboltz::nCharDescr];
  static char scrptn[Magboltz::nMaxNullTerms][Magboltz::nCharDescr];

  // Loop over the gases in the mixture.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    // Choose the energy range large enough to cover all relevant levels.
    const double emax = 400.;
    Magboltz::inpt_.efinal = emax;
    Magboltz::inpt_.estep = emax / Magboltz::nEnergySteps;
    char name[Magboltz::nCharName];
    // Number of inelastic, ionisation, attachment and null-collision levels.
    std::int64_t nIn = 0, nIon = 0, nAtt = 1, nNull = 0;
    // Virial coefficient (not used)
    double virial = 0.;
    // Thresholds/characteristic energies.
    static double e[6];
    // Energy losses and ionisation thresholds.
    static double eIn[Magboltz::nMaxInelasticTerms];
    static double eIon[Magboltz::nMaxIonisationTerms];
    // Scattering parameters.
    static std::int64_t kIn[Magboltz::nMaxInelasticTerms];
    static std::int64_t kEl[6];
    // Opal-Beaty parameters.
    static double eoby[Magboltz::nMaxIonisationTerms];
    // Scaling factor for "null-collision" terms
    static double scln[Magboltz::nMaxNullTerms];
    // Parameters for simulation of Auger and fluorescence processes.
    static std::int64_t nc0[Magboltz::nMaxIonisationTerms];
    static std::int64_t ng1[Magboltz::nMaxIonisationTerms];
    static std::int64_t ng2[Magboltz::nMaxIonisationTerms];
    static double ec0[Magboltz::nMaxIonisationTerms];
    static double wklm[Magboltz::nMaxIonisationTerms];
    static double efl[Magboltz::nMaxIonisationTerms];
    static double eg1[Magboltz::nMaxIonisationTerms];
    static double eg2[Magboltz::nMaxIonisationTerms];

    // Retrieve the cross-section data for this gas from Magboltz.
    std::int64_t ng = GetGasNumberMagboltz(m_gas[i]);
    if (ng <= 0) {
      std::cerr << m_className << "::GetExcitationIonisationLevels:\n\n"
                << "    Gas " << m_gas[i] << " not available in Magboltz.\n";
      continue;
    }
    Magboltz::gasmix_(
        &ng, q[0], qIn[0], &nIn, e, eIn, name, &virial, eoby, pEqEl[0], 
        pEqIn[0], penFra[0], kEl, kIn, qIon[0], pEqIon[0], eIon, &nIon, 
        qAtt[0], &nAtt, qNull[0], &nNull, scln, nc0, ec0, wklm, efl,
        ng1, eg1, ng2, eg2, scrpt, scrptn,
        Magboltz::nCharName, Magboltz::nCharDescr, Magboltz::nCharDescr);
    const double r = 1. + 0.5 * e[1];
    // Ionisation cross section(s).
    for (int j = 0; j < nIon; ++j) {
      const std::string descr = GetDescription(2 + j, scrpt);
      IonLevel ion;
      ion.label = m_gas[i] + descr;
      ion.energy = eIon[j] / r;
      m_ionLevels.push_back(std::move(ion));
    }
    // Excitation cross-sections.
    for (int j = 0; j < nIn; ++j) {
      const std::string descr = GetDescription(4 + nIon + nAtt + j, scrpt);
      if ((descr[1] == 'E' && descr[2] == 'X') ||
          (descr[0] == 'E' && descr[1] == 'X')) {
        // Excitation
        ExcLevel exc;
        exc.label = m_gas[i] + descr;
        exc.energy = eIn[j] / r;
        exc.prob = 0.;
        exc.rms = 0.;
        exc.dt = 0.;
        m_excLevels.push_back(std::move(exc));
      }
    }
  }
}

}
