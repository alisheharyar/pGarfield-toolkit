#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Garfield/Utilities.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackBichsel.hh"

namespace {
bool IsComment(const std::string& line) {
  if (line.empty()) return false;
  if (line[0] == '#') return true;
  if (line.size() > 1 && (line[0] == '/' && line[1] == '/')) return true;
  return false;
}

}

namespace Garfield {

TrackBichsel::TrackBichsel()
    : Track() {
  m_className = "TrackBichsel";
  Initialise();
}

bool TrackBichsel::Initialise() {

  std::cout << m_className << "::Initialize:\n";
  // Reset the tables.
  m_E.fill(0.);
  m_dfdE.fill(0.);
  m_eps1.fill(0.);
  m_eps2.fill(0.);
  m_int.fill(0.);
  m_k1.fill(0.);

  MediumSilicon si;
  m_density = si.GetNumberDensity();
  // Conversion from loss function to oscillator strength density.
  m_conv = ElectronMass / (2 * Pi2 * FineStructureConstant * pow(HbarC, 3) * m_density);
  // Number of bins for each factor of 2 in energy
  constexpr unsigned int n2 = 64;
  const double u = log(2.) / n2;
  const double um  = exp(u);
  constexpr double kEdge = 1839.;
  const double ken = log(kEdge / 1.5) / u;
  m_E[0] = kEdge / pow(2, ken / n2);
  if (m_debug) std::cout << "  Bin       Energy [eV]\n";
  for (size_t j = 0; j < NEnergyBins; ++j) {
    m_E[j + 1] = m_E[j] * um;
    if (!m_debug) continue;
    if ((j + 1) % 50 == 0) std::printf(" %4zu %16.8f\n", j + 1, m_E[j]);
  }

  // Read in the complex dielectric function (epsilon). This is used 
  // for the cross section of small momentum transfer excitations.
  std::string path = std::getenv("GARFIELD_INSTALL");
  if (!path.empty()) path += "/share/Garfield/Data/";
  std::ifstream infile;
  std::cout << "    Reading dielectric function.\n";
  infile.open(path + "heps.tab", std::ios::in);
  if (!infile) {
    std::cerr << "    Could not open heps.tab.\n";
    return false;
  }
  bool ok = true;
  for (std::string line; std::getline(infile, line);) {
    ltrim(line);
    // Skip comments and empty lines.
    if (IsComment(line)) continue;
    std::istringstream data(line);
    size_t j;
    double energy, eps1, eps2, lossFunction;
    data >> j >> energy >> eps1 >> eps2 >> lossFunction;
    if (j < 1 || j > NEnergyBins) {
      std::cerr << "    Index out of range.\n";
      ok = false;
      break;
    }
    --j;
    m_eps1[j] = eps1;
    m_eps2[j] = eps2;
    m_dfdE[j] = lossFunction * m_conv * m_E[j];
  }
  infile.close();
  if (!ok) {
    std::cerr << "    Error reading dielectric function.\n";
    return false;
  }

  // Read in a table of the generalized oscillator strength density 
  // integrated over the momentum transfer K, corresponding to the 
  // function A(E) in Equation (2.11) in (Bichsel, 1988).
  // These values are used for calculating the cross section for 
  // longitudinal excitations of K and L shell electrons 
  // with large momentum transfer.
  std::cout << "    Reading K and L shell calculations.\n";
  infile.open(path + "macom.tab", std::ios::in);
  if (!infile) {
    std::cerr << "    Could not open macom.tab.\n";
    return false;
  }
  for (std::string line; std::getline(infile, line);) {
    ltrim(line);
    if (IsComment(line)) continue;
    std::istringstream data(line);
    size_t j;
    double energy, val;
    data >> j >> energy >> val;
    if (j < 1 || j > NEnergyBins) {
      std::cerr << "    Index out of range.\n";
      ok = false;
      break;
    }
    --j;
    m_int[j] = val;
  }
  infile.close();
  if (!ok) {
    std::cerr << "    Error reading K/L shell data.\n";
    return false;
  }

  // Read in a table of the generalized oscillator strength density 
  // for M shell electrons integrated over the momentum transfer K.
  // Based on equations in the appendix in Emerson et al., 
  // Phys. Rev. B 7 (1973), 1798 (DOI 10.1103/PhysRevB.7.1798)
  std::cout << "    Reading M shell calculations.\n";
  infile.open(path + "emerc.tab", std::ios::in);
  if (!infile) {
    std::cerr << "    Could not open emerc.tab.\n";
    return false;
  }
  for (std::string line; std::getline(infile, line);) {
    ltrim(line);
    if (IsComment(line)) continue;
    std::istringstream data(line);
    size_t j;
    double energy, a, k1;
    data >> j >> energy >> a >> k1;
    if (j < 1 || j > NEnergyBins) {
      std::cerr << "    Index out of range.\n";
      ok = false;
      break;
    }
    --j;
    m_int[j] = a;
    m_k1[j] = k1;
    if (!m_debug) continue;
    if (j < 19 || j > 30) continue;
    std::printf(" %4zu %11.2f %11.2f %12.6f %12.6f\n", 
                j + 1, m_E[j], energy, m_int[j], m_k1[j]); 
  }
  infile.close();
  if (!ok) {
    std::cerr << "    Error reading M shell data.\n";
    return false;
  }

  double s0 = 0.;
  double s1 = 0.;
  double avI = 0.;
  for (size_t j = 0; j < NEnergyBins; ++j) {
    const auto logE = log(m_E[j]);
    const double dE = m_E[j + 1] - m_E[j];
    s0 += m_dfdE[j] * dE;
    s1 += m_dfdE[j] * dE * m_E[j];
    avI += m_dfdE[j] * dE * logE;
  }
  if (m_debug) {
    std::printf("    S0  = %15.5f\n", s0);
    std::printf("    S1  = %15.5f\n", s1);
    std::printf("    lnI = %15.5f\n", avI);
  }
  m_initialised = true;
  return true;
}

bool TrackBichsel::ComputeCrossSection() {

  if (!m_initialised) {
    std::cerr << m_className << "::ComputeCrossSection: Not initialised.\n";
    return false;
  }
  m_ready = false;
  const double bg = GetBetaGamma();
  if (m_debug) {
    std::cerr << m_className << "::ComputeCrossSection:\n"
              << "    Calculating differential cross-section for bg = "
              << bg << ".\n";
  }
  const double gamma = sqrt(bg * bg + 1.);
  // Prefactors in Bhabha formula.
  const double g1 = (gamma - 1.) * (gamma - 1.) / (gamma * gamma);
  const double g2 = (2. * gamma  - 1.) / (gamma * gamma);

  const double ek = GetKineticEnergy();
  const double rm = ElectronMass / m_mass;
  // Maximum energy transfer.
  double emax = 2 * ElectronMass * bg * bg;
  if (m_isElectron) {
    emax = 0.5 * ek;
  } else {
    emax = 2 * ElectronMass * bg * bg / (1. + 2 * gamma * rm + rm * rm);
  }
  if (m_debug) std::printf("    Max. energy transfer: %12.4f eV\n", emax);
  const double betaSq = bg * bg / (1. + bg * bg);
  m_speed = SpeedOfLight * sqrt(betaSq);
  constexpr size_t nTerms = 3;
  std::array<double, nTerms + 1> m0;
  std::array<double, nTerms + 1> m1;
  std::array<double, nTerms + 1> m2; 
  m0.fill(0.);
  m1.fill(0.);
  m2.fill(0.);

  std::array<std::array<double, NEnergyBins>, nTerms + 1> cs;
  std::vector<double> cdf(NEnergyBins, 0.);

  const auto betaSqOverEmax = betaSq / emax;
  const auto twoMeBetaSq = 2 * ElectronMass * betaSq;

  size_t jmax = 0; 
  for (size_t j = 0; j < NEnergyBins; ++j) {
    ++jmax;
    if (m_E[j] > emax) break;
    const auto e2 = m_E[j] * m_E[j];
    double q1 = RydbergEnergy;
    // CCS-33, 39 & 47
    if (m_E[j] < 11.9) {
      q1 = m_k1[j] * m_k1[j] * RydbergEnergy;
    } else if (m_E[j] < 100.) {
      constexpr double k1 = 0.025;
      q1 = k1 * k1 * RydbergEnergy;
    }
    const auto qmin = e2 / twoMeBetaSq;
    if (m_E[j] < 11.9 && q1 <= qmin) {
      cs[0][j] = 0.;
    } else {
      cs[0][j] = m_E[j] * m_dfdE[j] * log(q1 / qmin);
    }
    // Fano Eq. (47).
    auto b1 = 1. - betaSq * m_eps1[j];
    if (b1 == 0.) b1 = 1.e-20;
    const auto b2 = betaSq * m_eps2[j];
    const auto g = m_E[j] * m_dfdE[j] * (-0.5) * log(b1 * b1 + b2 * b2);
    auto theta = atan(b2 / b1);
    if (theta < 0.) theta += Pi;
    const auto epsSq = m_eps1[j] * m_eps1[j] + m_eps2[j] * m_eps2[j];
    const auto h = m_conv * e2 * (betaSq - m_eps1[j] / epsSq) * theta;
    cs[1][j] = g + h;
    // The integral was over d(lnK) rather than d(lnQ)
    cs[2][j] = 2 * m_int[j];
    if (m_isElectron) {
      // Uehling Eq. (9).
      const double u = m_E[j] / (ek - m_E[j]);
      const double v = m_E[j] / ek;
      cs[2][j] *= (1. + u * u + g1 * v * v - g2 * u);
    } else {
      // Uehling Eq. (2).
      cs[2][j] *= (1. - m_E[j] * betaSqOverEmax);
    } 
    cs[3][j] = 0.;
    cs[nTerms][j] = 0.;
    const double dE = m_E[j + 1] - m_E[j];
    for (size_t k = 0; k < nTerms; ++k) {
      m0[k] += cs[k][j] * dE / e2;
      m1[k] += cs[k][j] * dE / m_E[j];
      m2[k] += cs[k][j] * dE;
      cs[nTerms][j] += cs[k][j];
    } 
    m0[nTerms] += cs[nTerms][j] * dE / e2;
    m1[nTerms] += cs[nTerms][j] * dE / m_E[j]; 
    m2[nTerms] += cs[nTerms][j] * dE;
    cs[nTerms][j] /= e2;
    if (!m_debug) continue;
    if ((j + 1) % 10 != 0) continue;
    std::printf(" %4zu %9.1f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", 
                j + 1, m_E[j], m_dfdE[j], g, h, cs[0][j], cs[1][j], cs[2][j]);
  }
  if (jmax < NEnergyBins) cdf.resize(jmax);
  if (m_debug) {
    printf("  M0: %15.4f %15.4f %15.4f %15.4f\n", m0[0], m0[1], m0[2], m0[3]);
    printf("  M1: %15.4f %15.4f %15.4f %15.4f\n", m1[0], m1[1], m1[2], m1[3]);
    printf("  M2: %15.4f %15.4f %15.4f %15.4f\n", m2[0], m2[1], m2[2], m2[3]);
  }
  // Calculate the cumulative differential cross-section.
  cdf[0] = 0.5 * m_E[0] * cs[nTerms][0];
  for (size_t j = 1; j < jmax; ++j) {
    const double dE = m_E[j] - m_E[j - 1];
    cdf[j] = cdf[j - 1] + 0.5 * (cs[nTerms][j - 1] + cs[nTerms][j]) * dE;
  }

  constexpr double ary = BohrRadius * RydbergEnergy;
  constexpr double prefactor = 8 * Pi * ary * ary / ElectronMass;
  const double dec = m_q * m_q * m_density * prefactor / betaSq;
  m_imfp = m0.back() * dec;
  m_dEdx = m1.back() * dec;
  if (m_debug) {
    std::printf("    M0 = %12.4f cm-1      ... inverse mean free path\n", 
                m_imfp);
    std::printf("    M1 = %12.4f keV/cm    ... dE/dx\n", m_dEdx * 1.e-3);
    std::printf("    M2 = %12.4f keV2/cm\n", m2.back() * dec * 1.e-6);
  }
  // Calculate the residual cross-section.
  double integral = cdf.back();
  if (emax > m_E.back()) {
    const double e1 = m_E.back();
    const double rm0 = (1. - 720. * betaSqOverEmax) * (
      (1. / e1 - 1. / emax) + 
      2 * (1. / (e1 * e1) - 1. / (emax * emax))) -
      betaSqOverEmax * log(emax / e1);
    if (m_debug) std::printf("    Residual M0 = %15.5f\n", rm0 * 14. * dec);
    m_imfp += rm0 * 14. * dec;
    m0.back() += rm0;
    integral += 14. * rm0;
  }
  const double scale = 1. / integral;
  for (size_t j = 0; j < jmax; ++j) {
    cdf[j] *= scale;
    if (!m_debug) continue;
    if ((j + 1) % 20 != 0) continue;
    std::printf(" %4zu %9.1f %15.6f\n", j + 1, m_E[j], cdf[j]); 
  } 
  m_tab.fill(0.);
  for (size_t i = 0; i < NCdfBins; ++i) {
    constexpr double step = 1. / NCdfBins;
    const double x = (i + 1) * step;
    // Interpolate.
    m_tab[i] = 0.;
    const auto it1 = std::upper_bound(cdf.cbegin(), cdf.cend(), x);
    if (it1 == cdf.cbegin()) {
      m_tab[i] = m_E.front();
      continue;
    }
    const auto it0 = std::prev(it1);
    const double x0 = *it0;
    const double x1 = *it1;
    const double y0 = m_E[it0 - cdf.cbegin()];
    const double y1 = m_E[it1 - cdf.cbegin()];
    const double f0 = (x - x0) / (x1 - x0);
    const double f1 = 1. - f0;
    m_tab[i] = f0 * y0 + f1 * y1;
  } 
  m_ready = true;
  return true;
}

bool TrackBichsel::NewTrack(const double x0, const double y0, const double z0,
                            const double t0, const double dx0, const double dy0,
                            const double dz0) {
  // Make sure a sensor has been defined.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    m_isInMedium = false;
    return false;
  }

  // If not yet done, compute the cross-section table.
  if (!m_ready || m_isChanged) {
    if (!ComputeCrossSection()) {
      std::cerr << m_className << "::NewTrack:\n"
                << "    Could not calculate cross-section table.\n";
      return false;
    }
    m_isChanged = false;
  }

  // Make sure we are inside a medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::NewTrack: No medium at initial position.\n";
    m_isInMedium = false;
    return false;
  }

  // Check if the medium is silicon.
  if (medium->GetName() != "Si") {
    std::cerr << m_className << "::NewTrack: Medium is not silicon.\n";
    m_isInMedium = false;
    return false;
  }

  // Check if primary ionisation has been enabled.
  if (!medium->IsIonisable()) {
    std::cerr << m_className << "::NewTrack: Medium is not ionisable.\n";
    m_isInMedium = false;
    return false;
  }

  m_isInMedium = true;
  m_x = x0;
  m_y = y0;
  m_z = z0;
  m_t = t0;

  // Normalise the direction vector.
  const double d = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  if (d < Small) {
    // In case of a null vector, choose a random direction.
    RndmDirection(m_dx, m_dy, m_dz);
  } else {
    m_dx = dx0 / d;
    m_dy = dy0 / d;
    m_dz = dz0 / d;
  }
  return true;
}

bool TrackBichsel::GetCluster(double& xcls, double& ycls, double& zcls,
                              double& tcls, int& n, double& e, double& extra) {
  if (!m_ready || !m_isInMedium) return false;

  const double d = -log(RndmUniformPos()) / m_imfp;
  m_x += m_dx * d;
  m_y += m_dy * d;
  m_z += m_dz * d;
  m_t += d / m_speed;

  xcls = m_x;
  ycls = m_y;
  zcls = m_z;
  tcls = m_t;
  n = 0;
  e = 0.;
  extra = 0.;

  Medium* medium = m_sensor->GetMedium(m_x, m_y, m_z);
  if (!medium) {
    m_isInMedium = false;
    if (m_debug) {
      std::cout << m_className << "::GetCluster: Particle left the medium.\n";
    }
    return false;
  }

  if (medium->GetName() != "Si" || !medium->IsIonisable()) {
    m_isInMedium = false;
    if (m_debug) {
      std::cout << m_className << "::GetCluster: Particle left the medium.\n";
    }
    return false;
  }

  const double u = NCdfBins * RndmUniform();
  const size_t j = static_cast<size_t>(std::floor(u));
  if (j == 0) {
    e = 0. + u * m_tab.front();
  } else if (j >= NCdfBins) {
    e = m_tab.back();
  } else {
    e = m_tab[j - 1] + (u - j) * (m_tab[j] - m_tab[j - 1]);
  }

  return true;
}

double TrackBichsel::GetClusterDensity() {

  if (m_isChanged) {
    if (!ComputeCrossSection()) return 0.;
    m_isChanged = false;
  }
  return m_imfp;
}

double TrackBichsel::GetStoppingPower() {

  if (m_isChanged) {
    if (!ComputeCrossSection()) return 0.;
    m_isChanged = false;
  }
  return m_dEdx;
}

}
