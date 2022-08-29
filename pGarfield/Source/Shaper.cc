#include <iostream>
#include <algorithm>

#include <Math/SpecFuncMathCore.h>

#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Shaper.hh"

namespace {

double Heaviside(const double t, const double t0) {
  if (t < t0)
    return 0;
  else if (fabs(t - t0) < Garfield::Small)
    return 0.5;
  else
    return 1;
}

}

namespace Garfield {

Shaper::Shaper(const unsigned int n, const double tau, const double g,
               std::string shaperType) :
    m_n(n),
    m_tau(tau),
    m_g(g) {

  std::transform(shaperType.begin(), shaperType.end(), 
                 shaperType.begin(), toupper);
  if (shaperType == "UNIPOLAR") {
    m_type = ShaperType::Unipolar;
    m_tp = m_n * m_tau;
    m_prefactor = exp(m_n);
    m_transfer_func_sq = (exp(2 * m_n) / pow(2 * m_n, 2 * m_n)) * m_tp * 
                         ROOT::Math::tgamma(2 * m_n);
  } else if (shaperType == "BIPOLAR") {
    m_type = ShaperType::Bipolar;
    const double r = m_n - sqrt(m_n);
    m_tp = r * m_tau;
    m_prefactor = exp(r) / sqrt(m_n);
    m_transfer_func_sq = (exp(2 * r) / pow(2 * r, 2 * m_n)) * r * m_tp * 
                         ROOT::Math::tgamma(2 * m_n - 1);
  } else {
    std::cerr << m_className << ": Unknown shaper type.\n";
  } 
}

double Shaper::Shape(const double t) const {
  switch (m_type) {
    case ShaperType::Unipolar:
      return UnipolarShaper(t);
    case ShaperType::Bipolar:
      return BipolarShaper(t);
    default:
      break;
  }
  return 0;
}

double Shaper::UnipolarShaper(const double t) const {
  double f = m_prefactor * pow(t / m_tp, m_n) * exp(-t / m_tau) * Heaviside(t, 0.);
  return m_g * f;
}

double Shaper::BipolarShaper(const double t) const {
  double f = m_prefactor * (m_n - t / m_tau) * pow(t / m_tp, m_n - 1) * exp(-t / m_tau) * Heaviside(t, 0.);
  return m_g * f;
}

}
