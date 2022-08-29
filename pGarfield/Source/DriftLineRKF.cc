#include <iostream>
#include <cstdio>
#include <cmath>
#include <numeric>

#include "Garfield/DriftLineRKF.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/Random.hh"

namespace {

std::string PrintVec(const std::array<double, 3>& x) {

  return "(" + std::to_string(x[0]) + ", " + std::to_string(x[1]) + ", " +
         std::to_string(x[2]) + ")";
}

double Mag(const std::array<double, 3>& x) {
 
  return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

double Mag(const double x, const double y, const double z) {

  return sqrt(x * x + y * y + z * z);
}

double Mag(const double x, const double y) {

  return sqrt(x * x + y * y);
}

double Mag2(const double x, const double y) {

  return x * x + y * y;
}

double Dist(const std::array<double, 3>& x0,
            const std::array<double, 3>& x1) {

  return Mag(x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]);
}

std::array<double, 3> MidPoint(const std::array<double, 3>& x0,
                               const std::array<double, 3>& x1) {
  std::array<double, 3> xm;
  for (size_t k = 0; k < 3; ++k) xm[k] = 0.5 * (x0[k] + x1[k]);
  return xm;
}
}

namespace Garfield {

using Vec = std::array<double, 3>;

DriftLineRKF::DriftLineRKF() {
  m_t.reserve(1000);
  m_x.reserve(1000);
}

void DriftLineRKF::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }
  m_sensor = s;
}

void DriftLineRKF::SetIntegrationAccuracy(const double eps) {
  if (eps > 0.) {
    m_accuracy = eps;
  } else {
    std::cerr << m_className << "::SetIntegrationAccuracy:\n"
              << "    Accuracy must be greater than zero.\n";
  }
}

void DriftLineRKF::SetMaximumStepSize(const double ms) {
  if (ms > 0.) {
    m_maxStepSize = ms;
    m_useStepSizeLimit = true;
  } else {
    std::cerr << m_className << "::SetMaximumStepSize:\n"
              << "    Step size must be greater than zero.\n";
  }
}

void DriftLineRKF::EnablePlotting(ViewDrift* view) {
  if (!view) {
    std::cerr << m_className << "::EnablePlotting: Null pointer.\n";
    return;
  }
  m_view = view;
}

void DriftLineRKF::DisablePlotting() { m_view = nullptr; }

void DriftLineRKF::SetGainFluctuationsFixed(const double gain) {

  if (gain > 1.) {
    std::cout << m_className << "::SetGainFluctuationsFixed: "
              << "Avalanche size set to " << gain << ".\n";
  } else {
    std::cout << m_className << "::SetGainFluctuationsFixed:\n    "
              << "Avalanche size will be given by "
              << "the integrated Townsend coefficient.\n";
  }
  m_gain = gain; 
  m_gainFluctuations = GainFluctuations::None;
}

void DriftLineRKF::SetGainFluctuationsPolya(const double theta, 
                                            const double mean,
                                            const bool quiet) {

  if (theta < 0.) {
    std::cerr << m_className << "::SetGainFluctuationsPolya: "
              << "Shape parameter must be >= 0.\n";
    return;
  }  
  if (!quiet) {
    if (mean > 1.) {
      std::cout << m_className << "::SetGainFluctuationsPolya: "
                << "Mean avalanche size set to " << mean << ".\n";
    } else {
      std::cout << m_className << "::SetGainFluctuationsPolya:\n    "
                << "Mean avalanche size will be given by "
                << "the integrated Townsend coefficient.\n";
    }
  }
  m_gain = mean;
  m_theta = theta;
  m_gainFluctuations = GainFluctuations::Polya;
} 

bool DriftLineRKF::DriftElectron(const double x0, const double y0,
                                 const double z0, const double t0) {
  std::vector<std::array<double, 3> > x;
  std::vector<double> t;
  int status = 0;
  const bool ok = DriftLine({x0, y0, z0}, t0, Particle::Electron, 
                            t, x, status);
  std::vector<double> ne(t.size(), 1.);
  std::vector<double> ni(t.size(), 0.);
  std::vector<double> nn(t.size(), 0.);
  double scale = 1.;
  if (ok) {
    if (m_doAvalanche) Avalanche(Particle::Electron, x, ne, ni, nn, scale);
    if (m_doSignal) {
      ComputeSignal(Particle::Electron, scale * m_scaleE, t, x, ne);
    }
    if (m_doAvalanche) {
      if (m_doIonTail) AddIonTail(t, x, ni, scale);
      if (m_doNegativeIonTail) AddNegativeIonTail(t, x, nn, scale);
    }
  }
  m_nE = scale * ne.back();
  m_nI = scale * std::accumulate(ni.begin(), ni.end(), 0.);
  std::swap(m_x, x);
  std::swap(m_t, t);
  m_particle = Particle::Electron;
  m_status = status;
  return ok;
}

bool DriftLineRKF::AddIonTail(const std::vector<double>& te,
                              const std::vector<Vec>& xe,
                              const std::vector<double>& ni,
                              const double scale) const {
  // SIGETR, SIGIOR
  const size_t nPoints = te.size();
  if (nPoints < 2 || ni.size() != nPoints) return false;
  // Loop over the electron track.
  for (size_t i = 1; i < nPoints; ++i) {
    // Skip points at which there are no ions yet.
    if (scale * ni[i] < 1.) continue;
    // Skip also points with a negligible contribution.
    // if (scale * ni[i] < threshold * m_nI) continue;
    // Compute the ion drift line.
    std::vector<double> ti;
    std::vector<Vec> xi;
    int stat = 0;
    if (!DriftLine(xe[i], te[i], Particle::Ion, ti, xi, stat)) {
      std::cerr << m_className << "::AddIonTail:\n"
                << "    Unable to obtain an ion tail; tail not added.\n";
      return false;
    }
    if (m_debug) {
      std::cout << m_className << "::AddIonTail: Origin = " << PrintVec(xe[i])
                << ", n = " << ti.size() << ", status = " << stat << "\n";
    }
    // Compute the contribution of the drift line to the signal.
    ComputeSignal(Particle::Ion, scale * m_scaleI * ni[i], ti, xi, {});
  }
  return true;
}

bool DriftLineRKF::AddNegativeIonTail(
    const std::vector<double>& te, const std::vector<Vec>& xe,
    const std::vector<double>& nn, const double scale) const {
  const size_t nPoints = te.size();
  if (nPoints < 2 || nn.size() != nPoints) return false;
  // Loop over the electron track.
  for (size_t i = 1; i < nPoints; ++i) {
    // Skip points at which there are no negative ions yet.
    if (scale * nn[i] < Small) continue;
    // Compute the negative ion drift line.
    std::vector<double> tn;
    std::vector<Vec> xn;
    int stat = 0;
    if (!DriftLine(xe[i], te[i], Particle::NegativeIon, tn, xn, stat)) {
      std::cerr << m_className << "::AddNegativeIonTail:\n"
                << "    Unable to obtain a negative ion tail.\n";
      return false;
    }
    // Compute the contribution of the drift line to the signal.
    ComputeSignal(Particle::NegativeIon, scale * m_scaleI * nn[i], tn, xn, {});
  }
  return true;
}

bool DriftLineRKF::DriftPositron(const double x0, const double y0,
                                 const double z0, const double t0) {
  std::vector<std::array<double, 3> > x;
  std::vector<double> t;
  int status = 0;
  const bool ok = DriftLine({x0, y0, z0}, t0, Particle::Positron, 
                            t, x, status);
  if (ok && m_doSignal) {
    ComputeSignal(Particle::Positron, m_scaleE, t, x, {});
  }
  std::swap(m_x, x);
  std::swap(m_t, t);
  m_particle = Particle::Positron;
  m_status = status;
  return ok;
}

bool DriftLineRKF::DriftHole(const double x0, const double y0, const double z0,
                             const double t0) {
  std::vector<std::array<double, 3> > x;
  std::vector<double> t;
  int status = 0;
  const bool ok = DriftLine({x0, y0, z0}, t0, Particle::Hole, t, x, status);
  if (ok && m_doSignal) {
    ComputeSignal(Particle::Hole, m_scaleH, t, x, {});
  }
  std::swap(m_x, x);
  std::swap(m_t, t);
  m_particle = Particle::Hole;
  m_status = status;
  return ok;
}

bool DriftLineRKF::DriftIon(const double x0, const double y0, const double z0,
                            const double t0) {
  std::vector<std::array<double, 3> > x;
  std::vector<double> t;
  int status = 0;
  const bool ok = DriftLine({x0, y0, z0}, t0, Particle::Ion, t, x, status);
  if (ok && m_doSignal) {
    ComputeSignal(Particle::Ion, m_scaleI, t, x, {});
  }
  std::swap(m_x, x);
  std::swap(m_t, t);
  m_particle = Particle::Ion;
  m_status = status;
  return ok;
}

bool DriftLineRKF::DriftNegativeIon(const double x0, const double y0,
                                    const double z0, const double t0) {
  std::vector<std::array<double, 3> > x;
  std::vector<double> t;
  int status = 0;
  const bool ok = DriftLine({x0, y0, z0}, t0, Particle::NegativeIon, 
                            t, x, status);
  if (ok && m_doSignal) {
    ComputeSignal(Particle::NegativeIon, m_scaleI, t, x, {});
  }
  std::swap(m_x, x);
  std::swap(m_t, t);
  m_particle = Particle::NegativeIon;
  m_status = status;
  return ok;
}

bool DriftLineRKF::DriftLine(const Vec& xi, const double ti, 
                             const Particle particle,
                             std::vector<double>& ts,
                             std::vector<Vec>& xs, int& flag) const {

  // -----------------------------------------------------------------------
  //    DLCALC - Subroutine doing the actual drift line calculations. 
  //             The calculations are based on a Runge-Kutta-Fehlberg method
  //             which has the advantage of controlling the stepsize and the
  //             error while needing only relatively few calls to EFIELD.
  //             Full details are given in the reference quoted below.
  //    VARIABLES : H          : Current stepsize (it is in fact a delta t).
  //                HPREV      : Stores the previous value of H (comparison)
  //                INITCH     : Used for checking initial stepsize (1 = ok)
  //                Other variables such as F0, F1, F2, F3, PHII, PHIII,
  //                CI. ,CII. , BETA.. etc   are explained in the reference.
  //    REFERENCE : Stoer + Bulirsch, Einfuhrung in die Numerische
  //                Mathematic II, chapter 7, page 122, 1978, HTB, Springer.
  // -----------------------------------------------------------------------

  // Set the numerical constants for the RKF integration.
  constexpr double c10 = 214. / 891.;
  constexpr double c11 = 1. / 33.;
  constexpr double c12 = 650. / 891.;
  constexpr double c20 = 533. / 2106.;
  constexpr double c22 = 800. / 1053.;
  constexpr double c23 = -1. / 78.;

  constexpr double b10 = 1. / 4.;
  constexpr double b20 = -189. / 800.;
  constexpr double b21 = 729. / 800.;
  constexpr double b30 = 214. / 891.;
  constexpr double b31 = 1. / 33.;
  constexpr double b32 = 650. / 891.;

  // Reset the drift line.
  ts.clear();
  xs.clear();
  // Reset the status flag.
  flag = StatusAlive;

  // Check if the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::DriftLine: Sensor is not defined.\n";
    flag = StatusCalculationAbandoned;
    return false;
  }

  // Get the sensor's bounding box.
  double xmin = 0., xmax = 0.;
  double ymin = 0., ymax = 0.;
  double zmin = 0., zmax = 0.;
  bool bbox = m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax);

  // Set the charge of the drifting particle.
  const double charge = Charge(particle);

  // Initialise the current position and velocity.
  Vec x0 = xi;
  Vec v0 = GetVelocity(x0, particle, flag);
  if (flag != 0) {
    std::cerr << m_className << "::DriftLine:\n"
              << "    Cannot retrieve drift velocity at initial position "
              << PrintVec(x0) << ".\n";
    return false;
  }

  const double speed0 = Mag(v0);
  if (speed0 < Small) {
    std::cerr << m_className << "::DriftLine: "
              << "Zero velocity at initial position.\n";
    return false;
  }

  // Initialise time step and previous time step.
  double h = m_accuracy / speed0;
  double hprev = h;
  double t0 = ti;

  // Set the initial point.
  ts.push_back(t0);
  xs.push_back(x0);

  if (m_debug) {
    std::cout << m_className << "::DriftLine:\n"
              << "    Initial step size: " << h << " ns.\n";
  }
  int initCycle = 3;
  bool ok = true;
  while (ok) {
    // Get the velocity at the first probe point.
    Vec x1 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x1[i] += h * b10 * v0[i];
    }
    int stat = 0;
    const Vec v1 = GetVelocity(x1, particle, stat);
    if (stat == StatusCalculationAbandoned) {
      flag = stat;
      break;
    } else if (stat != 0) {
      if (m_debug) std::cout << "    Point 1 outside.\n";
      if (!Terminate(x0, x1, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      } else {
        flag = stat;
      }
      break;
    }
    // Get the velocity at the second probe point.
    Vec x2 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x2[i] += h * (b20 * v0[i] + b21 * v1[i]);
    }
    const Vec v2 = GetVelocity(x2, particle, stat);
    if (stat == StatusCalculationAbandoned) {
      flag = stat;
      break;
    } else if (stat != 0) {
      if (m_debug) std::cout << "    Point 2 outside.\n";
      if (!Terminate(x0, x2, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      } else {
        flag = stat;
      }
      break;
    }
    // Get the velocity at the third probe point.
    Vec x3 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x3[i] += h * (b30 * v0[i] + b31 * v1[i] + b32 * v2[i]);
    }
    const Vec v3 = GetVelocity(x3, particle, stat);
    if (stat == StatusCalculationAbandoned) {
      flag = stat;
      break;
    } else if (stat != 0) {
      if (m_debug) std::cout << "    Point 3 outside.\n";
      if (!Terminate(x0, x3, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      } else {
        flag = stat;
      }
      break;
    }
    // Check if we crossed a wire.
    double xw = 0., yw = 0., zw = 0., rw = 0.;
    if (m_sensor->CrossedWire(x0[0], x0[1], x0[2], 
                              x1[0], x1[1], x1[2], xw, yw, zw, true, rw) ||
        m_sensor->CrossedWire(x0[0], x0[1], x0[2], 
                              x2[0], x2[1], x2[2], xw, yw, zw, true, rw) ||
        m_sensor->CrossedWire(x0[0], x0[1], x0[2], 
                              x3[0], x3[1], x3[2], xw, yw, zw, true, rw)) {
      if (m_debug) std::cout << "    Crossed wire.\n";
      if (DriftToWire(xw, yw, rw, particle, ts, xs, stat)) {
        flag = stat;
      } else if (h > Small) {
        h *= 0.5;
        continue;
      } else {
        std::cerr << m_className << "::DriftLine: Step size too small. Stop.\n";
        flag = StatusCalculationAbandoned;
      }
      break;
    }
    // Check if we are inside the trap radius of a wire.
    if (particle != Particle::Ion) {
      if (m_sensor->InTrapRadius(charge, x1[0], x1[1], x1[2], xw, yw, rw)) {
        if (!DriftToWire(xw, yw, rw, particle, ts, xs, flag)) {
          flag = StatusCalculationAbandoned;
        }
        break;
      }
      if (m_sensor->InTrapRadius(charge, x2[0], x2[1], x2[2], xw, yw, rw)) {
        if (!DriftToWire(xw, yw, rw, particle, ts, xs, flag)) {
          flag = StatusCalculationAbandoned;
        }
        break;
      }
      if (m_sensor->InTrapRadius(charge, x3[0], x3[1], x3[2], xw, yw, rw)) {
        if (!DriftToWire(xw, yw, rw, particle, ts, xs, flag)) {
          flag = StatusCalculationAbandoned;
        }
        break;
      }
    } 
    // Check if we crossed a plane.
    Vec xp = {0., 0., 0.};
    if (m_sensor->CrossedPlane(x0[0], x0[1], x0[2], 
                               x1[0], x1[1], x1[2], xp[0], xp[1], xp[2]) ||
        m_sensor->CrossedPlane(x0[0], x0[1], x0[2], 
                               x2[0], x2[1], x2[2], xp[0], xp[1], xp[2]) ||
        m_sensor->CrossedPlane(x0[0], x0[1], x0[2], 
                               x3[0], x3[1], x3[2], xp[0], xp[1], xp[2])) {
      // DLCPLA
      ts.push_back(t0 + Dist(x0, xp) / Mag(v0));
      xs.push_back(xp);
      flag = StatusHitPlane;
      break; 
    }
    // Calculate the correction terms.
    Vec phi1 = {0., 0., 0.};
    Vec phi2 = {0., 0., 0.};
    for (size_t i = 0; i < 3; ++i) {
      phi1[i] = c10 * v0[i] + c11 * v1[i] + c12 * v2[i];
      phi2[i] = c20 * v0[i] + c22 * v2[i] + c23 * v3[i];
    }
    // Check if the step length is valid.
    const double phi1mag = Mag(phi1);
    if (phi1mag < Small) {
      std::cerr << m_className << "::DriftLine: Step has zero length. Stop.\n";
      flag = StatusCalculationAbandoned;
      break;
    } else if (m_useStepSizeLimit && h * phi1mag > m_maxStepSize) {
      if (m_debug) {
        std::cout << "    Step is considered too long. H is reduced.\n";
      }
      h = 0.5 * m_maxStepSize / phi1mag;
      continue;
    } else if (bbox) {
      // Don't allow h to become too large in view of the time resolution.
      if (h * fabs(phi1[0]) > 0.1 * fabs(xmax - xmin) ||
          h * fabs(phi1[1]) > 0.1 * fabs(ymax - ymin)) {
        h *= 0.5;
        if (m_debug) {
          std::cout << "    Step is considered too long. H is halved.\n";
        }
        continue;
      }
    } else if (m_rejectKinks && xs.size() > 1) {
      const unsigned int np = xs.size();
      const auto& x = xs[np - 1];
      const auto& xprev = xs[np - 2];
      if (phi1[0] * (x[0] - xprev[0]) + phi1[1] * (x[1] - xprev[1]) +
          phi1[2] * (x[2] - xprev[2]) < 0.) {
        std::cerr << m_className << "::DriftLine: Bending angle > 90 degree.\n";
        flag = StatusSharpKink;
        break;
      }
    }
    if (m_debug) std::cout << "    Step size ok.\n";
    // Update the position and time.
    for (size_t i = 0; i < 3; ++i) x0[i] += h * phi1[i];
    t0 += h;
    if (!m_sensor->IsInside(x0[0], x0[1], x0[2])) {
      // The new position is not inside a valid drift medium.
      // Terminate the drift line.
      if (m_debug) std::cout << "    New point is outside. Terminate.\n";
      if (!Terminate(xs.back(), x0, particle, ts, xs)) {
        flag = StatusCalculationAbandoned;
      }
      break;
    }
    // Add the new point to the drift line.
    ts.push_back(t0);
    xs.push_back(x0);
    // Adjust the step size according to the accuracy of the two estimates.
    hprev = h;
    const double dphi = fabs(phi1[0] - phi2[0]) + fabs(phi1[1] - phi2[1]) +
                        fabs(phi1[2] - phi2[2]);
    if (dphi > 0) {
      h = sqrt(h * m_accuracy / dphi);
      if (m_debug) std::cout << "    Adapting H to " << h << ".\n";
    } else {
      h *= 2;
      if (m_debug) std::cout << "    H increased by factor two.\n";
    }
    // Make sure that H is different from zero; this should always be ok.
    if (h < Small) {
      std::cerr << m_className << "::DriftLine: Step size is zero. Stop.\n";
      flag = StatusCalculationAbandoned;
      break;
    }
    // Check the initial step size.
    if (initCycle > 0 && h < 0.2 * hprev) {
      if (m_debug) std::cout << "    Reinitialise step size.\n";
      --initCycle;
      t0 = ti;
      x0 = xi;
      ts = {t0};
      xs = {x0};
      continue;
    }
    initCycle = 0;
    // Don't allow H to grow too quickly
    if (h > 10 * hprev) {
      h = 10 * hprev;
      if (m_debug) {
        std::cout << "    H restricted to 10 times the previous value.\n";
      }
    }
    // Stop in case H tends to become too small.
    if (h * (fabs(phi1[0]) + fabs(phi1[1]) + fabs(phi1[2])) < m_accuracy) {
      std::cerr << m_className << "::DriftLine: Step size has become smaller "
                << "than int. accuracy. Stop.\n";
      flag = StatusCalculationAbandoned;
      break;
    }
    // Update the velocity.
    v0 = v3;
  }
  if (m_view) {
    // If requested, add the drift line to a plot.
    size_t id = 0;
    const size_t nPoints = xs.size();
    m_view->NewDriftLine(particle, nPoints, id, xi[0], xi[1], xi[2]);
    for (size_t i = 0; i < nPoints; ++i) {
      const auto& x = xs[i];
      m_view->SetDriftLinePoint(id, i, x[0], x[1], x[2]);
    }
  }
  if (flag == StatusCalculationAbandoned) return false;
  return true;
}

bool DriftLineRKF::Avalanche(const Particle particle,
                             const std::vector<Vec>& xs,
                             std::vector<double>& ne,
                             std::vector<double>& ni, 
                             std::vector<double>& nn, double& scale) const {

  // SIGETR
  const size_t nPoints = xs.size();
  if (nPoints < 2) return true;
  // Locations and weights for 6-point Gaussian integration.
  const size_t nG = 6;
  auto tg = Numerics::GaussLegendreNodes6();
  auto wg = Numerics::GaussLegendreWeights6();

  ne.assign(nPoints, 1.);
  ni.assign(nPoints, 0.);
  nn.assign(nPoints, 0.);
  bool start = false;
  bool overflow = false;
  // Loop over the drift line.
  for (size_t i = 1; i < nPoints; ++i) {
    const auto& xp = xs[i - 1];
    const auto& x = xs[i];
    const Vec dx = {x[0] - xp[0], x[1] - xp[1], x[2] - xp[2]};
    // Calculate the integrated Townsend and attachment coefficients.
    double alpsum = 0.;
    double etasum = 0.;
    for (size_t j = 0; j < nG; ++j) {
      const double f = 0.5 * (1. + tg[j]);
      Vec xj = xp;
      for (size_t k = 0; k < 3; ++k) xj[k] += f * dx[k];
      const double alp = GetAlpha(xj, particle);
      if (alp < 0.) {
        std::cerr << m_className << "::Avalanche:\n    Cannot retrieve alpha at "
                  << "drift line point " << i  << ", segment " << j << ".\n";
        continue;
      }
      const double eta = GetEta(xj, particle);
      if (eta < 0.) {
        std::cerr << m_className << "::Avalanche:\n    Cannot retrieve eta at "
                  << "drift line point " << i  << ", segment " << j << ".\n";
        continue;
      }
      alpsum += wg[j] * alp;
      etasum += wg[j] * eta;
    }
    alpsum *= 0.5;
    etasum *= 0.5;
    if (alpsum > 1.e-6 && !start) {
      if (m_debug) {
        std::cout << m_className << "::Avalanche: Avalanche starts at step " 
                  << i << ".\n";
      }
      start = true;
    }
    const double d = Mag(dx);
    // Update the number of electrons.
    constexpr double expmax = 30.;
    const double logp = log(std::max(1., ne[i - 1]));
    if (logp + d * (alpsum - etasum) > expmax) {
      overflow = true;
      ne[i] = exp(expmax);
    } else { 
      ne[i] = ne[i - 1] * exp(d * (alpsum - etasum));
    }
    // Update the number of ions.
    if (logp + d * alpsum > expmax) {
      overflow = true;
      ni[i] = exp(expmax);
    } else {
      ni[i] = ne[i - 1] * (exp(d * alpsum) - 1);
    }
    nn[i] = std::max(ne[i - 1] + ni[i] - ne[i], 0.);
  } 
  if (overflow) {
    std::cerr << m_className << "::Avalanche:\n    "
              << "Warning: Integrating the Townsend coefficients "
              << "would lead to exponential overflow.\n    "
              << "Avalanche truncated.\n";
  }
  const double qe = ne.back();
  const double qi = std::accumulate(ni.begin(), ni.end(), 0.);
  scale = 1.;
  if (qi > 1. && 
      !(m_gainFluctuations == GainFluctuations::None && m_gain < 1.)) {
    constexpr double eps = 1.e-4;
    const double gain = m_gain > 1. ? m_gain : ComputeGain(xs, particle, eps);
    double q1 = gain;
    if (m_gainFluctuations == GainFluctuations::Polya) {
      for (unsigned int i = 0; i < 100; ++i) {
        q1 = gain * RndmPolya(m_theta);
        if (q1 >= 1.) break;
      }
      q1 = std::max(q1, 1.);
    }
    q1 *= ComputeLoss(xs, particle, eps);
    scale = (q1 + 1.) / (qi + 1.); 
  }
  if (m_debug) {
    const double qn = std::accumulate(nn.begin(), nn.end(), 0.);
    std::cout << m_className << "::Avalanche:\n    "
              << "Final number of electrons: " << qe << "\n    "
              << "Number of positive ions:   " << qi << "\n    "
              << "Number of negative ions:   " << qn << "\n    "
              << "Charge scaling factor:     " << scale << "\n    "
              << "Avalanche development:\n Step      Electrons     Ions\n";
    for (unsigned int i = 0; i < nPoints; ++i) {
      std::printf("%6d %15.7f %15.7f\n", i, scale * ne[i], scale * ni[i]);
    }
  }
  return true;
}

double DriftLineRKF::GetArrivalTimeSpread(const double eps) const {
  return ComputeSigma(m_x, m_particle, eps);
} 

double DriftLineRKF::ComputeSigma(const std::vector<Vec>& x,
                                  const Particle particle,
                                  const double eps) const {

  // -----------------------------------------------------------------------
  //    DLCDF1 - Routine returning the integrated diffusion coefficient of
  //             the current drift line. The routine uses an adaptive
  //             Simpson integration.
  // -----------------------------------------------------------------------

  const size_t nPoints = x.size();
  // Return straight away if there is only one point.
  if (nPoints < 2) return 0.;

  // First get a rough estimate of the result.
  double crude = 0.;
  double varPrev = 0.;
  for (size_t i = 0; i < nPoints; ++i) {
    // Get the variance at this point.
    const double var = GetVar(x[i], particle);
    if (var < 0.) {
      std::cerr << m_className << "::ComputeSigma:\n"
              << "    Cannot retrieve variance at point " << i << ".\n";
      continue;
    }
    if (i > 0) crude += 0.5 * Dist(x[i - 1], x[i]) * (var + varPrev);
    varPrev = var;
  }
  crude = sqrt(crude);

  const double tol = eps * crude; 
  double sum = 0.;
  for (size_t i = 0; i < nPoints - 1; ++i) {
    sum += IntegrateDiffusion(x[i], x[i + 1], particle, tol);
  }
  return sqrt(sum);
}

double DriftLineRKF::GetGain(const double eps) const {
  if (m_status == StatusCalculationAbandoned) return 1.;
  return ComputeGain(m_x, m_particle, eps);
}

double DriftLineRKF::ComputeGain(const std::vector<Vec>& x,
                                 const Particle particle, 
                                 const double eps) const {

  // -----------------------------------------------------------------------
  //    DLCTW1 - Routine returning the multiplication factor for the current
  //             drift line. The routine uses an adaptive Simpson style
  //             integration.
  // -----------------------------------------------------------------------

  const size_t nPoints = m_x.size();
  // Return straight away if there is only one point.
  if (nPoints < 2) return 1.;
  if (particle == Particle::Ion || particle == Particle::NegativeIon) {
    return 1.;
  }

  // First get a rough estimate of the result.
  double crude = 0.;
  double alphaPrev = 0.;
  for (size_t i = 0; i < nPoints; ++i) {
    // Get the Townsend coefficient at this point.
    const double alpha = GetAlpha(x[i], particle);
    if (alpha < 0.) {
      std::cerr << m_className << "::ComputeGain:\n"
                << "    Cannot retrieve alpha at point " << i << ".\n";
      continue;
    }
    if (i > 0) crude += 0.5 * Dist(x[i - 1], x[i]) * (alpha + alphaPrev);
    alphaPrev = alpha;
  }
  // Stop if the rough estimate is negligibly small.
  if (crude < Small) return 1.;

  // Calculate the integration tolerance based on the rough estimate.
  const double tol = eps * crude;
  double sum = 0.;
  for (size_t i = 0; i < nPoints - 1; ++i) {
    sum += IntegrateAlpha(x[i], x[i + 1], particle, tol);
  }
  return exp(sum);
}

double DriftLineRKF::GetLoss(const double eps) const {
  if (m_status == StatusCalculationAbandoned) return 1.;
  return ComputeLoss(m_x, m_particle, eps);
}

double DriftLineRKF::ComputeLoss(const std::vector<Vec>& x,
                                 const Particle particle,
                                 const double eps) const {

  // -----------------------------------------------------------------------
  //    DLCAT1 - Routine returning the attachment losses for the current
  //             drift line. The routine uses an adaptive Simpson style
  //             integration.
  // -----------------------------------------------------------------------

  const size_t nPoints = x.size();
  // Return straight away if there is only one point.
  if (nPoints < 2) return 1.;
  if (particle == Particle::Ion || particle == Particle::NegativeIon) {
    return 1.;
  }

  // First get a rough estimate of the result.
  double crude = 0.;
  double etaPrev = 0.;
  for (size_t i = 0; i < nPoints; ++i) {
    // Get the attachment coefficient at this point.
    const double eta = GetEta(x[i], particle);
    if (eta < 0.) {
      std::cerr << m_className << "::ComputeLoss:\n"
                << "    Cannot retrieve eta at point " << i << ".\n";
      continue;
    }
    if (i > 0) crude += 0.5 * Dist(x[i - 1], x[i]) * (eta + etaPrev);
    etaPrev = eta;
  }

  // Calculate the integration tolerance based on the rough estimate.
  const double tol = eps * crude;
  double sum = 0.;
  for (size_t i = 0; i < nPoints - 1; ++i) {
    sum += IntegrateEta(x[i], x[i + 1], particle, tol);
  }
  return exp(-sum);
}

double DriftLineRKF::GetPathLength() const {

  const size_t nPoints = m_x.size();
  if (nPoints < 2) return 0.;
  double path = 0.;
  for (size_t i = 1; i < nPoints; ++i) {
    path += Dist(m_x[i - 1], m_x[i]);
  }
  return path;
}

int DriftLineRKF::GetField(const std::array<double, 3>& x,
                           double& ex, double& ey, double& ez,
                           double& bx, double& by, double& bz,
                           Medium*& medium) const {
  int status = 0;
  m_sensor->MagneticField(x[0], x[1], x[2], bx, by, bz, status);
  m_sensor->ElectricField(x[0], x[1], x[2], ex, ey, ez, medium, status);
  return status;
}

Vec DriftLineRKF::GetVelocity(const std::array<double, 3>& x,
                              const Particle particle,
                              int& status) const {
  Vec v = {0., 0., 0.};
  status = 0;
  // Stop if we are outside the drift area.
  if (!m_sensor->IsInArea(x[0], x[1], x[2])) {
    status = StatusLeftDriftArea;
    return v;
  } 
  if (m_useVelocityMap && 
      particle != Particle::Ion && particle != Particle::NegativeIon) {
    // We assume there is only one component with a velocity map.
    const auto nComponents = m_sensor->GetNumberOfComponents();
    for (size_t i = 0; i < nComponents; ++i) {
      auto cmp = m_sensor->GetComponent(i);
      if (!cmp->HasVelocityMap()) continue;
      bool ok = false;
      if (particle == Particle::Electron || particle == Particle::Positron) {
        ok = cmp->ElectronVelocity(x[0], x[1], x[2], v[0], v[1], v[2]);
      } else if (particle == Particle::Hole) {
        ok = cmp->HoleVelocity(x[0], x[1], x[2], v[0], v[1], v[2]);
      }
      if (!ok) continue;
      // Seems to have worked.
      if (particle == Particle::Positron) {
        for (unsigned int k = 0; k < 3; ++k) v[k] *= -1;
      }
      return v;
    }
  }
  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;
  // Stop if we are outside a valid drift medium.
  status = GetField(x, ex, ey, ez, bx, by, bz, medium);
  if (status != 0) return v;
  bool ok = false;
  if (particle == Particle::Electron) {
    ok = medium->ElectronVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Ion) {
    ok = medium->IonVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Hole) {
    ok = medium->HoleVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } else if (particle == Particle::Positron) {
    ok = medium->ElectronVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
    for (unsigned int i = 0; i < 3; ++i) v[i] *= -1;
  } else if (particle == Particle::NegativeIon) {
    ok = medium->NegativeIonVelocity(ex, ey, ez, bx, by, bz, v[0], v[1], v[2]);
  } 
  if (!ok) {
    std::cerr << m_className << "::GetVelocity:\n"
              << "    Cannot retrieve drift velocity at " 
              << PrintVec(x) << ".\n";
    status = StatusCalculationAbandoned;
  }
  return v;
}

bool DriftLineRKF::GetDiffusion(const std::array<double, 3>& x,
                                const Particle particle,
                                double& dl, double& dt) const {
  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;
  if (GetField(x, ex, ey, ez, bx, by, bz, medium) != 0) return false;

  if (particle == Particle::Electron || particle == Particle::Positron) {
    return medium->ElectronDiffusion(ex, ey, ez, bx, by, bz, dl, dt);
  } else if (particle == Particle::Ion || particle == Particle::NegativeIon) {
    return medium->IonDiffusion(ex, ey, ez, bx, by, bz, dl, dt);
  } else if (particle == Particle::Hole) {
    return medium->HoleDiffusion(ex, ey, ez, bx, by, bz, dl, dt);
  }
  return false;
}

double DriftLineRKF::GetVar(const std::array<double, 3>& x,
                            const Particle particle) const {
  // Get the drift velocity.
  int stat = 0;
  const Vec v = GetVelocity(x, particle, stat);
  if (stat != 0) return -1.;

  const double speed = Mag(v);
  if (speed < Small) {
    std::cerr << m_className << "::GetVariance: Zero velocity.\n";
    return -1.;
  }
  // Get the diffusion coefficients.
  double dl = 0., dt = 0.;
  if (!GetDiffusion(x, particle, dl, dt)) return -1.;

  const double sigma = dl / speed;
  return sigma * sigma;
}

double DriftLineRKF::GetAlpha(const std::array<double, 3>& x,
                              const Particle particle) const {
  double alpha = 0.;
  if (m_useTownsendMap && (particle == Particle::Electron || 
      particle == Particle::Hole || particle == Particle::Positron)) {
    const auto nComponents = m_sensor->GetNumberOfComponents();
    for (size_t i = 0; i < nComponents; ++i) {
      auto cmp = m_sensor->GetComponent(i);
      if (!cmp->HasTownsendMap()) continue;
      if (particle == Particle::Electron || particle == Particle::Positron) {
        if (!cmp->ElectronTownsend(x[0], x[1], x[2], alpha)) continue;
      } else {
        if (!cmp->HoleTownsend(x[0], x[1], x[2], alpha)) continue;
      }
      return alpha;
    }
  }
  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;
  if (GetField(x, ex, ey, ez, bx, by, bz, medium) != 0) return -1.;

  if (particle == Particle::Electron || particle == Particle::Positron) {
    medium->ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  } else if (particle == Particle::Hole) {
    medium->HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
  return alpha;
}

double DriftLineRKF::GetEta(const std::array<double, 3>& x,
                            const Particle particle) const {

  double ex = 0., ey = 0., ez = 0.;
  double bx = 0., by = 0., bz = 0.;
  Medium* medium = nullptr;
  if (GetField(x, ex, ey, ez, bx, by, bz, medium) != 0) return -1.;
  double eta = 0.;
  if (particle == Particle::Electron) {
    medium->ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  } else if (particle == Particle::Hole) {
    medium->HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  return eta;
}

bool DriftLineRKF::Terminate(const std::array<double, 3>& xx0,
                             const std::array<double, 3>& xx1,
                             const Particle particle,
                             std::vector<double>& ts,
                             std::vector<Vec>& xs) const {

  // -----------------------------------------------------------------------
  //    DLCFMP - Terminates drift line calculation by making a last step
  //             to the boundary of the mesh or the drift medium.
  //    VARIABLES : XX0: Last point in drift medium.
  //                XX1: Estimated step, outside drift medium.
  // -----------------------------------------------------------------------

  // Check the validity of the initial point.
  int status = 0;
  const Vec vv0 = GetVelocity(xx0, particle, status);
  if (status != 0) {
    std::cerr << m_className << "::Terminate:\n"
              << "    Cannot retrieve drift velocity at initial point.\n";
    return false;
  }
  double speed = Mag(vv0);
  if (speed < Small) {
    std::cerr << m_className << "::Terminate:\n"
              << "    Zero velocity at initial position.\n";
    return false;
  }

  // Final point just inside the medium.
  Vec x0 = xx0;
  // Final point just outside the medium.
  Vec x1 = xx1;
  // Perform some bisections.
  constexpr unsigned int nBisections = 20;
  for (unsigned int i = 0; i < nBisections; ++i) {
    // Quit bisection when interval becomes too small.
    bool small = true;
    for (unsigned int j = 0; j < 3; ++j) {
      if (fabs(x1[j] - x0[j]) > 1.e-6 * (fabs(x0[j]) + fabs(x1[j]))) {
        small = false;
        break;
      }
    } 
    if (small) {
      if (m_debug) {
        std::cout << m_className << "::Terminate:\n"
                  << "    Bisection ended at cycle " << i << ".\n";
      }
      break; 
    }
    // Calculate the mid point.
    const Vec xm = MidPoint(x0, x1);
    // Halve the step.
    if (m_sensor->IsInside(xm[0], xm[1], xm[2]) &&
        m_sensor->IsInArea(xm[0], xm[1], xm[2])) {
      x0 = xm;
    } else {
      x1 = xm;
    }
  }

  // Compute drift velocity at the end of the step.
  Vec v0 = GetVelocity(x0, particle, status);
  if (status != 0) {
    std::cerr << m_className << "::Terminate:\n"
              << "    Warning: Unable to compute mean velocity at last step.\n";
  } else {
    speed = 0.5 * (speed + Mag(v0));
  }
  // Calculate the time step.
  const double dt = Dist(xx0, x0) / speed;
  // Add the last point, just inside the drift area.
  ts.push_back(ts.back() + dt);
  xs.push_back(x0);
  return true;
}

bool DriftLineRKF::DriftToWire(const double xw, const double yw,
                               const double rw, const Particle particle,
                               std::vector<double>& ts, 
                               std::vector<Vec>& xs, int& stat) const {

  // -----------------------------------------------------------------------
  //   DLCWIR - Terminates drift line calculation by making some last steps
  //            towards the surface of the wire on which it is supposed to
  //            end. The precision is controlled in order to obtain a
  //            good estimate of the total remaining drift-time.
  // -----------------------------------------------------------------------

  // Get the starting point.
  Vec x0 = xs.back();
  double t0 = ts.back() - ts.front();
  if (m_debug) {
    std::cout << m_className << "::DriftToWire:\n    Drifting from ("
              << x0[0] << ", " << x0[1] << ") to wire at ("
              << xw << ", " << yw << ") with radius " << rw << " cm.\n";
  }

  // Get the initial drift velocity.
  int status = 0;
  Vec v0 = GetVelocity(x0, particle, status);
  if (status != 0) {
    std::cerr << m_className << "::DriftToWire:\n"
              << "    Cannot retrieve initial drift velocity.\n";
    return false;
  } 

  // Estimate the time needed to reach the wire
  // assuming a straight-line trajectory and constant velocity.
  double dt = (Mag(xw - x0[0], yw - x0[1]) - rw) / Mag(v0[0], v0[1]);
  if (m_debug) {
    std::cout << "    Estimated time needed to reach the wire: " 
              << dt << " ns.\n";
  }

  constexpr unsigned int nMaxSplit = 10;
  unsigned int nSplit = 0;
  // Move towards the wire.
  bool onwire = false;
  const double r2 = rw * rw;
  while (!onwire && dt > 1.e-6 * t0) {
    // Calculate the estimated end point.
    Vec x1 = x0;
    for (unsigned int j = 0; j < 3; ++j) x1[j] += dt * v0[j];
    // Make sure we are not moving away from the wire.
    const double xinp0 = (x1[0] - x0[0]) * (xw - x0[0]) + 
                         (x1[1] - x0[1]) * (yw - x0[1]);
    if (xinp0 < 0.) {
      if (m_debug) {
        std::cerr << "    Particle moves away from the wire. Quit.\n";
      }
      return false;
    }
    // Check if the end point is inside the wire or the wire was crossed.
    const double xinp1 = (x0[0] - x1[0]) * (xw - x1[0]) + 
                         (x0[1] - x1[1]) * (yw - x1[1]);
    if (xinp1 < 0.) {
      if (Mag2(xw - x1[0], yw - x1[1]) <= r2) {
        onwire = true;
        if (m_debug) std::cout << "    Inside.\n";
      }
    } else {
      if (m_debug) std::cout << "    Wire crossed.\n";
      onwire = true;
    }
    if (onwire) {
      const double dw0 = Mag(xw - x0[0], yw - x0[1]);
      x1[0] = xw - (rw + BoundaryDistance) * (xw - x0[0]) / dw0;
      x1[1] = yw - (rw + BoundaryDistance) * (yw - x0[1]) / dw0;
      dt = Mag(x1[0] - x0[0], x1[1] - x0[1]) / Mag(v0[0], v0[1]);
      x1[2] = x0[2] + dt * v0[2];
    }
    // Calculate the drift velocity at the end point.
    Vec v1 = GetVelocity(x1, particle, status);
    if (status != 0) {
      std::cerr << m_className << "::DriftToWire:\n"
                << "    Cannot retrieve drift velocity at end point. Quit.\n";
      return false;
    } 
    // Get a point halfway between for an accuracy check.
    const Vec xm = MidPoint(x0, x1);
    // Calculate the drift velocity at the mid point.
    Vec vm = GetVelocity(xm, particle, status);
    if (status != 0) {
      std::cerr << m_className << "::DriftToWire:\n"
                << "    Cannot retrieve drift velocity at mid point. Quit.\n";
      return false;
    }
    // Make sure the velocities are non-zero.
    const double speed0 = Mag(v0[0], v0[1]);
    const double speed1 = Mag(v1[0], v1[1]);
    const double speedm = Mag(vm[0], vm[1]);
    if (speed0 < Small || speed1 < Small || speedm < Small) {
      std::cerr << m_className << "::DriftToWire: Zero velocity. Stop.\n";
      return false;
    }
    const double p0 = 1. / speed0;
    const double p1 = 1. / speed1;
    const double pm = 1. / speedm;
    // Compare first and second order estimates.
    const double d = Mag(x0[0] - x1[0], x0[1] - x1[1]);
    const double tol = 1.e-4 * (1. + fabs(t0));
    if (d * fabs(p0 - 2 * pm + p1) / 3. > tol && nSplit < nMaxSplit) {
      // Accuracy was not good enough so halve the step time.
      if (m_debug) std::cout << "    Reducing step size.\n";
      dt *= 0.5;
      onwire = false;
      ++nSplit;
      continue;
    }
    const double t1 = t0 + d * (p0 + 4 * pm + p1) / 6.;
    // Add a new point to the drift line.
    ts.push_back(ts[0] + t1);
    xs.push_back(x1);
    // Proceed to the next step.
    x0 = x1;
    t0 = t1;
    v0 = v1;
  }
  // Get the wire index (status code inside the wire).
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  m_sensor->ElectricField(xw, yw, 0., ex, ey, ez, medium, stat);
  return true;
}

void DriftLineRKF::PrintDriftLine() const {

  std::cout << m_className << "::PrintDriftLine:\n";
  if (m_x.empty()) {
    std::cout << "    No drift line present.\n";
    return;
  }
  if (m_particle == Particle::Electron) {
    std::cout << "    Particle: electron\n";
  } else if (m_particle == Particle::Ion) {
    std::cout << "    Particle: ion\n";
  } else if (m_particle == Particle::Hole) {
    std::cout << "    Particle: hole\n";
  } else if (m_particle == Particle::Positron) {
    std::cout << "    Particle: positive electron\n";
  } else if (m_particle == Particle::NegativeIon) {
    std::cout << "    Particle: negative ion\n";
  } else {
    std::cout << "    Particle: unknown\n";
  }
  std::cout << "    Status: " << m_status << "\n"
            << "  Step       time [ns]        "
            << "x [cm]          y [cm]          z [cm]\n";
  const unsigned int nPoints = m_x.size();
  for (unsigned int i = 0; i < nPoints; ++i) {
    std::printf("%6d %15.7f %15.7f %15.7f %15.7f\n", 
                i, m_t[i], m_x[i][0], m_x[i][1], m_x[i][2]);
  }
 
}

void DriftLineRKF::GetEndPoint(double& x, double& y, double& z, double& t,
                               int& stat) const {
  if (m_x.empty()) {
    x = y = z = t = 0.;
    stat = m_status;
    return;
  }
  const auto& p = m_x.back();
  x = p[0];
  y = p[1];
  z = p[2];
  t = m_t.back();
  stat = m_status;
}

void DriftLineRKF::GetDriftLinePoint(const size_t i, double& x, double& y,
                                     double& z, double& t) const {
  if (i >= m_x.size()) {
    std::cerr << m_className << "::GetDriftLinePoint: Index out of range.\n";
    return;
  }

  const auto& p = m_x[i];
  x = p[0];
  y = p[1];
  z = p[2];
  t = m_t[i];
}

double DriftLineRKF::IntegrateDiffusion(const std::array<double, 3>& xi,
                                        const std::array<double, 3>& xe,
                                        const Particle particle, 
                                        const double tol) const {
  // Make sure the starting and end points are valid.
  Vec x0 = xi;
  double var0 = GetVar(x0, particle);
  Vec x1 = xe;
  double var1 = GetVar(x1, particle);
  if (var0 < 0. || var1 < 0.) {
    std::cerr << m_className << "::IntegrateDiffusion:\n"
              << "    Cannot retrieve variance at initial point.\n";
    return 0.;
  }

  double integral = 0.;
  while (Dist(xe, x0) > 1.e-6) {
    const double d = Dist(x1, x0);
    if (d < 1.e-6) {
      // Step length has become very small.
      if (m_debug) {
        std::cout << m_className << "::IntegrateDiffusion: Small step.\n";
      }
      integral += var0 * d;
      // Proceed with the next step.
      x0 = x1;
      x1 = xe;
      continue;
    }
    // Determine the variance at the end point of the step.
    var1 = GetVar(x1, particle);
    // Determine the variance at the mid point of the step.
    const Vec xm = MidPoint(x0, x1);
    const double varm = GetVar(xm, particle);
    if (var1 < 0. || varm < 0.) {
      std::cerr << m_className << "::IntegrateDiffusion:\n"
                << "    Cannot retrieve variance at mid or end point.\n";
      break;
    }
    // Compare first and second order estimates 
    // (integrals calculated using trapezoidal and Simpson's rule).
    if (fabs(var0 - 2 * varm + var1) * sqrt(d * 2 / (var0 + var1)) / 6. < tol) {
      // Accuracy is good enough.
      integral += d * (var0 + 4 * varm + var1) / 6.;
      // Proceed to the next step.
      x0 = x1;
      x1 = xe;
      var0 = var1;
    } else {
      // Accuracy is not good enough, so halve the step.
      x1 = xm;
      var1 = varm;
    }
  }
  return integral;
}

double DriftLineRKF::IntegrateAlpha(const std::array<double, 3>& xi, 
                                    const std::array<double, 3>& xe,
                                    const Particle particle, 
                                    const double tol) const {

  // Determine the Townsend coefficient at the initial point.
  Vec x0 = xi;
  double alpha0 = GetAlpha(x0, particle);
  // Determine the Townsend coefficient at the end point.
  Vec x1 = xe;
  double alpha1 = GetAlpha(x1, particle);
  if (alpha0 < 0. || alpha1 < 0.) {
    std::cerr << m_className << "::IntegrateAlpha:\n"
              << "    Cannot retrieve alpha at start point or end point.\n";
    return 0.;
  }
  double integral = 0.;
  while (Dist(xe, x0) > 1.e-6) {
    const double d = Dist(x1, x0);
    if (d < 1.e-6) {
      // Step length has become very small.
      if (m_debug) {
        std::cout << m_className << "::IntegrateAlpha: Small step.\n";
      }
      integral += alpha0 * d;
      // Proceed with the next step.
      x0 = x1;
      x1 = xe;
      continue;
    }
    // Calculate the Townsend coefficient at the end point of the step.
    alpha1 = GetAlpha(x1, particle);
    // Calculate the Townsend coefficient at the mid point of the step.
    const Vec xm = MidPoint(x0, x1);
    const double alpham = GetAlpha(xm, particle);
    if (alpha1 < 0. || alpham < 0.) {
      std::cerr << m_className << "::IntegrateAlpha:\n"
                << "    Cannot retrieve alpha at mid point or end point.\n";
      break;
    }
    // Compare first and second order estimates.
    if (d * fabs(alpha0 - 2 * alpham + alpha1) / 3. < tol) {
      // Accuracy is good enough.
      integral += d * (alpha0 + 4 * alpham + alpha1) / 6.;
      // Proceed to the next step.
      x0 = x1;
      x1 = xe;
      alpha0 = alpha1;
    } else {
      // Accuracy is not good enough, so halve the step.
      x1 = xm;
      alpha1 = alpham;
    }
  }
  return integral;
}

double DriftLineRKF::IntegrateEta(const std::array<double, 3>& xi, 
                                  const std::array<double, 3>& xe,
                                  const Particle particle, 
                                  const double tol) const {
  // Determine the attachment coefficient at the initial point.
  Vec x0 = xi;
  double eta0 = GetEta(x0, particle);
  // Determine the attachment coefficient at the end point.
  Vec x1 = xe;
  double eta1 = GetEta(x1, particle);
  if (eta0 < 0. || eta1 < 0.) {
    std::cerr << m_className << "::IntegrateEta:\n"
              << "    Cannot retrieve eta at start point or end point.\n";
    return 0.;
  }
  double integral = 0.;
  while (Dist(xe, x0) > 1.e-6) {
    const double d = Dist(x1, x0);
    if (d < 1.e-6) {
      // Step length has become very small.
      if (m_debug) {
        std::cout << m_className << "::IntegrateEta: Small step.\n";
      }
      integral += eta0 * d;
      // Proceed with the next step.
      x0 = x1;
      x1 = xe;
      continue;
    }
    // Calculate the attachment coefficient at the end point of the step.
    eta1 = GetEta(x1, particle);
    // Calculate the attachment coefficient at the mid point of the step.
    const Vec xm = MidPoint(x0, x1);
    const double etam = GetEta(xm, particle);
    if (eta1 < 0. || etam < 0.) {
      std::cerr << m_className << "::IntegrateEta:\n"
                << "    Cannot retrieve eta at mid point or end point.\n";
      break;
    }
    // Compare first and second order estimates.
    if (d * fabs(eta0 - 2 * etam + eta1) / 3. < tol) {
      // Accuracy is good enough.
      integral += d * (eta0 + 4 * etam + eta1) / 6.;
      // Proceed to the next step.
      x0 = x1;
      x1 = xe;
      eta0 = eta1;
    } else {
      // Accuracy is not good enough, so halve the step.
      x1 = xm;
      eta1 = etam;
    }
  }
  return integral;
}

void DriftLineRKF::ComputeSignal(const Particle particle, const double scale,
                                 const std::vector<double>& ts,
                                 const std::vector<Vec>& xs,
                                 const std::vector<double>& ne) const {

  const auto nPoints = ts.size();
  if (nPoints < 2) return;
  const double q0 = Charge(particle) * scale;
  
  if (m_useWeightingPotential) {
    const bool aval = ne.size() == nPoints; 
    for (size_t i = 0; i < nPoints - 1; ++i) {
      const auto& x0 = xs[i];
      const auto& x1 = xs[i + 1];
      const double q = aval ? q0 * 0.5 * (ne[i] + ne[i + 1]) : q0; 
      m_sensor->AddSignal(q, ts[i], ts[i + 1], 
                          x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], 
                          false, true);
    }
    return;
  }
  // Get the drift velocity at each point.
  std::vector<std::array<double, 3> > vs;
  for (const auto& x : xs) {
    int stat = 0;
    Vec v = GetVelocity(x, particle, stat);
    if (stat != 0) {
      std::cerr << m_className << "::ComputeSignal:\n"
                << "    Cannot retrieve velocity at " << PrintVec(x) << "\n";
    }
    vs.push_back(std::move(v));
  }
  m_sensor->AddSignal(q0, ts, xs, vs, ne, m_navg);
}

bool DriftLineRKF::FieldLine(const double xi, const double yi, const double zi,
                             std::vector<std::array<float, 3> >& xl, 
                             const bool electron) const {

  xl.clear();
  // Is the sensor?
  if (!m_sensor) {
    std::cerr << m_className << "::FieldLine: Sensor is not defined.\n";
    return false;
  }

  // Get the sensor's bounding box.
  double xmin = 0., xmax = 0.;
  double ymin = 0., ymax = 0.;
  double zmin = 0., zmax = 0.;
  bool bbox = m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax);

  // Make sure the initial position is at a valid location.
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  int stat = 0;
  m_sensor->ElectricField(xi, yi, zi, ex, ey, ez, medium, stat);
  if (!medium || stat != 0) {
    std::cerr << m_className << "::FieldLine: "
              << "No valid field at initial position.\n";
    return false;
  }
  Vec x0 = {xi, yi, zi};
  Vec f0 = {ex, ey, ez};
  if (electron) for (auto& f : f0) f *= -1; 

  // Set the numerical constants for the RKF integration.
  constexpr double c10 = 214. / 891.;
  constexpr double c11 = 1. / 33.;
  constexpr double c12 = 650. / 891.;
  constexpr double c20 = 533. / 2106.;
  constexpr double c22 = 800. / 1053.;
  constexpr double c23 = -1. / 78.;

  constexpr double b10 = 1. / 4.;
  constexpr double b20 = -189. / 800.;
  constexpr double b21 = 729. / 800.;
  constexpr double b30 = 214. / 891.;
  constexpr double b31 = 1. / 33.;
  constexpr double b32 = 650. / 891.;

  const double fmag0 = Mag(f0);
  if (fmag0 < Small) {
    std::cerr << m_className << "::FieldLine:\n"
              << "    Zero field at initial position.\n";
    return false;
  }

  // Initialise time step and previous time step.
  double h = m_accuracy / fmag0;
  double hprev = h;
  if (m_debug) {
    std::cout << m_className << "::FieldLine:\n"
              << "    Initial step size: " << h << ".\n";
  }
  // Set the initial point.
  xl.push_back({float(x0[0]), float(x0[1]), float(x0[2])});

  int initCycle = 3;
  bool ok = true;
  while (ok) {
    // Get the field at the first probe point.
    Vec x1 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x1[i] += h * b10 * f0[i];
    }
    m_sensor->ElectricField(x1[0], x1[1], x1[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      if (m_debug) std::cout << "    Point 1 outside.\n";
      Terminate(x0, x1, xl);
      return true;
    }
    Vec f1 = {ex, ey, ez};
    if (electron) for (auto& f : f1) f *= -1;
    // Get the field at the second probe point.
    Vec x2 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x2[i] += h * (b20 * f0[i] + b21 * f1[i]);
    }
    m_sensor->ElectricField(x2[0], x2[1], x2[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      if (m_debug) std::cout << "    Point 2 outside.\n";
      Terminate(x0, x2, xl);
      return true;
    }
    Vec f2 = {ex, ey, ez};
    if (electron) for (auto& f : f2) f *= -1;
    // Get the field at the third probe point.
    Vec x3 = x0;
    for (unsigned int i = 0; i < 3; ++i) {
      x3[i] += h * (b30 * f0[i] + b31 * f1[i] + b32 * f2[i]);
    }
    m_sensor->ElectricField(x3[0], x3[1], x3[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      if (m_debug) std::cout << "    Point 3 outside.\n";
      Terminate(x0, x3, xl);
      return true;
    }
    Vec f3 = {ex, ey, ez};
    if (electron) for (auto& f : f3) f *= -1;
    // Check if we crossed a wire.
    double xw = 0., yw = 0., zw = 0., rw = 0.;
    if (m_sensor->CrossedWire(x0[0], x0[1], x0[2], 
                              x1[0], x1[1], x1[2], xw, yw, zw, false, rw) ||
        m_sensor->CrossedWire(x0[0], x0[1], x0[2], 
                              x2[0], x2[1], x2[2], xw, yw, zw, false, rw) ||
        m_sensor->CrossedWire(x0[0], x0[1], x0[2], 
                              x3[0], x3[1], x3[2], xw, yw, zw, false, rw)) {
      // TODO!
      xl.push_back({float(xw), float(yw), float(zw)});
      return true;
      // Drift to wire.
      if (h > Small) {
        h *= 0.5;
        continue;
      } else {
        std::cerr << m_className << "::FieldLine: Step size too small. Stop.\n";
        return false;
      }
    }
    // Calculate the correction terms.
    Vec phi1 = {0., 0., 0.};
    Vec phi2 = {0., 0., 0.};
    for (unsigned int i = 0; i < 3; ++i) {
      phi1[i] = c10 * f0[i] + c11 * f1[i] + c12 * f2[i];
      phi2[i] = c20 * f0[i] + c22 * f2[i] + c23 * f3[i];
    }
    // Check if the step length is valid.
    const double phi1mag = Mag(phi1);
    if (phi1mag < Small) {
      std::cerr << m_className << "::FieldLine: Step has zero length. Stop.\n";
      break;
    } else if (m_useStepSizeLimit && h * phi1mag > m_maxStepSize) {
      if (m_debug) {
        std::cout << "    Step is considered too long. H is reduced.\n";
      }
      h = 0.5 * m_maxStepSize / phi1mag;
      continue;
    } else if (bbox) {
      // Don't allow h to become too large.
      if (h * fabs(phi1[0]) > 0.1 * fabs(xmax - xmin) ||
          h * fabs(phi1[1]) > 0.1 * fabs(ymax - ymin)) {
        h *= 0.5;
        if (m_debug) {
          std::cout << "    Step is considered too long. H is halved.\n";
        }
        continue;
      }
    }
    if (m_debug) std::cout << "    Step size ok.\n";
    // Update the position.
    for (size_t i = 0; i < 3; ++i) x0[i] += h * phi1[i];
    // Check the new position.
    m_sensor->ElectricField(x0[0], x0[1], x0[2], ex, ey, ez, medium, stat);
    if (stat != 0) {
      // The new position is not inside a valid drift medium.
      // Terminate the drift line.
      if (m_debug) std::cout << "    Point outside. Terminate.\n";
      std::array<double, 3> xp = {xl.back()[0], xl.back()[1], xl.back()[2]};
      Terminate(xp, x0, xl);
      return true;
    }
    // Add the new point to the drift line.
    xl.push_back({float(x0[0]), float(x0[1]), float(x0[2])});
    // Adjust the step size according to the accuracy of the two estimates.
    hprev = h;
    const double dphi = fabs(phi1[0] - phi2[0]) + fabs(phi1[1] - phi2[1]) +
                        fabs(phi1[2] - phi2[2]);
    if (dphi > 0) {
      h = sqrt(h * m_accuracy / dphi);
      if (m_debug) std::cout << "    Adapting H to " << h << ".\n";
    } else {
      h *= 2;
      if (m_debug) std::cout << "    H increased by factor two.\n";
    }
    // Make sure that H is different from zero; this should always be ok.
    if (h < Small) {
      std::cerr << m_className << "::FieldLine: Step size is zero. Stop.\n";
      return false;
    }
    // Check the initial step size.
    if (initCycle > 0 && h < 0.2 * hprev) {
      if (m_debug) std::cout << "    Reinitialise step size.\n";
      --initCycle;
      x0 = {xi, yi, zi};
      xl.clear();
      xl.push_back({float(xi), float(yi), float(zi)});
      continue;
    }
    initCycle = 0;
    // Don't allow H to grow too quickly
    if (h > 10 * hprev) {
      h = 10 * hprev;
      if (m_debug) {
        std::cout << "    H restricted to 10 times the previous value.\n";
      }
    }
    // Stop in case H tends to become too small.
    if (h * (fabs(phi1[0]) + fabs(phi1[1]) + fabs(phi1[2])) < m_accuracy) {
      std::cerr << m_className << "::FieldLine: Step size has become smaller "
                << "than int. accuracy. Stop.\n";
      return false;
    }
    // Update the field.
    f0 = f3;
  }
  return true;
}

void DriftLineRKF::Terminate(const std::array<double, 3>& xx0,
                             const std::array<double, 3>& xx1,
                             std::vector<std::array<float, 3> >& xs) const {

  // Final point just inside the medium.
  Vec x0 = xx0;
  // Final point just outside the medium.
  Vec x1 = xx1;
  // Perform some bisections.
  constexpr unsigned int nBisections = 20;
  for (unsigned int i = 0; i < nBisections; ++i) {
    // Quit bisection when interval becomes too small.
    bool small = true;
    for (unsigned int j = 0; j < 3; ++j) {
      if (fabs(x1[j] - x0[j]) > 1.e-6 * (fabs(x0[j]) + fabs(x1[j]))) {
        small = false;
        break;
      }
    } 
    if (small) {
      if (m_debug) {
        std::cout << m_className << "::Terminate: Bisection ends at cycle "
                  << i << ".\n";
      }
      break; 
    }
    // Calculate the mid point.
    const Vec xm = MidPoint(x0, x1);
    // Halve the step.
    if (m_sensor->IsInside(xm[0], xm[1], xm[2]) &&
        m_sensor->IsInArea(xm[0], xm[1], xm[2])) {
      x0 = xm;
    } else {
      x1 = xm;
    }
  }

  xs.push_back({float(x0[0]), float(x0[1]), float(x0[2])});
}
}
