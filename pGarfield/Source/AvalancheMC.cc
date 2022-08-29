#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

#include "Garfield/AvalancheMC.hh"
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

double Dist(const std::array<double, 3>& x0, 
            const std::array<double, 3>& x1) {
  std::array<double, 3> d = x1;
  for (size_t i = 0; i < 3; ++i) d[i] -= x0[i];
  return Mag(d); 
}

double Slope(const double x0, const double x1) {
  return (x0 > 0. && x1 > 0.) ? fabs(x1 - x0) / (x0 + x1) : 0.;
}

std::array<double, 3> MidPoint(const std::array<double, 3>& x0,
                               const std::array<double, 3>& x1) {
  std::array<double, 3> xm;
  for (size_t k = 0; k < 3; ++k) xm[k] = 0.5 * (x0[k] + x1[k]);
  return xm;
}

}  // namespace

namespace Garfield {

AvalancheMC::AvalancheMC() { m_drift.reserve(10000); }

void AvalancheMC::SetSensor(Sensor* sensor) {
  if (!sensor) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }

  m_sensor = sensor;
}

void AvalancheMC::EnablePlotting(ViewDrift* view) {
  if (!view) {
    std::cerr << m_className << "::EnablePlotting: Null pointer.\n";
    return;
  }

  m_viewer = view;
}

void AvalancheMC::SetTimeSteps(const double d) {
  m_stepModel = StepModel::FixedTime;
  if (d < Small) {
    std::cerr << m_className << "::SetTimeSteps:\n    "
              << "Step size is too small. Using default (20 ps) instead.\n";
    m_tMc = 0.02;
    return;
  }
  if (m_debug) {
    std::cout << m_className << "::SetTimeSteps:\n"
              << "    Step size set to " << d << " ns.\n";
  }
  m_tMc = d;
}

void AvalancheMC::SetDistanceSteps(const double d) {
  m_stepModel = StepModel::FixedDistance;
  if (d < Small) {
    std::cerr << m_className << "::SetDistanceSteps:\n    "
              << "Step size is too small. Using default (10 um) instead.\n";
    m_dMc = 0.001;
    return;
  }
  if (m_debug) {
    std::cout << m_className << "::SetDistanceSteps:\n"
              << "    Step size set to " << d << " cm.\n";
  }
  m_dMc = d;
}

void AvalancheMC::SetCollisionSteps(const unsigned int n) {
  m_stepModel = StepModel::CollisionTime;
  if (n < 1) {
    std::cerr << m_className << "::SetCollisionSteps:\n    "
              << "Number of collisions set to default value (100).\n";
    m_nMc = 100;
    return;
  }
  if (m_debug) {
    std::cout << m_className << "::SetCollisionSteps:\n    "
              << "Number of collisions to be skipped set to " << n << ".\n";
  }
  m_nMc = n;
}

void AvalancheMC::SetStepDistanceFunction(double (*f)(double x, double y,
                                                      double z)) {
  if (!f) {
    std::cerr << m_className << "::SetStepDistanceFunction: Null pointer.\n";
    return;
  }
  m_fStep = f;
  m_stepModel = StepModel::UserDistance;
}

void AvalancheMC::SetTimeWindow(const double t0, const double t1) {
  if (fabs(t1 - t0) < Small) {
    std::cerr << m_className << "::SetTimeWindow:\n"
              << "    Time interval must be greater than zero.\n";
    return;
  }

  m_tMin = std::min(t0, t1);
  m_tMax = std::max(t0, t1);
  m_hasTimeWindow = true;
}

void AvalancheMC::GetDriftLinePoint(const size_t i, double& x, double& y,
                                    double& z, double& t) const {
  if (i >= m_drift.size()) {
    std::cerr << m_className << "::GetDriftLinePoint: Index out of range.\n";
    return;
  }

  x = m_drift[i].x[0];
  y = m_drift[i].x[1];
  z = m_drift[i].x[2];
  t = m_drift[i].t;
}

void AvalancheMC::GetHoleEndpoint(const size_t i, double& x0, double& y0,
                                  double& z0, double& t0, double& x1,
                                  double& y1, double& z1, double& t1,
                                  int& status) const {
  if (i >= m_endpointsHoles.size()) {
    std::cerr << m_className << "::GetHoleEndpoint: Index out of range.\n";
    return;
  }

  x0 = m_endpointsHoles[i].x0[0];
  y0 = m_endpointsHoles[i].x0[1];
  z0 = m_endpointsHoles[i].x0[2];
  t0 = m_endpointsHoles[i].t0;
  x1 = m_endpointsHoles[i].x1[0];
  y1 = m_endpointsHoles[i].x1[1];
  z1 = m_endpointsHoles[i].x1[2];
  t1 = m_endpointsHoles[i].t1;
  status = m_endpointsHoles[i].status;
}

void AvalancheMC::GetIonEndpoint(const size_t i, double& x0, double& y0,
                                 double& z0, double& t0, double& x1, double& y1,
                                 double& z1, double& t1, int& status) const {
  if (i >= m_endpointsIons.size()) {
    std::cerr << m_className << "::GetIonEndpoint: Index out of range.\n";
    return;
  }

  x0 = m_endpointsIons[i].x0[0];
  y0 = m_endpointsIons[i].x0[1];
  z0 = m_endpointsIons[i].x0[2];
  t0 = m_endpointsIons[i].t0;
  x1 = m_endpointsIons[i].x1[0];
  y1 = m_endpointsIons[i].x1[1];
  z1 = m_endpointsIons[i].x1[2];
  t1 = m_endpointsIons[i].t1;
  status = m_endpointsIons[i].status;
}

void AvalancheMC::GetElectronEndpoint(const size_t i, double& x0,
                                      double& y0, double& z0, double& t0,
                                      double& x1, double& y1, double& z1,
                                      double& t1, int& status) const {
  if (i >= m_endpointsElectrons.size()) {
    std::cerr << m_className << "::GetElectronEndpoint: Index out of range.\n";
    return;
  }

  x0 = m_endpointsElectrons[i].x0[0];
  y0 = m_endpointsElectrons[i].x0[1];
  z0 = m_endpointsElectrons[i].x0[2];
  t0 = m_endpointsElectrons[i].t0;
  x1 = m_endpointsElectrons[i].x1[0];
  y1 = m_endpointsElectrons[i].x1[1];
  z1 = m_endpointsElectrons[i].x1[2];
  t1 = m_endpointsElectrons[i].t1;
  status = m_endpointsElectrons[i].status;
}

bool AvalancheMC::DriftElectron(const double x0, const double y0,
                                const double z0, const double t0) {
  if (!m_sensor) {
    std::cerr << m_className << "::DriftElectron: Sensor is not defined.\n";
    return false;
  }

  m_endpointsElectrons.clear();
  m_endpointsHoles.clear();
  m_endpointsIons.clear();

  m_nElectrons = 1;
  m_nHoles = 0;
  m_nIons = 0;

  std::vector<DriftPoint> secondaries;
  return DriftLine({x0, y0, z0}, t0, Particle::Electron, secondaries);
}

bool AvalancheMC::DriftHole(const double x0, const double y0, const double z0,
                            const double t0) {
  if (!m_sensor) {
    std::cerr << m_className << "::DriftHole: Sensor is not defined.\n";
    return false;
  }

  m_endpointsElectrons.clear();
  m_endpointsHoles.clear();
  m_endpointsIons.clear();

  m_nElectrons = 0;
  m_nHoles = 1;
  m_nIons = 0;

  std::vector<DriftPoint> secondaries;
  return DriftLine({x0, y0, z0}, t0, Particle::Hole, secondaries);
}

bool AvalancheMC::DriftIon(const double x0, const double y0, const double z0,
                           const double t0) {
  if (!m_sensor) {
    std::cerr << m_className << "::DriftIon: Sensor is not defined.\n";
    return false;
  }

  m_endpointsElectrons.clear();
  m_endpointsHoles.clear();
  m_endpointsIons.clear();

  m_nElectrons = 0;
  m_nHoles = 0;
  m_nIons = 1;

  std::vector<DriftPoint> secondaries;
  return DriftLine({x0, y0, z0}, t0, Particle::Ion, secondaries);
}

bool AvalancheMC::DriftNegativeIon(const double x0, const double y0, 
                                   const double z0, const double t0) {
  if (!m_sensor) {
    std::cerr << m_className << "::DriftNegativeIon: Sensor is not defined.\n";
    return false;
  }

  m_endpointsElectrons.clear();
  m_endpointsHoles.clear();
  m_endpointsIons.clear();

  m_nElectrons = 0;
  m_nHoles = 0;
  m_nIons = 1;

  std::vector<DriftPoint> secondaries;
  return DriftLine({x0, y0, z0}, t0, Particle::NegativeIon, secondaries);
}

bool AvalancheMC::DriftLine(const std::array<double, 3>& xi, const double ti,
                            const Particle particle,
                            std::vector<DriftPoint>& secondaries,
                            const bool aval) {
  // Reset the drift line.
  m_drift.clear();
  // Check the initial position.
  std::array<double, 3> x0 = xi;
  std::array<double, 3> e0 = {0., 0., 0.};
  std::array<double, 3> b0 = {0., 0., 0.};
  Medium* medium = nullptr;
  int status = GetField(x0, e0, b0, medium);
  if (status != 0) {
    std::cerr << m_className + "::DriftLine: "
              << PrintVec(x0) + " is not in a valid drift region.\n";
  }
  // Check the initial time.
  double t0 = ti;
  if (m_hasTimeWindow && (t0 < m_tMin || t0 > m_tMax)) {
    status = StatusOutsideTimeWindow;
    std::cerr << m_className + "::DriftLine: " << t0
              << " is outside the time window.\n";
  }
  // Stop here if initial position or time are invalid.
  if (status != 0) return false;
  // Add the first point to the line.
  AddPoint(x0, t0, particle, 1, m_drift);
  if (m_debug) {
    std::cout << m_className + "::DriftLine: Starting at "
              << PrintVec(x0) + ".\n";
  }
  // Determine if the medium is a gas or semiconductor.
  const bool semiconductor = medium->IsSemiconductor();

  while (0 == status) {
    constexpr double tol = 1.e-10;
    // Make sure the electric field has a non-vanishing component.
    const double emag = Mag(e0);
    if (emag < tol && !m_useDiffusion) {
      std::cerr << m_className + "::DriftLine: Too small electric field at "
                << PrintVec(x0) + ".\n";
      status = StatusCalculationAbandoned;
      break;
    }
    // Compute the drift velocity at this point.
    std::array<double, 3> v0;
    if (!GetVelocity(particle, medium, x0, e0, b0, v0)) {
      status = StatusCalculationAbandoned;
      break;
    }

    // Make sure the drift velocity vector has a non-vanishing component.
    double vmag = Mag(v0);
    if (vmag < tol && !m_useDiffusion) {
      std::cerr << m_className + "::DriftLine: Too small drift velocity at "
                << PrintVec(x0) + ".\n";
      status = StatusCalculationAbandoned;
      break;
    }

    // Coordinates after the step.
    std::array<double, 3> x1 = x0;
    // Time after the step.
    double t1 = t0;
    if (vmag < tol || emag < tol) {
      // Diffusion only. Get the mobility.
      const double mu = GetMobility(particle, medium);
      if (mu < 0.) {
        std::cerr << m_className + "::DriftLine: Invalid mobility.\n";
        status = StatusCalculationAbandoned;
        break;
      }
      // Calculate the diffusion coefficient.
      const double dif = mu * BoltzmannConstant * medium->GetTemperature();
      double sigma = 0.;
      switch (m_stepModel) {
        case StepModel::FixedTime:
          sigma = sqrt(2. * dif * m_tMc);
          t1 += m_tMc;
          break;
        case StepModel::FixedDistance:
          sigma = m_dMc;
          break;
        case StepModel::CollisionTime: {
          // Thermal velocity.
          const double vth =
              SpeedOfLight * sqrt(2 * BoltzmannConstant *
                                  medium->GetTemperature() / ElectronMass);
          sigma = m_nMc * dif / vth;
        } break;
        case StepModel::UserDistance:
          sigma = m_fStep(x0[0], x0[1], x0[2]);
          break;
        default:
          std::cerr << m_className + "::DriftLine: Unknown stepping model.\n";
          status = StatusCalculationAbandoned;
          break;
      }
      if (status != 0) break;
      if (m_stepModel != StepModel::FixedTime) {
        t1 += sigma * sigma / (2 * dif);
      }
      for (size_t i = 0; i < 3; ++i) x1[i] += RndmGaussian(0., sigma);
    } else {
      // Drift and diffusion. Determine the time step.
      double dt = 0.;
      switch (m_stepModel) {
        case StepModel::FixedTime:
          dt = m_tMc;
          break;
        case StepModel::FixedDistance:
          dt = m_dMc / vmag;
          break;
        case StepModel::CollisionTime:
          if (particle == Particle::Ion || particle == Particle::NegativeIon) {
            constexpr double c1 = AtomicMassUnitElectronVolt / 
              (SpeedOfLight * SpeedOfLight);
            dt = -m_nMc * (c1 * vmag / emag) * log(RndmUniformPos());
          } else { 
            constexpr double c1 = ElectronMass / (SpeedOfLight * SpeedOfLight);
            dt = -m_nMc * (c1 * vmag / emag) * log(RndmUniformPos());
          }
          break;
        case StepModel::UserDistance:
          dt = m_fStep(x0[0], x0[1], x0[2]) / vmag;
          break;
        default:
          std::cerr << m_className + "::DriftLine: Unknown stepping model.\n";
          status = StatusCalculationAbandoned;
          break;
      }
      if (status != 0) break;

      double difl = 0., dift = 0.;
      if (m_useDiffusion) {
        // Get the diffusion coefficients.
        if (!GetDiffusion(particle, medium, e0, b0, difl, dift)) {
          PrintError("DriftLine", "diffusion", particle, x0);
          status = StatusCalculationAbandoned;
          break;
        }
        if (m_stepModel != StepModel::FixedTime) {
          const double ds = vmag * dt;
          const double dif = std::max(difl, dift);
          if (dif * dif > ds) {
            dt = ds * ds / (dif * dif * vmag);
          }
        }
      }
      // Compute the proposed end point of this step.
      for (size_t k = 0; k < 3; ++k) x1[k] += dt * v0[k];
      std::array<double, 3> v1 = v0;
      constexpr unsigned int nMaxIter = 3;
      for (unsigned int i = 0; i < nMaxIter; ++i) {
        status = GetField(x1, e0, b0, medium);
        if (status != 0) {
          // Point is outside the active region. Reduce the step size.
          x1 = MidPoint(x0, x1);
          dt *= 0.5;
          continue;
        }
        // Compute the velocity at the proposed end point.
        if (!GetVelocity(particle, medium, x1, e0, b0, v1)) {
          status = StatusCalculationAbandoned;
          break;
        }
        if (Slope(vmag, Mag(v1)) < 0.05) break;
        // Halve the step.
        x1 = MidPoint(x0, x1);
        dt *= 0.5;
      }
      if (status == StatusCalculationAbandoned) break;
      if (m_doRKF) {
        StepRKF(particle, x0, v0, dt, x1, v1, status);
        vmag = Mag(v1);
      }
      if (m_useDiffusion) AddDiffusion(sqrt(vmag * dt), difl, dift, x1, v1);
      t1 += dt;
    }
    if (m_debug) {
      std::cout << m_className + "::DriftLine: Next point: "
                << PrintVec(x1) + ".\n";
    }

    // Get the electric and magnetic field at the new position.
    status = GetField(x1, e0, b0, medium);
    if (status == StatusLeftDriftMedium || status == StatusLeftDriftArea) {
      // Point is not inside a "driftable" medium or outside the drift area.
      // Try terminating the drift line close to the boundary.
      Terminate(x0, t0, x1, t1);
      if (m_debug) {
        std::cout << m_className + "::DriftLine: Left drift region at "
                  << PrintVec(x1) + ".\n";
      }
      // Add the point to the drift line.
      AddPoint(x1, t1, particle, 1, m_drift);
      break;
    }
    // Check if the particle has crossed a wire.
    std::array<double, 3> xc = x0;
    double rc = 0.;
    if (m_sensor->CrossedWire(x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], 
                              xc[0], xc[1], xc[2], false, rc)) {
      if (m_debug) {
        std::cout << m_className + "::DriftLine: Hit a wire at "
                  << PrintVec(xc) + ".\n";
      }
      status = StatusLeftDriftMedium;
      // Adjust the time step.
      const double tc = t0 + (t1 - t0) * Dist(x0, xc) / Dist(x0, x1);
      // Add the point to the drift line.
      AddPoint(xc, tc, particle, 1, m_drift);
      break;
    }
    if (m_sensor->CrossedPlane(x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], 
                               xc[0], xc[1], xc[2])) {
      if (m_debug) {
        std::cout << m_className + "::DriftLine: Hit a plane at "
                  << PrintVec(xc) + ".\n";
      }
      status = StatusHitPlane;
      // Adjust the time step.
      const double tc = t0 + (t1 - t0) * Dist(x0, xc) / Dist(x0, x1);
      // Add the point to the drift line.
      AddPoint(xc, tc, particle, 1, m_drift);
      break;
    }

    // Make sure the time is still within the specified interval.
    if (m_hasTimeWindow && (t1 < m_tMin || t1 > m_tMax)) {
      status = StatusOutsideTimeWindow;
    }
    // Add the point to the drift line.
    AddPoint(x1, t1, particle, 1, m_drift);
    // Update the current position and time.
    x0 = x1;
    t0 = t1;
  }
 
  if (status == StatusCalculationAbandoned) {
    std::cerr << m_className + "::DriftLine: Abandoned the calculation.\n";
  }

  // Compute Townsend and attachment coefficients for each drift step.
  unsigned int nElectronsOld = m_nElectrons;
  unsigned int nHolesOld = m_nHoles;
  unsigned int nIonsOld = m_nIons;

  if ((particle == Particle::Electron || particle == Particle::Hole) &&
      (aval || m_useAttachment) &&
      (m_sizeCut == 0 || m_nElectrons < m_sizeCut)) {
    ComputeGainLoss(particle, m_drift, status, secondaries, semiconductor);
    if (status == StatusAttached && m_debug) {
      std::cout << m_className + "::DriftLine: Attached at "
                << PrintVec(m_drift.back().x) + ".\n";
    }
  }

  if (m_debug) {
    std::cout << m_className << "::DriftLine: Stopped at "
              << PrintVec(m_drift.back().x) + ".\n";
  }
  // Create an "endpoint".
  AddEndPoint(xi, ti, m_drift.back().x, m_drift.back().t,
              status, particle);

  if (m_debug) {
    const int nNewElectrons = m_nElectrons - nElectronsOld;
    const int nNewHoles = m_nHoles - nHolesOld;
    const int nNewIons = m_nIons - nIonsOld;
    std::cout << m_className << "::DriftLine: Produced\n"
              << "      " << nNewElectrons << " electrons,\n"
              << "      " << nNewHoles << " holes, and\n"
              << "      " << nNewIons << " ions.\n";
  }

  // Compute the induced signal and induced charge if requested.
  const double scale = particle == Particle::Electron
                           ? -m_scaleE
                           : particle == Particle::Hole ? m_scaleH : m_scaleI;
  if (m_doSignal) ComputeSignal(particle, scale, m_drift);
  if (m_doInducedCharge) ComputeInducedCharge(scale, m_drift);

  // Plot the drift line if requested.
  if (m_viewer && !m_drift.empty()) {
    const size_t nPoints = m_drift.size();
    // Register the new drift line and get its ID.
    size_t id;
    m_viewer->NewDriftLine(particle, nPoints, id, xi[0], xi[1], xi[2]);
    // Set the points along the trajectory.
    for (size_t i = 0; i < nPoints; ++i) {
      const auto& x = m_drift[i].x;
      m_viewer->SetDriftLinePoint(id, i, x[0], x[1], x[2]);
    }
  }

  if (status == StatusCalculationAbandoned) return false;
  return true;
}

bool AvalancheMC::AvalancheElectron(const double x0, const double y0,
                                    const double z0, const double t0,
                                    const bool holes) {
  std::vector<DriftPoint> aval;
  AddPoint({x0, y0, z0}, t0, Particle::Electron, 1, aval);
  return Avalanche(aval, true, holes);
}

bool AvalancheMC::AvalancheHole(const double x0, const double y0,
                                const double z0, const double t0,
                                const bool electrons) {
  std::vector<DriftPoint> aval;
  AddPoint({x0, y0, z0}, t0, Particle::Hole, 1, aval);
  return Avalanche(aval, electrons, true);
}

bool AvalancheMC::AvalancheElectronHole(const double x0, const double y0,
                                        const double z0, const double t0) {
  std::vector<DriftPoint> aval;
  AddPoint({x0, y0, z0}, t0, Particle::Electron, 1, aval);
  AddPoint({x0, y0, z0}, t0, Particle::Hole, 1, aval);
  return Avalanche(aval, true, true);
}

void AvalancheMC::AddElectron(const double x, const double y, const double z,
                              const double t) {
  AddEndPoint({x, y, z}, t, {x, y, z}, t, StatusAlive, Particle::Electron);
  ++m_nElectrons;
}

void AvalancheMC::AddHole(const double x, const double y, const double z,
                          const double t) {
  AddEndPoint({x, y, z}, t, {x, y, z}, t, StatusAlive, Particle::Hole);
  ++m_nHoles;
}
void AvalancheMC::AddIon(const double x, const double y, const double z,
                         const double t) {
  AddEndPoint({x, y, z}, t, {x, y, z}, t, StatusAlive, Particle::Ion);
  ++m_nIons;
}

bool AvalancheMC::ResumeAvalanche(const bool electrons, const bool holes) {

  std::vector<DriftPoint> aval;
  for (const auto& p: m_endpointsElectrons) {
    if (p.status == StatusAlive || p.status == StatusOutsideTimeWindow) {
      AddPoint(p.x1, p.t1, Particle::Electron, 1, aval);
    } 
  }
  for (const auto& p: m_endpointsHoles) {
    if (p.status == StatusAlive || p.status == StatusOutsideTimeWindow) {
      AddPoint(p.x1, p.t1, Particle::Hole, 1, aval);
    } 
  }
  for (const auto& p: m_endpointsIons) {
    if (p.status == StatusAlive || p.status == StatusOutsideTimeWindow) {
      AddPoint(p.x1, p.t1, Particle::Ion, 1, aval);
    } 
  }
  return Avalanche(aval, electrons, holes);
}

bool AvalancheMC::Avalanche(std::vector<DriftPoint>& aval, 
                            const bool withE, const bool withH) {
  // -----------------------------------------------------------------------
  //   DLCMCA - Subroutine that computes a drift line using a Monte-Carlo
  //            technique to take account of diffusion and of avalanche
  //            formation.
  // -----------------------------------------------------------------------

  m_endpointsElectrons.clear();
  m_endpointsHoles.clear();
  m_endpointsIons.clear();

  // Make sure the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::Avalanche: Sensor is not defined.\n";
    return false;
  }

  m_nElectrons = 0;
  m_nHoles = 0;
  m_nIons = 0;
  for (const auto& point : aval) {
    if (point.particle == Particle::Electron) {
      ++m_nElectrons;
    } else if (point.particle == Particle::Hole) {
      ++m_nHoles;
    } else {
      ++m_nIons;
    }
  }

  if (!withH && !withE) {
    std::cerr << m_className + "::Avalanche: "
              << "Neither electron nor hole/ion component requested.\n";
  }

  std::vector<DriftPoint> secondaries;
  while (!aval.empty()) {
    for (const auto& point : aval) {
      if (!withE && point.particle == Particle::Electron) continue;
      if (!withH && point.particle != Particle::Electron) continue; 
      // Skip points outside the time window.
      if (m_hasTimeWindow && (point.t < m_tMin || point.t > m_tMax)) {
        for (unsigned int i = 0; i < point.n; ++i) {
          AddEndPoint(point.x, point.t, point.x, point.t,
                      StatusOutsideTimeWindow, point.particle);
        }
        continue;
      }
      for (unsigned int i = 0; i < point.n; ++i) {
        // Compute a drift line.
        DriftLine(point.x, point.t, point.particle, secondaries, true);
      }
    }
    aval.swap(secondaries);
    secondaries.clear();
  }
  return true;
}

int AvalancheMC::GetField(const std::array<double, 3>& x,
                          std::array<double, 3>& e, std::array<double, 3>& b,
                          Medium*& medium) const {
  e.fill(0.);
  b.fill(0.);
  int status = 0;
  // Get the magnetic field.
  m_sensor->MagneticField(x[0], x[1], x[2], b[0], b[1], b[2], status);
  // Get the electric field.
  m_sensor->ElectricField(x[0], x[1], x[2], e[0], e[1], e[2], medium, status);
  // Make sure the point is inside a drift medium.
  if (status != 0 || !medium) return StatusLeftDriftMedium;
  // Make sure the point is inside the drift area.
  if (!m_sensor->IsInArea(x[0], x[1], x[2])) return StatusLeftDriftArea;

  return 0;
}

double AvalancheMC::GetMobility(const Particle particle, Medium* medium) const {
  if (particle == Particle::Electron) {
    return medium->ElectronMobility();
  } else if (particle == Particle::Hole) {
    return medium->HoleMobility();
  } else if (particle == Particle::Ion) {
    return medium->IonMobility();
  } else if (particle == Particle::NegativeIon) {
    return medium->NegativeIonMobility();
  }
  return -1.;
}

bool AvalancheMC::GetVelocity(const Particle particle, Medium* medium,
                              const std::array<double, 3>& x,
                              const std::array<double, 3>& e,
                              const std::array<double, 3>& b,
                              std::array<double, 3>& v) const {
  v.fill(0.);
  bool ok = false;
  if (m_useVelocityMap && 
      particle != Particle::Ion && particle != Particle::NegativeIon) {
    // We assume there is only one component with a velocity map.
    const auto nComponents = m_sensor->GetNumberOfComponents();
    for (size_t i = 0; i < nComponents; ++i) {
      auto cmp = m_sensor->GetComponent(i);
      if (!cmp->HasVelocityMap()) continue;
      if (particle == Particle::Electron) {
        ok = cmp->ElectronVelocity(x[0], x[1], x[2], v[0], v[1], v[2]);
      } else if (particle == Particle::Hole) {
        ok = cmp->HoleVelocity(x[0], x[1], x[2], v[0], v[1], v[2]);
      }
      if (!ok) continue;
      // Seems to have worked.
      if (m_debug) {
        std::cout << m_className << "::GetVelocity: Velocity at "
                  << PrintVec(x) << " = " << PrintVec(v) << "\n";
      }
      return true;
    }
  }
  if (particle == Particle::Electron) {
    ok = medium->ElectronVelocity(e[0], e[1], e[2], b[0], b[1], b[2], 
                                  v[0], v[1], v[2]);
  } else if (particle == Particle::Hole) {
    ok = medium->HoleVelocity(e[0], e[1], e[2], b[0], b[1], b[2], 
                              v[0], v[1], v[2]);
  } else if (particle == Particle::Ion) {
    ok = medium->IonVelocity(e[0], e[1], e[2], b[0], b[1], b[2], 
                             v[0], v[1], v[2]);
  } else if (particle == Particle::NegativeIon) {
    ok = medium->NegativeIonVelocity(e[0], e[1], e[2], b[0], b[1], b[2], 
                                     v[0], v[1], v[2]);
  }
  if (!ok) {
    PrintError("GetVelocity", "velocity", particle, x);
    return false;
  }
  if (m_debug) {
    std::cout << m_className << "::GetVelocity: Velocity at " << PrintVec(x)
              << " = " << PrintVec(v) << "\n";
  }
  return true;
}

bool AvalancheMC::GetDiffusion(const Particle particle, Medium* medium,
                               const std::array<double, 3>& e,
                               const std::array<double, 3>& b, double& dl,
                               double& dt) const {
  dl = 0.;
  dt = 0.;
  bool ok = false;
  if (particle == Particle::Electron) {
    ok = medium->ElectronDiffusion(e[0], e[1], e[2], b[0], b[1], b[2], dl, dt);
  } else if (particle == Particle::Hole) {
    ok = medium->HoleDiffusion(e[0], e[1], e[2], b[0], b[1], b[2], dl, dt);
  } else if (particle == Particle::Ion || particle == Particle::NegativeIon) {
    ok = medium->IonDiffusion(e[0], e[1], e[2], b[0], b[1], b[2], dl, dt);
  }
  return ok;
}

double AvalancheMC::GetAttachment(const Particle particle, Medium* medium,
                                  const std::array<double, 3>& x,
                                  const std::array<double, 3>& e,
                                  const std::array<double, 3>& b) const {
  double eta = 0.;
  if (m_useAttachmentMap) {
    const auto nComponents = m_sensor->GetNumberOfComponents();
    for (size_t i = 0; i < nComponents; ++i) {
      auto cmp = m_sensor->GetComponent(i);
      if (!cmp->HasAttachmentMap()) continue;
      if (particle == Particle::Electron) {
        if (!cmp->ElectronAttachment(x[0], x[1], x[2], eta)) continue;
      } else {
        if (!cmp->HoleAttachment(x[0], x[1], x[2], eta)) continue;
      }
      return eta;
    }
  }
  if (particle == Particle::Electron) {
    medium->ElectronAttachment(e[0], e[1], e[2], b[0], b[1], b[2], eta);
  } else {
    medium->HoleAttachment(e[0], e[1], e[2], b[0], b[1], b[2], eta);
  }
  return eta;
}

double AvalancheMC::GetTownsend(const Particle particle, Medium* medium,
                                const std::array<double, 3>& x,
                                const std::array<double, 3>& e,
                                const std::array<double, 3>& b) const {
  double alpha = 0.;
  if (m_useTownsendMap) {
    const auto nComponents = m_sensor->GetNumberOfComponents();
    for (size_t i = 0; i < nComponents; ++i) {
      auto cmp = m_sensor->GetComponent(i);
      if (!cmp->HasTownsendMap()) continue;
      if (particle == Particle::Electron) {
        if (!cmp->ElectronTownsend(x[0], x[1], x[2], alpha)) continue;
      } else {
        if (!cmp->HoleTownsend(x[0], x[1], x[2], alpha)) continue;
      }
      return alpha;
    }
  }
  if (particle == Particle::Electron) {
    medium->ElectronTownsend(e[0], e[1], e[2], b[0], b[1], b[2], alpha);
  } else {
    medium->HoleTownsend(e[0], e[1], e[2], b[0], b[1], b[2], alpha);
  }
  return alpha;
}

void AvalancheMC::StepRKF(const Particle particle,
                          const std::array<double, 3>& x0,
                          const std::array<double, 3>& v0, const double dt,
                          std::array<double, 3>& xf, std::array<double, 3>& vf,
                          int& status) const {
  // Constants appearing in the RKF formulas.
  constexpr double ci0 = 214. / 891.;
  constexpr double ci1 = 1. / 33.;
  constexpr double ci2 = 650. / 891.;
  constexpr double beta10 = 1. / 4.;
  constexpr double beta20 = -189. / 800.;
  constexpr double beta21 = 729. / 800.;

  vf = v0;
  // First probe point.
  for (size_t k = 0; k < 3; ++k) {
    xf[k] = x0[k] + dt * beta10 * v0[k];
  }
  std::array<double, 3> e;
  std::array<double, 3> b;
  Medium* medium = nullptr;
  status = GetField(xf, e, b, medium);
  if (status != 0) return;

  // Get the velocity at the first point.
  std::array<double, 3> v1;
  if (!GetVelocity(particle, medium, xf, e, b, v1)) {
    status = StatusCalculationAbandoned;
    return;
  }

  // Second point.
  for (size_t k = 0; k < 3; ++k) {
    xf[k] = x0[k] + dt * (beta20 * v0[k] + beta21 * v1[k]);
  }
  status = GetField(xf, e, b, medium);
  if (status != 0) return;

  // Get the velocity at the second point.
  std::array<double, 3> v2;
  if (!GetVelocity(particle, medium, xf, e, b, v2)) {
    status = StatusCalculationAbandoned;
    return;
  }

  // Compute the mean velocity and endpoint of the step.
  for (size_t k = 0; k < 3; ++k) {
    vf[k] = ci0 * v0[k] + ci1 * v1[k] + ci2 * v2[k];
    xf[k] = x0[k] + dt * vf[k];
  }
}

void AvalancheMC::AddDiffusion(const double step, const double dl,
                               const double dt, std::array<double, 3>& x,
                               const std::array<double, 3>& v) const {
  // Draw a random diffusion direction in the particle frame.
  const std::array<double, 3> d = {step * RndmGaussian(0., dl),
                                   step * RndmGaussian(0., dt),
                                   step * RndmGaussian(0., dt)};
  if (m_debug) {
    std::cout << m_className << "::AddDiffusion: Adding diffusion step "
              << PrintVec(d) << "\n";
  }
  // Compute the rotation angles to align diffusion and drift velocity vectors.
  const double vt = sqrt(v[0] * v[0] + v[1] * v[1]);
  const double phi = vt > Small ? atan2(v[1], v[0]) : 0.;
  const double theta =
      vt > Small ? atan2(v[2], vt) : v[2] < 0. ? -HalfPi : HalfPi;
  const double cphi = cos(phi);
  const double sphi = sin(phi);
  const double ctheta = cos(theta);
  const double stheta = sin(theta);

  x[0] += cphi * ctheta * d[0] - sphi * d[1] - cphi * stheta * d[2];
  x[1] += sphi * ctheta * d[0] + cphi * d[1] - sphi * stheta * d[2];
  x[2] += stheta * d[0] + ctheta * d[2];
}

void AvalancheMC::Terminate(const std::array<double, 3>& x0, const double t0,
                            std::array<double, 3>& x1, double& t1) const {
  double dt = t1 - t0;
  // Calculate the normalised direction vector.
  std::array<double, 3> dx = {x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]};
  double ds = Mag(dx);
  if (ds > 0.) {
    const double scale = 1. / ds;
    for (unsigned int k = 0; k < 3; ++k) dx[k] *= scale;
  }
  x1 = x0;
  t1 = t0;
  while (ds > BoundaryDistance) {
    dt *= 0.5;
    ds *= 0.5;
    std::array<double, 3> xm = x1;
    for (unsigned int k = 0; k < 3; ++k) xm[k] += dx[k] * ds;
    // Check if the mid-point is inside the drift medium and the drift area.
    double ex = 0., ey = 0., ez = 0.;
    int status = 0;
    Medium* medium = nullptr;
    m_sensor->ElectricField(xm[0], xm[1], xm[2], ex, ey, ez, medium, status);
    if (status == 0 && m_sensor->IsInArea(xm[0], xm[1], xm[2])) {
      x1 = xm;
      t1 += dt;
    }
  }
}

bool AvalancheMC::ComputeGainLoss(const Particle particle,
                                  std::vector<DriftPoint>& driftLine,
                                  int& status, std::vector<DriftPoint>& secondaries, 
                                  const bool semiconductor) {

  std::vector<double> alps;
  std::vector<double> etas;
  // Compute the integrated Townsend and attachment coefficients.
  if (!ComputeAlphaEta(particle, driftLine, alps, etas)) return false;

  // Opposite-charge particle produced in the avalanche.
  Particle other = Particle::Electron;
  if (particle == Particle::Electron) {
    other = semiconductor ? Particle::Hole : Particle::Ion;
  } 
  const size_t nPoints = driftLine.size();
  // Loop over the drift line.
  for (size_t i = 0; i < nPoints - 1; ++i) {
    // Start with the initial electron (or hole).
    int ne = 1;
    int ni = 0;
    if (etas[i] < Small) {
      ne = RndmYuleFurry(std::exp(alps[i]));
      ni = ne - 1;
    } else {
      // Subdivision of a step.
      constexpr double probth = 0.01;
      // Compute the number of subdivisions.
      const int nDiv = std::max(int((alps[i] + etas[i]) / probth), 1);
      // Compute the probabilities for gain and loss.
      const double p = std::max(alps[i] / nDiv, 0.);
      const double q = std::max(etas[i] / nDiv, 0.);
      // Loop over the subdivisions.
      for (int j = 0; j < nDiv; ++j) {
        if (ne > 100) {
          // Gaussian approximation.
          const int gain = int(ne * p + RndmGaussian() * sqrt(ne * p * (1. - p)));
          const int loss = int(ne * q + RndmGaussian() * sqrt(ne * q * (1. - q)));
          ne += gain - loss;
          ni += gain;
        } else {
          // Binomial approximation
          for (int k = ne; k--;) {
            if (RndmUniform() < p) {
              ++ne;
              ++ni;
            }
            if (RndmUniform() < q) --ne;
          }
        }
        // Check if the electron (or hole) has survived.
        if (ne <= 0) {
          status = StatusAttached;
          if (particle == Particle::Electron) {
            --m_nElectrons;
          } else if (particle == Particle::Hole) {
            --m_nHoles;
          } else {
            --m_nIons;
          }
          const double f0 = (j + 0.5) / nDiv;
          const double f1 = 1. - f0;
          const auto x0 = driftLine[i].x;
          const auto x1 = driftLine[i + 1].x;
          driftLine.resize(i + 2);
          for (size_t k = 0; k < 3; ++k) {
            driftLine[i + 1].x[k] = f0 * x0[k] + f1 * x1[k];
          }
          driftLine[i + 1].t = f0 * driftLine[i].t + f1 * driftLine[i + 1].t;
          break;
        }
      }
    }
    // Add the new electrons to the table.
    if (ne > 1) {
      DriftPoint point = driftLine[i + 1];
      point.n = ne - 1; 
      secondaries.push_back(std::move(point));
      if (particle == Particle::Electron) {
        m_nElectrons += ne - 1;
      } else if (particle == Particle::Hole) {
        m_nHoles += ne - 1;
      }
    }
    // Add the new holes/ions to the table.
    if (ni > 0) {
      if (other == Particle::Hole) {
        m_nHoles += ni;
      } else if (other == Particle::Ion) {
        m_nIons += ni;
      } else {
        m_nElectrons += ni;
      } 
      const auto x0 = driftLine[i].x;
      const auto x1 = driftLine[i + 1].x;
      const double n1 = std::exp(alps[i]) - 1;
      const double a1 = n1 > 0. ? 1. / std::log1p(n1) : 0.;
      for (int j = 0; j < ni; ++j) {
        const double f1 = n1 > 0. ? a1 * std::log1p(RndmUniform() * n1) : 0.5;
        const double f0 = 1. - f1;
        DriftPoint point;
        for (size_t k = 0; k < 3; ++k) {
          point.x[k] = f0 * x0[k] + f1 * x1[k]; 
        }
        point.t = f0 * driftLine[i].t + f1 * driftLine[i + 1].t;
        point.particle = other;
        point.n = 1;
        secondaries.push_back(std::move(point));
      }
    }
    // If trapped, exit the loop over the drift line.
    if (status == StatusAttached) return true;
  }
  return true;
}

bool AvalancheMC::ComputeAlphaEta(const Particle particle,
                                  std::vector<DriftPoint>& driftLine,
                                  std::vector<double>& alps,
                                  std::vector<double>& etas) const {
  // -----------------------------------------------------------------------
  //    DLCEQU - Computes equilibrated alphas and etas.
  // -----------------------------------------------------------------------

  // Loop a first time over the drift line and get alpha at each point.
  size_t nPoints = driftLine.size();
  alps.assign(nPoints, 0.);
  etas.assign(nPoints, 0.);
  for (size_t i = 0; i < nPoints; ++i) {
    const auto& x0 = driftLine[i].x;
    std::array<double, 3> e0, b0;
    Medium* medium = nullptr;
    if (GetField(x0, e0, b0, medium) != 0) continue;
    alps[i] = GetTownsend(particle, medium, x0, e0, b0);
  }

  std::vector<DriftPoint> driftLineExt;
  for (size_t i = 0; i < nPoints - 1; ++i) {
    const auto& p0 = driftLine[i];
    const auto& p1 = driftLine[i + 1];
    driftLineExt.push_back(p0);
    if (Slope(alps[i], alps[i + 1]) < 0.5) continue;
    auto xm = MidPoint(p0.x, p1.x);
    std::array<double, 3> em, bm;
    Medium* medium = nullptr;
    if (GetField(xm, em, bm, medium) != 0) continue;
    auto pm = p0;
    pm.x = xm;
    pm.t = 0.5 * (p0.t + p1.t);
    driftLineExt.push_back(std::move(pm));
  }
  driftLineExt.push_back(driftLine.back());
  driftLine.swap(driftLineExt);

  nPoints = driftLine.size();
  alps.assign(nPoints, 0.);
  etas.assign(nPoints, 0.);
  if (nPoints < 2) return true;

  // Locations and weights for Gaussian integration.
  constexpr size_t nG = 6;
  auto tg = Numerics::GaussLegendreNodes6();
  auto wg = Numerics::GaussLegendreWeights6();

  bool equilibrate = m_doEquilibration;
  // Loop over the drift line.
  for (size_t i = 0; i < nPoints - 1; ++i) {
    const auto& x0 = driftLine[i].x;
    const auto& x1 = driftLine[i + 1].x;
    // Compute the step length.
    const std::array<double, 3> del = {x1[0] - x0[0], x1[1] - x0[1],
                                       x1[2] - x0[2]};
    const double dmag = Mag(del);
    if (dmag < Small) continue;
    const double veff = dmag / (driftLine[i + 1].t - driftLine[i].t);
    // Integrate drift velocity and Townsend and attachment coefficients.
    std::array<double, 3> vd = {0., 0., 0.};
    for (size_t j = 0; j < nG; ++j) {
      const double f = 0.5 * (1. + tg[j]);
      std::array<double, 3> x = x0;
      for (size_t k = 0; k < 3; ++k) x[k] += f * del[k];
      // Get the field.
      std::array<double, 3> e;
      std::array<double, 3> b;
      Medium* medium = nullptr;
      const int status = GetField(x, e, b, medium);
      // Make sure we are in a drift medium.
      if (status != 0) {
        // Check if this point is the last but one.
        if (i < nPoints - 2) {
          std::cerr << m_className << "::ComputeAlphaEta: Got status " << status
                    << " at segment " << j + 1 << "/" << nG 
                    << ", drift point " << i + 1 << "/" << nPoints << ".\n";
          return false;
        }
        continue;
      }
      // Get the drift velocity.
      std::array<double, 3> v;
      if (!GetVelocity(particle, medium, x, e, b, v)) continue;
      // Get Townsend and attachment coefficients.
      double alpha = GetTownsend(particle, medium, x, e, b);
      double eta = GetAttachment(particle, medium, x, e, b);
      if (eta < 0.) {
        eta = std::abs(eta) * Mag(v) / veff;
        equilibrate = false;
      }
      for (size_t k = 0; k < 3; ++k) vd[k] += wg[j] * v[k];
      alps[i] += wg[j] * alpha;
      etas[i] += wg[j] * eta;
    }

    // Compute the scaling factor for the projected length.
    double scale = 1.;
    if (equilibrate) {
      const double vdmag = Mag(vd);
      if (vdmag * dmag <= 0.) {
        scale = 0.;
      } else {
        const double dinv = del[0] * vd[0] + del[1] * vd[1] + del[2] * vd[2];
        scale = dinv < 0. ? 0. : dinv / (vdmag * dmag);
      }
    }
    alps[i] *= 0.5 * dmag * scale;
    etas[i] *= 0.5 * dmag * scale;
  }

  // Skip equilibration if projection has not been requested.
  if (!equilibrate) return true;
  if (!Equilibrate(alps)) {
    if (m_debug) {
      std::cerr << m_className << "::ComputeAlphaEta:\n"
                << "    Unable to even out alpha steps.\n"
                << "    Calculation is probably inaccurate.\n";
    }
    return false;
  }
  if (!Equilibrate(etas)) {
    if (m_debug) {
      std::cerr << m_className << "::ComputeAlphaEta:\n"
                << "    Unable to even out eta steps.\n"
                << "    Calculation is probably inaccurate.\n";
    }
    return false;
  }
  // Seems to have worked.
  return true;
}

bool AvalancheMC::Equilibrate(std::vector<double>& alphas) const {
  const size_t nPoints = alphas.size();
  // Try to equilibrate the returning parts.
  for (size_t i = 0; i < nPoints - 1; ++i) {
    // Skip non-negative points.
    if (alphas[i] >= 0.) continue;
    // Targets for subtracting
    double sub1 = -0.5 * alphas[i];
    double sub2 = sub1;
    bool try1 = false;
    bool try2 = false;
    // Try to subtract half in earlier points.
    for (size_t j = 0; j < i - 1; ++j) {
      if (alphas[i - j] > sub1) {
        alphas[i - j] -= sub1;
        alphas[i] += sub1;
        sub1 = 0.;
        try1 = true;
        break;
      } else if (alphas[i - j] > 0.) {
        alphas[i] += alphas[i - j];
        sub1 -= alphas[i - j];
        alphas[i - j] = 0.;
      }
    }
    // Try to subtract the other half in later points.
    for (size_t j = 0; j < nPoints - i - 1; ++j) {
      if (alphas[i + j] > sub2) {
        alphas[i + j] -= sub2;
        alphas[i] += sub2;
        sub2 = 0.;
        try2 = true;
        break;
      } else if (alphas[i + j] > 0.) {
        alphas[i] += alphas[i + j];
        sub2 -= alphas[i + j];
        alphas[i + j] = 0.;
      }
    }

    // Done if both sides have margin left.
    bool done = false;
    if (try1 && try2) {
      done = true;
    } else if (try1) {
      // Try earlier points again.
      sub1 = -alphas[i];
      for (size_t j = 0; j < i - 1; ++j) {
        if (alphas[i - j] > sub1) {
          alphas[i - j] -= sub1;
          alphas[i] += sub1;
          sub1 = 0.;
          done = true;
          break;
        } else if (alphas[i - j] > 0.) {
          alphas[i] += alphas[i - j];
          sub1 -= alphas[i - j];
          alphas[i - j] = 0.;
        }
      }
    } else if (try2) {
      // Try later points again.
      sub2 = -alphas[i];
      for (size_t j = 0; j < nPoints - i - 1; ++j) {
        if (alphas[i + j] > sub2) {
          alphas[i + j] -= sub2;
          alphas[i] += sub2;
          sub2 = 0.;
          done = true;
          break;
        } else if (alphas[i + j] > 0.) {
          alphas[i] += alphas[i + j];
          sub2 -= alphas[i + j];
          alphas[i + j] = 0.;
        }
      }
    }
    // See whether we succeeded.
    if (!done) return false;
  }
  return true;
}

void AvalancheMC::ComputeSignal(
    const Particle particle, const double q,
    const std::vector<DriftPoint>& driftLine) const {
  const size_t nPoints = driftLine.size();
  if (nPoints < 2) return;

  if (m_useWeightingPotential) {
    for (size_t i = 0; i < nPoints - 1; ++i) {
      const auto& x0 = driftLine[i].x;
      const auto& x1 = driftLine[i + 1].x;
      const double t0 = driftLine[i].t;
      const double t1 = driftLine[i + 1].t;
      m_sensor->AddSignal(q, t0, t1, x0[0], x0[1], x0[2], x1[0], x1[1], x1[2],
                          false, true);
    }
    return;
  }
  // Get the drift velocity at each point.
  std::vector<double> ts;
  std::vector<std::array<double, 3> > xs;
  std::vector<std::array<double, 3> > vs;
  for (const auto& p : driftLine) {
    std::array<double, 3> e;
    std::array<double, 3> b;
    Medium* medium = nullptr;
    int status = GetField(p.x, e, b, medium);
    if (status != 0) continue;
    std::array<double, 3> v;
    if (!GetVelocity(particle, medium, p.x, e, b, v)) continue;
    ts.push_back(p.t);
    xs.push_back(p.x);
    vs.push_back(std::move(v));
  }
  m_sensor->AddSignal(q, ts, xs, vs, {}, m_navg);
}

void AvalancheMC::ComputeInducedCharge(
    const double q, const std::vector<DriftPoint>& driftLine) const {
  if (driftLine.size() < 2) return;
  const auto& x0 = driftLine.front().x;
  const auto& x1 = driftLine.back().x;
  m_sensor->AddInducedCharge(q, x0[0], x0[1], x0[2], x1[0], x1[1], x1[2]);
}

void AvalancheMC::PrintError(const std::string& fcn, const std::string& par,
                             const Particle particle,
                             const std::array<double, 3>& x) const {
  const std::string ehi = particle == Particle::Electron
                              ? "electron"
                              : particle == Particle::Hole ? "hole" : "ion";
  std::cerr << m_className + "::" + fcn + ": Error calculating " + ehi + " "
            << par + " at " + PrintVec(x) << ".\n";
}
}  // namespace Garfield
