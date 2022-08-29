#include <iostream>
#include <algorithm>
#include <cctype>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/Track.hh"
#include "Garfield/ViewDrift.hh"

namespace Garfield {

Track::Track() : m_mass(MuonMass) { SetBetaGamma(3.); }

void Track::SetParticle(const std::string& part) {
  std::string id = part;
  std::transform(id.begin(), id.end(), id.begin(), 
                 [](unsigned char c) -> unsigned char { 
                   return std::toupper(c);
                 });
  m_isElectron = false;
  if (id == "ELECTRON" || id == "E-") {
    m_q = -1;
    m_mass = ElectronMass;
    m_spin = 1;
    m_isElectron = true;
    m_particleName = "e-";
  } else if (id == "POSITRON" || id == "E+") {
    m_q = 1;
    m_mass = ElectronMass;
    m_spin = 1;
    m_particleName = "e+";
  } else if (id == "MUON" || id == "MU" || id == "MU-") {
    m_q = -1;
    m_mass = MuonMass;
    m_spin = 1;
    m_particleName = "mu-";
  } else if (id == "MU+") {
    m_q = 1;
    m_mass = MuonMass;
    m_spin = 1;
    m_particleName = "mu+";
  } else if (id == "PION" || id == "PI" || id == "PI-") {
    m_q = -1;
    m_mass = 139.57018e6;
    m_spin = 0;
    m_particleName = "pi-";
  } else if (id == "PI+") {
    m_q = 1;
    m_mass = 139.57018e6;
    m_spin = 0;
    m_particleName = "pi+";
  } else if (id == "KAON" || id == "K" || id == "K-") {
    m_q = -1;
    m_mass = 493.677e6;
    m_spin = 0;
    m_particleName = "K-";
  } else if (id == "K+") {
    m_q = 1;
    m_mass = 493.677e6;
    m_spin = 0;
    m_particleName = "K+";
  } else if (id == "PROTON" || id == "P") {
    m_q = 1;
    m_mass = ProtonMass;
    m_spin = 1;
    m_particleName = "p";
  } else if (id == "ANTI-PROTON" || id == "ANTIPROTON" || 
             id == "P-BAR" || id == "PBAR") {
    m_q = -1;
    m_mass = ProtonMass;
    m_spin = 1;
    m_particleName = "pbar";
  } else if (id == "DEUTERON" || id == "D") {
    m_q = 1;
    m_mass = 1875.612793e6;
    m_spin = 2;
    m_particleName = "d";
  } else if (id == "ALPHA") {
    m_q = 2;
    m_mass = 3.727379240e9;
    m_spin = 0;
    m_particleName = "alpha";
  } else {
    std::cerr << m_className << "::SetParticle:\n"
              << "    Particle " << part << " is not defined.\n";
  }
}

void Track::SetEnergy(const double e) {
  if (e <= m_mass) {
    std::cerr << m_className << "::SetEnergy:\n"
              << "    Particle energy must be greater than the mass.\n";
    return;
  }

  m_energy = e;
  const double gamma = m_energy / m_mass;
  m_beta2 = 1. - 1. / (gamma * gamma);
  m_isChanged = true;
}

void Track::SetBetaGamma(const double bg) {
  if (bg <= 0.) {
    std::cerr << m_className << "::SetBetaGamma:\n"
              << "    Particle speed must be greater than zero.\n";
    return;
  }

  const double bg2 = bg * bg;
  m_energy = m_mass * sqrt(1. + bg2);
  m_beta2 = bg2 / (1. + bg2);
  m_isChanged = true;
}

void Track::SetBeta(const double beta) {
  if (beta <= 0. || beta >= 1.) {
    std::cerr << m_className << "::SetBeta:\n"
              << "    Beta must be between zero and one.\n";
    return;
  }

  m_beta2 = beta * beta;
  m_energy = m_mass * sqrt(1. / (1. - m_beta2));
  m_isChanged = true;
}

void Track::SetGamma(const double gamma) {
  if (gamma <= 1.) {
    std::cerr << m_className << "::SetGamma:\n"
              << "    Gamma must be greater than one.\n";
    return;
  }

  m_energy = m_mass * gamma;
  m_beta2 = 1. - 1. / (gamma * gamma);
  m_isChanged = true;
}

void Track::SetMomentum(const double p) {
  if (p <= 0.) {
    std::cerr << m_className << "::SetMomentum:\n"
              << "    Particle momentum must be greater than zero.\n";
    return;
  }

  m_energy = sqrt(m_mass * m_mass + p * p);
  const double bg = p / m_mass;
  m_beta2 = bg * bg / (1. + bg * bg);
  m_isChanged = true;
}

void Track::SetKineticEnergy(const double ekin) {
  if (ekin <= 0.) {
    std::cerr << m_className << "::SetKineticEnergy:\n"
              << "    Kinetic energy must be greater than zero.\n";
    return;
  }

  m_energy = m_mass + ekin;
  const double gamma = 1. + ekin / m_mass;
  m_beta2 = 1. - 1. / (gamma * gamma);
  m_isChanged = true;
}

void Track::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }

  m_sensor = s;
}

void Track::EnablePlotting(ViewDrift* view) {
  if (!view) {
    std::cerr << m_className << "::EnablePlotting: Null pointer.\n";
    return;
  }
  m_viewer = view;
}

void Track::DisablePlotting() {
  m_viewer = nullptr;
}

void Track::PlotNewTrack(const double x0, const double y0, const double z0) {
  if (!m_viewer) return;
  m_viewer->NewChargedParticleTrack(1, m_plotId, x0, y0, z0);
}

void Track::PlotCluster(const double x0, const double y0, const double z0) {
  if (m_viewer) m_viewer->AddTrackPoint(m_plotId, x0, y0, z0);
}

std::array<double, 3> Track::StepBfield(const double dt, 
    const double qoverm, const double vmag, double bx, double by, double bz,
    std::array<double, 3>& dir) {

  double bmag = sqrt(bx * bx + by * by + bz * bz);
  if (bmag < Garfield::Small) {
    const double step = vmag * dt;
    return {step * dir[0], step * dir[1], step * dir[2]};
  }
  std::array<std::array<double, 3>, 3> rot = {{{1, 0, 0}, {0, 1, 0}, {0, 0,      1}}};

  bx /= bmag;
  by /= bmag;
  bz /= bmag;
  const double bt = by * by + bz * bz;
  if (bt > Garfield::Small) {
    const double btInv = 1. / bt;
    rot[0][0] = bx;
    rot[0][1] = by;
    rot[0][2] = bz;
    rot[1][0] = -by;
    rot[2][0] = -bz;
    rot[1][1] = (bx * by * by + bz * bz) * btInv;
    rot[2][2] = (bx * bz * bz + by * by) * btInv;
    rot[1][2] = rot[2][1] = (bx - 1.) * by * bz * btInv;
  } else if (bx < 0.) {
    // B field is anti-parallel to x.
    rot[0][0] = -1.;
    rot[1][1] = -1.;
  }
  bmag *= Garfield::Tesla2Internal;
  const double omega = qoverm * Garfield::OmegaCyclotronOverB * bmag * 
                       Garfield::ElectronMass;
  const double cphi = cos(omega * dt);
  const double sphi = sin(omega * dt);

  // Rotate the initial direction vector to the local frame.
  std::array<double, 3> v0;
  v0[0] = rot[0][0] * dir[0] + rot[0][1] * dir[1] + rot[0][2] * dir[2];
  v0[1] = rot[1][0] * dir[0] + rot[1][1] * dir[1] + rot[1][2] * dir[2];
  v0[2] = rot[2][0] * dir[0] + rot[2][1] * dir[1] + rot[2][2] * dir[2];

  // Calculate the new direction in the local frame. 
  std::array<double, 3> v1;
  v1[0] = v0[0];
  v1[1] = v0[1] * cphi + v0[2] * sphi;
  v1[2] = v0[2] * cphi - v0[1] * sphi;

  // Rotate the direction vector back to the global frame.
  dir[0] = rot[0][0] * v1[0] + rot[1][0] * v1[1] + rot[2][0] * v1[2];
  dir[1] = rot[0][1] * v1[0] + rot[1][1] * v1[1] + rot[2][1] * v1[2];
  dir[2] = rot[0][2] * v1[0] + rot[1][2] * v1[1] + rot[2][2] * v1[2];

  // Calculate the new position in the local frame...
  const double rho = vmag / omega;
  const double u = vmag * v0[0] * dt;
  const double v = rho * (v0[1] * sphi + v0[2] * (1. - cphi));
  const double w = rho * (v0[2] * sphi - v0[1] * (1. - cphi));
  // .... and in the global frame.
  std::array<double, 3> pos;
  pos[0] = rot[0][0] * u + rot[1][0] * v + rot[2][0] * w; 
  pos[1] = rot[0][1] * u + rot[1][1] * v + rot[2][1] * w;
  pos[2] = rot[0][2] * u + rot[1][2] * v + rot[2][2] * w;
  return pos;
}

}
