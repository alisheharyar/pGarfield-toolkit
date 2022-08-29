#include <iostream>

#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/matter/MatterDef.h"

#include "heed++/code/ElElasticScat.h"
#include "heed++/code/EnTransfCS.h"
#include "heed++/code/HeedCluster.h"
#include "heed++/code/HeedCondElectron.h"
#include "heed++/code/HeedDeltaElectron.h"
#include "heed++/code/HeedDeltaElectronCS.h"
#include "heed++/code/HeedMatterDef.h"
#include "heed++/code/HeedParticle.h"
#include "heed++/code/HeedPhoton.h"
#include "heed++/code/PhotoAbsCSLib.h"

#include "HeedChamber.hh"
#include "HeedFieldMap.h"

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"

#include "Garfield/TrackHeed.hh"

namespace {

void ClearBank(std::vector<Heed::gparticle*>& bank) {
  for (auto particle : bank)
    if (particle) delete particle;
  bank.clear();
}

Heed::MolecPhotoAbsCS makeMPACS(const std::string& atom, const int n,
                                const double w = 0.) {
  return Heed::MolecPhotoAbsCS(Heed::PhotoAbsCSLib::getAPACS(atom), n, w);
}

Heed::MolecPhotoAbsCS makeMPACS(const std::string& atom1, const int n1,
                                const std::string& atom2, const int n2,
                                const double w = 0.) {
  return Heed::MolecPhotoAbsCS(Heed::PhotoAbsCSLib::getAPACS(atom1), n1, 
                               Heed::PhotoAbsCSLib::getAPACS(atom2), n2, w);
}

Heed::MolecPhotoAbsCS makeMPACS(const std::string& atom1, const int n1,
                                const std::string& atom2, const int n2,
                                const std::string& atom3, const int n3,
                                const double w = 0.) {
  return Heed::MolecPhotoAbsCS(Heed::PhotoAbsCSLib::getAPACS(atom1), n1, 
                               Heed::PhotoAbsCSLib::getAPACS(atom2), n2,
                               Heed::PhotoAbsCSLib::getAPACS(atom3), n3, w);
}

}

// Actual class implementation

namespace Garfield {

TrackHeed::TrackHeed() : Track() {
  m_className = "TrackHeed";
  m_conductionElectrons.reserve(1000);
  m_conductionIons.reserve(1000);

  m_fieldMap.reset(new Heed::HeedFieldMap());
}

TrackHeed::~TrackHeed() {}

bool TrackHeed::NewTrack(const double x0, const double y0, const double z0,
                         const double t0, const double dx0, const double dy0,
                         const double dz0) {
  m_hasActiveTrack = false;
  m_ready = false;

  // Make sure the sensor has been set.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    return false;
  }

  bool update = false;
  if (!UpdateBoundingBox(update)) return false;

  // Make sure the initial position is inside an ionisable medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::NewTrack:\n"
              << "    No medium at initial position.\n";
    return false;
  } else if (!medium->IsIonisable()) {
    std::cerr << m_className << "::NewTrack:\n"
              << "    Medium at initial position is not ionisable.\n";
    return false;
  }

  // Check if the medium has changed since the last call.
  if (medium->GetName() != m_mediumName ||
      fabs(medium->GetMassDensity() - m_mediumDensity) > 1.e-9) {
    m_isChanged = true;
  }

  // If medium, particle or bounding box have changed,
  // update the cross-sections.
  if (m_isChanged) {
    if (!Initialise(medium)) return false;
    m_isChanged = false;
    m_mediumName = medium->GetName();
    m_mediumDensity = medium->GetMassDensity();
  }

  ClearParticleBank();
  m_photons.clear();
  m_deltaElectrons.clear();
  m_conductionElectrons.clear();
  m_conductionIons.clear();

  // Check the direction vector.
  double dx = dx0, dy = dy0, dz = dz0;
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d < Small) {
    if (m_debug) {
      std::cout << m_className << "::NewTrack:\n"
                << "    Direction vector has zero norm.\n"
                << "    Initial direction is randomized.\n";
    }
    // Null vector. Sample the direction isotropically.
    RndmDirection(dx, dy, dz);
  } else {
    // Normalise the direction vector.
    dx /= d;
    dy /= d;
    dz /= d;
  }
  Heed::vec velocity(dx, dy, dz);
  velocity = velocity * Heed::CLHEP::c_light * GetBeta();

  if (m_debug) {
    std::cout << m_className << "::NewTrack:\n    Track starts at (" << x0
              << ", " << y0 << ", " << z0 << ") at time " << t0 << "\n"
              << "    Direction: (" << dx << ", " << dy << ", " << dz << ")\n";
  }

  // Initial position (shift with respect to bounding box center and
  // convert from cm to mm).
  Heed::point p0((x0 - m_cX) * 10., (y0 - m_cY) * 10., (z0 - m_cZ) * 10.);

  // Setup the particle.
  Heed::particle_def* particleType = &Heed::muon_minus_def;
  if (m_particleName == "e-") {
    particleType = &Heed::electron_def;
  } else if (m_particleName == "e+") {
    particleType = &Heed::positron_def;
  } else if (m_particleName == "mu-") {
    particleType = &Heed::muon_minus_def;
  } else if (m_particleName == "mu+") {
    particleType = &Heed::muon_plus_def;
  } else if (m_particleName == "pi-") {
    particleType = &Heed::pi_minus_meson_def;
  } else if (m_particleName == "pi+") {
    particleType = &Heed::pi_plus_meson_def;
  } else if (m_particleName == "K-") {
    particleType = &Heed::K_minus_meson_def;
  } else if (m_particleName == "K+") {
    particleType = &Heed::K_plus_meson_def;
  } else if (m_particleName == "p") {
    particleType = &Heed::proton_def;
  } else if (m_particleName == "pbar") {
    particleType = &Heed::anti_proton_def;
  } else if (m_particleName == "d") {
    particleType = &Heed::deuteron_def;
  } else if (m_particleName == "alpha") {
    particleType = &Heed::alpha_particle_def;
  } else if (m_particleName == "exotic") {
    // User defined particle
    m_particle_def.reset(new Heed::particle_def(Heed::pi_plus_meson_def));
    m_particle_def->set_mass(m_mass * 1.e-6);
    m_particle_def->set_charge(m_q);
    particleType = m_particle_def.get();
  } else {
    // Not a predefined particle, use muon definition.
    if (m_q > 0.) {
      particleType = &Heed::muon_minus_def;
    } else {
      particleType = &Heed::muon_plus_def;
    }
  }

  Heed::HeedParticle particle(m_chamber.get(), p0, velocity, t0, particleType,
                              m_fieldMap.get(), m_coulombScattering);
  // Set the step limits.
  particle.set_step_limits(m_maxStep * Heed::CLHEP::cm,
                           m_radStraight * Heed::CLHEP::cm,
                           m_stepAngleStraight * Heed::CLHEP::rad,
                           m_stepAngleCurved * Heed::CLHEP::rad);
  // Transport the particle.
  if (m_oneStepFly) {
    particle.fly(m_particleBank, true);
  } else {
    particle.fly(m_particleBank);
  }

  m_bankIterator = m_particleBank.begin();
  m_hasActiveTrack = true;
  m_ready = true;

  // Plot the new track.
  if (m_viewer) PlotNewTrack(x0, y0, z0);
  return true;
}

double TrackHeed::GetClusterDensity() {

  if (!m_transferCs) {
    std::cerr << m_className << "::GetClusterDensity:\n"
              << "    Ionisation cross-section is not available.\n";
    return 0.;
  }

  return m_transferCs->quanC;
}

double TrackHeed::GetStoppingPower() {

  if (!m_transferCs) {
    std::cerr << m_className << "::GetStoppingPower:\n"
              << "    Ionisation cross-section is not available.\n";
    return 0.;
  }

  return m_transferCs->meanC1 * 1.e6;
}

bool TrackHeed::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& n, double& e, double& extra) {
  int ni = 0, np = 0;
  return GetCluster(xcls, ycls, zcls, tcls, n, ni, np, e, extra);
}

bool TrackHeed::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& ne, int& ni, double& e, 
                           double& extra) {
  int np = 0;
  return GetCluster(xcls, ycls, zcls, tcls, ne, ni, np, e, extra);
}

bool TrackHeed::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& ne, int& ni, int& np, 
                           double& e, double& extra) {
  // Initialise and reset.
  xcls = ycls = zcls = tcls = 0.;
  extra = 0.;
  ne = ni = np = 0;
  e = 0.;
  m_photons.clear();
  m_deltaElectrons.clear();
  m_conductionElectrons.clear();
  m_conductionIons.clear();

  // Make sure NewTrack has been called successfully.
  if (!m_ready) {
    std::cerr << m_className << "::GetCluster:\n"
              << "    Track has not been initialized. Call NewTrack first.\n";
    return false;
  }

  if (m_particleBank.empty()) return false;
  std::vector<Heed::gparticle*>::const_iterator end = m_particleBank.end();
  if (m_bankIterator == end) return false;

  // Look for the next cluster (i. e. virtual photon) in the list.
  Heed::HeedPhoton* virtualPhoton = nullptr;
  for (; m_bankIterator != end; ++m_bankIterator) {
    // Convert the particle to a (virtual) photon.
    virtualPhoton = dynamic_cast<Heed::HeedPhoton*>(*m_bankIterator);
    if (!virtualPhoton) {
      std::cerr << m_className << "::GetCluster:\n"
                << "    Particle is not a virtual photon. Program bug!\n";
      // Try the next element.
      continue;
    }
    // Get the location of the interaction (convert from mm to cm
    // and shift with respect to bounding box center).
    xcls = virtualPhoton->position().x * 0.1 + m_cX;
    ycls = virtualPhoton->position().y * 0.1 + m_cY;
    zcls = virtualPhoton->position().z * 0.1 + m_cZ;
    tcls = virtualPhoton->time();
    // Skip clusters outside the drift area or outside the active medium.
    if (!IsInside(xcls, ycls, zcls)) continue;
    // Add the first ion (at the position of the cluster).
    m_conductionIons.emplace_back(
        Heed::HeedCondElectron(Heed::point(virtualPhoton->position()), tcls));
    ++m_bankIterator;
    break;
  }

  // Stop if we did not find a virtual photon.
  if (!virtualPhoton) return false;
  // Plot the cluster, if requested.
  if (m_viewer) PlotCluster(xcls, ycls, zcls);

  std::vector<Heed::gparticle*> secondaries;
  // Transport the virtual photon.
  virtualPhoton->fly(secondaries);
  // Get the transferred energy (convert from MeV to eV).
  e = virtualPhoton->m_energy * 1.e6;

  while (!secondaries.empty()) {
    std::vector<Heed::gparticle*> newSecondaries;
    // Loop over the secondaries.
    for (auto secondary : secondaries) {
      // Check if it is a delta electron.
      auto delta = dynamic_cast<Heed::HeedDeltaElectron*>(secondary);
      if (delta) {
        extra += delta->kinetic_energy() * 1.e6;
        const double x = delta->position().x * 0.1 + m_cX;
        const double y = delta->position().y * 0.1 + m_cY;
        const double z = delta->position().z * 0.1 + m_cZ;
        if (!IsInside(x, y, z)) continue;
        if (m_doDeltaTransport) {
          // Transport the delta electron.
          delta->fly(newSecondaries);
          // Add the conduction electrons and ions to the list.
          m_conductionElectrons.insert(m_conductionElectrons.end(),
                                       delta->conduction_electrons.begin(),
                                       delta->conduction_electrons.end());
          m_conductionIons.insert(m_conductionIons.end(),
                                  delta->conduction_ions.begin(),
                                  delta->conduction_ions.end());
        } else {
          // Add the delta electron to the list, for later use.
          DeltaElectron deltaElectron;
          deltaElectron.x = delta->position().x * 0.1 + m_cX;
          deltaElectron.y = delta->position().y * 0.1 + m_cY;
          deltaElectron.z = delta->position().z * 0.1 + m_cZ;
          deltaElectron.t = delta->time();
          deltaElectron.e = delta->kinetic_energy() * 1.e6;
          deltaElectron.dx = delta->direction().x;
          deltaElectron.dy = delta->direction().y;
          deltaElectron.dz = delta->direction().z;
          m_deltaElectrons.push_back(std::move(deltaElectron));
        }
        continue;
      }
      // Check if it is a real photon.
      auto photon = dynamic_cast<Heed::HeedPhoton*>(secondary);
      if (!photon) {
        std::cerr << m_className << "::GetCluster:\n"
                  << "    Particle is neither an electron nor a photon.\n";
        continue;
      }
      extra += photon->m_energy * 1.e6;
      const double x = photon->position().x * 0.1 + m_cX;
      const double y = photon->position().y * 0.1 + m_cY;
      const double z = photon->position().z * 0.1 + m_cZ;
      if (!IsInside(x, y, z)) continue;
      // Transport the photon.
      if (m_doPhotonReabsorption) {
        photon->fly(newSecondaries);
      } else {
        Photon unabsorbedPhoton;
        unabsorbedPhoton.x = photon->position().x * 0.1 + m_cX;
        unabsorbedPhoton.y = photon->position().y * 0.1 + m_cY;
        unabsorbedPhoton.z = photon->position().z * 0.1 + m_cZ;
        unabsorbedPhoton.t = photon->time();
        unabsorbedPhoton.e = photon->m_energy * 1.e6;
        unabsorbedPhoton.dx = photon->direction().x;
        unabsorbedPhoton.dy = photon->direction().y;
        unabsorbedPhoton.dz = photon->direction().z;
        m_photons.push_back(std::move(unabsorbedPhoton));
      }
    }
    for (auto secondary : secondaries)
      if (secondary) delete secondary;
    secondaries.clear();
    secondaries.swap(newSecondaries);
  }
  // Get the total number of electrons produced in this step.
  ne = m_doDeltaTransport ? m_conductionElectrons.size()
                          : m_deltaElectrons.size();
  ni = m_conductionIons.size();
  np = m_photons.size();
  return true;
}

bool TrackHeed::GetElectron(const unsigned int i, double& x, double& y,
                            double& z, double& t, double& e, double& dx,
                            double& dy, double& dz) {
  // Make sure NewTrack has successfully been called.
  if (!m_ready) {
    std::cerr << m_className << "::GetElectron:\n"
              << "    Track has not been initialized. Call NewTrack first.\n";
    return false;
  }

  if (m_doDeltaTransport) {
    // Make sure an electron with this number exists.
    if (i >= m_conductionElectrons.size()) {
      std::cerr << m_className << "::GetElectron: Index out of range.\n";
      return false;
    }

    x = m_conductionElectrons[i].x * 0.1 + m_cX;
    y = m_conductionElectrons[i].y * 0.1 + m_cY;
    z = m_conductionElectrons[i].z * 0.1 + m_cZ;
    t = m_conductionElectrons[i].time;
    e = 0.;
    dx = dy = dz = 0.;

  } else {
    // Make sure a delta electron with this number exists.
    if (i >= m_deltaElectrons.size()) {
      std::cerr << m_className << "::GetElectron:\n"
                << "    Delta electron number out of range.\n";
      return false;
    }

    x = m_deltaElectrons[i].x;
    y = m_deltaElectrons[i].y;
    z = m_deltaElectrons[i].z;
    t = m_deltaElectrons[i].t;
    e = m_deltaElectrons[i].e;
    dx = m_deltaElectrons[i].dx;
    dy = m_deltaElectrons[i].dy;
    dz = m_deltaElectrons[i].dz;
  }
  return true;
}

bool TrackHeed::GetIon(const unsigned int i, double& x, double& y, double& z,
                       double& t) const {
  // Make sure a "conduction" ion with this number exists.
  if (i >= m_conductionIons.size()) {
    std::cerr << m_className << "::GetIon: Index out of range.\n";
    return false;
  }

  x = m_conductionIons[i].x * 0.1 + m_cX;
  y = m_conductionIons[i].y * 0.1 + m_cY;
  z = m_conductionIons[i].z * 0.1 + m_cZ;
  t = m_conductionIons[i].time;
  return true;
}

bool TrackHeed::GetPhoton(const unsigned int i, double& x, double& y,
                          double& z, double& t, double& e, double& dx,
                          double& dy, double& dz) const {
  // Make sure a photon with this number exists.
  if (i >= m_photons.size()) {
    std::cerr << m_className << "::GetPhoton: Index out of range.\n";
    return false;
  }

  x = m_photons[i].x;
  y = m_photons[i].y;
  z = m_photons[i].z;
  t = m_photons[i].t;
  e = m_photons[i].e;
  dx = m_photons[i].dx;
  dy = m_photons[i].dy;
  dz = m_photons[i].dz;
  return true;
}

void TrackHeed::TransportDeltaElectron(const double x0, const double y0,
                                       const double z0, const double t0,
                                       const double e0, const double dx0,
                                       const double dy0, const double dz0,
                                       int& nel) {
  int ni = 0;
  return TransportDeltaElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0, nel, ni);
}

void TrackHeed::TransportDeltaElectron(const double x0, const double y0,
                                       const double z0, const double t0,
                                       const double e0, const double dx0,
                                       const double dy0, const double dz0,
                                       int& nel, int& ni) {
  nel = 0;
  ni = 0;

  // Check if delta electron transport was disabled.
  if (!m_doDeltaTransport) {
    std::cerr << m_className << "::TransportDeltaElectron:\n"
              << "    Delta electron transport has been switched off.\n";
    return;
  }

  // Make sure the sensor has been set.
  if (!m_sensor) {
    std::cerr << m_className << "::TransportDeltaElectron:\n"
              << "    Sensor is not defined.\n";
    m_ready = false;
    return;
  }

  bool update = false;
  if (!UpdateBoundingBox(update)) return;

  // Make sure the initial position is inside an ionisable medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::TransportDeltaElectron:\n"
              << "    No medium at initial position.\n";
    return;
  } else if (!medium->IsIonisable()) {
    std::cerr << "TrackHeed:TransportDeltaElectron:\n"
              << "    Medium at initial position is not ionisable.\n";
    m_ready = false;
    return;
  }

  // Check if the medium has changed since the last call.
  if (medium->GetName() != m_mediumName ||
      fabs(medium->GetMassDensity() - m_mediumDensity) > 1.e-9) {
    m_isChanged = true;
    update = true;
    m_ready = false;
    m_hasActiveTrack = false;
  }

  // If medium or bounding box have changed, update the "chamber".
  if (update) {
    if (!Initialise(medium)) return;
    m_ready = true;
    m_mediumName = medium->GetName();
    m_mediumDensity = medium->GetMassDensity();
  }
  m_photons.clear();
  m_deltaElectrons.clear();
  m_conductionElectrons.clear();
  m_conductionIons.clear();

  // Initial position (shift with respect to bounding box center and
  // convert from cm to mm).
  Heed::point p0((x0 - m_cX) * 10., (y0 - m_cY) * 10., (z0 - m_cZ) * 10.);

  // Make sure the kinetic energy is positive.
  if (e0 <= 0.) {
    // Just create a conduction electron on the spot.
    m_conductionElectrons.emplace_back(Heed::HeedCondElectron(p0, t0));
    nel = 1;
    return;
  }

  // Check the direction vector.
  double dx = dx0, dy = dy0, dz = dz0;
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d <= 0.) {
    // Null vector. Sample the direction isotropically.
    RndmDirection(dx, dy, dz);
  } else {
    // Normalise the direction vector.
    dx /= d;
    dy /= d;
    dz /= d;
  }
  Heed::vec velocity(dx, dy, dz);

  // Calculate the speed for the given kinetic energy.
  const double gamma = 1. + e0 / ElectronMass;
  const double beta = sqrt(1. - 1. / (gamma * gamma));
  double speed = Heed::CLHEP::c_light * beta;
  velocity = velocity * speed;

  // Transport the electron.
  std::vector<Heed::gparticle*> secondaries;
  Heed::HeedDeltaElectron delta(m_chamber.get(), p0, velocity, t0, 0,
                                m_fieldMap.get());
  delta.fly(secondaries);
  ClearBank(secondaries);

  m_conductionElectrons.swap(delta.conduction_electrons);
  m_conductionIons.swap(delta.conduction_ions);
  nel = m_conductionElectrons.size();
  ni = m_conductionIons.size();
}

void TrackHeed::TransportPhoton(const double x0, const double y0,
                                const double z0, const double t0,
                                const double e0, const double dx0,
                                const double dy0, const double dz0, int& ne) {
  int ni = 0, np = 0;
  TransportPhoton(x0, y0, z0, t0, e0, dx0, dy0, dz0, ne, ni, np);
}

void TrackHeed::TransportPhoton(const double x0, const double y0,
                                const double z0, const double t0,
                                const double e0, const double dx0,
                                const double dy0, const double dz0, int& ne,
                                int& ni) {
  int np = 0;
  TransportPhoton(x0, y0, z0, t0, e0, dx0, dy0, dz0, ne, ni, np);
} 

void TrackHeed::TransportPhoton(const double x0, const double y0,
                                const double z0, const double t0,
                                const double e0, const double dx0,
                                const double dy0, const double dz0, int& ne,
                                int& ni, int& np) {
  ne = 0;
  ni = 0;
  np = 0;

  // Make sure the energy is positive.
  if (e0 <= 0.) {
    std::cerr << m_className << "::TransportPhoton:\n"
              << "    Photon energy must be positive.\n";
    return;
  }

  // Make sure the sensor has been set.
  if (!m_sensor) {
    std::cerr << m_className << "::TransportPhoton: Sensor is not defined.\n";
    m_ready = false;
    return;
  }

  bool update = false;
  if (!UpdateBoundingBox(update)) return;

  // Make sure the initial position is inside an ionisable medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::TransportPhoton:\n"
              << "    No medium at initial position.\n";
    return;
  } else if (!medium->IsIonisable()) {
    std::cerr << "TrackHeed:TransportPhoton:\n"
              << "    Medium at initial position is not ionisable.\n";
    m_ready = false;
    return;
  }

  // Check if the medium has changed since the last call.
  if (medium->GetName() != m_mediumName ||
      fabs(medium->GetMassDensity() - m_mediumDensity) > 1.e-9) {
    m_isChanged = true;
    update = true;
    m_ready = false;
  }

  // If medium or bounding box have changed, update the "chamber".
  if (update) {
    if (!Initialise(medium)) return;
    m_ready = true;
    m_mediumName = medium->GetName();
    m_mediumDensity = medium->GetMassDensity();
  }

  // Delete the particle bank.
  // Clusters from the current track will be lost.
  m_hasActiveTrack = false;
  ClearParticleBank();
  m_photons.clear();
  m_deltaElectrons.clear();
  m_conductionElectrons.clear();
  m_conductionIons.clear();

  // Check the direction vector.
  double dx = dx0, dy = dy0, dz = dz0;
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d <= 0.) {
    // Null vector. Sample the direction isotropically.
    RndmDirection(dx, dy, dz);
  } else {
    // Normalise the direction vector.
    dx /= d;
    dy /= d;
    dz /= d;
  }
  Heed::vec velocity(dx, dy, dz);
  velocity = velocity * Heed::CLHEP::c_light;

  // Initial position (shift with respect to bounding box center and
  // convert from cm to mm).
  Heed::point p0((x0 - m_cX) * 10., (y0 - m_cY) * 10., (z0 - m_cZ) * 10.);

  // Create and transport the photon.
  Heed::HeedPhoton photon(m_chamber.get(), p0, velocity, t0, 0, e0 * 1.e-6,
                          m_fieldMap.get());
  std::vector<Heed::gparticle*> secondaries;
  photon.fly(secondaries);
  if (secondaries.empty()) {
    Photon unabsorbedPhoton;
    unabsorbedPhoton.x = photon.position().x * 0.1 + m_cX;
    unabsorbedPhoton.y = photon.position().y * 0.1 + m_cY;
    unabsorbedPhoton.z = photon.position().z * 0.1 + m_cZ;
    unabsorbedPhoton.t = photon.time();
    unabsorbedPhoton.e = photon.m_energy * 1.e6;
    unabsorbedPhoton.dx = photon.direction().x;
    unabsorbedPhoton.dy = photon.direction().y;
    unabsorbedPhoton.dz = photon.direction().z;
    m_photons.push_back(std::move(unabsorbedPhoton));
  }

  while (!secondaries.empty()) {
    std::vector<Heed::gparticle*> newSecondaries;
    // Loop over the particle bank and look for daughter particles.
    std::vector<Heed::gparticle*>::iterator it;
    for (it = secondaries.begin(); it != secondaries.end(); ++it) {
      // Check if it is a delta electron.
      auto delta = dynamic_cast<Heed::HeedDeltaElectron*>(*it);
      if (delta) {
        if (m_doDeltaTransport) {
          // Transport the delta electron.
          delta->fly(newSecondaries);
          // Add the conduction electrons to the list.
          m_conductionElectrons.insert(m_conductionElectrons.end(),
                                       delta->conduction_electrons.begin(),
                                       delta->conduction_electrons.end());
          m_conductionIons.insert(m_conductionIons.end(),
                                  delta->conduction_ions.begin(),
                                  delta->conduction_ions.end());
        } else {
          // Add the delta electron to the list, for later use.
          DeltaElectron deltaElectron;
          deltaElectron.x = delta->position().x * 0.1 + m_cX;
          deltaElectron.y = delta->position().y * 0.1 + m_cY;
          deltaElectron.z = delta->position().z * 0.1 + m_cZ;
          deltaElectron.t = delta->time();
          deltaElectron.e = delta->kinetic_energy() * 1.e6;
          deltaElectron.dx = delta->direction().x;
          deltaElectron.dy = delta->direction().y;
          deltaElectron.dz = delta->direction().z;
          m_deltaElectrons.push_back(std::move(deltaElectron));
        }
        continue;
      }
      // Check if it is a fluorescence photon.
      auto fluorescencePhoton = dynamic_cast<Heed::HeedPhoton*>(*it);
      if (!fluorescencePhoton) {
        std::cerr << m_className << "::TransportPhoton:\n"
                  << "    Unknown secondary particle.\n";
        ClearBank(secondaries);
        ClearBank(newSecondaries);
        return;
      }
      if (m_doPhotonReabsorption) {
        fluorescencePhoton->fly(newSecondaries);
      } else {
        Photon unabsorbedPhoton;
        unabsorbedPhoton.x = fluorescencePhoton->position().x * 0.1 + m_cX;
        unabsorbedPhoton.y = fluorescencePhoton->position().y * 0.1 + m_cY;
        unabsorbedPhoton.z = fluorescencePhoton->position().z * 0.1 + m_cZ;
        unabsorbedPhoton.t = fluorescencePhoton->time();
        unabsorbedPhoton.e = fluorescencePhoton->m_energy * 1.e6;
        unabsorbedPhoton.dx = fluorescencePhoton->direction().x;
        unabsorbedPhoton.dy = fluorescencePhoton->direction().y;
        unabsorbedPhoton.dz = fluorescencePhoton->direction().z;
        m_photons.push_back(std::move(unabsorbedPhoton));
      }
    }
    secondaries.swap(newSecondaries);
    ClearBank(newSecondaries);
  }
  ClearBank(secondaries);
  // Get the total number of electrons produced in this step.
  ne = m_doDeltaTransport ? m_conductionElectrons.size()
                          : m_deltaElectrons.size();
  ni = m_conductionIons.size();
  np = m_photons.size();
}

void TrackHeed::EnableElectricField() { m_fieldMap->UseEfield(true); }
void TrackHeed::DisableElectricField() { m_fieldMap->UseEfield(false); }
void TrackHeed::EnableMagneticField() { m_fieldMap->UseBfield(true); }
void TrackHeed::DisableMagneticField() { m_fieldMap->UseBfield(false); }

void TrackHeed::SetEnergyMesh(const double e0, const double e1,
                              const int nsteps) {
  if (fabs(e1 - e0) < Small) {
    std::cerr << m_className << "::SetEnergyMesh:\n"
              << "    Invalid energy range:\n"
              << "    " << e0 << " < E [eV] < " << e1 << "\n";
    return;
  }

  if (nsteps <= 0) {
    std::cerr << m_className << "::SetEnergyMesh:\n"
              << "    Number of intervals must be > 0.\n";
    return;
  }

  m_emin = std::min(e0, e1);
  m_emax = std::max(e0, e1);
  m_emin *= 1.e-6;
  m_emax *= 1.e-6;
  m_nEnergyIntervals = nsteps;
}

void TrackHeed::SetParticleUser(const double m, const double z) {
  if (fabs(z) < Small) {
    std::cerr << m_className << "::SetParticleUser:\n"
              << "    Particle cannot have zero charge.\n";
    return;
  }
  if (m < Small) {
    std::cerr << m_className << "::SetParticleUser:\n"
              << "    Particle mass must be greater than zero.\n";
  }
  m_q = z;
  m_mass = m;
  m_isElectron = false;
  m_spin = 0;
  m_particleName = "exotic";
}

bool TrackHeed::Initialise(Medium* medium, const bool verbose) {
  // Make sure the path to the Heed database is known.
  std::string databasePath;
  char* dbPath = std::getenv("HEED_DATABASE");
  if (dbPath) {
    databasePath = dbPath;
  } else {
    // Try GARFIELD_INSTALL.
    dbPath = std::getenv("GARFIELD_INSTALL");
    if (dbPath) {
      databasePath = std::string(dbPath) + "/share/Heed/database";
    } else {
      // Try GARFIELD_HOME.
      dbPath = std::getenv("GARFIELD_HOME");
      if (dbPath) {
        databasePath = std::string(dbPath) + "/Heed/heed++/database";
      } else {
      }
    }
  }
  if (databasePath.empty()) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Cannot retrieve database path (none of the"
              << " environment variables HEED_DATABASE, GARFIELD_INSTALL," 
              << " GARFIELD_HOME is defined).\n"
              << "    Cannot proceed.\n";
    return false;
  }
  if (databasePath.back() != '/') databasePath.append("/");

  if (m_debug || verbose) {
    std::cout << m_className << "::Initialise:\n"
              << "    Database path: " << databasePath << "\n";
  }

  // Check once more that the medium exists.
  if (!medium) {
    std::cerr << m_className << "::Initialise: Null pointer.\n";
    return false;
  }

  // Setup the energy mesh.
  m_energyMesh.reset(new Heed::EnergyMesh(m_emin, m_emax, m_nEnergyIntervals));

  if (medium->IsGas()) {
    if (!SetupGas(medium)) return false;
  } else {
    if (!SetupMaterial(medium)) return false;
  }

  // Energy transfer cross-section
  // Set a flag indicating whether the primary particle is an electron.
  m_transferCs.reset(new Heed::EnTransfCS(1.e-6 * m_mass, GetGamma() - 1.,
                                          m_isElectron, m_matter.get(),
                                          long(m_q)));
  if (!m_transferCs->m_ok) {
    std::cerr << m_className << "::Initialise:\n"
              << "    Problems occured when calculating the differential"
              << " cross-section table.\n"
              << "    Results will be inaccurate.\n";
  }
  if (!SetupDelta(databasePath)) return false;

  if (m_debug || verbose) {
    const double nc = m_transferCs->quanC;
    const double dedx = m_transferCs->meanC * 1.e3;
    const double dedx1 = m_transferCs->meanC1 * 1.e3;
    const double w = m_matter->W * 1.e6;
    const double f = m_matter->F;
    const double minI = m_matter->min_ioniz_pot * 1.e6;
    std::cout << m_className << "::Initialise:\n";
    std::cout << "    Cluster density:             " << nc << " cm-1\n";
    std::cout << "    Stopping power (restricted): " << dedx << " keV/cm\n";
    std::cout << "    Stopping power (incl. tail): " << dedx1 << " keV/cm\n";
    std::cout << "    W value:                     " << w << " eV\n";
    std::cout << "    Fano factor:                 " << f << "\n";
    std::cout << "    Min. ionization potential:   " << minI << " eV\n";
  }

  Heed::fixsyscoor primSys(Heed::point(0., 0., 0.), Heed::basis("primary"),
                           "primary");
  m_chamber.reset(new HeedChamber(primSys, m_lX, m_lY, m_lZ,
                                  *m_transferCs.get(), *m_deltaCs.get()));
  m_fieldMap->SetSensor(m_sensor);
  return true;
}

bool TrackHeed::SetupGas(Medium* medium) {
  // Get temperature and pressure.
  double pressure = medium->GetPressure();
  pressure = (pressure / AtmosphericPressure) * Heed::CLHEP::atmosphere;
  double temperature = medium->GetTemperature();

  const unsigned int nComponents = medium->GetNumberOfComponents();
  if (nComponents < 1) {
    std::cerr << m_className << "::SetupGas:\n"
              << "    Gas " << medium->GetName() << " has zero constituents.\n";
    return false;
  }

  std::vector<Heed::MolecPhotoAbsCS> mpacs;
  std::vector<std::string> notations;
  std::vector<double> fractions;
  for (unsigned int i = 0; i < nComponents; ++i) {
    std::string gasname;
    double frac;
    medium->GetComponent(i, gasname, frac);
    if (gasname == "paraH2" || gasname == "orthoD2" ||
        gasname == "D2") {
      gasname = "H2";
    } else if (gasname == "He-3") {
      gasname = "He";
    } else if (gasname == "CD4") {
      gasname = "CH4";
    } else if (gasname == "iC4H10" || gasname == "nC4H10") {
      gasname = "C4H10";
    } else if (gasname == "neoC5H12" || gasname == "nC5H12") {
      gasname = "C5H12";
    } else if (gasname == "cC3H6") {
      gasname = "C3H6";
    } else if (gasname == "nC3H7OH") {
      gasname = "C3H7OH";
    }
    // Assemble the molecular photoabsorption cross-section.
    if (gasname == "CF4") {
      mpacs.emplace_back(makeMPACS("C for CF4", 1, "F", 4));
    } else if (gasname == "Ar") {
      mpacs.emplace_back(makeMPACS("Ar", 1, 26.4e-6));
    } else if (gasname == "He" || gasname == "He-3") {
      mpacs.emplace_back(makeMPACS("He", 1, 41.3e-6));
    } else if (gasname == "Ne") {
      mpacs.emplace_back(makeMPACS("Ne", 1, 35.4e-6));
    } else if (gasname == "Kr") {
      mpacs.emplace_back(makeMPACS("Kr", 1, 24.4e-6));
    } else if (gasname == "Xe") {
      mpacs.emplace_back(makeMPACS("Xe", 1, 22.1e-6));
    } else if (gasname == "CH4" || gasname == "CD4") {
      mpacs.emplace_back(makeMPACS("C for CH4", 1, "H for CH4", 4, 27.3e-6));
    } else if (gasname == "C2H6") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 2, "H for H2", 6, 25.0e-6));
    } else if (gasname == "C3H8") {
      mpacs.emplace_back(makeMPACS("C for CH4", 3, "H for H2", 8, 24.0e-6));
    } else if (gasname == "C4H10") {
      mpacs.emplace_back(makeMPACS("C for C4H10", 4, "H for H2", 10, 23.4e-6));
    } else if (gasname == "CO2") {
      mpacs.emplace_back(makeMPACS("C for CO2", 1, "O for CO2", 2, 33.0e-6));
    } else if (gasname == "C5H12") {
      mpacs.emplace_back(makeMPACS("C for C4H10", 5, "H for H2", 12, 23.2e-6));
    } else if (gasname == "Water" || gasname == "H2O") {
      mpacs.emplace_back(makeMPACS("H for H2", 2, "O", 1, 29.6e-6));
    } else if (gasname == "O2") {
      mpacs.emplace_back(makeMPACS("O", 2, 30.8e-6));
    } else if (gasname == "N2" || gasname == "N2 (Phelps)") {
      mpacs.emplace_back(makeMPACS("N", 2, 34.8e-6));
    } else if (gasname == "NO") {
      mpacs.emplace_back(makeMPACS("N", 1, "O", 1));
    } else if (gasname == "N2O") {
      mpacs.emplace_back(makeMPACS("N", 2, "O", 1, 34.8e-6));
    } else if (gasname == "C2H4") {
      mpacs.emplace_back(makeMPACS("C for C2H4", 2, "H for H2", 4, 25.8e-6));
    } else if (gasname == "C2H2") {
      mpacs.emplace_back(makeMPACS("C for CH4", 2, "H for H2", 2, 25.8e-6));
    } else if (gasname == "H2" || gasname == "D2") {
      mpacs.emplace_back(makeMPACS("H for H2", 2));
    } else if (gasname == "CO") {
      mpacs.emplace_back(makeMPACS("C for CO2", 1, "O", 1));
    } else if (gasname == "Methylal") {
      // W similar to C4H10
      mpacs.emplace_back(makeMPACS("O", 2, "C for Methylal", 3,
                                   "H for H2", 8, 10.0e-6 * 23.4 / 10.55));
    } else if (gasname == "DME") {
      mpacs.emplace_back(makeMPACS("C for Methylal", 2, "H for H2", 6, "O", 1));
    } else if (gasname == "C2F6") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 2, "F", 6));
    } else if (gasname == "SF6") {
      mpacs.emplace_back(makeMPACS("S", 1, "F", 6));
    } else if (gasname == "NH3") {
      mpacs.emplace_back(makeMPACS("N", 1, "H for NH4", 3, 26.6e-6));
    } else if (gasname == "C3H6" || gasname == "cC3H6") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 3, "H for H2", 6));
    } else if (gasname == "CH3OH") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 1, "H for H2", 4,
                                    "O", 1, 24.7e-6));
    } else if (gasname == "C2H5OH") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 2, "H for H2", 6,
                                    "O", 1, 24.8e-6));
    } else if (gasname == "C3H7OH") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 3, "H for H2", 8, "O", 1));
    } else if (gasname == "Cs") {
      mpacs.emplace_back(makeMPACS("Cs", 1));
    } else if (gasname == "F2") {
      mpacs.emplace_back(makeMPACS("F", 2));
    } else if (gasname == "CS2") {
      mpacs.emplace_back(makeMPACS("C for CO2", 1, "S", 2));
    } else if (gasname == "COS") {
      mpacs.emplace_back(makeMPACS("C for CO2", 1, "O", 1, "S", 1));
    } else if (gasname == "BF3") {
      mpacs.emplace_back(makeMPACS("B", 1, "F", 3));
    } else if (gasname == "C2HF5") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 2, "H for H2", 1, "F", 5));
    } else if (gasname == "C2H2F4") {
      mpacs.emplace_back(makeMPACS("C for C2H6", 2, "F", 4, "H for H2", 2));
    } else if (gasname == "CHF3") {
      mpacs.emplace_back(makeMPACS("C for CF4", 1, "H for H2", 1, "F", 3));
    } else if (gasname == "CF3Br") {
      mpacs.emplace_back(makeMPACS("C for CF4", 1, "F", 3, "Br", 1));
    } else if (gasname == "C3F8") {
      mpacs.emplace_back(makeMPACS("C for CF4", 3, "F", 8));
    } else if (gasname == "O3") {
      mpacs.emplace_back(makeMPACS("O", 3));
    } else if (gasname == "Hg") {
      mpacs.emplace_back(makeMPACS("Hg", 1));
    } else if (gasname == "H2S") {
      mpacs.emplace_back(makeMPACS("H for H2", 2, "S", 1));
    } else if (gasname == "GeH4") {
      mpacs.emplace_back(makeMPACS("Ge", 1, "H for H2", 4));
    } else if (gasname == "SiH4") {
      mpacs.emplace_back(makeMPACS("Si", 1, "H for H2", 4));
    } else {
      std::cerr << m_className << "::SetupGas:\n"
                << "    Photoabsorption cross-section for " 
                << gasname << " is not implemented.\n";
      return false;
    }
    notations.push_back(gasname);
    fractions.push_back(frac);
  }
  if (m_usePacsOutput) {
    std::ofstream pacsfile;
    pacsfile.open("heed_pacs.txt", std::ios::out);
    const int nValues = m_energyMesh->get_q();
    if (nValues > 0) {
      for (int i = 0; i < nValues; ++i) {
        double e = m_energyMesh->get_e(i);
        pacsfile << 1.e6 * e << "  ";
        for (unsigned int j = 0; j < nComponents; ++j) {
          pacsfile << mpacs[j].get_ACS(e) << "  " << mpacs[j].get_ICS(e) << "  ";
        }
        pacsfile << "\n";
      }
    }
    pacsfile.close();
  }

  const std::string gasname = medium->GetName();
  m_gas.reset(new Heed::GasDef(gasname, gasname, nComponents, notations,
                               fractions, pressure, temperature, -1.));

  const double w = std::max(medium->GetW() * 1.e-6, 0.);
  double f = medium->GetFanoFactor();
  if (f <= 0.) f = Heed::standard_factor_Fano;

  m_matter.reset(
      new Heed::HeedMatterDef(m_energyMesh.get(), m_gas.get(), mpacs, w, f));
  return true;
}

bool TrackHeed::SetupMaterial(Medium* medium) {
  // Get temperature and density.
  double temperature = medium->GetTemperature();
  const double density =
      medium->GetMassDensity() * Heed::CLHEP::gram / Heed::CLHEP::cm3;

  const unsigned int nComponents = medium->GetNumberOfComponents();
  std::vector<Heed::AtomPhotoAbsCS*> atPacs(nComponents, nullptr);
  std::vector<std::string> notations;
  std::vector<double> fractions;
  for (unsigned int i = 0; i < nComponents; ++i) {
    std::string materialName;
    double frac;
    medium->GetComponent(i, materialName, frac);
    if (materialName == "C") {
      if (medium->GetName() == "Diamond") {
        atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Diamond");
      } else {
        atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("C");
      }
    } else if (materialName == "Si") {
      atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Si crystal");
      // atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Si G4");
    } else if (materialName == "Ga") {
      atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Ga for GaAs");
    } else if (materialName == "Ge") {
      atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Ge crystal");
    } else if (materialName == "As") {
      atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("As for GaAs");
    } else if (materialName == "Cd") {
      atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Cd for CdTe");
    } else if (materialName == "Te") {
      atPacs[i] = Heed::PhotoAbsCSLib::getAPACS("Te for CdTe");
    } else {
      std::cerr << m_className << "::SetupMaterial:\n"
                << "    Photoabsorption cross-section for " << materialName
                << " is not implemented.\n";
      return false;
    }
    notations.push_back(materialName);
    fractions.push_back(frac);
  }
  if (m_usePacsOutput) {
    std::ofstream pacsfile;
    pacsfile.open("heed_pacs.txt", std::ios::out);
    const int nValues = m_energyMesh->get_q();
    if (nValues > 0) {
      for (int i = 0; i < nValues; ++i) {
        double e = m_energyMesh->get_e(i);
        pacsfile << 1.e6 * e << "  ";
        for (unsigned int j = 0; j < nComponents; ++j) {
          pacsfile << atPacs[j]->get_ACS(e) << "  " << atPacs[j]->get_ICS(e)
                   << "  ";
        }
        pacsfile << "\n";
      }
    }
    pacsfile.close();
  }
  const std::string materialName = medium->GetName();
  m_material.reset(new Heed::MatterDef(materialName, materialName, nComponents,
                                       notations, fractions, density,
                                       temperature));

  double w = medium->GetW() * 1.e-6;
  if (w < 0.) w = 0.;
  double f = medium->GetFanoFactor();
  if (f <= 0.) f = Heed::standard_factor_Fano;

  m_matter.reset(new Heed::HeedMatterDef(m_energyMesh.get(), m_material.get(),
                                         atPacs, w, f));

  return true;
}

bool TrackHeed::SetupDelta(const std::string& databasePath) {
  // Load elastic scattering data.
  std::string filename = databasePath + "cbdel.dat";
  m_elScat.reset(new Heed::ElElasticScat(filename));

  filename = databasePath + "elastic_disp.dat";
  m_lowSigma.reset(new Heed::ElElasticScatLowSigma(m_elScat.get(), filename));

  // Load data for calculation of ionization.
  // Get W value and Fano factor.
  const double w = m_matter->W * 1.e6;
  const double f = m_matter->F;
  filename = databasePath + "delta_path.dat";
  m_pairProd.reset(new Heed::PairProd(filename, w, f));

  m_deltaCs.reset(new Heed::HeedDeltaElectronCS(
      m_matter.get(), m_elScat.get(), m_lowSigma.get(), m_pairProd.get()));
  return true;
}

double TrackHeed::GetW() const { return m_matter->W * 1.e6; }
double TrackHeed::GetFanoFactor() const { return m_matter->F; }

double TrackHeed::GetPhotoAbsorptionCrossSection(const double en) const {

  if (!m_matter) return 0.;
  // Convert eV to MeV.
  const double e = 1.e-6 * en;
  double cs = 0.;
  const auto n = m_matter->apacs.size();
  for (size_t i = 0; i < n; ++i) {
    const double w = m_matter->matter->weight_quan(i);
    cs += m_matter->apacs[i]->get_ACS(e) * w;
  }
  // Convert Mbarn to cm-2.
  return cs * 1.e-18;
}

void TrackHeed::ClearParticleBank() {
  Heed::gparticle::reset_counter();
  ClearBank(m_particleBank);
  m_bankIterator = m_particleBank.end();
}

bool TrackHeed::IsInside(const double x, const double y, const double z) {
  // Check if the point is inside the drift area.
  if (!m_sensor->IsInArea(x, y, z)) return false;
  // Check if the point is inside a medium.
  Medium* medium = m_sensor->GetMedium(x, y, z);
  if (!medium) return false;
  // Make sure the medium has not changed.
  if (medium->GetName() != m_mediumName ||
      fabs(medium->GetMassDensity() - m_mediumDensity) > 1.e-9 ||
      !medium->IsIonisable()) {
    return false;
  }
  return true;
}

bool TrackHeed::UpdateBoundingBox(bool& update) {
  // Get the bounding box.
  double xmin = 0., ymin = 0., zmin = 0.;
  double xmax = 0., ymax = 0., zmax = 0.;
  if (!m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cerr << m_className << "::UpdateBoundingBox: Drift area is not set.\n";
    m_ready = false;
    return false;
  }
  // Check if the bounding box has changed.
  const double lx = fabs(xmax - xmin);
  const double ly = fabs(ymax - ymin);
  const double lz = fabs(zmax - zmin);
  if (m_debug) {
    std::cout << m_className << "::UpdateBoundingBox:\n"
              << "    Bounding box dimensions:\n"
              << "      x: " << lx << " cm\n"
              << "      y: " << ly << " cm\n"
              << "      z: " << lz << " cm\n";
  }
  if (fabs(lx - m_lX) > Small || fabs(ly - m_lY) > Small ||
      fabs(lz - m_lZ) > Small) {
    m_lX = lx;
    m_lY = ly;
    m_lZ = lz;
    m_isChanged = true;
    update = true;
    m_hasActiveTrack = false;
  }
  // Update the center of the bounding box.
  m_cX = (std::isinf(xmin) || std::isinf(xmax)) ? 0. : 0.5 * (xmin + xmax);
  m_cY = (std::isinf(ymin) || std::isinf(ymax)) ? 0. : 0.5 * (ymin + ymax);
  m_cZ = (std::isinf(zmin) || std::isinf(zmax)) ? 0. : 0.5 * (zmin + zmax);
  if (m_debug) {
    std::cout << m_className << "::UpdateBoundingBox:\n"
              << "    Center of bounding box:\n"
              << "      x: " << m_cX << " cm\n"
              << "      y: " << m_cY << " cm\n"
              << "      z: " << m_cZ << " cm\n";
  }

  m_fieldMap->SetSensor(m_sensor);
  m_fieldMap->SetCentre(m_cX, m_cY, m_cZ);

  return true;
}
}
