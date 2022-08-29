#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/Utilities.hh"
#include "Garfield/TrackTrim.hh"

namespace {

std::tuple<double, double, double, double> GetRotation(
    const double x, const double y, const double z) {
  double phi = 0., theta = 0.;
  const double t = sqrt(x * x + y * y);
  if (t > 0.) {
    phi = atan2(y, x);
    theta = atan2(z, t);
  } else {
    theta = z < 0. ? -Garfield::HalfPi : Garfield::HalfPi;
  }
  return std::make_tuple(cos(theta), sin(theta), cos(phi), sin(phi));
}

std::array<double, 3> Step(const std::array<float, 5>& p0, 
                           const std::array<float, 5>& p1,
                           const double ctheta, const double stheta,
                           const double cphi, const double sphi) { 
  const double x = p1[0] - p0[0];
  const double y = p1[1] - p0[1];
  const double z = p1[2] - p0[2];
  return {cphi * ctheta * x - sphi * y - cphi * stheta * z,
          sphi * ctheta * x + cphi * y - sphi * stheta * z,
          stheta * x + ctheta * z};
}

double Speed(const double ekin, const double mass) {
  const double rk = ekin / mass;
  const double gamma = 1. + rk;
  const double beta2 = rk > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rk;
  return sqrt(beta2) * Garfield::SpeedOfLight;
}

}

namespace Garfield {

TrackTrim::TrackTrim() : Track() { 
  m_className = "TrackTrim";
  m_q = 1.;
}

void TrackTrim::SetParticle(const std::string& /*particle*/) {
  std::cerr << m_className << "::SetParticle: Not applicable.\n";
}

bool TrackTrim::ReadFile(const std::string& filename, 
                         const unsigned int nIons, const unsigned int nSkip) {

  // TRMREE - Reads the TRIM EXYZ file.

  // Reset.
  m_ekin = 0.;
  m_ions.clear();
  m_ion = 0;
  m_clusters.clear();
  m_cluster = 0;

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << m_className << "::ReadFile:\n"
              << "    Unable to open the EXYZ file (" << filename << ").\n";
    return false;
  }

  constexpr double Angstrom = 1.e-8;
  unsigned int nRead = 0;
  unsigned int ionNumber = 0;
  bool header = true;
  double mass = 0.;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> dedx;
  std::vector<float> ekin;
  for (std::string line; std::getline(infile, line);) {
    if (line.find("------- ") != std::string::npos) {
      // Reached the end of the header.
      header = false;
      continue;
    } else if (header) {
      if (line.find("Ion Data: ") != std::string::npos) {
        // Read the next line.
        std::getline(infile, line);
        auto words = tokenize(line);
        if (words.size() >= 3) {
          m_particleName = words[0];
          mass = std::stod(words[1]);
          auto pos = words[2].find("keV");
          if (pos != std::string::npos) {
            m_ekin = 1.e3 * std::stod(words[2].substr(0, pos));
          }
        } 
      } 
      // Otherwise, skip the header.
      continue;
    } 
    auto words = tokenize(line);
    if (words.size() < 6) {
      std::cerr << m_className << "::ReadFile: Unexpected line:\n"
                << line << "\n";
      continue;
    }
    if (ionNumber != std::stoul(words[0])) {
      // New ion.
      if (ionNumber > 0) {
        if (nRead >= nSkip) AddIon(x, y, z, dedx, ekin);
        x.clear();
        y.clear();
        z.clear();
        dedx.clear(); 
        ekin.clear();
        ++nRead;
        // Stop if we are done reading the requested number of ions.
        if (nIons > 0 && m_ions.size() >= nIons) break;
      }
      ionNumber = std::stoi(words[0]);
    }
    if (nRead < nSkip) continue;
    // Convert coordinates to cm.
    x.push_back(std::stof(words[2]) * Angstrom);
    y.push_back(std::stof(words[3]) * Angstrom);
    z.push_back(std::stof(words[4]) * Angstrom);
    // Convert stopping power from eV/A to eV/cm.
    dedx.push_back(std::stof(words[5]) / Angstrom);
    // Convert ion energy from keV to eV.
    ekin.push_back(std::stof(words[1]) * 1.e3);
  }
  infile.close();
  AddIon(x, y, z, dedx, ekin);
  std::cout << m_className << "::ReadFile: Read energy vs position for " 
            << m_ions.size() << " ions.\n";
  if (m_ekin > 0. && mass > 0.) {
    std::cout << "    Initial kinetic energy set to " 
              << m_ekin * 1.e-3 << " keV. Mass number: " << mass << ".\n";
    m_mass = AtomicMassUnitElectronVolt * mass;
    SetKineticEnergy(m_ekin);
  }
  return true;
}

void TrackTrim::AddIon(const std::vector<float>& x,
                       const std::vector<float>& y,
                       const std::vector<float>& z,
                       const std::vector<float>& dedx, 
                       const std::vector<float>& ekin) {

  const size_t nPoints = x.size();
  if (nPoints < 2) return;
  std::vector<std::array<float, 5> > path;
  for (size_t i = 0; i < nPoints; ++i) {
    float eloss = 0.;
    if (i < nPoints - 1) {
      const float dx = x[i + 1] - x[i];
      const float dy = y[i + 1] - y[i];
      const float dz = z[i + 1] - z[i];
      const float step = sqrt(dx * dx + dy * dy + dz * dz);
      if (i == 0 && dedx[i] > 10. * dedx[i + 1]) {
        eloss = step * dedx[i + 1];
      } else { 
        eloss = step * dedx[i];
      }
      const float dekin = ekin[i] - ekin[i + 1];
      if (dekin > 0.) eloss = std::min(eloss, dekin); 
    }
    path.push_back({x[i], y[i], z[i], eloss, ekin[i]});
  }
  m_ions.push_back(std::move(path));
}

void TrackTrim::Print() {
  std::cout << m_className << "::Print:\n";
  if (m_ions.empty()) {
    std::cerr << "    No TRIM data present.\n";
    return;
  }
  std::cout << "    Projectile: " << m_particleName << ", "
            << m_ekin * 1.e-3 << " keV\n"
            << "    Number of tracks: " << m_ions.size() << "\n";
  if (m_work > 0.) {
    std::cout << "    Work function: " << m_work << " eV\n";
  } else {
    std::cout << "    Work function: Automatic\n";
  }
  if (m_fset) {
    std::cout << "    Fano factor: " << m_fano << "\n";
  } else {
    std::cout << "    Fano factor: Automatic\n";
  }
}

bool TrackTrim::NewTrack(const double x0, const double y0, const double z0,
    const double t0, const double dx0, const double dy0, const double dz0) {

  // TRMGEN - Generates TRIM clusters

  if (m_ions.empty()) {
    std::cerr << m_className << "::NewTrack: No TRIM data present.\n";
    return false;
  }

  if (m_ion >= m_ions.size()) {
    // Rewind.
    std::cout << m_className << "::NewTrack: Rewinding.\n";
    m_ion = 0;
  }

  // Verify that a sensor has been set.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    return false;
  }

  // Get the rotation parameters.
  double ctheta = 1.;
  double stheta = 0.;
  double cphi = 1.;
  double sphi = 0.;
  if (sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0) < Small) {
    if (m_debug) {
      std::cout << m_className << "::NewTrack: Randomizing initial direction.\n";
    }
    // Null vector. Sample the direction isotropically.
    ctheta = 2 * RndmUniform() - 1.;
    stheta = sqrt(1. - ctheta * ctheta);
    const double phi = TwoPi * RndmUniform();
    cphi = cos(phi);
    sphi = sin(phi);
  } else {
    std::tie(ctheta, stheta, cphi, sphi) = GetRotation(dx0, dy0, dz0); 
  }

  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::NewTrack: No medium at initial point.\n";
  }
  double w = m_work > 0. ? m_work : medium->GetW();
  // Warn if the W value is not defined.
  if (w < Small) {
    std::cerr << m_className << "::NewTrack: "
              << "Work function at initial point is not defined.\n";
  }
  double fano = m_fset ? m_fano : medium->GetFanoFactor();
  fano = std::max(fano, 0.);

  // Charge-over-mass ratio.
  const double qoverm = m_mass > 0. ? m_q / m_mass : 0.;

  // Plot.
  if (m_viewer) PlotNewTrack(x0, y0, z0);
 
  // Reset the cluster count.
  m_cluster = 0;
  m_clusters.clear();

  // Pool of unused energy
  double epool = 0.0;

  const double ekin0 = GetKineticEnergy();
  const auto& path = m_ions[m_ion];
  const size_t nPoints = path.size();
  std::array<double, 3> xp = {x0, y0, z0};
  double tp = t0;
  for (size_t i = 1; i < nPoints; ++i) {
    // Skip points with kinetic energy below the initial one set by the user.
    const double ekin = path[i][4];
    if (ekin > ekin0) continue;
    double eloss = path[i - 1][3];
    std::array<double, 3> d = Step(path[i - 1], path[i], ctheta, stheta,
                                   cphi, sphi);
    double dmag = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    unsigned int nSteps = 1;
    if (m_maxStepSize > 0. && dmag > m_maxStepSize) {
      nSteps = std::ceil(dmag / m_maxStepSize);
    }
    if (m_maxLossPerStep > 0. && eloss > nSteps * m_maxLossPerStep) {
      nSteps = std::ceil(eloss / m_maxLossPerStep);
    }
    if (nSteps > 1) {
      const double scale = 1. / nSteps;
      for (size_t j = 0; j < 3; ++j) d[j] *= scale;
      dmag *= scale;
      eloss *= scale;
    }
    // Compute the particle velocity.
    const double vmag = Speed(ekin, m_mass);
    // Compute the timestep.
    const double dt = vmag > 0. ? dmag / vmag : 0.;
    for (unsigned int j = 0; j < nSteps; ++j) {
      Cluster cluster;
      if (m_sensor->HasMagneticField()) {
        double bx = 0., by = 0., bz = 0.;
        int status = 0;
        m_sensor->MagneticField(xp[0], xp[1], xp[2], bx, by, bz, status);
        // Direction vector.
        std::array<double, 3> v = {d[0] / dmag, d[1] / dmag, d[2] / dmag};
        d = StepBfield(dt, qoverm, vmag, bx, by, bz, v);
        // Update the rotation angles.
        std::tie(ctheta, stheta, cphi, sphi) = GetRotation(v[0], v[1], v[2]);
      }
      cluster.x = {xp[0] + d[0], xp[1] + d[1], xp[2] + d[2]};
      cluster.t = tp + dt;
      xp = cluster.x;
      tp = cluster.t;
      // Is this point inside an ionisable medium?
      medium = m_sensor->GetMedium(cluster.x[0], cluster.x[1], cluster.x[2]);
      if (!medium || !medium->IsIonisable()) continue;
      if (m_work < Small) w = medium->GetW();
      if (w < Small) continue;
      if (fano < Small) {
        // No fluctuations.
        cluster.ne = int((eloss + epool) / w);
        cluster.ec = w * cluster.ne;
      } else {
        double ec = eloss + epool;
        cluster.ne = 0;
        cluster.ec = 0.0;
        while (true) {
          const double er = RndmHeedWF(w, fano);
          if (er > ec) break;
          cluster.ne++;
          cluster.ec += er;
          ec -= er;
        }
      }
      cluster.ekin = path[i][4];
      epool += eloss - cluster.ec;
      if (cluster.ne == 0) continue;
      m_clusters.push_back(std::move(cluster));
      if (m_viewer) PlotCluster(cluster.x[0], cluster.x[1], cluster.x[2]);
    } 
  }
  // Move to the next ion in the list.
  ++m_ion;
  return true;
}

bool TrackTrim::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& n, double& e, double& extra) {
  if (m_debug) {
    std::cout << m_className << "::GetCluster: Cluster " << m_cluster
              << " of " << m_clusters.size() << "\n";
  }
  // Stop if we have exhausted the list of clusters.
  if (m_cluster >= m_clusters.size()) return false;

  const auto& cluster = m_clusters[m_cluster];
  xcls = cluster.x[0];
  ycls = cluster.x[1];
  zcls = cluster.x[2];
  tcls = cluster.t;

  n = cluster.ne;
  e = cluster.ec;
  extra = cluster.ekin;
  // Move to the next cluster.
  ++m_cluster;
  return true;
}
}
