#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackSimple.hh"

namespace Garfield {

TrackSimple::TrackSimple() : Track() { m_className = "TrackSimple"; }

void TrackSimple::SetClusterDensity(const double d) {
  if (d < Small) {
    std::cerr << m_className << "::SetClusterDensity:\n"
              << "    Cluster density (number of clusters per cm)"
              << " must be positive.\n";
    return;
  }

  m_mfp = 1. / d;
}

double TrackSimple::GetClusterDensity() { return 1. / m_mfp; }

void TrackSimple::SetStoppingPower(const double dedx) {
  if (dedx < Small) {
    std::cerr << m_className << "::SetStoppingPower:\n"
              << "    Stopping power (average energy loss [eV] per cm)"
              << " must be positive.\n";
    return;
  }

  m_eloss = dedx;
}

double TrackSimple::GetStoppingPower() { return m_eloss; }

bool TrackSimple::NewTrack(const double x0, const double y0, const double z0,
                           const double t0, const double dx0, const double dy0,
                           const double dz0) {
  // Check if a sensor has been defined
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    m_isReady = false;
    return false;
  }

  // Make sure we are inside a medium
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::NewTrack: No medium at initial position.\n";
    m_isReady = false;
    return false;
  }

  m_isReady = true;

  m_x = x0;
  m_y = y0;
  m_z = z0;
  m_t = t0;

  // Normalise the direction.
  const double d = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  if (d < Small) {
    // Choose random direction.
    RndmDirection(m_dx, m_dy, m_dz);
  } else {
    m_dx = dx0 / d;
    m_dy = dy0 / d;
    m_dz = dz0 / d;
  }
  return true;
}

bool TrackSimple::GetCluster(double& xcls, double& ycls, double& zcls,
                             double& tcls, int& n, double& e, double& extra) {
  extra = 0.;
  if (!m_isReady) return false;

  if (m_useEqualSpacing) {
    m_x += m_dx * m_mfp;
    m_y += m_dy * m_mfp;
    m_z += m_dz * m_mfp;
  } else {
    const double d = -m_mfp * log(RndmUniformPos());
    m_x += m_dx * d;
    m_y += m_dy * d;
    m_z += m_dz * d;
  }

  xcls = m_x;
  ycls = m_y;
  zcls = m_z;
  tcls = m_t;

  n = 1;
  e = m_eloss * m_mfp;

  Medium* medium = m_sensor->GetMedium(m_x, m_y, m_z);
  if (!medium) {
    m_isReady = false;
    if (m_debug) {
      std::cout << m_className << "::GetCluster: Particle left the medium.\n";
    }
    return false;
  }

  return true;
}
}
