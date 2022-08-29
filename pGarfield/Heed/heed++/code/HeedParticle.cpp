#include <iomanip>
#include <numeric>
#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/random/ranluxint.h"
#include "wcpplib/random/pois.h"
#include "wcpplib/random/rnorm.h"
#include "wcpplib/math/kinem.h"
#include "wcpplib/math/tline.h"
#include "heed++/code/HeedParticle.h"
#include "heed++/code/HeedCluster.h"
#include "heed++/code/HeedPhoton.h"
#include "heed++/code/EnTransfCS.h"

// 2003-2008, I. Smirnov

namespace Heed {

using CLHEP::c_light;
using CLHEP::c_squared;
using CLHEP::cm;
using CLHEP::MeV;
using CLHEP::electron_mass_c2;

HeedParticle::HeedParticle(manip_absvol* primvol, const point& pt,
                           const vec& vel, vfloat ftime, particle_def* fpardef,
                           HeedFieldMap* fieldmap, 
                           const bool fcoulomb_scattering,
                           const bool floss_only,
                           const bool fprint_listing)
    : eparticle(primvol, pt, vel, ftime, fpardef, fieldmap),
      m_coulomb_scattering(fcoulomb_scattering),
      m_loss_only(floss_only),
      m_print_listing(fprint_listing),
      m_particle_number(s_counter++) {}

void HeedParticle::physics(std::vector<gparticle*>& secondaries) {
  mfunname("void HeedParticle::physics()");
  if (m_print_listing) {
    mcout << "HeedParticle::physics is started\n";
    Iprintn(mcout, m_currpos.prange);
  }
  // Get the step.
  if (m_currpos.prange <= 0.0) return;
  const double stp = m_currpos.prange / cm;
  const vec dir = unit_vec(m_currpos.pt - m_prevpos.pt);
  const double range = (m_currpos.pt - m_prevpos.pt).length();
  if (m_print_listing) Iprint3n(mcout, m_prevpos.pt, dir, range);
  // Get local volume.
  const absvol* av = m_currpos.volume();
  auto etcs = dynamic_cast<const EnTransfCS*>(av);
  if (!etcs) return;
  HeedMatterDef* hmd = etcs->hmd;
  MatterDef* matter = hmd->matter;
  EnergyMesh* emesh = hmd->energy_mesh;
  const double* aetemp = emesh->get_ae();
  PointCoorMesh<double, const double*> pcm(emesh->get_q() + 1, &(aetemp));
  basis tempbas(m_currpos.dir, "tempbas");
  // Particle mass and energy.
  const double mp = m_mass * c_squared;
  const double ep = mp + m_curr_ekin;
  // Electron mass.
  const double mt = electron_mass_c2;
  // Particle velocity.
  const double invSpeed = 1. / m_prevpos.speed;
  // Shorthand.
  const auto sampleTransfer = t_hisran_step_ar<double, std::vector<double>,
                                               PointCoorMesh<double, const double*> >;
  const long qa = matter->qatom();
  if (m_print_listing) Iprintn(mcout, qa);
  for (long na = 0; na < qa; ++na) {
    if (m_print_listing) Iprintn(mcout, na);
    const long qs = hmd->apacs[na]->get_qshell();
    for (long ns = 0; ns < qs; ++ns) {
      if (m_print_listing) Iprintn(mcout, ns);
      if (etcs->quan[na][ns] <= 0.0) continue;
      // Sample the number of collisions for this shell.
      const long qt = pois(etcs->quan[na][ns] * stp);
      if (m_print_listing) Iprintn(mcout, qt);
      if (qt <= 0) continue;
      for (long nt = 0; nt < qt; ++nt) {
        // Sample the energy transfer in this collision.
        const double r = sampleTransfer(pcm, etcs->fadda[na][ns], SRANLUX());
        // Convert to internal units.
        const double et = r * MeV;
        m_edep += et;
        if (m_print_listing) Iprint2n(mcout, nt, et);
        // Sample the position of the collision.
        const double arange = SRANLUX() * range;
        point pt = m_prevpos.pt + dir * arange;
        if (m_loss_only) continue;
        if (m_print_listing) mcout << "generating new cluster\n";
        if (m_store_clusters) {
          m_clusterBank.emplace_back(HeedCluster(et, pt, na, ns));
        }
        // Generate a virtual photon.
        double theta_p, theta_t;
        theta_two_part(ep, ep - et, mp, mt, theta_p, theta_t);
        vec vel;
        vel.random_conic_vec(fabs(theta_t));
        vel.down(&tempbas);  // direction is OK
        vel *= c_light;
        const double t = m_prevpos.time + arange * invSpeed;
        if (m_print_listing) mcout << "generating new virtual photon\n";
        HeedPhoton* hp = new HeedPhoton(m_currpos.tid.eid[0], pt, vel, t,
                                        m_particle_number, et, m_fieldMap);
        if (!hp->alive()) {
          delete hp;
          continue;
        }
        hp->m_photon_absorbed = true;
        hp->m_delta_generated = false;
        hp->m_na_absorbing = na;
        hp->m_ns_absorbing = ns;
        secondaries.push_back(hp);
      }
    }
  }

  if (m_edep >= m_curr_ekin) {
    // Accumulated energy loss exceeds the particle's kinetic energy.
    m_alive = false;
  } 

  if (m_coulomb_scattering) {
    if (hmd->radiation_length > 0.) {
      const double x = range / hmd->radiation_length;
      const double sigma = etcs->sigma_ms * sqrt(x);
      double theta = sigma * rnorm_improved();
      turn(cos(theta), sin(theta));
    }
  }
  if (m_print_listing) {
    Iprintn(mcout, m_edep);
    mcout << "Exiting HeedParticle::physics\n";
  }
}

void HeedParticle::physics_mrange(double& fmrange) {
  if (!m_coulomb_scattering) return;
  // Get local volume and convert it to a cross-section object.
  const absvol* av = m_currpos.volume();
  auto etcs = dynamic_cast<const EnTransfCS*>(av);
  if (!etcs) return;
  if (etcs->quanC > 0.) {
    // Make sure the step is smaller than the mean free path between 
    // ionising collisions.
    fmrange = std::min(fmrange, 0.1 / etcs->quanC);
  }
} 

void HeedParticle::print(std::ostream& file, int l) const {
  if (l < 0) return;
  file << "HeedParticle: particle_number=" << m_particle_number << " type=";
  if (!m_pardef) {
    file << "none";
  } else {
    file << m_pardef->notation;
  }
  file << "\n  edep=" << m_edep << "\n";
  if (l <= 1) return;
  mparticle::print(file, l - 1);
}
}
