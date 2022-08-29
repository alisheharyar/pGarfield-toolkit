#include <iomanip>
#include <numeric>
#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/random/ranluxint.h"
#include "wcpplib/random/pois.h"
#include "wcpplib/math/kinem.h"
#include "wcpplib/math/tline.h"
#include "heed++/code/HeedParticle_BGM.h"
#include "heed++/code/HeedCluster.h"
#include "heed++/code/HeedPhoton.h"
#include "heed++/code/EnTransfCS_BGM.h"

// 2003-2008, I. Smirnov

namespace Heed {

using CLHEP::c_light;
using CLHEP::c_squared;
using CLHEP::cm;
using CLHEP::MeV;
using CLHEP::electron_mass_c2;

HeedParticle_BGM::HeedParticle_BGM(manip_absvol* primvol, const point& pt,
                                   const vec& vel, vfloat ftime,
                                   particle_def* fpardef,
                                   HeedFieldMap* fieldmap,
                                   const bool floss_only,
                                   const bool fprint_listing)
    : eparticle(primvol, pt, vel, ftime, fpardef, fieldmap),
      m_print_listing(fprint_listing),
      m_loss_only(floss_only),
      m_particle_number(s_counter++) {}

void HeedParticle_BGM::physics(std::vector<gparticle*>& secondaries) {
  mfunname("void HeedParticle_BGM::physics()");
  if (m_print_listing) {
    mcout << "HeedParticle_BGM::physics is started\n";
    Iprintn(mcout, m_currpos.prange);
  }
  // Get the step.
  if (m_currpos.prange <= 0.0) return;
  const double stp = m_currpos.prange / cm;
  const vec dir = unit_vec(m_currpos.pt - m_prevpos.pt);
  // This approximation ignores curvature
  const double range = (m_currpos.pt - m_prevpos.pt).length();
  if (m_print_listing) Iprint3n(mcout, m_prevpos.pt, dir, range);
  // Get local volume.
  const absvol* av = m_currpos.volume();
  auto etcs = dynamic_cast<const EnTransfCS_BGM*>(av);
  if (!etcs) return;
  HeedMatterDef* hmd = etcs->hmd;
  MatterDef* matter = hmd->matter;
  EnergyMesh* emesh = hmd->energy_mesh;
  const double* aetemp = hmd->energy_mesh->get_ae();
  PointCoorMesh<double, const double*> pcm_e(emesh->get_q() + 1, &(aetemp));
  // Particle mass, energy and momentum.
  const double mp = m_mass * c_squared;
  const double ep = mp + m_curr_ekin;
  const double bg = sqrt(m_curr_gamma_1 * (m_curr_gamma_1 + 2.0));
  // Electron mass.
  const double mt = electron_mass_c2;
  // Particle velocity.
  const double invSpeed = 1. / m_prevpos.speed;
  PointCoorMesh<double, std::vector<double> > pcm(etcs->mesh->q,
                                                  &(etcs->mesh->x));
  long n1, n2;
  double b1, b2;
  int s_ret = pcm.get_interval(bg, n1, b1, n2, b2);
  if (s_ret != 1) {
    mcerr << "ERROR in void HeedParticle_BGM::physics()\n";
    mcerr << "beta*gamma is outside range of cross-section table\n";
    std::streamsize old_prec = mcerr.precision(15);
    Iprint2n(mcerr, m_curr_gamma_1, bg);
    mcerr.precision(old_prec);
    Iprint2n(mcerr, n1, n2);
    Iprint2n(mcerr, b1, b2);
    Iprintn(mcerr, etcs->mesh);
    mcerr << "This particle is:\n";
    print(mcerr, 2);
    mcerr << "This volume is:\n";
    av->print(mcerr, 2);
    spexit(mcerr);
    return;
  }

  const double f2 = (bg - b1) * (b2 - b1);
  const double f1 = 1. - f2;
  const long qa = matter->qatom();
  if (m_print_listing) Iprintn(mcout, qa);
  basis tempbas(m_currpos.dir, "tempbas");
  // Shorthand.
  const auto sampleTransfer = t_hisran_step_ar<double, std::vector<double>,
                                               PointCoorMesh<double, const double*> >;
  for (long na = 0; na < qa; ++na) {
    if (m_print_listing) Iprintn(mcout, na);
    long qs = hmd->apacs[na]->get_qshell();
    for (long ns = 0; ns < qs; ++ns) {
      if (m_print_listing) Iprintn(mcout, ns);
      const double y1 = etcs->etcs_bgm[n1].quan[na][ns];
      const double y2 = etcs->etcs_bgm[n2].quan[na][ns];
      const double mean_pois = f1 * y1 + f2 * y2;
      if (mean_pois <= 0.) continue;
      const long qt = pois(mean_pois * stp);
      if (m_print_listing) Iprintn(mcout, qt);
      if (qt <= 0) continue;
      for (long nt = 0; nt < qt; nt++) {
        // Sample the energy transfer in this collision.
        const double rn = SRANLUX();
        const double r1 = sampleTransfer(pcm_e, etcs->etcs_bgm[n1].fadda[na][ns], rn);
        const double r2 = sampleTransfer(pcm_e, etcs->etcs_bgm[n2].fadda[na][ns], rn);
        const double r = f1 * r1 + f2 * r2;
        if (m_print_listing) {
          Iprintn(mcout, rn);
          Iprint3n(mcout, r1, r2, r);
        }
        // Convert to internal units.
        const double et = r * MeV;
        m_edep += et;
        if (m_print_listing) Iprint2n(mcout, nt, et);
        // Sample the position of the collision.
        const double arange = SRANLUX() * range;
        point pt = m_prevpos.pt + dir * arange;
        if (m_loss_only) continue;
        if (m_print_listing) mcout << "generating new cluster\n";
        m_clusterBank.push_back(HeedCluster(et, pt, na, ns));
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

  if (m_print_listing) {
    Iprintn(mcout, m_edep);
    mcout << "Exiting HeedParticle_BGM::physics\n";
  }
}

void HeedParticle_BGM::print(std::ostream& file, int l) const {
  if (l < 0) return;
  Ifile << "HeedParticle_BGM (l=" << l
        << "): particle_number=" << m_particle_number << " type=";
  if (!m_pardef) {
    file << "none";
  } else {
    file << m_pardef->notation;
  }
  file << std::endl;
  if (l == 1) return;
  mparticle::print(file, l - 1);
  Iprintn(mcout, m_edep);
}
}
