#include <iomanip>
#include <fstream>
#include "wcpplib/clhep_units/WSystemOfUnits.h"
#include "wcpplib/math/lorgamma.h"
#include "wcpplib/math/tline.h"
#include "heed++/code/EnTransfCS.h"
#include "heed++/code/HeedMatterDef.h"
/*
2003, I. Smirnov
*/

namespace Heed {

EnTransfCS::EnTransfCS(double fparticle_mass, double fgamma_1,
                       int fs_primary_electron, HeedMatterDef* fhmd,
                       long fparticle_charge)
    : particle_mass(fparticle_mass),
      particle_charge(fparticle_charge),
      beta(lorbeta(fgamma_1)),
      beta2(0.0),
      beta12(0.0),
      gamma(fgamma_1 + 1.0),
      gamma_1(fgamma_1),
      s_simple_form(1),
      s_primary_electron(fs_primary_electron),
      hmd(fhmd),
      quanC(0.0)
#ifndef EXCLUDE_MEAN
      ,
      meanC(0.0),
      meanC1(0.0),
      meaneleC(0.0),
      meaneleC1(0.0)
#endif
      {
  mfunnamep("EnTransfCS::EnTransfCS(double fparticle_mass, double fbeta, "
            "HeedMatterDef* fhmd)");
  beta2 = beta * beta;
  beta12 = 1.0 - beta2;
  particle_tkin = particle_mass * gamma_1;
  particle_ener = particle_mass * gamma;
  if (s_primary_electron == 1) {
    maximal_energy_trans = 0.5 * particle_tkin;
  } else {
    double rm2 = particle_mass * particle_mass;
    double rme = ELMAS;
    if (beta12 > 1.0e-10) {
      maximal_energy_trans =
          2.0 * rm2 * ELMAS * beta2 /
          ((rm2 + rme * rme + 2.0 * rme * gamma * particle_mass) * (beta12));
      if (maximal_energy_trans > particle_tkin) {
        maximal_energy_trans = particle_tkin;
      }
    } else {
      maximal_energy_trans = particle_tkin;
    }
  }
  long qe = hmd->energy_mesh->get_q();
  long ne;
  log1C.put_qel(qe, 0.0);
  log2C.put_qel(qe, 0.0);
  chereC.put_qel(qe, 0.0);
  chereCangle.put_qel(qe, 0.0);
  Rreser.put_qel(qe, 0.0);
#ifdef DEBUG_EnTransfCS
  treser.put_qel(qe, 0.0);
#endif
  addaC.put_qel(qe, 0.0);
#ifndef EXCLUDE_A_VALUES
  addaC_a.put_qel(qe, 0.0);
#endif

  long qa = hmd->matter->qatom();
  cher.put_qel(qa);
  frezer.put_qel(qa);
  adda.put_qel(qa);
  fadda.put_qel(qa);
  quan.put_qel(qa);
#ifndef EXCLUDE_A_VALUES
  cher_a.put_qel(qa);
  adda_a.put_qel(qa);
  fadda_a.put_qel(qa);
  quan_a.put_qel(qa);
#endif

#ifndef EXCLUDE_VAL_FADDA
  val_fadda.put_qel(qa);
#ifndef EXCLUDE_A_VALUES
  val_fadda_a.put_qel(qa);
#endif
#endif

#ifndef EXCLUDE_MEAN
  mean.put_qel(qa);
#ifndef EXCLUDE_A_VALUES
  mean_a.put_qel(qa);
#endif
#endif

  for (long na = 0; na < qa; na++) {
    long qs = hmd->apacs[na]->get_qshell();
    cher[na].put_qel(qs);
    frezer[na].put_qel(qs);
    adda[na].put_qel(qs);
    fadda[na].put_qel(qs);
    quan[na].put_qel(qs);
#ifndef EXCLUDE_A_VALUES
    cher_a[na].put_qel(qs);
    adda_a[na].put_qel(qs);
    fadda_a[na].put_qel(qs);
    quan_a[na].put_qel(qs);
#endif
#ifndef EXCLUDE_VAL_FADDA
    val_fadda[na].put_qel(qs);
#ifndef EXCLUDE_A_VALUES
    val_fadda_a[na].put_qel(qs);
#endif
#endif

#ifndef EXCLUDE_MEAN
    mean[na].put_qel(qs);
#ifndef EXCLUDE_A_VALUES
    mean_a[na].put_qel(qs);
#endif
#endif
    for (long ns = 0; ns < qs; ns++) {
      cher[na][ns].put_qel(qe, 0.0);
      frezer[na][ns].put_qel(qe, 0.0);
      adda[na][ns].put_qel(qe, 0.0);
      fadda[na][ns].put_qel(qe, 0.0);
#ifndef EXCLUDE_A_VALUES
      cher_a[na][ns].put_qel(qe, 0.0);
      adda_a[na][ns].put_qel(qe, 0.0);
      fadda_a[na][ns].put_qel(qe, 0.0);
#endif
    }
  }
  double r;
  for (ne = 0; ne < qe; ne++) {
    r = -hmd->epsi1[ne] + (1.0 + hmd->epsi1[ne]) * beta12;
    r = r * r + beta2 * beta2 * hmd->epsi2[ne] * hmd->epsi2[ne];
    r = 1.0 / sqrt(r);
    r = log(r);
    log1C[ne] = r;
  }
  for (ne = 0; ne < qe; ne++) {
    r = 2.0 * 0.511 * beta2 / hmd->energy_mesh->get_ec(ne);
    if (r > 0.0) {
      r = log(r);
    } else {
      r = 0.0;
    }
    log2C[ne] = r;
  }
  double r0, rr12, rr22, r1, r2, r3;
  double coefpa =
      double(particle_charge * particle_charge) / (FSCON * beta2 * M_PI);
  for (ne = 0; ne < qe; ne++) {
    r0 = 1.0 + hmd->epsi1[ne];
    r = -hmd->epsi1[ne] + r0 * beta12;
    rr12 = r0 * r0;
    rr22 = hmd->epsi2[ne] * hmd->epsi2[ne];
    r1 = (-r0 * r + beta2 * rr22) / (rr12 + rr22);
    r2 = hmd->epsi2[ne] * beta2 / r;
    r3 = atan(r2);
    if (r < 0) r3 = M_PI + r3;
    chereCangle[ne] = r3;
    chereC[ne] = (coefpa / hmd->eldens) * r1 * r3;
  }

  for (ne = 0; ne < qe; ne++) {
    double ec = hmd->energy_mesh->get_ec(ne);
    if (s_simple_form == 1) {
      if (s_primary_electron == 0) {
        Rreser[ne] =
            1.0 / (ec * ec) * (1.0 - beta2 * ec / maximal_energy_trans);
      } else {
        Rreser[ne] = 1.0 / (ec * ec);
      }
    } else {
      if (s_primary_electron == 0) {
        Rreser[ne] =
            1.0 / (ec * ec) * (1.0 - beta2 * ec / maximal_energy_trans +
                               ec * ec / (2.0 * particle_ener * particle_ener));
      } else {
        double delta = ec / particle_mass;
        double pg2 = gamma * gamma;
        Rreser[ne] =
            beta2 / (particle_mass * particle_mass) * 1.0 / (pg2 - 1.0) *
            (gamma_1 * gamma_1 * pg2 / (pow(delta * (gamma_1 - delta), 2.0)) -
             (2.0 * pg2 + 2.0 * gamma - 1.0) / (delta * (gamma_1 - delta)) +
             1.0);
      }
    }
  }

  double Z_mean = hmd->matter->Z_mean();
  for (long na = 0; na < qa; na++) {
    long qs = hmd->apacs[na]->get_qshell();
    double at_weight_quan = hmd->matter->weight_quan(na);
    for (long ns = 0; ns < qs; ns++) {
      DynLinArr<double>& acher = cher[na][ns];
#ifndef EXCLUDE_A_VALUES
      DynLinArr<double>& acher_a = cher_a[na][ns];
#endif
      DynLinArr<double>& afrezer = frezer[na][ns];

      for (ne = 0; ne < qe; ne++) {
        double e1 = hmd->energy_mesh->get_e(ne);
        double e2 = hmd->energy_mesh->get_e(ne + 1);
        double ics = 0.;
        if (s_use_mixture_thresholds == 1) {
          ics = hmd->apacs[na]->get_integral_TICS(
              ns, e1, e2, hmd->min_ioniz_pot) / (e2 - e1) * C1_MEV2_MBN;
        } else {
          ics = hmd->apacs[na]->get_integral_ICS(ns, e1, e2) / (e2 - e1) *
                C1_MEV2_MBN;
        }

#ifndef EXCLUDE_A_VALUES
        double acs = hmd->apacs[na]->get_integral_ACS(ns, e1, e2) / (e2 - e1) *
                     C1_MEV2_MBN;
#endif
        check_econd11a(ics, < 0,
                       "na=" << na << " ns=" << ns << " ne=" << ne << '\n',
                       mcerr);
        if (hmd->ACS[ne] > 0.0) {
          acher[ne] =
              chereC[ne] * at_weight_quan * ics / (hmd->ACS[ne] * C1_MEV2_MBN);
#ifndef EXCLUDE_A_VALUES
          acher_a[ne] =
              chereC[ne] * at_weight_quan * acs / (hmd->ACS[ne] * C1_MEV2_MBN);
#endif
        } else {
          acher[ne] = 0.0;
#ifndef EXCLUDE_A_VALUES
          acher_a[ne] = 0.0;
#endif
        }
      }
      // Calculate the integral.
      double s = 0.;
      for (ne = 0; ne < qe; ne++) {
        double e1 = hmd->energy_mesh->get_e(ne);
        double ec = hmd->energy_mesh->get_ec(ne);
        double e2 = hmd->energy_mesh->get_e(ne + 1);
        r = hmd->apacs[na]->get_integral_ACS(ns, e1, e2) * C1_MEV2_MBN *
            at_weight_quan;
        // Here it must be ACS to satisfy sum rule for rezerford
        check_econd11a(r, < 0.0, "na=" << na << " ns=" << ns << " na=" << na,
                       mcerr);
        if (ec > hmd->min_ioniz_pot && ec < maximal_energy_trans) {
          afrezer[ne] = (s + 0.5 * r) * coefpa * Rreser[ne] / Z_mean;
          check_econd11a(afrezer[ne], < 0,
                         "na=" << na << " ns=" << ns << " na=" << na, mcerr);
        } else {
          afrezer[ne] = 0.0;
        }
        s += r;
#ifdef DEBUG_EnTransfCS
        treser[ne] += afrezer[ne];
#endif
      }
    }
  }
  for (ne = 0; ne < qe; ++ne) {
    double s = 0.0;
#ifndef EXCLUDE_A_VALUES
    double s_a = 0.0;
#endif
    double e1 = hmd->energy_mesh->get_e(ne);
    double ec = hmd->energy_mesh->get_ec(ne);
    double e2 = hmd->energy_mesh->get_e(ne + 1);
    double sqepsi = pow((1 + hmd->epsi1[ne]), 2.0) + pow(hmd->epsi2[ne], 2.0);
    for (long na = 0; na < qa; na++) {
      double at_weight_quan = hmd->matter->weight_quan(na);
      long qs = hmd->apacs[na]->get_qshell();
      for (long ns = 0; ns < qs; ns++) {
        double ics;
        if (s_use_mixture_thresholds == 1) {
          ics = hmd->apacs[na]->get_integral_TICS(
              ns, e1, e2, hmd->min_ioniz_pot) / (e2 - e1) * C1_MEV2_MBN;
        } else {
          ics = hmd->apacs[na]->get_integral_ICS(ns, e1, e2) / (e2 - e1) *
                C1_MEV2_MBN;
        }
#ifndef EXCLUDE_A_VALUES
        double acs = hmd->apacs[na]->get_integral_ACS(ns, e1, e2) / (e2 - e1) *
                     C1_MEV2_MBN;
#endif
        double r1 =
            at_weight_quan * log1C[ne] * coefpa * ics / (ec * Z_mean * sqepsi);
        double r2 =
            at_weight_quan * log2C[ne] * coefpa * ics / (ec * Z_mean * sqepsi);
        double& r_adda = adda[na][ns][ne];
        double& r_frezer = frezer[na][ns][ne];
        r_adda = r1 + r2 + r_frezer;
        if (r_adda < 0.0) r_adda = 0.0;

#ifndef EXCLUDE_A_VALUES
        double r1_a =
            at_weight_quan * log1C[ne] * coefpa * acs / (ec * Z_mean * sqepsi);
        double r2_a =
            at_weight_quan * log2C[ne] * coefpa * acs / (ec * Z_mean * sqepsi);
        double& r_adda_a = adda_a[na][ns][ne];
        r_adda_a = r1_a + r2_a + frezer;
        if (r_adda_a < 0.0) r_adda_a = 0.0;
#endif
        if (ec > hmd->min_ioniz_pot) {
          r_adda += cher[na][ns][ne];
          if (r_adda < 0.0) {
            funnw.whdr(mcout);
            mcout << "negative adda\n";
            mcout << "na=" << na << " ns=" << ns << " ne=" << ne
                  << " adda[na][ns][ne] = " << adda[na][ns][ne] << '\n';
            adda[na][ns][ne] = 0.0;
          }
        }
#ifndef EXCLUDE_A_VALUES
        adda_a[na][ns][ne] += cher[na][ns][ne];
        check_econd11a(adda_a[na][ns][ne], < 0,
                       "na=" << na << " ns=" << ns << " na=" << na, mcerr);
#endif
        s += r_adda;
#ifndef EXCLUDE_A_VALUES
        s_a += r_adda_a;
#endif
      }
    }
    addaC[ne] = s;
#ifndef EXCLUDE_A_VALUES
    addaC_a[ne] = s_a;
#endif
  }

  const double* aetemp = hmd->energy_mesh->get_ae();
  PointCoorMesh<double, const double*> pcm_e(qe + 1, &(aetemp));
  double emin = hmd->energy_mesh->get_emin();
  double emax = hmd->energy_mesh->get_emax();

  quanC = t_integ_step_ar<double, DynLinArr<double>,
                          PointCoorMesh<double, const double*> >(
      pcm_e, addaC, emin, emax, 0) * hmd->xeldens;

#ifndef EXCLUDE_A_VALUES
  quanC_a = t_integ_step_ar<double, DynLinArr<double>,
                            PointCoorMesh<double, const double*> >(
      pcm_e, addaC_a, emin, emax, 0) * hmd->xeldens;
#endif

#ifndef EXCLUDE_MEAN
  meanC = t_integ_step_ar<double, DynLinArr<double>,
                          PointCoorMesh<double, const double*> >(
      pcm_e, addaC, emin, emax, 1) * hmd->xeldens;

#ifndef EXCLUDE_A_VALUES
  meanC_a = t_integ_step_ar<double, DynLinArr<double>,
                            PointCoorMesh<double, const double*> >(
      pcm_e, addaC_a, emin, emax, 1) * hmd->xeldens;
#endif

  meanCleft = meanC;  // meanCleft does not have sense  in this approach

  if (s_simple_form == 1) {
    if (s_primary_electron == 0) {
      meanC1 = meanC;
      if (maximal_energy_trans > hmd->energy_mesh->get_e(qe)) {
        double e1 = hmd->energy_mesh->get_e(qe);
        double e2 = maximal_energy_trans;
        meanC1 += double(particle_charge * particle_charge) * 2.0 * M_PI /
                  (FSCON * FSCON * ELMAS * beta2) * hmd->xeldens *
                  (log(e2 / e1) - beta2 / maximal_energy_trans * (e2 - e1));
      }
    } else {
      meanC1 = meanC;
      if (maximal_energy_trans > hmd->energy_mesh->get_e(qe)) {
        double e1 = hmd->energy_mesh->get_e(qe);
        double e2 = maximal_energy_trans;
        meanC1 +=
            double(particle_charge * particle_charge) * 2.0 * M_PI /
            (pow(FSCON, 2.0) * ELMAS * beta2) * hmd->xeldens * log(e2 / e1);
      }
    }
  } else {
    if (s_primary_electron == 0) {
      meanC1 = meanC;
      if (maximal_energy_trans > hmd->energy_mesh->get_e(qe)) {
        double e1 = hmd->energy_mesh->get_e(qe);
        double e2 = maximal_energy_trans;
        meanC1 += double(particle_charge * particle_charge) * 2.0 * M_PI /
                  (FSCON * FSCON * ELMAS * beta2) * hmd->xeldens *
                  (log(e2 / e1) - beta2 / maximal_energy_trans * (e2 - e1) +
                   (e2 * e2 - e1 * e1) / (4.0 * particle_ener * particle_ener));
      }
#ifndef EXCLUDE_A_VALUES
      meanC1_a = meanC_a;
      if (maximal_energy_trans > hmd->energy_mesh->get_e(qe)) {
        double e1 = hmd->energy_mesh->get_e(qe);
        double e2 = maximal_energy_trans;
        meanC1_a +=
            double(particle_charge * particle_charge) * 2.0 * M_PI /
            (pow(FSCON, 2.0) * ELMAS * beta2) * hmd->xeldens *
            (log(e2 / e1) - beta2 / maximal_energy_trans * (e2 - e1) +
             (e2 * e2 - e1 * e1) / (4.0 * particle_ener * particle_ener));
      }
#endif
    }
  }

  meaneleC = meanC / hmd->W;
  meaneleC1 = meanC1 / hmd->W;
#endif

  for (long na = 0; na < qa; na++) {
    long qs = hmd->apacs[na]->get_qshell();
    for (long ns = 0; ns < qs; ns++) {
      quan[na][ns] = t_integ_step_ar<double, DynLinArr<double>,
                                     PointCoorMesh<double, const double*> >(
          pcm_e, adda[na][ns], emin, emax, 0) * hmd->xeldens;
#ifndef EXCLUDE_A_VALUES
      quan_a[na][ns] = t_integ_step_ar<double, DynLinArr<double>,
                                       PointCoorMesh<double, const double*> >(
          pcm_e, adda_a[na][ns], emin, emax, 0) * hmd->xeldens;
#endif
#ifndef EXCLUDE_MEAN
      mean[na][ns] = t_integ_step_ar<double, DynLinArr<double>,
                                     PointCoorMesh<double, const double*> >(
          pcm_e, adda[na][ns], emin, emax, 1) * hmd->xeldens;
#ifndef EXCLUDE_A_VALUES
      mean_a[na][ns] = t_integ_step_ar<double, DynLinArr<double>,
                                       PointCoorMesh<double, const double*> >(
          pcm_e, adda_a[na][ns], emin, emax, 1) * hmd->xeldens;
#endif
#endif
    }
  }

  for (long na = 0; na < qa; na++) {
    long qs = hmd->apacs[na]->get_qshell();
    for (long ns = 0; ns < qs; ns++) {
      if (quan[na][ns] > 0.0)
#ifndef EXCLUDE_VAL_FADDA
        val_fadda[na][ns] =
#endif
            t_hispre_step_ar<double, DynLinArr<double>,
                             PointCoorMesh<double, const double*> >(
                pcm_e, adda[na][ns], fadda[na][ns]);

#ifndef EXCLUDE_A_VALUES
      if (quan_a[na][ns] > 0.0)
#ifndef EXCLUDE_VAL_FADDA
        val_fadda_a[na][ns] =
#endif
            t_hispre_step_ar<double, DynLinArr<double>,
                             PointCoorMesh<double, const double*> >(
                pcm_e, adda_a[na][ns], fadda_a[na][ns]);
#endif
    }
  }

  length_y0 = DynLinArr<double>(qe, 0.0);
  for (ne = 0; ne < qe; ne++) {
    double k0 = hmd->energy_mesh->get_ec(ne) / PLANKCLIGHT;
    double det_value = 1.0 / (gamma * gamma) - hmd->epsi1[ne] * beta2;
    if (det_value <= 0.0) {
      length_y0[ne] = 0.0;
    } else {
      length_y0[ne] = beta / k0 * 1.0 / sqrt(det_value);
    }
  }

#ifdef CLEAR_ARRAYS
  log1C.clear();
  log2C.clear();
  chereC.clear();
  chereCangle.clear();
  Rreser.clear();
#ifdef DEBUG_EnTransfCS
  treser.clear();
#endif
  std::ofstream dcsfile;
  dcsfile.open("dcs.txt", std::ios::out);
  dcsfile << "# energy [MeV] vs. diff. cs per electron [Mbarn / MeV]\n";
  for (int i = 0; i < qe; ++i) {
    dcsfile << hmd->energy_mesh->get_ec(i) << "  " << addaC[i] / C1_MEV2_MBN
            << "\n";
  }
  dcsfile.close();

  addaC.clear();
#ifndef EXCLUDE_A_VALUES
  addaC_a.clear();
#endif
  cher.clear();
#ifndef EXCLUDE_A_VALUES
  cher_a.clear();
#endif
  frezer.clear();
  adda.clear();
#ifndef EXCLUDE_A_VALUES
  adda_a.clear();
#endif
#ifndef EXCLUDE_A_VALUES
  fadda_a.clear();
#endif
#ifndef EXCLUDE_VAL_FADDA
  val_fadda.clear();
#ifndef EXCLUDE_A_VALUES
  val_fadda_a.clear();
#endif
#endif
#ifndef EXCLUDE_MEAN
  mean.clear();
#ifndef EXCLUDE_A_VALUES
  mean_a.clear();
#endif
#endif

#endif

}

void EnTransfCS::print(std::ostream& file, int l) const {
  if (l <= 0) return;
  Ifile << "EnTransfCS(l=" << l << "):\n";
  indn.n += 2;
  Ifile << "particle_mass=" << particle_mass
        << " particle_tkin=" << particle_tkin
        << " particle_ener=" << particle_ener
        << " particle_charge=" << particle_charge << std::endl;
  Ifile << "beta=" << beta << " beta2=" << beta2 << " beta12=" << beta12
        << " gamma=" << gamma << std::endl;
  Ifile << "maximal_energy_trans=" << maximal_energy_trans << std::endl;
  Ifile << "s_primary_electron=" << s_primary_electron << std::endl;
  Ifile << "hmd:\n";
  hmd->print(file, 1);
#ifndef EXCLUDE_MEAN
#ifndef EXCLUDE_A_VALUES
  Ifile << "quanC=" << quanC << " quanC_a=" << quanC_a << '\n';
  Ifile << "meanC=" << meanC << " meanC_a=" << meanC_a << '\n';
  Ifile << " meanCleft=" << meanCleft << " meaneleC=" << meaneleC << '\n';
  Ifile << "meanC1=" << meanC1 << " meanC1_a=" << meanC1_a << '\n';
#else
  Ifile << "quanC=" << quanC << '\n';
  Ifile << "meanC=" << meanC << '\n';
  Ifile << " meanCleft=" << meanCleft << " meaneleC=" << meaneleC << '\n';
  Ifile << "meanC1=" << meanC1 << '\n';
#endif
  Ifile << " meaneleC1=" << meaneleC1 << '\n';
#else
#ifndef EXCLUDE_A_VALUES
  Ifile << "quanC=" << quanC << " quanC_a=" << quanC_a << '\n';
#else
  Ifile << "quanC=" << quanC << '\n';
#endif
#endif
  if (l > 2) {
    long qe = hmd->energy_mesh->get_q();
    long ne;
    if (l > 4) {
#ifdef DEBUG_EnTransfCS
      Ifile << "       enerc,      log1C,      log2C,      chereC,     addaC, "
               "chereCangle     Rreser      treser    length_y0\n";
#else
      Ifile << "       enerc,      log1C,      log2C,      chereC,     addaC, "
               "chereCangle   Rreser   length_y0\n";
#endif
      for (ne = 0; ne < qe; ne++) {
        Ifile << std::setw(12) << hmd->energy_mesh->get_ec(ne) << std::setw(12)
              << log1C[ne] << std::setw(12) << log2C[ne] << std::setw(12)
              << chereC[ne] << std::setw(12) << addaC[ne] << std::setw(12)
              << chereCangle[ne] << std::setw(12) << Rreser[ne]
#ifdef DEBUG_EnTransfCS
              << std::setw(12) << treser[ne]
#endif
              << std::setw(12) << length_y0[ne] << '\n';
      }
    }
    if (l > 3) {
      long qa = hmd->matter->qatom();
      long na;
      Iprintn(file, hmd->matter->qatom());
      for (na = 0; na < qa; na++) {
        Iprintn(file, na);
        long qs = hmd->apacs[na]->get_qshell();
        long ns;
        Iprintn(file, hmd->apacs[na]->get_qshell());
        for (ns = 0; ns < qs; ns++) {
          Iprintn(file, ns);
#ifndef EXCLUDE_MEAN
          Ifile << "quan      =" << std::setw(13) << quan[na][ns]
                << " mean  =" << std::setw(13) << mean[na][ns] << '\n';
#ifndef EXCLUDE_A_VALUES
          Ifile << "quan_a    =" << std::setw(13) << quan_a[na][ns]
                << " mean_a=" << std::setw(13) << mean_a[na][ns] << '\n';
#endif
#else
          Ifile << "quan      =" << std::setw(13) << quan[na][ns] << '\n';
#ifndef EXCLUDE_A_VALUES
          Ifile << "quan_a    =" << std::setw(13) << quan_a[na][ns] << '\n';
#endif
#endif
#ifndef EXCLUDE_VAL_FADDA
          Ifile << "val_fadda=" << std::setw(13) << val_fadda[na][ns]
#ifndef EXCLUDE_A_VALUES
                << " val_fadda_a=" << std::setw(13) << val_fadda_a[na][ns]
#endif
                << '\n';
#endif
          if (l > 5) {
            Ifile << "   enerc,        cher,       cher_a,     frezer,   adda, "
                     "  adda_a,  fadda,   fadda_a\n";
            for (ne = 0; ne < qe; ne++) {
              Ifile << std::setw(12) << hmd->energy_mesh->get_ec(ne)
                  //    << std::setw(12) << flog1[na][ns][ne]
                  //    << std::setw(12) << flog2[na][ns][ne]
                    << std::setw(12) << cher[na][ns][ne]
#ifndef EXCLUDE_A_VALUES
                    << std::setw(12) << cher_a[na][ns][ne]
#endif
                  //    << std::setw(12) << rezer[na][ns][ne]
                    << std::setw(12) << frezer[na][ns][ne] << std::setw(12)
                    << adda[na][ns][ne]
#ifndef EXCLUDE_A_VALUES
                    << std::setw(12) << adda_a[na][ns][ne]
#endif
                    << std::setw(12) << fadda[na][ns][ne]
#ifndef EXCLUDE_A_VALUES
                    << std::setw(12) << fadda_a[na][ns][ne]
#endif
                    << '\n';
            }
          }
        }
      }
    }
  }
  indn.n -= 2;

}

std::ostream& operator<<(std::ostream& file, const EnTransfCSType& f) {
  mfunname("std::ostream& operator << (std::ostream& file, const "
           "EnTransfCSType& f)");
  if (f.etcs.get() == NULL) {
    Ifile << "EnTransfCSType: type is not initialized\n";
  } else {
    Ifile << "EnTransfCSType: =";
    f.etcs->print(file, 1);
  }
  return file;
}

}
