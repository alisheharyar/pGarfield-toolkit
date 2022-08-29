#include <fstream>
#include <iostream>
#include <algorithm>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLatex.h>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewBase.hh"
#include "Garfield/TrackSrim.hh"

namespace {

double StepVavilov(const double rkappa) {
  double xlmin = -3.7;
  if (rkappa < 0.1) {
    xlmin = -2.7;
  } else if (rkappa < 1) {
    xlmin = -2.9;
  } else if (rkappa < 2) {
    xlmin = -3.0;
  } else if (rkappa < 3) {
    xlmin = -3.1;
  } else if (rkappa < 4) {
    xlmin = -3.2;
  } else if (rkappa < 5) {
    xlmin = -3.3;
  } else if (rkappa < 6) {
    xlmin = -3.4;
  } else if (rkappa < 7) {
    xlmin = -3.5;
  } else if (rkappa < 8) {
    xlmin = -3.6;
  }
  return xlmin;
}

double Interpolate(const double x, const std::vector<double>& xtab,
                   const std::vector<double>& ytab) {
  if (x < xtab[0]) {
    return ytab[0];
  } else if (x > xtab.back()) {
    return ytab.back();
  }
  return Garfield::Numerics::Divdif(ytab, xtab, xtab.size(), x, 2);
}

void PrintSettings(const std::string& hdr, const double de, const double step,
                   const double ekin, const double beta2, const double gamma,
                   const double edensity, const double qp, const double mp,
                   const double emax, const double xi, const double kappa) {
  std::cout << hdr << "Settings:\n"
            << "    dE = " << de << " MeV,\n"
            << "    step = " << step << " cm.\n"
            << "    Ekin = " << ekin << " MeV,\n"
            << "    beta2 = " << beta2 << ",\n"
            << "    gamma = " << gamma << ".\n"
            << "    electron density = " << edensity << " / cm3.\n"
            << "    Qpart = " << qp << ", mpart = " << 1.e-6 * mp << " MeV.\n"
            << "    Emax = " << emax << " MeV,\n"
            << "    xi = " << xi << " MeV,\n"
            << "    kappa = " << kappa << ".\n";
}
}

namespace Garfield {

TrackSrim::TrackSrim() : Track() { m_className = "TrackSrim"; }

bool TrackSrim::ReadFile(const std::string& file) {
  // SRMREA

  const std::string hdr = m_className + "::ReadFile:\n    ";
  // Open the material list.
  std::ifstream fsrim(file);
  if (!fsrim) {
    std::cerr << hdr << "Could not open SRIM file " << file
              << " for reading.\n    The file perhaps does not exist.\n";
    return false;
  }
  unsigned int nread = 0;

  // Read the header
  if (m_debug) {
    std::cout << hdr << "SRIM header records from file " << file << "\n";
  }
  constexpr size_t size = 100;
  char line[size];
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "SRIM version") != NULL) {
      if (m_debug) std::cout << "\t" << line << "\n";
    } else if (strstr(line, "Calc. date") != NULL) {
      if (m_debug) std::cout << "\t" << line << "\n";
    } else if (strstr(line, "Ion =") != NULL) {
      break;
    }
  }

  // Identify the ion
  char* token = NULL;
  token = strtok(line, " []=");
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  // Set the ion charge.
  m_qion = std::atof(token);
  m_chargeset = true;
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  // Set the ion mass (convert amu to eV).
  m_mion = std::atof(token) * AtomicMassUnitElectronVolt;

  // Find the target density
  if (!fsrim.getline(line, size, '\n')) {
    std::cerr << hdr << "Premature EOF looking for target density (line "
              << nread << ").\n";
    return false;
  }
  nread++;
  if (!fsrim.getline(line, size, '\n')) {
    std::cerr << hdr << "Premature EOF looking for target density (line "
              << nread << ").\n";
    return false;
  }
  nread++;
  const bool pre2013 = (strstr(line, "Target Density") != NULL); 
  token = strtok(line, " ");
  token = strtok(NULL, " ");
  token = strtok(NULL, " ");
  if (pre2013) token = strtok(NULL, " ");
  SetDensity(std::atof(token));

  // Check the stopping units
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "Stopping Units") == NULL) continue;
    if (strstr(line, "Stopping Units =  MeV / (mg/cm2)") != NULL ||
        strstr(line, "Stopping Units =  MeV/(mg/cm2)") != NULL) {
      if (m_debug) {
        std::cout << hdr << "Stopping units: MeV / (mg/cm2) as expected.\n";
      }
      break;
    }
    std::cerr << hdr << "Unknown stopping units. Aborting (line " << nread
              << ").\n";
    return false;
  }

  // Skip to the table
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "-----------") != NULL) break;
  }

  // Read the table line by line
  m_ekin.clear();
  m_emloss.clear();
  m_hdloss.clear();
  m_range.clear();
  m_transstraggle.clear();
  m_longstraggle.clear();
  unsigned int ntable = 0;
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "-----------") != NULL) break;
    // Energy
    token = strtok(line, " ");
    m_ekin.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "eV") == 0) {
      m_ekin[ntable] *= 1.0e-6;
    } else if (strcmp(token, "keV") == 0) {
      m_ekin[ntable] *= 1.0e-3;
    } else if (strcmp(token, "GeV") == 0) {
      m_ekin[ntable] *= 1.0e3;
    } else if (strcmp(token, "MeV") != 0) {
      std::cerr << hdr << "Unknown energy unit " << token << "; aborting\n";
      return false;
    }
    // EM loss
    token = strtok(NULL, " ");
    m_emloss.push_back(atof(token));
    // HD loss
    token = strtok(NULL, " ");
    m_hdloss.push_back(atof(token));
    // Projected range
    token = strtok(NULL, " ");
    m_range.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      m_range[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      m_range[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      m_range[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      m_range[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      m_range[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }
    // Longitudinal straggling
    token = strtok(NULL, " ");
    m_longstraggle.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      m_longstraggle[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      m_longstraggle[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      m_longstraggle[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      m_longstraggle[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      m_longstraggle[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }
    // Transverse straggling
    token = strtok(NULL, " ");
    m_transstraggle.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      m_transstraggle[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      m_transstraggle[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      m_transstraggle[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      m_transstraggle[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      m_transstraggle[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }

    // Increment table line counter
    ++ntable;
  }

  // Find the scaling factor and convert to MeV/cm
  double scale = -1.;
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "=============") != NULL) {
      break;
    } else if (strstr(line, "MeV / (mg/cm2)") != NULL ||
               strstr(line, "MeV/(mg/cm2)") != NULL) {
      token = strtok(line, " ");
      scale = std::atof(token);
    }
  }
  if (scale < 0) {
    std::cerr << hdr << "Did not find stopping unit scaling; aborting.\n";
    return false;
  }
  scale *= 1.e3;
  for (unsigned int i = 0; i < ntable; ++i) {
    m_emloss[i] *= scale;
    m_hdloss[i] *= scale;
  }

  // Seems to have worked
  if (m_debug) {
    std::cout << hdr << "Successfully read " << file << "(" << nread
              << " lines).\n";
  }
  return true;
}

void TrackSrim::Print() {
  std::cout << "TrackSrim::Print:\n    SRIM energy loss table\n\n"
            << "    Energy     EM Loss     HD loss       Range  "
            << "l straggle  t straggle\n"
            << "     [MeV]    [MeV/cm]    [MeV/cm]        [cm] "
            << "      [cm]        [cm]\n\n";
  const unsigned int nPoints = m_emloss.size();
  for (unsigned int i = 0; i < nPoints; ++i) {
    printf("%10g  %10g  %10g  %10g  %10g  %10g\n", m_ekin[i],
           m_emloss[i] * m_rho, m_hdloss[i] * m_rho, m_range[i],
           m_longstraggle[i], m_transstraggle[i]);
  }
  std::cout << "\n";
  if (m_work > 0.) {
    std::printf("    Work function:  %g eV\n", m_work);
  } else {
    std::cout << "    Work function:  automatic\n";
  }
  if (m_fset) {
    std::printf("    Fano factor:    %g\n", m_fano);
  } else {
    std::cout << "    Fano factor:    automatic\n";
  }
  printf("    Ion charge:     %g\n", m_qion);
  printf("    Mass:           %g MeV\n", 1.e-6 * m_mion);
  printf("    Density:        %g g/cm3\n", m_rho);
  if (m_a > 0. && m_z > 0.) {
    std::printf("    A, Z:           %g, %g\n", m_a, m_z);
  } else {
    std::cout << "    A, Z:           automatic\n";
  }
}

void TrackSrim::PlotEnergyLoss() {

  const unsigned int nPoints = m_ekin.size();
  std::vector<double> yE;
  std::vector<double> yH;
  std::vector<double> yT;
  for (unsigned int i = 0; i < nPoints; ++i) {
    const double em = m_emloss[i] * m_rho;
    const double hd = m_hdloss[i] * m_rho;
    yE.push_back(em);
    yH.push_back(hd);
    yT.push_back(em + hd);
  }
  const double xmin = *std::min_element(std::begin(m_ekin), std::end(m_ekin));
  const double xmax = *std::max_element(std::begin(m_ekin), std::end(m_ekin));
  const double ymax = *std::max_element(std::begin(yT), std::end(yT));  
  // Prepare a plot frame.
  const std::string name = ViewBase::FindUnusedCanvasName("cSRIM"); 
  TCanvas* celoss = new TCanvas(name.c_str(), "Energy loss");
  celoss->SetLogx();
  celoss->SetGridx();
  celoss->SetGridy();
  celoss->DrawFrame(xmin, 0., xmax, 1.05 * ymax, ";Ion energy [MeV];Energy loss [MeV/cm]");

  // Make a graph for the 3 curves to plot.
  TGraph gr;
  gr.SetLineStyle(kSolid);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(21);
  gr.SetLineColor(kBlue + 1);
  gr.SetMarkerColor(kBlue + 1);
  gr.DrawGraph(nPoints, m_ekin.data(), yE.data(), "plsame");

  gr.SetLineColor(kGreen + 2);
  gr.SetMarkerColor(kGreen + 2);
  gr.DrawGraph(nPoints, m_ekin.data(), yH.data(), "plsame");

  gr.SetLineColor(kOrange - 3);
  gr.SetMarkerColor(kOrange - 3);
  gr.DrawGraph(nPoints, m_ekin.data(), yT.data(), "plsame");

  TLatex label;
  double xLabel = 0.4 * xmax;
  double yLabel = 0.9 * ymax;
  label.SetTextColor(kBlue + 1);
  label.SetText(xLabel, yLabel, "EM energy loss");
  label.DrawLatex(xLabel, yLabel, "EM energy loss");
  yLabel -= 1.5 * label.GetYsize();
  label.SetTextColor(kGreen + 2);
  label.DrawLatex(xLabel, yLabel, "HD energy loss");
  yLabel -= 1.5 * label.GetYsize();
  label.SetTextColor(kOrange - 3);
  label.DrawLatex(xLabel, yLabel, "Total energy loss");
  celoss->Update();
}

void TrackSrim::PlotRange() {

  const double xmin = *std::min_element(std::begin(m_ekin), std::end(m_ekin));
  const double xmax = *std::max_element(std::begin(m_ekin), std::end(m_ekin));
  const double ymax = *std::max_element(std::begin(m_range), std::end(m_range));

  // Prepare a plot frame.
  const std::string name = ViewBase::FindUnusedCanvasName("cSRIM");
  TCanvas* crange = new TCanvas(name.c_str(), "Range");
  crange->SetLogx();
  crange->SetGridx();
  crange->SetGridy();
  crange->DrawFrame(xmin, 0., xmax, 1.05 * ymax, ";Ion energy [MeV];Projected range [cm]");
  // Make a graph.
  TGraph gr;
  gr.SetLineColor(kOrange - 3);
  gr.SetMarkerColor(kOrange - 3);
  gr.SetLineStyle(kSolid);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(21);
  gr.DrawGraph(m_ekin.size(), m_ekin.data(), m_range.data(), "plsame");
  crange->Update();
}

void TrackSrim::PlotStraggling() {

  const double xmin = *std::min_element(std::begin(m_ekin), std::end(m_ekin));
  const double xmax = *std::max_element(std::begin(m_ekin), std::end(m_ekin));
  const double ymax = std::max(*std::max_element(std::begin(m_longstraggle),
                                                 std::end(m_longstraggle)),
                                *std::max_element(std::begin(m_transstraggle),
                                                  std::end(m_transstraggle)));
  // Prepare a plot frame.
  const std::string name = ViewBase::FindUnusedCanvasName("cSRIM");
  TCanvas* cstraggle = new TCanvas(name.c_str(), "Straggling");
  cstraggle->SetLogx();
  cstraggle->SetGridx();
  cstraggle->SetGridy();
  cstraggle->DrawFrame(xmin, 0., xmax, 1.05 * ymax, ";Ion energy [MeV];Straggling [cm]");

  // Make a graph for the 2 curves to plot.
  const unsigned int nPoints = m_ekin.size();
  TGraph gr;
  gr.SetLineStyle(kSolid);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(21);

  gr.SetLineColor(kOrange - 3);
  gr.SetMarkerColor(kOrange - 3);
  gr.DrawGraph(nPoints, m_ekin.data(), m_longstraggle.data(), "plsame");

  gr.SetLineColor(kGreen + 2);
  gr.SetMarkerColor(kGreen + 2);
  gr.DrawGraph(nPoints, m_ekin.data(), m_transstraggle.data(), "plsame");

  TLatex label;
  double xLabel = 1.2 * xmin;
  double yLabel = 0.9 * ymax;
  label.SetTextColor(kOrange - 3);
  label.SetText(xLabel, yLabel, "Longitudinal");
  label.DrawLatex(xLabel, yLabel, "Longitudinal");
  yLabel -= 1.5 * label.GetYsize();
  label.SetTextColor(kGreen + 2);
  label.DrawLatex(xLabel, yLabel, "Transverse");
  cstraggle->Update();
}

double TrackSrim::DedxEM(const double e) const {
  return Interpolate(e, m_ekin, m_emloss);
}

double TrackSrim::DedxHD(const double e) const {
  return Interpolate(e, m_ekin, m_hdloss);
}

double TrackSrim::Xi(const double x, const double beta2, 
                     const double edens) const {

  constexpr double fconst = 1.e-6 * TwoPi * (
    FineStructureConstant * FineStructureConstant * HbarC * HbarC) / 
    ElectronMass;
  return fconst * m_qion * m_qion * edens * x / beta2;
}

bool TrackSrim::PreciseLoss(const double step, const double estart,
                            double& deem, double& dehd) const {
  // SRMRKS

  // Debugging
  if (m_debug) printf("    Integrating energy losses over %g cm.\n", step);
  // Precision aimed for.
  const double eps = 1.0e-2;
  // Number of intervals.
  unsigned int ndiv = 1;
  // Loop until precision achieved
  const unsigned int nMaxIter = 10;
  bool converged = false;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    double e4 = estart;
    double e2 = estart;
    deem = 0.;
    dehd = 0.;
    // Compute rk2 and rk4 over the number of sub-divisions
    const double s = m_rho * step / ndiv;
    for (unsigned int i = 0; i < ndiv; i++) {
      // rk2: initial point
      const double de21 = s * (DedxEM(e2) + DedxHD(e2));
      // Mid-way point
      const double em22 = s * DedxEM(e2 - 0.5 * de21);
      const double hd22 = s * DedxHD(e2 - 0.5 * de21);
      // Trace the rk2 energy
      e2 -= em22 + hd22;
      // rk4: initial point
      const double em41 = s * DedxEM(e4);
      const double hd41 = s * DedxHD(e4);
      const double de41 = em41 + hd41;
      // Mid-way point
      const double em42 = s * DedxEM(e4 - 0.5 * de41);
      const double hd42 = s * DedxHD(e4 - 0.5 * de41);
      const double de42 = em42 + hd42;
      // Second mid-point estimate
      const double em43 = s * DedxEM(e4 - 0.5 * de42);
      const double hd43 = s * DedxHD(e4 - 0.5 * de42);
      const double de43 = em43 + hd43;
      // End point estimate
      const double em44 = s * DedxEM(e4 - de43);
      const double hd44 = s * DedxHD(e4 - de43);
      const double de44 = em44 + hd44;
      // Store the energy loss terms (according to rk4)
      deem += (em41 + em44) / 6. + (em42 + em43) / 3.;
      dehd += (hd41 + hd44) / 6. + (hd42 + hd43) / 3.;
      // Store the new energy computed with rk4
      e4 -= (de41 + de44) / 6. + (de42 + de43) / 3.;
    }
    if (m_debug) {
      printf("    Iteration %u has %d division(s). Losses:\n", iter, ndiv);
      printf("      de4 = %12g, de2 = %12g MeV\n", estart - e2, estart - e4);
      printf("      em4 = %12g, hd4 = %12g MeV\n", deem, dehd);
    }
    // Compare the two estimates
    if (fabs(e2 - e4) > eps * (fabs(e2) + fabs(e4) + fabs(estart))) {
      // Repeat with twice the number of steps.
      ndiv *= 2;
    } else {
      converged = true;
      break;
    }
  }

  if (!converged) {
    std::cerr << m_className << "::PreciseLoss: "
              << "No convergence achieved integrating energy loss.\n";
  } else if (m_debug) {
    std::cout << "    Convergence at eps = " << eps << ".\n";
  }
  return converged;
}

bool TrackSrim::EstimateRange(const double ekin, const double step,
                              double& stpmax) const {
  // Find distance over which the ion just does not lose all its energy
  // ekin       : Kinetic energy [MeV]
  // step       : Step length as guessed [cm]
  // stpmax     : Maximum step
  // SRMDEZ

  const std::string hdr = m_className + "::EstimateRange: ";
  // Initial estimate
  stpmax = step;

  // Find the energy loss expected for the present step length.
  double st1 = step;
  double deem = 0., dehd = 0.;
  PreciseLoss(st1, ekin, deem, dehd);
  double de1 = deem + dehd;
  // Do nothing if this is ok
  if (de1 < ekin) {
    if (m_debug) std::cout << "    Initial step OK.\n";
    return true;
  }
  // Find a smaller step for which the energy loss is less than EKIN.
  double st2 = 0.5 * step;
  double de2 = de1;
  const unsigned int nMaxIter = 20;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    // See where we stand
    PreciseLoss(st2, ekin, deem, dehd);
    de2 = deem + dehd;
    // Below the kinetic energy: done
    if (de2 < ekin) break;
    // Not yet below the kinetic energy: new iteration.
    st1 = st2;
    de1 = de2;
    st2 *= 0.5;
  }
  if (de2 >= ekin) {
    std::cerr << hdr << "\n    Did not find a smaller step in " << nMaxIter
              << " iterations. Abandoned.\n";
    stpmax = 0.5 * (st1 + st2);
    return false;
  }
  if (m_debug) {
    printf("    Step 1 = %g cm, dE 1 = %g MeV\n", st1, de1 - ekin);
    printf("    Step 2 = %g cm, dE 2 = %g MeV\n", st2, de2 - ekin);
  }
  // Now perform a bisection
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    // Avoid division by zero.
    if (de2 == de1) {
      if (m_debug) {
        std::cerr << "    Bisection failed due to equal energy loss for "
                  << "two step sizes. Abandoned.\n";
      }
      stpmax = 0.5 * (st1 + st2);
      return false;
    }
    // Estimate step to give total energy loss.
    double st3;
    if ((fabs(de1 - ekin) < 0.01 * fabs(de2 - de1)) ||
        (fabs(de1 - ekin) > 0.99 * fabs(de2 - de1))) {
      st3 = 0.5 * (st1 + st2);
    } else {
      st3 = st1 - (st2 - st1) * (de1 - ekin) / (de2 - de1);
    }
    // See how well we are doing.
    PreciseLoss(st3, ekin, deem, dehd);
    const double de3 = deem + dehd;
    if (m_debug) {
      std::printf("    Step 1 = %g cm, dE 1 = %g MeV\n", st1, de1 - ekin);
      std::printf("    Step 2 = %g cm, dE 2 = %g MeV\n", st2, de2 - ekin);
      std::printf("    Step 3 = %g cm, dE 3 = %g MeV\n", st3, de3 - ekin);
    }
    //  Update the estimates above and below.
    if (de3 > ekin) {
      st1 = st3;
      de1 = de3;
    } else {
      st2 = st3;
      de2 = de3;
    }
    // See whether we've converged.
    if (fabs(de3 - ekin) < 1e-3 * (fabs(de3) + fabs(ekin)) ||
        fabs(st1 - st2) < 1e-3 * (fabs(st1) + fabs(st2))) {
      stpmax = st1 - (st2 - st1) * (de1 - ekin) / (de2 - de1);
      return true;
    }
  }
  if (m_debug) {
    std::cout << "    Bisection did not converge in " << nMaxIter
              << " steps. Abandoned.\n";
  }
  stpmax = st1 - (st2 - st1) * (de1 - ekin) / (de2 - de1);
  return false;
}

Medium* TrackSrim::GetMedium(const std::array<double, 3>& x) const {
  Medium* medium = m_sensor->GetMedium(x[0], x[1], x[2]);
  if (medium && medium->IsIonisable() && m_sensor->IsInArea(x[0], x[1], x[2])) {
    return medium;
  } 
  return nullptr;
}

bool TrackSrim::NewTrack(const double x0, const double y0, const double z0,
                         const double t0, const double dx0, const double dy0,
                         const double dz0) {
  // Generates electrons for a SRIM track
  // SRMGEN
  const std::string hdr = m_className + "::NewTrack: ";

  // Verify that a sensor has been set.
  if (!m_sensor) {
    std::cerr << hdr << "Sensor is not defined.\n";
    return false;
  }

  // Make sure the initial position is inside an ionisable medium.
  Medium* medium = GetMedium({x0, y0, z0});
  if (!medium) {
    std::cerr << hdr << "No valid medium at initial position.\n";
    return false;
  } 
  // Get the W value of the medium (unless it has been set explicitly
  // by the user).
  double w = m_work < Small ? medium->GetW() : m_work;
  if (w < Small) {
    std::cerr << hdr << "Work function not defined.\n";
    return false;
  }
  // Get the Fano factor
  double fano = m_fset ? m_fano : medium->GetFanoFactor();
  fano = std::max(fano, 0.); 

  // Compute the electron (number) density of the target.
  double edens = 0.;
  if (m_a > 0. && m_z > 0. && m_rho > 0.) {
    edens = m_z * m_rho / (AtomicMassUnit * m_a);
  } else {
    edens = medium->GetNumberDensity() * medium->GetAtomicNumber();
    if (m_rho > 0.) {
      edens *= m_rho / medium->GetMassDensity();
    }
  }
  if (edens < Small) {
    std::cerr << hdr << "Invalid target density.\n";
    return false;
  } 

  // Normalise and store the direction.
  const double normdir = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  std::array<double, 3> v = {dx0, dy0, dz0};
  if (normdir < Small) {
    if (m_debug) {
      std::cout << hdr << "Direction vector has zero norm.\n"
                << "    Initial direction is randomized.\n";
    }
    // Null vector. Sample the direction isotropically.
    RndmDirection(v[0], v[1], v[2]);
  } else {
    // Normalise the direction vector.
    v[0] /= normdir;
    v[1] /= normdir;
    v[2] /= normdir;
  }

  // Make sure all necessary parameters have been set.
  if (m_mion < Small) {
    std::cerr << hdr << "Particle mass not set.\n";
    return false;
  } else if (!m_chargeset) {
    std::cerr << hdr << "Particle charge not set.\n";
    return false;
  } else if (m_energy < Small) {
    std::cerr << hdr << "Initial particle energy not set.\n";
    return false;
  } 

  // Check the initial energy (in MeV).
  const double ekin0 = 1.e-6 * GetKineticEnergy();
  if (ekin0 < 1.e-14 * m_mion || ekin0 < 1.e-6 * w) {
    std::cerr << hdr << "Initial kinetic energy E = " << ekin0
              << " MeV such that beta2 = 0 or E << W; particle stopped.\n";
    return true;
  }

  // Get an upper limit for the track length.
  const double tracklength = 10 * Interpolate(ekin0, m_ekin, m_range);

  // Header of debugging output.
  if (m_debug) {
    std::cout << hdr << "Track generation with the following parameters:\n";
    const unsigned int nTable = m_ekin.size();
    printf("      Table size           %u\n", nTable);
    printf("      Particle kin. energy %g MeV\n", ekin0);
    printf("      Particle mass        %g MeV\n", 1.e-6 * m_mion);
    printf("      Particle charge      %g\n", m_qion);
    printf("      Work function        %g eV\n", w);
    printf("      Fano factor          %g\n", fano);
    printf("      Long. straggling:    %d\n", m_useLongStraggle);
    printf("      Trans. straggling:   %d\n", m_useTransStraggle);
    printf("      Cluster size         %d\n", m_nsize);
  }

  // Plot.
  if (m_viewer) PlotNewTrack(x0, y0, z0);
 
  // Reset the cluster count.
  m_currcluster = 0;
  m_clusters.clear();

  // Initial situation: starting position
  std::array<double, 3> x = {x0, y0, z0};
  double t = t0;
  // Store the energy [MeV].
  double ekin = ekin0;
  // Total distance covered
  double dsum = 0.0;
  // Pool of unused energy
  double epool = 0.0;

  // Loop generating clusters
  int iter = 0;
  while (iter < m_maxclusters || m_maxclusters < 0) {
    if (m_debug) printf("    Iteration %d. Ekin = %g MeV.\n", iter, ekin);
    // Work out the energy loss per cm at the start of the step.
    const double dedxem = DedxEM(ekin) * m_rho;
    const double dedxhd = DedxHD(ekin) * m_rho;
    // Find the step size for which we get approximately the 
    // requested number of clusters or cluster size.
    double step = 0.;
    double eloss = 0.;
    if (m_nsize > 0) {
      step = m_nsize * 1.e-6 * w / dedxem;
    } else {
      const double ncls = m_maxclusters > 0 ? 0.5 * m_maxclusters : 100;
      step = ekin0 / (ncls * (dedxem + dedxhd));
    }
    if (m_debug) printf("    Estimated step size: %g cm\n", step);
    bool finish = false;
    // Truncate if this step exceeds the length.
    if (dsum + step > tracklength) {
      step = tracklength - dsum;
      if (m_debug) printf("    Truncating step to %g cm.\n", step);
      finish = true;
    }
    // Make an accurate integration of the energy loss over the step.
    double deem = 0., dehd = 0.;
    PreciseLoss(step, ekin, deem, dehd);
    double stpmax = tracklength - dsum;
    // If the energy loss exceeds the particle energy, truncate the step.
    if (deem + dehd > ekin) {
      EstimateRange(ekin, step, stpmax);
      step = stpmax;
      PreciseLoss(step, ekin, deem, dehd);
      deem = ekin * deem / (dehd + deem);
      dehd = ekin - deem;
      finish = true;
      if (m_debug) std::cout << "    Track length reached.\n";
    }
    if (m_debug) {
      std::cout << "    Maximum step size set to " << stpmax << " cm.\n";
    }
    // Ensure that this is larger than the minimum modelable step size.
    double stpmin;
    if (!SmallestStep(ekin, edens, deem, step, stpmin)) {
      std::cerr << hdr << "Failure computing the minimum step size.\n"
                << "    Clustering abandoned.\n";
      return false;
    }
    if (stpmin > stpmax) {
      // No way to find a suitable step size: use fixed energy loss.
      if (m_debug) {
        std::cout << "    Min. step > max. step; depositing all energy.\n";
      }
      if (deem + dehd > ekin) {
        eloss = ekin - dehd;
      } else {
        eloss = deem;
      }
      finish = true;
    } else if (step < stpmin) {
      // If needed enlarge the step size.
      if (m_debug) std::cout << "    Enlarging step size.\n";
      step = stpmin;
      PreciseLoss(step, ekin, deem, dehd);
      if (deem + dehd > ekin) {
        if (m_debug) std::cout << "    Excess loss. Recomputing max. step.\n";
        EstimateRange(ekin, step, stpmax);
        step = stpmax;
        PreciseLoss(step, ekin, deem, dehd);
        deem = ekin * deem / (dehd + deem);
        dehd = ekin - deem;
        eloss = deem;
      } else {
        eloss = RndmEnergyLoss(ekin, deem, step, edens);
      }
    } else {
      // Draw an actual energy loss for such a step.
      if (m_debug) std::cout << "    Step size ok.\n";
      eloss = RndmEnergyLoss(ekin, deem, step, edens);
    }
    // Ensure we are neither below 0 nor above the total energy.
    if (eloss < 0) {
      if (m_debug) std::cout << "    Truncating negative energy loss.\n";
      eloss = 0;
    } else if (eloss > ekin - dehd) {
      if (m_debug) std::cout << "    Excess energy loss, using mean.\n";
      if (deem + dehd > ekin) {
        eloss = ekin - dehd;
        finish = true;
      } else {
        eloss = deem;
      }
    }
    if (m_debug) {
      std::cout << "    Step length = " << step << " cm.\n"
                << "    Mean loss =   " << deem << " MeV.\n"
                << "    Actual loss = " << eloss << " MeV.\n";
    }
    // Get the particle velocity.
    const double rk = 1.e6 * ekin / m_mion;
    const double gamma = 1. + rk;
    const double beta2 = rk > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rk;
    const double vmag = sqrt(beta2) * SpeedOfLight; 
    // Compute the timestep.
    const double dt = vmag > 0. ? step / vmag : 0.;
    // Compute the endpoint of this step.
    std::array<double, 3> x1 = x;
    std::array<double, 3> v1 = v;
    if (m_sensor->HasMagneticField()) {
      double bx = 0., by = 0., bz = 0.;
      int status = 0;
      m_sensor->MagneticField(x[0], x[1], x[2], bx, by, bz, status);
      const double qoverm = m_qion / m_mion;
      std::array<double, 3> d = StepBfield(dt, qoverm, vmag, bx, by, bz, v1);
      x1[0] = x[0] + d[0];
      x1[1] = x[1] + d[1]; 
      x1[2] = x[2] + d[2]; 
    } else {
      x1[0] = x[0] + step * v[0];
      x1[1] = x[1] + step * v[1];
      x1[2] = x[2] + step * v[2];
    }
    // Is the endpoint in an ionisable medium and inside the active area?
    if (!GetMedium(x1)) {
      if (!m_sensor->HasMagneticField()) {
        step = Terminate(x, v, step);
      } else {
        step = TerminateBfield(x, v, dt, vmag);
      }
      PreciseLoss(step, ekin, deem, dehd);
      if (deem + dehd > ekin) {
        if (m_debug) std::cout << "    Excess loss.\n";
        deem = ekin * deem / (dehd + deem);
        dehd = ekin - deem;
      }
      eloss = RndmEnergyLoss(ekin, deem, step, edens);
      if (eloss > ekin - dehd) eloss = deem;
      finish = true;
    }
    // Add a cluster.
    Cluster cluster;
    cluster.x = x[0];
    cluster.y = x[1];
    cluster.z = x[2];
    cluster.t = t;
    if (fano < Small) {
      // No fluctuations.
      cluster.electrons = int((eloss + epool) / (1.e-6 * w));
      cluster.ec = w * cluster.electrons;
    } else {
      double ecl = 1.e6 * (eloss + epool);
      cluster.electrons = 0;
      cluster.ec = 0.0;
      while (true) {
        // if (cluster.ec < 100) printf("ec = %g\n", cluster.ec);
        const double ernd1 = RndmHeedWF(w, fano);
        if (ernd1 > ecl) break;
        cluster.electrons++;
        cluster.ec += ernd1;
        ecl -= ernd1;
      }
      if (m_debug) {
        std::cout << "    EM + pool: " << 1.e6 * (eloss + epool)
                  << " eV, W: " << w
                  << " eV, E/w: " << (eloss + epool) / (1.e-6 * w)
                  << ", n: " << cluster.electrons << ".\n";
      }
    }
    cluster.kinetic = ekin;
    epool += eloss - 1.e-6 * cluster.ec;
    if (m_debug) {
      std::cout << "    Adding cluster " << m_clusters.size() << " at ("
                << cluster.x << ", " << cluster.y << ", " << cluster.z
                << "), e = " << cluster.ec << ", n = " << cluster.electrons
                << ".\n"
                << "    Pool = " << epool << " MeV.\n";
    }
    m_clusters.push_back(std::move(cluster));
    if (m_viewer) PlotCluster(x[0], x[1], x[2]);

    // Keep track of the length and energy
    dsum += step;
    ekin -= eloss + dehd;
    if (finish) {
      // Stop if the flag is raised.
      if (m_debug) std::cout << "    Finishing flag raised.\n";
      break;
    } else if (ekin < ekin0 * 1.e-9) {
      // No energy left.
      if (m_debug) std::cout << "    Energy exhausted.\n";
      break;
    }
    // Update time, position and direction.
    t += dt;
    x = x1;
    v = v1;
    // Get the projected range and straggling.
    const double prange = Interpolate(ekin, m_ekin, m_range);
    const double strlat = m_useTransStraggle ? 
        Interpolate(ekin, m_ekin, m_transstraggle) : 0.;
    const double strlon = m_useLongStraggle ? 
        Interpolate(ekin, m_ekin, m_longstraggle) : 0.;
    // Draw scattering distances
    const double scale = sqrt(step / prange);
    const double sigt1 = RndmGaussian(0., scale * strlat);
    const double sigt2 = RndmGaussian(0., scale * strlat);
    const double sigl = RndmGaussian(0., scale * strlon);
    if (m_debug) {
      printf("    Sigma longitudinal: %g cm\n", sigl);
      printf("    Sigma transverse: %g, %g cm\n", sigt1, sigt2);
    }
    // Rotation angles to bring z-axis in line
    double theta, phi;
    if (v[0] * v[0] + v[2] * v[2] <= 0) {
      if (v[1] < 0) {
        theta = -HalfPi;
      } else if (v[1] > 0) {
        theta = +HalfPi;
      } else {
        std::cerr << hdr << "Zero step length; clustering abandoned.\n";
        return false;
      }
      phi = 0;
    } else {
      phi = atan2(v[0], v[2]);
      theta = atan2(v[1], sqrt(v[0] * v[0] + v[2] * v[2]));
    }
    // Update the position.
    const double cp = cos(phi);
    const double ct = cos(theta);
    const double sp = sin(phi);
    const double st = sin(theta);
    x[0] += +cp * sigt1 - sp * st * sigt2 + sp * ct * sigl;
    x[1] += +ct * sigt2 + st * sigl;
    x[2] += -sp * sigt1 - cp * st * sigt2 + cp * ct * sigl;
    medium = GetMedium(x);
    if (!medium) break;
    // Update the W value (unless it is set explicitly).
    if (m_work < Small) w = medium->GetW();
    if (w < Small) {
      std::cerr << hdr << "W value in medium " << medium->GetName() 
                << " is not defined.\n";
      break;
    }
    // Next cluster
    iter++;
  }
  if (iter == m_maxclusters) {
    std::cerr << hdr << "Exceeded maximum number of clusters.\n";
  }
  return true;
  // finished generating
}

double TrackSrim::Terminate(const std::array<double, 3>& x0, 
                            const std::array<double, 3>& v0,
                            const double step0) const {
  double step = 0.;
  std::array<double, 3> x1 = x0;
  double s = step0;
  const double tol = std::max(1.e-6 * step0, 1.e-6); 
  while (s > tol) {
    s *= 0.5;
    const std::array<double, 3> x2 = {x1[0] + s * v0[0], x1[1] + s * v0[1],
                                      x1[2] + s * v0[2]};
    if (GetMedium(x2)) {
      x1 = x2;
      step += s;
    }
  }
  return step;
}

double TrackSrim::TerminateBfield(const std::array<double, 3>& x0,
                                  const std::array<double, 3>& v0,
                                  const double dt0, const double vmag) const {

  const double qoverm = m_qion / m_mion;
  double dt = 0.;
  std::array<double, 3> x1 = x0;
  std::array<double, 3> v1 = v0;
  double s = dt0;
  const double tol = std::max(1.e-6 * dt0, 1.e-6);
  while (s > tol) {
    s *= 0.5;
    double bx = 0., by = 0., bz = 0.;
    int status = 0;
    m_sensor->MagneticField(x1[0], x1[1], x1[2], bx, by, bz, status);
    std::array<double, 3> v2 = v1;
    std::array<double, 3> d = StepBfield(s, qoverm, vmag, bx, by, bz, v2);
    std::array<double, 3> x2 = {x1[0] + d[0], x1[1] + d[1], x1[2] + d[2]};
    if (GetMedium(x2)) {
      x1 = x2;
      v1 = v2;
      dt += s;
    }
  }
  return dt * vmag;
}

bool TrackSrim::SmallestStep(const double ekin, const double edens,
                             double de, double step, double& stpmin) {
  // Determines the smallest step size for which there is little
  // or no risk of finding negative energy fluctuations.
  // SRMMST

  const std::string hdr = m_className + "::SmallestStep: ";
  constexpr double expmax = 30;

  // By default, assume the step is right.
  stpmin = step;
  // Check correctness.
  if (ekin <= 0 || de <= 0 || step <= 0) {
    std::cerr << hdr << "Input parameters not valid.\n    Ekin = " << ekin
              << " MeV, dE = " << de << " MeV, step length = " << step
              << " cm.\n";
    return false;
  } else if (m_mion <= 0 || fabs(m_qion) <= 0) {
    std::cerr << hdr
              << "Track parameters not valid.\n    Mass = " << 1.e-6 * m_mion
              << " MeV, charge = " << m_qion << ".\n";
    return false;
  } else if (edens <= 0) {
    std::cerr << hdr << "Target parameters not valid.\n"
              << "    electron density = " << edens << " / cm3.\n";
    return false;
  }

  // Basic kinematic parameters
  const double rkin = 1.e6 * ekin / m_mion;
  const double gamma = 1. + rkin;
  const double beta2 = rkin > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rkin;

  // Compute the maximum energy transfer [MeV]
  const double rm = ElectronMass / m_mion;
  const double emax = 2 * ElectronMass * 1.e-6 * beta2 * gamma * gamma /
                      (1. + 2 * gamma * rm + rm * rm);
  // Compute the Rutherford term.
  double xi = Xi(step, beta2, edens);
  // Compute the scaling parameter
  double rkappa = xi / emax;
  // Step size and energy loss
  double denow = de;
  double stpnow = step;
  constexpr unsigned int nMaxIter = 10;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    bool retry = false;
    // Debugging output.
    if (m_debug) {
      PrintSettings(hdr, denow, stpnow, ekin, beta2, gamma, edens,
                    m_qion, m_mion, emax, xi, rkappa);
    }
    double xinew = xi;
    double rknew = rkappa;
    if (m_model <= 0 || m_model > 4) {
      // No fluctuations: any step is permitted
      stpmin = stpnow;
    } else if (m_model == 1) {
      // Landau distribution
      constexpr double xlmin = -3.;
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      if (m_debug) {
        std::cout << "    Landau distribution is imposed (kappa_min = "
                  << rklim << ", d_min = " << stpmin << " cm).\n";
      }
    } else if (m_model == 2) {
      // Vavilov distribution, ensure we're in range.
      const double xlmin = StepVavilov(rkappa);
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      xinew = Xi(stpmin, beta2, edens);
      rknew = xinew / emax;
      if (m_debug) {
        std::cout << "    Vavilov distribution is imposed (kappa_min = "
                  << rklim << ", d_min = " << stpmin << " cm, kappa_new = " 
                  << rknew << ", xi_new = " << xinew << " MeV).\n";
      }
      if (stpmin > stpnow * 1.1) {
        if (m_debug) std::cout << "    Step size increase. New pass.\n";
        retry = true;
      }
    } else if (m_model == 3) {
      // Gaussian model
      const double sigma2 = xi * emax * (1 - 0.5 * beta2);
      stpmin = stpnow * 16 * sigma2 / (denow * denow);
      if (m_debug) {
        const double sigmaMin2 = Xi(stpmin, beta2, edens) * emax * (1 - 0.5 * beta2);
        std::cout << "    Gaussian distribution is imposed, d_min = " 
                  << stpmin << " cm, sigma/mu (old) = "
                  << sqrt(sigma2) / de << ", sigma/mu (min) = " 
                  << sqrt(sigmaMin2) / (stpmin * denow / stpnow) << ".\n";
      }
    } else if (rkappa < 0.05) {
      // Combined model: for low kappa, use the Landau distribution.
      constexpr double xlmin = -3.;
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      xinew = Xi(stpmin, beta2, edens);
      rknew = xinew / emax;
      if (m_debug) {
        std::cout << "    Landau distribution automatic (kappa_min = " 
                  << rklim << ", d_min = " << stpmin << " cm).\n";
      }
      if (rknew > 0.05 || stpmin > stpnow * 1.1) {
        retry = true;
        if (m_debug) {
          std::cout << "    Model change or step increase. New pass.\n";
        }
      }
    } else if (rkappa < 5) {
      // For medium kappa, use the Vavilov distribution
      const double xlmin = StepVavilov(rkappa);
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      xinew = Xi(stpmin, beta2, edens);
      rknew = xinew / emax;
      if (m_debug) {
        std::cout << "    Vavilov distribution automatic (kappa_min = "
                  << rklim << ", d_min = " << stpmin << " cm, kappa_new = " 
                  << rknew << ", xi_new = " << xinew << " MeV).\n";
      }
      if (rknew > 5 || stpmin > stpnow * 1.1) {
        retry = true;
        if (m_debug) {
          std::cout << "    Model change or step increase. New pass.\n";
        }
      }
    } else {
      // And for large kappa, use the Gaussian values.
      const double sigma2 = xi * emax * (1 - 0.5 * beta2);
      stpmin = stpnow * 16 * sigma2 / (denow * denow);
      if (m_debug) {
        const double sigmaMin2 = Xi(stpmin, beta2, edens) * emax * (1 - 0.5 * beta2);
        std::cout << "    Gaussian distribution automatic, d_min = " 
                  << stpmin << " cm, sigma/mu (old) = "
                  << sqrt(sigma2) / de << ", sigma/mu (min) = " 
                  << sqrt(sigmaMin2) / (stpmin * denow / stpnow) << ".\n";
      }
    }
    // See whether we should do another pass.
    if (stpnow > stpmin) {
      if (m_debug) {
        std::cout << "    Step size ok, minimum: " << stpmin << " cm\n";
      }
      break;
    }
    if (!retry) {
      if (m_debug) {
        std::cerr << "    Step size must be increased to " << stpmin
                  << "cm.\n";
      }
      break;
    }
    // New iteration
    rkappa = rknew;
    xi = xinew;
    denow *= stpmin / stpnow;
    stpnow = stpmin;
    if (m_debug) std::cout << "    Iteration " << iter << "\n";
    if (iter == nMaxIter - 1) {
      // Need interation, but ran out of tries
      std::cerr << hdr << "No convergence reached on step size.\n";
    }
  }
  return true;
}

double TrackSrim::RndmEnergyLoss(const double ekin, const double de,
                                 const double step, const double edens) const {
  //   RNDDE  - Generates a random energy loss.
  //   VARIABLES : EKIN       : Kinetic energy [MeV]
  //            DE         : Mean energy loss over the step [MeV]
  //            STEP       : Step length [cm]
  //            BETA2      : Velocity-squared
  //            GAMMA      : Projectile gamma
  //            EMAX       : Maximum energy transfer per collision [MeV]
  //            XI         : Rutherford term [MeV]
  //            FCONST     : Proportionality constant
  //            EMASS      : Electron mass [MeV]

  const std::string hdr = "TrackSrim::RndmEnergyLoss: ";
  // Check correctness.
  if (ekin <= 0 || de <= 0 || step <= 0) {
    std::cerr << hdr << "Input parameters not valid.\n    Ekin = " << ekin
              << " MeV, dE = " << de << " MeV, step length = " << step
              << " cm.\n";
    return 0.;
  } else if (m_mion <= 0 || fabs(m_qion) <= 0) {
    std::cerr << hdr << "Track parameters not valid.\n    Mass = " 
              << m_mion << " MeV, charge = " << m_qion << ".\n";
    return 0.;
  } else if (edens <= 0.) {
    std::cerr << hdr << "Target parameters not valid.\n"
              << "    electron density = " << edens << " / cm3.\n";
    return 0.;
  }
  // Basic kinematic parameters
  const double rkin = 1.e6 * ekin / m_mion;
  const double gamma = 1. + rkin;
  const double beta2 = rkin > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rkin;

  // Compute maximum energy transfer
  const double rm = ElectronMass / m_mion;
  const double emax = 2 * ElectronMass * 1.e-6 * beta2 * gamma * gamma /
                      (1. + 2 * gamma * rm + rm * rm);
  // Compute the Rutherford term.
  const double xi = Xi(step, beta2, edens);
  // Compute the scaling parameter
  const double rkappa = xi / emax;
  // Debugging output.
  if (m_debug) {
    PrintSettings(hdr, de, step, ekin, beta2, gamma, edens, 
                  m_qion, m_mion, emax, xi, rkappa);
  }
  double rndde = de;
  if (m_model <= 0 || m_model > 4) {
    // No fluctuations.
    if (m_debug) std::cout << "    Fixed energy loss.\n";
  } else if (m_model == 1) {
    // Landau distribution
    if (m_debug) std::cout << "    Landau imposed.\n";
    const double xlmean = -(log(rkappa) + beta2 + 1. - Gamma);
    rndde += xi * (RndmLandau() - xlmean);
  } else if (m_model == 2) {
    // Vavilov distribution, ensure we are in range.
    if (m_debug) std::cout << "    Vavilov imposed.\n";
    if (rkappa > 0.01 && rkappa < 12) {
      const double xvav = RndmVavilov(rkappa, beta2);
      rndde += xi * (xvav + log(rkappa) + beta2 + (1 - Gamma));
    }
  } else if (m_model == 3) {
    // Gaussian model
    if (m_debug) std::cout << "    Gaussian imposed.\n";
    rndde += RndmGaussian(0., sqrt(xi * emax * (1 - 0.5 * beta2)));
  } else if (rkappa < 0.05) {
    // Combined model: for low kappa, use the landau distribution.
    if (m_debug) std::cout << "    Landau automatic.\n";
    const double xlmean = -(log(rkappa) + beta2 + (1 - Gamma));
    const double par[] = {0.50884,    1.26116, 0.0346688,  1.46314,
                          0.15088e-2, 1.00324, -0.13049e-3};
    const double xlmax = par[0] + par[1] * xlmean + par[2] * xlmean * xlmean +
                         par[6] * xlmean * xlmean * xlmean +
                         (par[3] + xlmean * par[4]) * exp(par[5] * xlmean);
    double xlan = RndmLandau();
    for (unsigned int iter = 0; iter < 100; ++iter) {
      if (xlan < xlmax) break;
      xlan = RndmLandau();
    }
    rndde += xi * (xlan - xlmean);
  } else if (rkappa < 5) {
    // For medium kappa, use the Vavilov distribution.
    if (m_debug) std::cout << "    Vavilov fast automatic.\n";
    const double xvav = RndmVavilov(rkappa, beta2);
    rndde += xi * (xvav + log(rkappa) + beta2 + (1 - Gamma));
  } else {
    // And for large kappa, use the Gaussian values.
    if (m_debug) std::cout << "    Gaussian automatic.\n";
    rndde = RndmGaussian(de, sqrt(xi * emax * (1 - 0.5 * beta2)));
  }
  // Debugging output
  if (m_debug) {
    std::cout << "    Energy loss generated = " << rndde << " MeV.\n";
  }
  return rndde;
}

bool TrackSrim::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& n, double& e, double& extra) {
  if (m_debug) {
    std::cout << m_className << "::GetCluster: Cluster " << m_currcluster
              << " of " << m_clusters.size() << "\n";
  }
  // Stop if we have exhausted the list of clusters.
  if (m_currcluster >= m_clusters.size()) return false;

  const auto& cluster = m_clusters[m_currcluster];
  xcls = cluster.x;
  ycls = cluster.y;
  zcls = cluster.z;
  tcls = cluster.t;

  n = cluster.electrons;
  e = cluster.ec;
  extra = cluster.kinetic;
  // Move to next cluster
  ++m_currcluster;
  return true;
}
}
