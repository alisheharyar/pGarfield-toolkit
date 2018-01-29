
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include <TAxis.h>

#include "Plotting.hh"
#include "Medium.hh"
#include "ViewMedium.hh"

namespace Garfield {

ViewMedium::ViewMedium()
    : m_debug(false),
      m_canvas(0),
      m_hasExternalCanvas(false),
      m_medium(0),
      m_eMin(0.),
      m_eMax(1000.),
      m_bMin(0.),
      m_bMax(5.),
      m_aMin(0.),
      m_aMax(3.14),
      m_vMin(0.),
      m_vMax(0.),
      m_efield(500.),
      m_bfield(1.e2),
      m_angle(0.),
      m_etolerance(1.),
      m_btolerance(0.01),
      m_atolerance(0.05),
      m_nFunctions(0),
      m_nGraphs(0) {

  m_className = "ViewMedium";

  m_functions.clear();
  m_graphs.clear();
  plottingEngine.SetDefaultStyle();
}

ViewMedium::~ViewMedium() {

  if (!m_hasExternalCanvas && m_canvas != 0) delete m_canvas;
}

void ViewMedium::SetCanvas(TCanvas* c) {

  if (c == 0) return;
  if (!m_hasExternalCanvas && m_canvas != 0) {
    delete m_canvas;
    m_canvas = 0;
  }
  m_canvas = c;
  m_hasExternalCanvas = true;
}

void ViewMedium::SetMedium(Medium* m) {

  if (m == 0) {
    std::cerr << m_className << "::SetMedium:\n";
    std::cerr << "    Medium pointer is null.\n";
    return;
  }

  m_medium = m;
}

void ViewMedium::SetElectricFieldRange(const double emin, const double emax) {

  if (emin >= emax || emin < 0.) {
    std::cerr << m_className << "::SetElectricFieldRange:\n";
    std::cerr << "    Incorrect field range.\n";
    return;
  }

  m_eMin = emin;
  m_eMax = emax;
}

void ViewMedium::SetMagneticFieldRange(const double bmin, const double bmax) {

  if (bmin >= bmax || bmin < 0.) {
    std::cerr << m_className << "::SetMagneticFieldRange:\n";
    std::cerr << "    Incorrect field range.\n";
    return;
  }

  m_bMin = bmin;
  m_bMax = bmax;
}

void ViewMedium::SetBAngleRange(const double amin, const double amax) {

  if (amin >= amax || amin < 0.) {
    std::cerr << m_className << "::SetBAngleRange:\n";
    std::cerr << "    Incorrect field range.\n";
    return;
  }

  m_aMin = amin;
  m_aMax = amax;
}

void ViewMedium::SetFunctionRange(const double vmin, const double vmax) {

  if (vmin >= vmax || vmin < 0.) {
    std::cerr << m_className << "::SetFunctionRange:\n";
    std::cerr << "    Incorrect range.\n";
    return;
  }

  m_vMin = vmin;
  m_vMax = vmax;
}

void ViewMedium::SetFunctionRange() { m_vMin = m_vMax = 0.; }

void ViewMedium::PlotElectronVelocity(const char xaxis, const double e,
                                      const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  double min = 0., max = 0.;
  std::string title = "";
  if (xaxis == 'e') {
    title = "electric field [V/cm]";
    min = m_eMin;
    max = m_eMax;
  } else if (xaxis == 'b') {
    title = "magnetic field [T]";
    min = m_bMin;
    max = m_bMax;
  } else if (xaxis == 'a') {
    title = "magnetic field angle [rad]";
    min = m_aMin;
    max = m_aMax;
  }
  AddFunction(min, max, m_vMin, m_vMax, keep, title, "drift velocity [cm/ns]",
              0, xaxis, e, b, a);
  keep = true;
  AddFunction(min, max, m_vMin, m_vMax, keep, title, "drift velocity [cm/ns]",
              23, xaxis, e, b, a);
  keep = true;
  AddFunction(min, max, m_vMin, m_vMax, keep, title, "drift velocity [cm/ns]",
              24, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotHoleVelocity(const char xaxis, const double e,
                                  const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  double min = 0., max = 0.;
  std::string title = "";
  if (xaxis == 'e') {
    title = "electric field [V/cm]";
    min = m_eMin;
    max = m_eMax;
  } else if (xaxis == 'b') {
    title = "magnetic field [T]";
    min = m_bMin;
    max = m_bMax;
  } else if (xaxis == 'a') {
    title = "magnetic field angle [rad]";
    min = m_aMin;
    max = m_aMax;
  }
  AddFunction(min, max, m_vMin, m_vMax, keep, title,
              "drift velocity [cm/ns]", 10, xaxis, e, b, a);
  keep = true;
  AddFunction(min, max, m_vMin, m_vMax, keep, title,
              "drift velocity [cm/ns]", 25, xaxis, e, b, a);
  keep = true;
  AddFunction(min, max, m_vMin, m_vMax, keep, title,
              "drift velocity [cm/ns]", 26, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotIonVelocity(const char xaxis, const double e,
                                 const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "drift velocity [cm/ns]", 20, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotElectronDiffusion(const char xaxis, const double e,
                                       const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "diffusion coefficient [#sqrt{cm}]", 1, xaxis, e, b, a);
  keep = true;
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "diffusion coefficient [#sqrt{cm}]", 2, xaxis, e, b, a);

  m_canvas->Update();
}

void ViewMedium::PlotHoleDiffusion(const char xaxis, const double e,
                                   const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "diffusion coefficient [#sqrt{cm}]", 11, xaxis, e, b, a);
  keep = true;
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "diffusion coefficient [#sqrt{cm}]", 12, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotIonDiffusion(const char xaxis, const double e,
                                  const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "diffusion coefficient [#sqrt{cm}]", 21, xaxis, e, b, a);
  keep = true;
  AddFunction(m_eMin, m_eMax, m_vMin, m_vMax, keep, "electric field [V/cm]",
              "diffusion coefficient [#sqrt{cm}]", 22, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotElectronTownsend(const char xaxis, const double e,
                                      const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, 0., 0., keep, "electric field [V/cm]",
              "Townsend coefficient [1/cm]", 3, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotHoleTownsend(const char xaxis, const double e,
                                  const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, 0., 0., keep, "electric field [V/cm]",
              "Townsend coefficient [1/cm]", 13, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotElectronAttachment(const char xaxis, const double e,
                                        const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, 0., 0., keep, "electric field [V/cm]",
              "Attachment coefficient [1/cm]", 4, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotHoleAttachment(const char xaxis, const double e,
                                    const double b, const double a) {

  bool keep = false;
  SetupCanvas();
  AddFunction(m_eMin, m_eMax, 0., 0., keep, "electric field [V/cm]",
              "Attachment coefficient [1/cm]", 14, xaxis, e, b, a);
  m_canvas->Update();
}

void ViewMedium::PlotElectronCrossSections() {

  std::cerr << m_className << "::PlotElectronCrossSections:\n";
  std::cerr << "    Function not yet implemented.\n";
  SetupCanvas();
}

void ViewMedium::SetupCanvas() {

  if (m_canvas == 0) {
    m_canvas = new TCanvas();
    m_canvas->SetTitle("Medium View");
    if (m_hasExternalCanvas) m_hasExternalCanvas = false;
  }
  m_canvas->cd();
  gPad->SetLeftMargin(0.15);
}

void ViewMedium::AddFunction(const double xmin, const double xmax,
                             const double ymin, const double ymax,
                             const bool keep, const std::string xlabel,
                             const std::string ylabel, const int type,
                             const char xaxis, const double e, const double b,
                             const double a) {

  // Make sure the medium pointer is set.
  if (m_medium == 0) {
    std::cerr << m_className << "::AddFunction:\n";
    std::cerr << "    Medium is not defined.\n";
    return;
  }

  // Look for an unused function name.
  int idx = 0;
  std::string fname = "fMediumView_0";
  while (gROOT->GetListOfFunctions()->FindObject(fname.c_str())) {
    ++idx;
    std::stringstream ss;
    ss << "fMediumView_";
    ss << idx;
    fname = ss.str();
  }
  if (m_debug) {
    std::cout << m_className << "::AddFunction:\n";
    std::cout << "    Adding function " << fname << "\n";
  }
  if (!keep) {
    m_functions.clear();
    m_nFunctions = 0;
    m_graphs.clear();
    m_nGraphs = 0;
  }

  // Set global variables used later on (effectively we will only use two)
  m_efield = e;
  m_bfield = b;
  m_angle = a;

  // Create a TF1 and add it to the list of functions.
  TF1 fNew(fname.c_str(), this, &ViewMedium::EvaluateFunction, xmin, xmax, 2,
           "ViewMedium", "EvaluateFunction");
  fNew.SetNpx(1000);
  m_functions.push_back(fNew);
  ++m_nFunctions;

  const std::string title = m_medium->GetName() + ";" + xlabel + ";" + ylabel;
  m_functions.back().SetRange(xmin, xmax);
  if ((fabs(ymax - ymin) > 0.)) {
    m_functions.back().SetMinimum(ymin);
    m_functions.back().SetMaximum(ymax);
  }
  m_functions.back().GetXaxis()->SetTitle(xlabel.c_str());
  m_functions.back().GetXaxis()->SetTitleOffset(1.2);
  m_functions.back().GetYaxis()->SetTitle(ylabel.c_str());
  m_functions.back().SetTitle(title.c_str());
  m_functions.back().SetParameter(0, type);
  m_functions.back().SetParameter(1, xaxis);
  // Set the color and marker style.
  int color;
  int marker = 20;
  if (type == 2 || type == 4) {
    color = plottingEngine.GetRootColorLine1();
  } else if (type == 12 || type == 14 || type == 22) {
    color = plottingEngine.GetRootColorLine2();
  } else if (type < 10) {
    color = plottingEngine.GetRootColorElectron();
  } else if (type == 23) {
    color = kGreen;
  } else if (type == 24) {
    color = kRed;
  } else if (type < 20 || type == 25 || type == 26) {
    color = plottingEngine.GetRootColorHole();
  } else {
    color = plottingEngine.GetRootColorIon();
  }
  m_functions.back().SetLineColor(color);
  // Get the field grid.
  std::vector<double> efields;
  std::vector<double> bfields;
  std::vector<double> bangles;
  m_medium->GetFieldGrid(efields, bfields, bangles);
  const int nEfields = efields.size();
  const int nBfields = bfields.size();
  const int nBangles = bangles.size();
  bool withGraph = true;

  if (nEfields <= 0 || nBfields <= 0 || nBangles <= 0) {
    withGraph = false;
  }

  if (withGraph) {
    int n = 0;
    TGraph graph;
    if (xaxis == 'e') {  // plot with respect to E field
      graph.Set(nEfields);
      n = nEfields;
    } else if (xaxis == 'b') {  // plot wrt B field
      graph.Set(nBfields);
      n = nBfields;
    } else if (xaxis == 'a') {  // plot wrt angle
      graph.Set(nBangles);
      n = nBangles;
    } else {
      std::cerr << m_className << "::AddFunction:\n";
      std::cerr << "    Error specifying X-axis.\n";
      std::cerr << "    Invalid parameter type.\n";
      return;
    }
    graph.SetMarkerStyle(marker);
    graph.SetMarkerColor(color);
    bool ok = true;

    int epoint = -1;
    int bpoint = -1;
    int apoint = -1;
    // if one of these stays -1, this means there is no corresponding point in
    // the table

    // search for matching point in table with a certain accuracy
    for (int i = 0; i < n; ++i) {
      if (fabs(efields[i] - m_efield) <= m_etolerance) {
        epoint = i;
        break;
      }
    }
    for (int i = 0; i < n; ++i) {
      if (fabs(bfields[i] - m_bfield) <= m_btolerance) {
        bpoint = i;
        break;
      }
    }
    for (int i = 0; i < n; ++i) {
      if (fabs(bangles[i] - m_angle) <= m_atolerance) {
        apoint = i;
        break;
      }
    }
    if ((xaxis == 'e' && (bpoint == -1 || apoint == -1)) ||
        (xaxis == 'b' && (epoint == -1 || apoint == -1)) ||
        (xaxis == 'a' && (epoint == -1 || bpoint == -1)))
      ok = false;

    for (int i = 0; i < n; ++i) {
      double value = 0.;
      // auxiliary variables
      double alongx = 0, alongy = 0., alongz = 0.;
      if (!ok) {
        withGraph = false;
        break;
      }
      switch (type) {
        case 0:
          // Electron drift velocity along E
          if (xaxis == 'e')
            ok = m_medium->GetElectronVelocityE(i, bpoint, apoint, value);
          else if (xaxis == 'b')
            ok = m_medium->GetElectronVelocityE(epoint, i, apoint, value);
          else if (xaxis == 'a')
            ok = m_medium->GetElectronVelocityE(epoint, bpoint, i, value);
          else
            value = 0.;
          value = m_medium->ScaleVelocity(value);
          break;
        case 1:
          // Electron transverse diffusion
          ok = m_medium->GetElectronTransverseDiffusion(i, bpoint, apoint,
                                                        value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case 2:
          // Electron longitudinal diffusion
          ok = m_medium->GetElectronLongitudinalDiffusion(i, bpoint, apoint,
                                                          value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case 3:
          // Electron Townsend coefficient
          ok = m_medium->GetElectronTownsend(i, bpoint, apoint, value);
          value = m_medium->ScaleTownsend(exp(value));
          break;
        case 4:
          // Electron attachment coefficient
          ok = m_medium->GetElectronAttachment(i, bpoint, apoint, value);
          value = m_medium->ScaleAttachment(exp(value));
          break;
        case 10:
          // Hole drift velocity along E
          ok = m_medium->GetHoleVelocityE(i, bpoint, apoint, value);
          value = m_medium->ScaleVelocity(value);
          break;
        case 11:
          // Hole transverse diffusion
          ok = m_medium->GetHoleTransverseDiffusion(i, bpoint, apoint, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case 12:
          // Hole longitudinal diffusion
          ok = m_medium->GetHoleLongitudinalDiffusion(i, bpoint, apoint, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case 13:
          // Hole Townsend coefficient
          ok = m_medium->GetHoleTownsend(i, bpoint, apoint, value);
          value = m_medium->ScaleTownsend(exp(value));
          break;
        case 14:
          // Hole attachment coefficient
          ok = m_medium->GetHoleAttachment(i, bpoint, apoint, value);
          value = m_medium->ScaleAttachment(exp(value));
          break;
        case 20:
          // Ion drift velocity
          ok = m_medium->GetIonMobility(i, bpoint, apoint, value);
          value *= m_medium->UnScaleElectricField(efields[i]);
          break;
        case 21:
          // Ion transverse diffusion
          ok = m_medium->GetIonTransverseDiffusion(i, bpoint, apoint, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case 22:
          // Ion longitudinal diffusion
          ok = m_medium->GetIonLongitudinalDiffusion(i, bpoint, apoint, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case 23:
          // Electron drift velocity along B
          if (xaxis == 'e')
            ok = m_medium->ElectronVelocity(efields[i] * cos(bangles[apoint]),
                                            efields[i] * sin(bangles[apoint]),
                                            0, bfields[bpoint], 0, 0, value,
                                            alongy, alongz);
          else if (xaxis == 'b')
            ok = m_medium->ElectronVelocity(
                efields[epoint] * cos(bangles[apoint]),
                efields[epoint] * sin(bangles[apoint]), 0, bfields[i], 0, 0,
                value, alongy, alongz);
          else if (xaxis == 'a') {
            ok = m_medium->ElectronVelocity(efields[epoint] * cos(bangles[i]),
                                            efields[epoint] * sin(bangles[i]),
                                            0, bfields[bpoint], 0, 0, value,
                                            alongy, alongz);
          } else
            value = 0.;
          value = fabs(m_medium->ScaleVelocity(value));
          break;
        case 24:
          // Electron drift velocity along ExB
          if (xaxis == 'e')
            ok = m_medium->ElectronVelocity(efields[i] * cos(bangles[apoint]),
                                            efields[i] * sin(bangles[apoint]),
                                            0, bfields[bpoint], 0, 0, alongx,
                                            alongy, value);
          else if (xaxis == 'b')
            ok = m_medium->ElectronVelocity(
                efields[epoint] * cos(bangles[apoint]),
                efields[epoint] * sin(bangles[apoint]), 0, bfields[i], 0, 0,
                alongx, alongy, value);
          else if (xaxis == 'a')
            ok = m_medium->ElectronVelocity(
                m_efield * cos(bangles[i]), m_efield * sin(bangles[i]), 0,
                m_bfield, 0, 0, alongx, alongy, value);
          else
            value = 0.;
          value = fabs(m_medium->ScaleVelocity(value));
          break;
        case 25:
          // Hole drift velocity along B
          if (xaxis == 'e')
            ok = m_medium->HoleVelocity(efields[i] * cos(bangles[apoint]),
                                        efields[i] * sin(bangles[apoint]), 0,
                                        bfields[bpoint], 0, 0, value, alongy,
                                        alongz);
          else if (xaxis == 'b')
            ok = m_medium->HoleVelocity(efields[epoint] * cos(bangles[apoint]),
                                        efields[epoint] * sin(bangles[apoint]),
                                        0, bfields[i], 0, 0, value, alongy,
                                        alongz);
          else if (xaxis == 'a')
            ok = m_medium->HoleVelocity(efields[epoint] * cos(bangles[i]),
                                        efields[epoint] * sin(bangles[i]), 0,
                                        bfields[bpoint], 0, 0, value, alongy,
                                        alongz);
          break;
        case 26:
          // Hole velocity along ExB
          if (xaxis == 'e')
            ok = m_medium->HoleVelocity(efields[i] * cos(bangles[apoint]),
                                        efields[i] * sin(bangles[apoint]), 0,
                                        bfields[bpoint], 0, 0, alongx, alongy,
                                        value);
          else if (xaxis == 'b')
            ok = m_medium->HoleVelocity(efields[epoint] * cos(bangles[apoint]),
                                        efields[epoint] * sin(bangles[apoint]),
                                        0, bfields[i], 0, 0, alongx, alongy,
                                        value);
          else if (xaxis == 'a')
            ok = m_medium->HoleVelocity(m_efield * cos(bangles[i]),
                                        m_efield * sin(bangles[i]), 0, m_bfield,
                                        0, 0, alongx, alongy, value);
          else
            value = 0.;
          value = fabs(m_medium->ScaleVelocity(value));
          break;
      }
      if (xaxis == 'e')
        graph.SetPoint(i, m_medium->UnScaleElectricField(efields[i]), value);
      else if (xaxis == 'b')
        graph.SetPoint(i, bfields[i], value);
      else if (xaxis == 'a')
        graph.SetPoint(i, bangles[i], value);
    }
    if (ok) {
      m_graphs.push_back(graph);
      ++m_nGraphs;
    } else {
      std::cerr << m_className << "::AddFunction:\n";
      std::cerr << "    Error retrieving data table.\n";
      std::cerr << "    Suppress plotting of graph.\n";
    }
  }

  if (keep && m_nFunctions > 1) {
    m_functions[0].GetYaxis()->SetTitleOffset(1.5);
    m_functions[0].Draw("");
    for (int i = 1; i < m_nFunctions; ++i) {
      m_functions[i].Draw("lsame");
    }
  } else {
    m_functions.back().GetYaxis()->SetTitleOffset(1.5);
    m_functions.back().Draw("");
  }
  if (m_nGraphs > 0) {
    for (int i = 0; i < m_nGraphs; ++i) {
      m_graphs[i].Draw("p");
    }
  }
}

double ViewMedium::EvaluateFunction(double* pos, double* par) {
  // to be modified to include B and angle

  if (m_medium == 0) return 0.;
  int type = int(par[0]);
  char xaxis = char(par[1]);
  const double x = pos[0];
  double y = 0.;

  // Auxiliary variables
  double value = 0., a = 0., b = 0., c = 0., alongx = 0., alongy = 0.,
         alongz = 0.;

  switch (type) {
    case 0:
      // Electron drift velocity
      if (xaxis == 'e') {  // plot with respect to E field
        if (!m_medium->ElectronVelocity(x, 0, 0, m_bfield * cos(m_angle),
                                        m_bfield * sin(m_angle), 0, a, b, c))
          return 0.;
      } else if (xaxis == 'b') {  // plot wrt B field
        if (!m_medium->ElectronVelocity(m_efield, 0, 0, x * cos(m_angle),
                                        x * sin(m_angle), 0, a, b, c))
          return 0.;
      } else if (xaxis == 'a') {  // plot wrt angle
        if (!m_medium->ElectronVelocity(m_efield, 0, 0, m_bfield * cos(x),
                                        m_bfield * sin(x), 0, a, b, c))
          return 0.;
      }
      y = fabs(a);
      break;
    case 1:
      // Electron transverse diffusion
      if (!m_medium->ElectronDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = b;
      break;
    case 2:
      // Electron longitudinal diffusion
      if (!m_medium->ElectronDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = a;
      break;
    case 3:
      // Electron Townsend coefficient
      if (!m_medium->ElectronTownsend(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 4:
      // Electron attachment coefficient
      if (!m_medium->ElectronAttachment(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 10:
      // Hole drift velocity
      if (!m_medium->HoleVelocity(x, 0, 0, 0, 0, 0, a, b, c)) return 0.;
      y = a;
      break;
    case 11:
      // Hole transverse diffusion
      if (!m_medium->HoleDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = b;
      break;
    case 12:
      // Hole longitudinal diffusion
      if (!m_medium->HoleDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = a;
      break;
    case 13:
      // Hole Townsend coefficient
      if (!m_medium->HoleTownsend(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 14:
      // Hole attachment coefficient
      if (!m_medium->HoleAttachment(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 20:
      // Ion drift velocity
      if (!m_medium->IonVelocity(x, 0, 0, 0, 0, 0, a, b, c)) return 0.;
      y = fabs(a);
      break;
    case 21:
      // Ion transverse diffusion
      if (!m_medium->IonDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = b;
      break;
    case 22:
      // Ion longitudinal diffusion
      if (!m_medium->IonDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = a;
      break;
    case 23:
      // Electron drift velocity along B
      if (xaxis == 'e') {  // plot with respect to E field
        if (!m_medium->ElectronVelocity(x * cos(m_angle), x * sin(m_angle), 0,
                                        m_bfield, 0, 0, value, alongy, alongz))
          return 0.;
      } else if (xaxis == 'b') {  // plot wrt B field
        if (!m_medium->ElectronVelocity(m_efield * cos(m_angle),
                                        m_efield * sin(m_angle), 0, x, 0, 0,
                                        value, alongy, alongz))
          return 0.;
      } else if (xaxis == 'a') {  // plot wrt angle
        if (!m_medium->ElectronVelocity(m_efield * cos(x), m_efield * sin(x), 0,
                                        m_bfield, 0, 0, value, alongy, alongz))
          return 0.;
      }
      y = fabs(value);
      break;
    case 24:
      // Electron drift velocity along ExB
      if (xaxis == 'e') {  // plot with respect to E field
        if (!m_medium->ElectronVelocity(x * cos(m_angle), x * sin(m_angle), 0,
                                        m_bfield, 0, 0, alongx, alongy, value))
          return 0.;
      } else if (xaxis == 'b') {  // plot wrt B field
        if (!m_medium->ElectronVelocity(m_efield * cos(m_angle),
                                        m_efield * sin(m_angle), 0, x, 0, 0,
                                        alongx, alongy, value))
          return 0.;
      } else if (xaxis == 'a') {  // plot wrt angle
        if (!m_medium->ElectronVelocity(m_efield * cos(x), m_efield * sin(x), 0,
                                        m_bfield, 0, 0, alongx, alongy, value))
          return 0.;
      }
      y = fabs(value);
      break;
    case 25:
      // Hole drift velocity along B
      if (xaxis == 'e') {  // plot with respect to E field
        if (!m_medium->HoleVelocity(x * cos(m_angle), x * sin(m_angle), 0,
                                    m_bfield, 0, 0, value, alongy, alongz))
          return 0.;
      } else if (xaxis == 'b') {  // plot wrt B field
        if (!m_medium->HoleVelocity(m_efield * cos(m_angle),
                                    m_efield * sin(m_angle), 0, x, 0, 0, value,
                                    alongy, alongz))
          return 0.;
      } else if (xaxis == 'a') {  // plot wrt angle
        if (!m_medium->HoleVelocity(m_efield * cos(x), m_efield * sin(x), 0,
                                    m_bfield, 0, 0, value, alongy, alongz))
          return 0.;
      }
      y = fabs(value);
      break;
    case 26:
      // Hole drift velocity along ExB
      if (xaxis == 'e') {  // plot with respect to E field
        if (!m_medium->HoleVelocity(x * cos(m_angle), x * sin(m_angle), 0,
                                    m_bfield, 0, 0, alongx, alongy, value))
          return 0.;
      } else if (xaxis == 'b') {  // plot wrt B field
        if (!m_medium->HoleVelocity(m_efield * cos(m_angle),
                                    m_efield * sin(m_angle), 0, x, 0, 0, alongx,
                                    alongy, value))
          return 0.;
      } else if (xaxis == 'a') {  // plot wrt angle
        if (!m_medium->HoleVelocity(m_efield * cos(x), m_efield * sin(x), 0,
                                    m_bfield, 0, 0, alongx, alongy, value))
          return 0.;
      }
      y = fabs(value);
      break;
    default:
      std::cerr << m_className << "::EvaluateFunction:\n";
      std::cerr << "    Unknown type of transport coefficient requested.\n";
      std::cerr << "    Program bug!\n";
      return 0.;
  }

  return y;
}
}
