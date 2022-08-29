#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <limits>

#include <TH1F.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TLatex.h>

#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Medium.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewMedium.hh"

namespace {

int FindIndex(const std::vector<double>& fields, const double field,
              const double eps) {

  if (fields.empty()) return -1;
  const int n = fields.size();
  for (int i = 0; i < n; ++i) {
    const double sum = fabs(fields[i]) + fabs(field);
    const double tol = std::max(eps * sum, Garfield::Small);
    if (fabs(fields[i] - field) < tol) return i;
  }
  return -1;
} 

bool NonZero(const std::vector<double>& v) {
 
  constexpr double tol = 1.e-10;
  return std::any_of(v.cbegin(), v.cend(), [](double x){ return fabs(x) > tol; });
}

}

namespace Garfield {

ViewMedium::ViewMedium() : ViewBase("ViewMedium") {}

void ViewMedium::SetMedium(Medium* m) {
  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
    return;
  }

  m_medium = m;
}

void ViewMedium::SetRangeE(const double emin, const double emax, 
                           const bool logscale) {
  if (emin >= emax || emin < 0.) {
    std::cerr << m_className << "::SetRangeE: Incorrect range.\n";
    return;
  }

  m_eMin = emin;
  m_eMax = emax;
  m_logE = logscale;
}

void ViewMedium::SetRangeB(const double bmin, const double bmax, const bool logscale) {
  if (bmin >= bmax || bmin < 0.) {
    std::cerr << m_className << "::SetRangeB: Incorrect range.\n";
    return;
  }

  m_bMin = bmin;
  m_bMax = bmax;
  m_logB = logscale;
}

void ViewMedium::SetRangeA(const double amin, const double amax, const bool logscale) {
  if (amin >= amax || amin < 0.) {
    std::cerr << m_className << "::SetRangeA: Incorrect range.\n";
    return;
  }

  m_aMin = amin;
  m_aMax = amax;
  m_logA = logscale;
}

void ViewMedium::SetRangeY(const double ymin, const double ymax,
                           const bool logscale) {
  if (ymin >= ymax || ymin < 0.) {
    std::cerr << m_className << "::SetRangeY: Incorrect range.\n";
    return;
  }

  m_yMin = ymin;
  m_yMax = ymax;
  m_logY = logscale;
}

void ViewMedium::Draw() {

  if (m_yPlot.empty()) return;
  auto canvas = GetCanvas();
  canvas->cd();

  // Find the x and y ranges.
  const double xmin = m_xPlot.front();
  const double xmax = m_xPlot.back();
  double ymin = m_yMin;
  double ymax = m_yMax;
  if (m_autoRangeY) {
    ymin =  std::numeric_limits<double>::max();
    ymax = -std::numeric_limits<double>::max();
    for (const auto& plot : m_yPlot) {
      ymin = std::min(ymin, *std::min_element(plot.cbegin(), plot.cend()));
      ymax = std::max(ymax, *std::max_element(plot.cbegin(), plot.cend()));
    }
    const double dy = 0.05 * fabs(ymax - ymin);
    ymin = std::max(0., ymin - dy);
    ymax += dy;
  }
  auto frame = canvas->DrawFrame(xmin, ymin, xmax, ymax);
  // Set the axis labels.
  if (m_xaxis == Axis::E) {
    frame->GetXaxis()->SetTitle("electric field [V/cm]");
  } else if (m_xaxis == Axis::B) {
    frame->GetXaxis()->SetTitle("magnetic field [T]");
  } else if (m_xaxis == Axis::Angle) {
    frame->GetXaxis()->SetTitle("angle between #bf{E} and #bf{B} [rad]");
  } else {
    std::cerr << m_className << "::Draw: Invalid x axis.\n";
    return;
  }
  auto yaxis = frame->GetYaxis();
  switch (m_par[0]) {
    case Parameter::VelocityE:
    case Parameter::VelocityB:
    case Parameter::VelocityExB: 
      yaxis->SetTitle("drift velocity [cm/ns]");
      canvas->SetTitle("Drift velocity");
      break;
    case Parameter::LongitudinalDiffusion:
    case Parameter::TransverseDiffusion:
      yaxis->SetTitle("diffusion coefficient [#kern[-0.1]{#sqrt{cm}}#kern[0.1]{]}");
      canvas->SetTitle("Diffusion");
      break;
    case Parameter::Townsend:
    case Parameter::Attachment:
      if (std::all_of(m_par.cbegin(), m_par.cend(), [](Parameter p) { 
                      return p == Parameter::Townsend; })) {
        yaxis->SetTitle("#it{#alpha} [1/cm]");
        canvas->SetTitle("Multiplication");
      } else if (std::all_of(m_par.cbegin(), m_par.cend(), [](Parameter p) {
                             return p == Parameter::Attachment; })) {
        yaxis->SetTitle("#it{#eta} [1/cm]");
        canvas->SetTitle("Attachment");
      } else {
        yaxis->SetTitle("#it{#alpha}, #it{#eta} [1/cm]");
        canvas->SetTitle("Multiplication and attachment");
      }
      break;
    case Parameter::LorentzAngle:
      yaxis->SetTitle("Angle between #bf{v} and #bf{E} [rad]");
      canvas->SetTitle("Lorentz angle");
      break;
    default:
      break;
  }
  yaxis->SetTitleOffset(0);
  gPad->SetLeftMargin(0.15);
  canvas->Update();

  const size_t nPlots = m_yPlot.size();
  // Set colours.
  std::vector<short> cols(nPlots, 0);
  if (m_colours.empty()) {
    auto it = std::find(m_q.cbegin(), m_q.cend(), Charge::Electron);
    if (it != m_q.cend()) {
      cols[std::distance(m_q.cbegin(), it)] = kOrange - 3; 
    } 
    it = std::find(std::next(it), m_q.cend(), Charge::Electron);
    if (it != m_q.cend()) {
      cols[std::distance(m_q.cbegin(), it)] = kGreen + 3; 
    } 
    it = std::find(m_q.cbegin(), m_q.cend(), Charge::Hole); 
    if (it != m_q.cend()) {
      cols[std::distance(m_q.cbegin(), it)] = kRed + 1; 
    } 
    it = std::find(m_q.cbegin(), m_q.cend(), Charge::Ion); 
    if (it != m_q.cend()) {
      cols[std::distance(m_q.cbegin(), it)] = kRed + 1; 
    }
  } else {
    for (size_t i = 0; i < nPlots; ++i) {
      cols[i] = m_colours[i % m_colours.size()];  
    }
  } 
  // Set legend.
  std::vector<std::string> labels = m_labels;
  labels.resize(nPlots, "");
  if (nPlots > 1 && m_labels.empty()) {
    bool allEqual = std::equal(m_q.begin() + 1, m_q.end(), m_q.begin());
    if (!allEqual) {
      for (size_t i = 0; i < nPlots; ++i) {
        if (m_q[i] == Charge::Electron) {
          labels[i] = "electrons";
        } else if (m_q[i] == Charge::Hole) {
          labels[i] = "holes";
        } else {
          labels[i] = "ions";
        }
      } 
    } 
    allEqual = std::equal(m_par.begin() + 1, m_par.end(), m_par.begin());
    if (!allEqual) {
      for (size_t i = 0; i < nPlots; ++i) {
        if (!labels[i].empty()) labels[i] += ", ";
        switch (m_par[i]) { 
          case Parameter::VelocityE:
            labels[i] += "#it{v}_{#it{E}}";
            break;
          case Parameter::VelocityB:
            labels[i] += "#it{v}_{#it{B}#kern[0.1]{t}}";
            break;
          case Parameter::VelocityExB:
            labels[i] += "#it{v}_{#it{E}#kern[0.05]{#times}#it{B}}";
            break;
          case Parameter::LongitudinalDiffusion:
            labels[i] += "#it{D}_{L}";
            break; 
          case Parameter::TransverseDiffusion:
            labels[i] += "#it{D}_{T}";
            break; 
          case Parameter::Townsend:
            labels[i] += "#alpha";
            break;
          case Parameter::Attachment:
            labels[i] += "#eta";
            break;
          case Parameter::LorentzAngle:
            labels[i] += "Lorentz angle";
        }
      }
    }
  }
  TGraph graph;
  TLatex latex;
  double xsize = 0., ysize = 0.;
  for (const auto& label : labels) {
    latex.SetText(0, 0, label.c_str());
    xsize = std::max(xsize, latex.GetXsize());
    ysize = std::max(ysize, latex.GetYsize());
  }
  // Convert to NDC.
  xsize /= (gPad->GetX2() - gPad->GetX1());
  ysize /= (gPad->GetY2() - gPad->GetY1());

  const double lm = gPad->GetLeftMargin();
  const double rm = 1. - gPad->GetRightMargin();
  const double tm = 1. - gPad->GetTopMargin();
  const double bm = gPad->GetBottomMargin();
  double xLabel = lm + 0.1 * (rm - lm);
  if (m_par[0] == Parameter::LongitudinalDiffusion ||
      m_par[0] == Parameter::TransverseDiffusion) {
    xLabel = rm - 0.1 * (rm - lm) - xsize; 
  }
  double yLabel = tm - 0.1 * (tm - bm);
  const double colrange = gStyle->GetNumberOfColors() / double(nPlots);
  for (size_t i = 0; i < nPlots; ++i) {
    int col = cols[i] > 0 ? cols[i] : gStyle->GetColorPalette(i * colrange);
    graph.SetLineColor(col);
    graph.SetMarkerColor(col);
    graph.DrawGraph(m_xPlot.size(), 
                    m_xPlot.data(), m_yPlot[i].data(), "L");
    graph.SetMarkerStyle(20 + i);
    if (!m_xGraph[i].empty()) {
      graph.DrawGraph(m_xGraph[i].size(), 
                      m_xGraph[i].data(), m_yGraph[i].data(), "P");

    }
    if (labels[i].empty()) continue;
    latex.SetTextColor(col);
    latex.DrawLatexNDC(xLabel, yLabel, labels[i].c_str());
    yLabel -= 1.5 * ysize;
  }
  if (m_logX && xmin > 0.) {
    canvas->SetLogx(1);
  } else {
    canvas->SetLogx(0);
  }
  if (m_logY && ymin > 0. && !m_autoRangeY) {
    canvas->SetLogy(1);
  } else {
    canvas->SetLogy(0);
  }
  gPad->Update();
  if (!m_outfile.empty()) Export();
}

void ViewMedium::Export() {

  if (m_yPlot.empty()) return;
  const size_t nPlots = m_yPlot.size();
  std::vector<std::string> ylabel = m_labels;
  ylabel.resize(nPlots, "");
  for (size_t i = 0; i < nPlots; ++i) {
    if (!ylabel[i].empty()) continue;
    switch (m_q[i]) {
      case Charge::Electron:
        ylabel[i] = "electron ";
        break;
      case Charge::Hole:
        ylabel[i] = "hole ";
        break;
      case Charge::Ion:
        ylabel[i] = "ion ";
        break;
      default:
        break;
    }
    switch (m_par[i]) { 
      case Parameter::VelocityE:
        ylabel[i] += "drift velocity along E [cm/ns]";
        break;
      case Parameter::VelocityB:
        ylabel[i] += "drift velocity along Bt [cm/ns]";
        break;
      case Parameter::VelocityExB:
        ylabel[i] += "drift velocity along ExB [cm/ns]";
        break;
      case Parameter::LongitudinalDiffusion:
        ylabel[i] += "longitudinal diffusion [cm1/2]";
        break; 
      case Parameter::TransverseDiffusion:
        ylabel[i] += "transverse diffusion [cm1/2]";
        break; 
      case Parameter::Townsend:
        ylabel[i] += "Townsend coefficient [1/cm]";
        break;
      case Parameter::Attachment:
        ylabel[i] += "attachment coefficient [1/cm]";
        break;
      case Parameter::LorentzAngle:
        ylabel[i] += "Lorentz angle [rad]";
    }
  }

  const std::string sep = " ";
  std::ofstream outfile(m_outfile, std::ios_base::app);
  if (!outfile) return;
  outfile << "# x-axis: ";
  if (m_xaxis == Axis::E) {
    outfile << "E [V/cm]";
  } else if (m_xaxis == Axis::B) {
    outfile << "B [T]";
  } else if (m_xaxis == Axis::Angle) {
    outfile << "Theta [rad]";
  }
  outfile << "\n";
  outfile << "# " << nPlots << " plots:\n";
  for (size_t i = 0; i < nPlots; ++i) {
    outfile << "# " << ylabel[i] << "\n";
  }
  const size_t nX = m_xPlot.size();
  for (size_t i = 0; i < nX; ++i) {
    outfile << m_xPlot[i];
    for (size_t j = 0; j < nPlots; ++j) {
      outfile << sep << m_yPlot[j][i];
    }
    outfile << "\n";
  }
  for (size_t i = 0; i < nPlots; ++i) {
    if (m_xGraph[i].empty()) continue;
    outfile << "# " << ylabel[i] << "\n";
    for (size_t j = 0; j < m_xGraph[i].size(); ++j) {
      outfile << m_xGraph[i][j] << sep << m_yGraph[i][j] << "\n";
    }
  }
  outfile.close();
}

void ViewMedium::ResetY() {
  m_yPlot.clear();
  m_par.clear();
  m_q.clear();
  m_xGraph.clear();
  m_yGraph.clear();
}

void ViewMedium::ResetX(const Axis xaxis) {

  double xmin = 0., xmax = 0.;
  bool logx = false;
  if (xaxis == Axis::E) {
    xmin = m_eMin;
    xmax = m_eMax;
    logx = m_logE;
  } else if (xaxis == Axis::B) {
    xmin = m_bMin;
    xmax = m_bMax;
    logx = m_logB;
  } else if (xaxis == Axis::Angle) {
    xmin = m_aMin;
    xmax = m_aMax;
    logx = m_logA;
  } else {
    m_xaxis = Axis::None;
    return;
  }

  if (m_autoRangeX) {
    // Get the field grid.
    std::vector<double> efields;
    std::vector<double> bfields;
    std::vector<double> bangles;
    m_medium->GetFieldGrid(efields, bfields, bangles);
    if (xaxis == Axis::E && !efields.empty()) {
      if (efields.front() > 0. && efields.back() > 10. * efields.front()) {
        logx = true;
        xmin = efields.front() / 1.5;
        xmax = efields.back() * 1.5;
      } else {
        logx = false;
        const double dx = 0.05 * fabs(efields.back() - efields.front());
        xmin = std::max(0., efields.front() - dx);
        xmax = efields.back() + dx;
      }
    } else if (xaxis == Axis::B && !bfields.empty()) {
      logx = false;
      const double dx = 0.05 * fabs(bfields.back() - bfields.front());
      xmin = std::max(0., bfields.front() - dx);
      xmax = bfields.back() + dx;
    } else if (xaxis == Axis::Angle && !bangles.empty()) {
      logx = false;
      const double dx = 0.05 * fabs(bangles.back() - bangles.front());
      xmin = std::max(0., bangles.front() - dx);
      xmax = std::min(TwoPi, bangles.back() + dx);
    }
  }

  constexpr size_t nX = 1000;
  m_xPlot.assign(nX, 0.);
  if (logx) {
    m_xPlot[0] = xmin;
    const double r = pow(xmax / xmin, 1. / (nX - 1));
    for (size_t i = 1; i < nX; ++i) {
      m_xPlot[i] = m_xPlot[i - 1] * r;
    }
  } else {
    const double dx = (xmax - xmin) / (nX - 1);
    for (size_t i = 0; i < nX; ++i) {
      m_xPlot[i] = xmin + i * dx;
    }
  }
  m_logX = logx;
  m_xaxis = xaxis;
}

void ViewMedium::PlotDiffusion(const Axis xaxis, const Charge charge,
                               const bool same) {

  if (!m_medium) {
    std::cerr << m_className << "::PlotDiffusion: Medium is not defined.\n";
    return;
  }
  if (xaxis != m_xaxis) {
    ResetX(xaxis);
    ResetY();
  } else if (!same) {
    ResetY();
  } else if (!m_par.empty()) {
    if (m_par[0] != Parameter::TransverseDiffusion && 
        m_par[0] != Parameter::LongitudinalDiffusion) {
      ResetY();
    }
  } 
  const size_t nX = m_xPlot.size();
  std::array<std::vector<double>, 2> ypl;
  for (size_t i = 0; i < 2; ++i) ypl[i].assign(nX, 0.);

  double ex = m_efield;
  double ctheta = cos(m_angle);
  double stheta = sin(m_angle);
  double bx = m_bfield * ctheta;
  double by = m_bfield * stheta; 
  for (size_t i = 0; i < nX; ++i) {
    if (xaxis == Axis::E) {
      ex = m_xPlot[i];
    } else if (xaxis == Axis::B) {
      bx = m_xPlot[i] * ctheta;
      by = m_xPlot[i] * stheta;
    } else {
      bx = m_bfield * cos(m_xPlot[i]);
      by = m_bfield * sin(m_xPlot[i]);
    } 
    double dl = 0., dt = 0.;
    if (charge == Charge::Electron) {
      if (!m_medium->ElectronDiffusion(ex, 0, 0, bx, by, 0, dl, dt)) continue;
    } else if (charge == Charge::Hole) {
      if (!m_medium->HoleDiffusion(ex, 0, 0, bx, by, 0, dl, dt)) continue;
    } else {
      if (!m_medium->IonDiffusion(ex, 0, 0, bx, by, 0, dl, dt)) continue;
    }
    ypl[0][i] = dl;
    ypl[1][i] = dt;
  }

  std::array<std::vector<double>, 2> xgr;
  std::array<std::vector<double>, 2> ygr;

  std::array<std::vector<double>, 3> grid;
  int ie = 0, ib = 0, ia = 0;
  if (GetGrid(grid, ie, ib, ia, xaxis)) {
    const auto nPoints = xaxis == Axis::E ? grid[0].size() : 
                         xaxis == Axis::B ? grid[1].size() : grid[2].size();
    for (size_t j = 0; j < nPoints; ++j) {
      double x = 0., y = 0.;
      if (xaxis == Axis::E) {
        ie = j;
        x = m_medium->UnScaleElectricField(grid[0][j]);
      } else if (xaxis == Axis::B) {
        ib = j;
        x = grid[1][j];
      } else if (xaxis == Axis::Angle) {
        ia = j;
        x = grid[2][j];
      }
      if (charge == Charge::Electron) {
        if (m_medium->GetElectronTransverseDiffusion(ie, ib, ia, y)) {
          xgr[1].push_back(x);
          ygr[1].push_back(m_medium->ScaleDiffusion(y));
        }
        if (m_medium->GetElectronLongitudinalDiffusion(ie, ib, ia, y)) {
          xgr[0].push_back(x);
          ygr[0].push_back(m_medium->ScaleDiffusion(y));
        }
      } else if (charge == Charge::Hole) {
        if (m_medium->GetHoleTransverseDiffusion(ie, ib, ia, y)) {
          xgr[1].push_back(x);
          ygr[1].push_back(m_medium->ScaleDiffusion(y));
        }
        if (m_medium->GetHoleLongitudinalDiffusion(ie, ib, ia, y)) {
          xgr[0].push_back(x);
          ygr[0].push_back(m_medium->ScaleDiffusion(y));
        }
      } else {
        if (m_medium->GetIonTransverseDiffusion(ie, ib, ia, y)) {
          xgr[1].push_back(x);
          ygr[1].push_back(m_medium->ScaleDiffusion(y));
        }
        if (m_medium->GetIonLongitudinalDiffusion(ie, ib, ia, y)) {
          xgr[0].push_back(x);
          ygr[0].push_back(m_medium->ScaleDiffusion(y));
        }
      }
    }
  }

  const std::array<Parameter, 2> pars = {Parameter::LongitudinalDiffusion,
                                         Parameter::TransverseDiffusion};
  for (size_t i = 0; i < 2; ++i) {
    if (!NonZero(ypl[i])) continue;
    m_yPlot.push_back(std::move(ypl[i]));
    m_par.push_back(pars[i]);
    m_q.push_back(charge);
    m_xGraph.push_back(std::move(xgr[i]));
    m_yGraph.push_back(std::move(ygr[i]));
  }
  Draw();
}

void ViewMedium::PlotVelocity(const Axis xaxis, const Charge charge,
                              const bool same) {

  if (!m_medium) {
    std::cerr << m_className << "::PlotVelocity: Medium is not defined.\n";
    return;
  }
  if (xaxis != m_xaxis) {
    ResetX(xaxis);
    ResetY();
  } else if (!same) {
    ResetY();
  } else if (!m_par.empty()) {
    if (m_par[0] != Parameter::VelocityE && 
        m_par[0] != Parameter::VelocityB && 
        m_par[0] != Parameter::VelocityExB) {
      ResetY();
    }
  } 
  const size_t nX = m_xPlot.size();
  std::array<std::vector<double>, 3> ypl;
  for (size_t i = 0; i < 3; ++i) ypl[i].assign(nX, 0.);

  double e0 = m_efield;
  double b0 = m_bfield;
  double ctheta = cos(m_angle);
  double stheta = sin(m_angle);
  
  for (size_t i = 0; i < nX; ++i) {
    if (xaxis == Axis::E) {
      e0 = m_xPlot[i];
    } else if (xaxis == Axis::B) {
      b0 = m_xPlot[i];
    } else { 
      ctheta = cos(m_xPlot[i]);
      stheta = sin(m_xPlot[i]);
    }
    double vx = 0., vy = 0., vz = 0.;
    if (charge == Charge::Electron) {
      if (m_medium->ElectronVelocity(e0, 0, 0, b0 * ctheta, b0 * stheta, 0, 
                                     vx, vy, vz)) {
        ypl[0][i] = fabs(vx);
        ypl[1][i] = fabs(vy);
        ypl[2][i] = fabs(vz);
      }
    } else if (charge == Charge::Hole) {
      if (m_medium->HoleVelocity(e0, 0, 0, b0 * ctheta, b0 * stheta, 0, 
                                 vx, vy, vz)) {
        ypl[0][i] = fabs(vx);
        ypl[1][i] = fabs(vy);
        ypl[2][i] = fabs(vz);
      }
    } else {
      if (m_medium->IonVelocity(e0, 0, 0, b0 * ctheta, b0 * stheta, 0, 
                                vx, vy, vz)) {
        ypl[0][i] = fabs(vx);
      }
    }
  }
  std::array<std::vector<double>, 3> xgr;
  std::array<std::vector<double>, 3> ygr;

  std::array<std::vector<double>, 3> grid;
  int ie = 0, ib = 0, ia = 0;
  if (GetGrid(grid, ie, ib, ia, xaxis)) {
    const auto nPoints = xaxis == Axis::E ? grid[0].size() : 
                         xaxis == Axis::B ? grid[1].size() : grid[2].size();
    for (size_t j = 0; j < nPoints; ++j) {
      double x = 0., y = 0.;
      if (xaxis == Axis::E) {
        ie = j;
        x = m_medium->UnScaleElectricField(grid[0][j]);
      } else if (xaxis == Axis::B) {
        ib = j;
        x = grid[1][j];
      } else if (xaxis == Axis::Angle) {
        ia = j;
        x = grid[2][j];
      }
      if (charge == Charge::Electron) {
        if (m_medium->GetElectronVelocityE(ie, ib, ia, y)) {
          xgr[0].push_back(x);
          ygr[0].push_back(m_medium->ScaleVelocity(y));
        }
        if (m_medium->GetElectronVelocityB(ie, ib, ia, y)) {
          xgr[1].push_back(x);
          ygr[1].push_back(fabs(m_medium->ScaleVelocity(y)));
        }
        if (m_medium->GetElectronVelocityExB(ie, ib, ia, y)) {
          xgr[2].push_back(x);
          ygr[2].push_back(fabs(m_medium->ScaleVelocity(y)));
        }
      } else if (charge == Charge::Hole) {
        if (m_medium->GetHoleVelocityE(ie, ib, ia, y)) {
          xgr[0].push_back(x);
          ygr[0].push_back(m_medium->ScaleVelocity(y));
        }
        if (m_medium->GetHoleVelocityB(ie, ib, ia, y)) {
          xgr[1].push_back(x);
          ygr[1].push_back(fabs(m_medium->ScaleVelocity(y)));
        }
        if (m_medium->GetHoleVelocityExB(ie, ib, ia, y)) {
          xgr[2].push_back(x);
          ygr[2].push_back(fabs(m_medium->ScaleVelocity(y)));
        }
      } else {
        if (m_medium->GetIonMobility(ie, ib, ia, y)) {
          xgr[0].push_back(x);
          ygr[0].push_back(y * m_medium->UnScaleElectricField(grid[0][ie]));
        }
      }
    }
  }

  const std::array<Parameter, 3> pars = {Parameter::VelocityE, 
                                         Parameter::VelocityB,
                                         Parameter::VelocityExB};
  for (size_t i = 0; i < 3; ++i) {
    if (!NonZero(ypl[i])) continue;
    m_yPlot.push_back(std::move(ypl[i]));
    m_par.push_back(pars[i]);
    m_q.push_back(charge);
    m_xGraph.push_back(std::move(xgr[i]));
    m_yGraph.push_back(std::move(ygr[i]));
  }
  Draw();
}

void ViewMedium::Plot(const Axis xaxis, const Charge charge,
                      const Parameter par, const bool same) {

  // Make sure the medium is set.
  if (!m_medium) {
    std::cerr << m_className << "::Plot: Medium is not defined.\n";
    return;
  }
  if (xaxis != m_xaxis) {
    ResetX(xaxis);
    ResetY();
  } else if (!same) {
    ResetY();
  } else if (!m_par.empty()) {
    if (m_par[0] != Parameter::Townsend &&
        m_par[0] != Parameter::Attachment) {
      ResetY();
    }
  } 

  const size_t nX = m_xPlot.size();
  std::vector<double> ypl(nX, 0.);  
  double ex = m_efield;
  double ctheta = cos(m_angle);
  double stheta = sin(m_angle);
  double bx = m_bfield * ctheta;
  double by = m_bfield * stheta;
  for (size_t i = 0; i < nX; ++i) {
    if (xaxis == Axis::E) {
      ex = m_xPlot[i];
    } else if (xaxis == Axis::B) {
      bx = m_xPlot[i] * ctheta;
      by = m_xPlot[i] * stheta;
    } else {
      bx = m_bfield * cos(m_xPlot[i]);
      by = m_bfield * sin(m_xPlot[i]);
    }
    double y = 0.;
    if (charge == Charge::Electron) {
      if (par == Parameter::Townsend) {
        if (!m_medium->ElectronTownsend(ex, 0, 0, bx, by, 0, y)) continue;
      } else {
        if (!m_medium->ElectronAttachment(ex, 0, 0, bx, by, 0, y)) continue;
        y = std::abs(y);
      }
    } else {
      if (par == Parameter::Townsend) {
        if (!m_medium->HoleTownsend(ex, 0, 0, bx, by, 0, y)) continue;
      } else {
        if (!m_medium->HoleAttachment(ex, 0, 0, bx, by, 0, y)) continue;
        y = std::abs(y);
      }
    }
    ypl[i] = y;
  }

  std::vector<double> xgr;
  std::vector<double> ygr;

  std::array<std::vector<double>, 3> grid;
  int ie = 0, ib = 0, ia = 0;
  if (GetGrid(grid, ie, ib, ia, xaxis)) {
    const auto nPoints = xaxis == Axis::E ? grid[0].size() : 
                         xaxis == Axis::B ? grid[1].size() : grid[2].size();
    for (size_t j = 0; j < nPoints; ++j) {
      double x = 0., y = 0.;
      if (xaxis == Axis::E) {
        ie = j;
        x = m_medium->UnScaleElectricField(grid[0][j]);
      } else if (xaxis == Axis::B) {
        ib = j;
        x = grid[1][j];
      } else if (xaxis == Axis::Angle) {
        ia = j;
        x = grid[2][j];
      }
      if (charge == Charge::Electron) {
        if (par == Parameter::Townsend) {
          if (m_medium->GetElectronTownsend(ie, ib, ia, y)) {
            xgr.push_back(x);
            ygr.push_back(m_medium->ScaleTownsend(exp(y)));
          }
        } else {
          if (m_medium->GetElectronAttachment(ie, ib, ia, y)) {
            xgr.push_back(x);
            ygr.push_back(m_medium->ScaleAttachment(exp(y)));
          }
        }
      } else if (charge == Charge::Hole) {
        if (par == Parameter::Townsend) {
          if (m_medium->GetHoleTownsend(ie, ib, ia, y)) {
            xgr.push_back(x);
            ygr.push_back(m_medium->ScaleTownsend(exp(y)));
          }
        } else {
          if (m_medium->GetHoleAttachment(ie, ib, ia, y)) {
            xgr.push_back(x);
            ygr.push_back(m_medium->ScaleAttachment(exp(y)));
          }
        }
      }
    }
  }

  if (!NonZero(ypl)) return;
  m_yPlot.push_back(std::move(ypl));
  m_par.push_back(par);
  m_q.push_back(charge);
  m_xGraph.push_back(std::move(xgr));
  m_yGraph.push_back(std::move(ygr));
  Draw();
}

void ViewMedium::PlotLorentzAngle(const Axis xaxis, const Charge charge,
                                  const bool same) {

  // Make sure the medium is set.
  if (!m_medium) {
    std::cerr << m_className << "::PlotLorentzAngle: Medium is not defined.\n";
    return;
  }
  if (xaxis != m_xaxis) {
    ResetX(xaxis);
    ResetY();
  } else if (!same) {
    ResetY();
  } else if (!m_par.empty()) {
    if (m_par[0] != Parameter::LorentzAngle) ResetY();
  } 

  const size_t nX = m_xPlot.size();
  std::vector<double> ypl(nX, 0.);  
  double ex = m_efield;
  double ctheta = cos(m_angle);
  double stheta = sin(m_angle);
  double bx = m_bfield * ctheta;
  double by = m_bfield * stheta;
  for (size_t i = 0; i < nX; ++i) {
    if (xaxis == Axis::E) {
      ex = m_xPlot[i];
    } else if (xaxis == Axis::B) {
      bx = m_xPlot[i] * ctheta;
      by = m_xPlot[i] * stheta;
    } else {
      bx = m_bfield * cos(m_xPlot[i]);
      by = m_bfield * sin(m_xPlot[i]);
    }
    double y = 0.;
    if (!m_medium->ElectronLorentzAngle(ex, 0, 0, bx, by, 0, y)) continue;
    ypl[i] = y;
  }

  std::vector<double> xgr;
  std::vector<double> ygr;

  std::array<std::vector<double>, 3> grid;
  int ie = 0, ib = 0, ia = 0;
  if (GetGrid(grid, ie, ib, ia, xaxis)) {
    const auto nPoints = xaxis == Axis::E ? grid[0].size() : 
                         xaxis == Axis::B ? grid[1].size() : grid[2].size();
    for (size_t j = 0; j < nPoints; ++j) {
      double x = 0., y = 0.;
      if (xaxis == Axis::E) {
        ie = j;
        x = m_medium->UnScaleElectricField(grid[0][j]);
      } else if (xaxis == Axis::B) {
        ib = j;
        x = grid[1][j];
      } else if (xaxis == Axis::Angle) {
        ia = j;
        x = grid[2][j];
      }
      if (m_medium->GetElectronLorentzAngle(ie, ib, ia, y)) {
        xgr.push_back(x);
        ygr.push_back(y);
      }
    }
  }

  if (!NonZero(ypl)) return;
  m_yPlot.push_back(std::move(ypl));
  m_par.push_back(Parameter::LorentzAngle);
  m_q.push_back(charge);
  m_xGraph.push_back(std::move(xgr));
  m_yGraph.push_back(std::move(ygr));
  Draw();
}

ViewMedium::Axis ViewMedium::GetAxis(const char xaxis) const {

  if (std::toupper(xaxis) == 'E') {
    return Axis::E;
  } else if (std::toupper(xaxis) == 'B') {
    return Axis::B;
  } else if (std::toupper(xaxis) == 'A') {
    return Axis::Angle;
  }
  return Axis::None; 
}

bool ViewMedium::GetGrid(std::array<std::vector<double>, 3>& grid,
                         int& ie, int& ib, int& ia,
                         const Axis xaxis) const { 

  if (!m_medium) return false;
  m_medium->GetFieldGrid(grid[0], grid[1], grid[2]);
  if (grid[0].empty() || grid[1].empty() || grid[2].empty()) return false;
  constexpr double eps = 1.e-3;
  ie = FindIndex(grid[0], m_efield, eps);
  ib = FindIndex(grid[1], m_bfield, eps);
  ia = FindIndex(grid[2], m_angle, eps);
  if (xaxis == Axis::E) {
    if (ib < 0 || ia < 0) return false;
  } else if (xaxis == Axis::B) {
    if (ie < 0 || ia < 0) return false;
  } else if (xaxis == Axis::Angle) {
    if (ie < 0 || ib < 0) return false;
  } else {
    return false;
  }
  return true;
}

}
