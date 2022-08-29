#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <algorithm>

#include <TAxis.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>

#include "Garfield/Component.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewField.hh"

namespace {

void SampleRange(const double xmin, const double ymin, const double xmax,
                 const double ymax, TF2* f, double& zmin, double& zmax) {
  constexpr unsigned int n = 1000;
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  zmin = std::numeric_limits<double>::max();
  zmax = -zmin;
  for (unsigned int i = 0; i < n; ++i) {
    const double z = f->Eval(xmin + Garfield::RndmUniform() * dx,
                             ymin + Garfield::RndmUniform() * dy);
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }
}

void SampleRange(TF1* f, double& ymin, double& ymax) {
  constexpr unsigned int n = 1000;
  ymin = std::numeric_limits<double>::max();
  ymax = -ymin;
  double xmin = 0.;
  double xmax = 1.;
  f->GetRange(xmin, xmax);
  const double dx = xmax - xmin;
  for (unsigned int i = 0; i < n; ++i) {
    const double y = f->Eval(xmin + dx * Garfield::RndmUniform());
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }
}

double Interpolate(const std::array<double, 1000>& y,
                   const std::array<double, 1000>& x, const double xx) {

  const double tol = 1.e-6 * fabs(x.back() - x.front());
  if (xx < x[0]) return y[0];
  const auto it1 = std::upper_bound(x.cbegin(), x.cend(), xx);
  if (it1 == x.cend()) return y.back();
  const auto it0 = std::prev(it1);
  const double dx = (*it1 - *it0);
  if (dx < tol) return y[it0 - x.cbegin()];
  const double f = (xx - *it0) / dx;
  return y[it0 - x.cbegin()] * (1. - f) + f * y[it1 - x.cbegin()];
}

}

namespace Garfield {

ViewField::ViewField() : ViewBase("ViewField") { }

void ViewField::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }
  m_sensor = s;
  m_component = nullptr;
}

void ViewField::SetComponent(Component* c) {
  if (!c) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }
  m_component = c;
  m_sensor = nullptr;
}

void ViewField::SetVoltageRange(const double vmin, const double vmax) {
  m_vmin = std::min(vmin, vmax);
  m_vmax = std::max(vmin, vmax);
  m_useAutoRange = false;
}

void ViewField::SetElectricFieldRange(const double emin, const double emax) {
  m_emin = std::min(emin, emax);
  m_emax = std::max(emin, emax);
  m_useAutoRange = false;
}

void ViewField::SetWeightingFieldRange(const double wmin, const double wmax) {
  m_wmin = std::min(wmin, wmax);
  m_wmax = std::max(wmin, wmax);
  m_useAutoRange = false;
}

void ViewField::SetMagneticFieldRange(const double bmin, const double bmax) {
  m_bmin = std::min(bmin, bmax);
  m_bmax = std::max(bmin, bmax);
  m_useAutoRange = false;
}

void ViewField::SetNumberOfContours(const unsigned int n) {
  if (n > 0) m_nContours = n;
}

void ViewField::SetNumberOfSamples1d(const unsigned int n) {
  m_nSamples1d = std::max(4u, n);
}

void ViewField::SetNumberOfSamples2d(const unsigned int nx,
                                     const unsigned int ny) {
  m_nSamples2dX = std::max(4u, nx);
  m_nSamples2dY = std::max(4u, ny);
}

void ViewField::PlotContour(const std::string& option) {
  Draw2d(option, true, false, "", "CONT1Z");
}

void ViewField::Plot(const std::string& option, const std::string& drawopt) {
  std::string opt1;
  std::transform(drawopt.begin(), drawopt.end(), 
                 std::back_inserter(opt1), toupper);
  if (opt1.find("CONT") != std::string::npos) {
    Draw2d(option, true, false, "", drawopt);
  } else {
    Draw2d(option, false, false, "", drawopt);
  }
}

void ViewField::PlotProfile(const double x0, const double y0, const double z0,
                            const double x1, const double y1, const double z1,
                            const std::string& option, const bool normalised) {
  DrawProfile(x0, y0, z0, x1, y1, z1, option, false, "", normalised);
}

void ViewField::PlotWeightingField(const std::string& label,
                                   const std::string& option,
                                   const std::string& drawopt, const double t) {
  Draw2d(option, false, true, label, drawopt, t);
}

void ViewField::PlotContourWeightingField(const std::string& label,
                                          const std::string& option) {
  Draw2d(option, true, true, label, "CONT1Z");
}

void ViewField::PlotProfileWeightingField(const std::string& label,
                                          const double x0, const double y0,
                                          const double z0, const double x1,
                                          const double y1, const double z1,
                                          const std::string& option,
                                          const bool normalised) {
  DrawProfile(x0, y0, z0, x1, y1, z1, option, true, label, normalised);
}


ViewField::Parameter ViewField::GetPar(const std::string& option,
                                       std::string& title, bool& bfield) const {

  bfield = false;
  std::string opt;
  std::transform(option.begin(), option.end(), 
                 std::back_inserter(opt), toupper);
  if (opt == "BMAG") {
    title = "field";
    bfield = true;
    return Parameter::Bmag;
  } else if (opt == "BX") {
    title = "field (x-component)";
    bfield = true;
    return Parameter::Bx;
  } else if (opt == "BY") {
    title = "field (y-component)";
    bfield = true;
    return Parameter::By;
  } else if (opt == "BZ") {
    title = "field (z-component)";
    bfield = true;
    return Parameter::Bz;
  } else if (opt == "V" || opt == "P" || opt == "PHI" || 
      opt.find("VOLT") != std::string::npos ||
      opt.find("POT") != std::string::npos) {
    title = "potential";
    return Parameter::Potential;
  } else if (opt == "E" || opt == "FIELD" || opt == "NORM" ||
             opt.find("MAG") != std::string::npos) {
    title = "field";
    return Parameter::Emag;
  } else if (opt.find("X") != std::string::npos) {
    title = "field (x-component)";
    return Parameter::Ex;
  } else if (opt.find("Y") != std::string::npos) {
    title = "field (y-component)";
    return Parameter::Ey;
  } else if (opt.find("Z") != std::string::npos) {
    title = "field (z-component)";
    return Parameter::Ez;
  }
  std::cerr << m_className << "::GetPar: Unknown option (" << option << ").\n";
  title = "potential";
  return Parameter::Potential;
}

void ViewField::Draw2d(const std::string& option, const bool contour,
                       const bool wfield, const std::string& electrode,
                       const std::string& drawopt, const double t) {
  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::Draw2d:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }

  // Determine the x-y range.
  if (!SetPlotLimits()) return;

  // Determine the quantity to be plotted.
  std::string title;
  bool bfield = false;
  const Parameter par = GetPar(option, title, bfield);

  auto eval = [this, par, wfield, bfield, electrode, t](double* u, double* /*p*/) {
      // Transform to global coordinates.
      const double x = m_proj[0][0] * u[0] + m_proj[1][0] * u[1] + m_proj[2][0];
      const double y = m_proj[0][1] * u[0] + m_proj[1][1] * u[1] + m_proj[2][1];
      const double z = m_proj[0][2] * u[0] + m_proj[1][2] * u[1] + m_proj[2][2];
      return wfield ? Wfield(x, y, z, par, electrode, t) :
             bfield ? Bfield(x, y, z, par) : Efield(x, y, z, par);
  };
  const std::string fname = FindUnusedFunctionName("f2D");
  TF2 f2(fname.c_str(), eval, m_xMinPlot, m_xMaxPlot, m_yMinPlot, m_yMaxPlot, 0);

  // Set the x-y range.
  f2.SetRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);

  // Set the z-range.
  double zmin = m_vmin;
  double zmax = m_vmax;
  if (wfield) {
    if (contour) {
      title = "Contours of the weighting " + title;
    } else {
      title = "Weighting " + title;
    }
    if (m_useAutoRange) {
      SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, &f2, 
                  zmin, zmax);
    } else if (par == Parameter::Potential) {
      zmin = 0.;
      zmax = 1.;
    } else {
      zmin = m_wmin;
      zmax = m_wmax;
    }
  } else if (bfield) {
    if (contour) {
      title = "Contours of the magnetic " + title;
    } else {
      title = "Magnetic " + title;
    }
    if (m_useAutoRange) {
      SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                  &f2, zmin, zmax);
    } else {
      zmin = m_bmin;
      zmax = m_bmax;
    }
  } else {
    if (contour) {
      title = "Contours of the electric " + title;
    } else {
      title = "Electric " + title;
    }
    if (par == Parameter::Potential) {
      if (m_useAutoRange) {
        if (m_component) {
          if (m_samplePotential || !m_component->GetVoltageRange(zmin, zmax)) {
            SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                        &f2, zmin, zmax);
          }
        } else if (m_sensor) {
          if (m_samplePotential || !m_sensor->GetVoltageRange(zmin, zmax)) {
            SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                        &f2, zmin, zmax);
          }
        }
      } else {
        zmin = m_vmin;
        zmax = m_vmax;
      }
    } else {
      if (m_useAutoRange) {
        SampleRange(m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot, 
                    &f2, zmin, zmax);
      } else {
        zmin = m_emin;
        zmax = m_emax;
      }
    }
  }
  f2.SetMinimum(zmin);
  f2.SetMaximum(zmax);

  // Set the contours if requested.
  if (contour) {
    std::vector<double> level(m_nContours, 0.);
    if (m_nContours > 1) {
      const double step = (zmax - zmin) / (m_nContours - 1.);
      for (unsigned int i = 0; i < m_nContours; ++i) {
        level[i] = zmin + i * step;
      }
    } else {
      level[0] = 0.5 * (zmax + zmin);
    }
    if (m_debug) {
      std::cout << m_className << "::Draw2d:\n"
                << "    Number of contours: " << m_nContours << "\n";
      for (unsigned int i = 0; i < m_nContours; ++i) {
        std::cout << "        Level " << i << " = " << level[i] << "\n";
      }
    }
    f2.SetContour(m_nContours, level.data());
  }

  // Set the resolution.
  f2.SetNpx(m_nSamples2dX);
  f2.SetNpy(m_nSamples2dY);

  // Set the labels.
  std::string labels = ";" + LabelX() + ";" + LabelY();
  f2.SetTitle(labels.c_str());

  auto pad = GetCanvas();
  pad->cd();
  pad->SetTitle(title.c_str());
  f2.DrawCopy(drawopt.c_str());
  gPad->SetRightMargin(0.15);
  gPad->Update();
}

void ViewField::DrawProfile(const double x0, const double y0, const double z0,
                            const double x1, const double y1, const double z1,
                            const std::string& option, 
                            const bool wfield, const std::string& electrode,
                            const bool normalised) {
  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::DrawProfile:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }

  // Check the distance between the two points.
  double dx = x1 - x0;
  double dy = y1 - y0;
  double dz = z1 - z0;
  if (dx * dx + dy * dy + dz * dz <= 0.) {
    std::cerr << m_className << "::DrawProfile:\n"
              << "    Start and end points coincide.\n";
    return;
  }

  // Determine the quantity to be plotted.
  std::string title;
  bool bfield = false;
  const Parameter par = GetPar(option, title, bfield);

  double t0 = 0.;
  double t1 = 1.;
  unsigned int dir = 3;
  if (fabs(dy) + fabs(dz) < 1.e-6 * fabs(dx)) {
    t0 = x0;
    t1 = x1;
    dir = 0;
  } else if (fabs(dx) + fabs(dz) < 1.e-6 * fabs(dy)) {
    t0 = y0;
    t1 = y1;
    dir = 1;
  } else if (fabs(dx) + fabs(dy) < 1.e-6 * fabs(dz)) {
    t0 = z0;
    t1 = z1;
    dir = 2;
  } else if (!normalised) {
    t1 = sqrt(dx * dx + dy * dy + dz * dz);
    dx /= t1;
    dy /= t1;
    dz /= t1; 
  }

  auto eval = [this, par, wfield, bfield, electrode, dir, 
               x0, y0, z0, dx, dy, dz](double* u, double* /*p*/) {
    // Get the position.
    const double t = u[0];
    double x = dir == 0 ? t : x0;
    double y = dir == 1 ? t : y0;
    double z = dir == 2 ? t : z0;
    if (dir > 2) {
      x += t * dx;
      y += t * dy;
      z += t * dz;
    }
    return wfield ? Wfield(x, y, z, par, electrode) : 
           bfield ? Bfield(x, y, z, par) : Efield(x, y, z, par);
  };

  const std::string fname = FindUnusedFunctionName("fProfile");
  TF1 f1(fname.c_str(), eval, t0, t1, 0);

  double fmin = m_vmin;
  double fmax = m_vmax;
  if (wfield) {
    title = "weighting " + title;
    if (par == Parameter::Potential) {
      if (m_useAutoRange && m_samplePotential) {
        SampleRange(&f1, fmin, fmax);
      } else {
        fmin = 0.;
        fmax = 1.;
      }
    } else {
      if (m_useAutoRange) {
        SampleRange(&f1, fmin, fmax);
      } else {
        fmin = m_wmin;
        fmax = m_wmax;
      }
    }
  } else if (bfield) {
    title = "magnetic " + title;
    if (m_useAutoRange) {
      SampleRange(&f1, fmin, fmax);
    } else {
      fmin = m_bmin;
      fmax = m_bmax;
    }
  } else {
    title = "electric " + title;
    if (par == Parameter::Potential) {
      if (m_useAutoRange) {
        if (m_component) {
          if (m_samplePotential || !m_component->GetVoltageRange(fmin, fmax)) {
            SampleRange(&f1, fmin, fmax);
          }
        } else if (m_sensor) {
          if (m_samplePotential || !m_sensor->GetVoltageRange(fmin, fmax)) {
            SampleRange(&f1, fmin, fmax);
          }
        }
      } else {
        fmin = m_vmin;
        fmax = m_vmax;
      }
    } else {
      if (m_useAutoRange) {
        SampleRange(&f1, fmin, fmax);
      } else {
        fmin = m_emin;
        fmax = m_emax;
      }
    }
  }
  f1.SetMinimum(fmin);
  f1.SetMaximum(fmax);

  std::string labels = ";normalised distance;";
  if (dir == 0) {
    labels = ";#it{x} [cm];";
  } else if (dir == 1) {
    labels = ";#it{y} [cm];";
  } else if (dir == 2) {
    labels = ";#it{z} [cm];";
  } else if (!normalised) {
    labels = ";distance [cm];";
  }
  if (bfield) {
    labels += "#it{B}";
    if (par == Parameter::Bx) {
      labels += "_{x}";
    } else if (par == Parameter::By) {
      labels += "_{y}";
    } else if (par == Parameter::Bz) {
      labels += "_{z}";
    }
    labels += " [T]";
  } else if (par == Parameter::Potential) {
    labels += "#phi";
    if (wfield) {
      labels += "_w";
    } else {
      labels += " [V]";
    }
  } else {
    labels += "#it{E}";
    if (wfield) {
      labels += "_{w";
      if (par != Parameter::Emag) labels += ",";
    } else if (par != Parameter::Emag) {
      labels += "_{";
    }
    if (par == Parameter::Ex) {
      labels += "x";
    } else if (par == Parameter::Ey) {
      labels += "y";
    } else if (par == Parameter::Ez) {
      labels += "z";
    }
    if (wfield || par != Parameter::Emag) labels += "}";
    if (wfield) {
      labels += " [1/cm]";
    } else {
      labels += " [V/cm]";
    }
  }
  f1.SetTitle(labels.c_str());
  f1.SetNpx(m_nSamples1d);

  auto pad = GetCanvas();
  pad->cd();
  title = "Profile plot of the " + title;
  pad->SetTitle(title.c_str());
  f1.DrawCopy();
  gPad->Update();
}

bool ViewField::SetPlotLimits() {

  if (m_userPlotLimits) return true;
  double xmin = 0., ymin = 0., xmax = 0., ymax = 0.;
  if (m_userBox) {
    if (PlotLimitsFromUserBox(xmin, ymin, xmax, ymax)) {
      m_xMinPlot = xmin;
      m_xMaxPlot = xmax;
      m_yMinPlot = ymin;
      m_yMaxPlot = ymax;
      return true;
    } 
  }
  // Try to get the area/bounding box from the sensor/component.
  bool ok = false;

  if (m_sensor) {
    ok = PlotLimits(m_sensor, xmin, ymin, xmax, ymax);
  } else {
    ok = PlotLimits(m_component, xmin, ymin, xmax, ymax);
  } 
  if (ok) {
    m_xMinPlot = xmin;
    m_xMaxPlot = xmax;
    m_yMinPlot = ymin;
    m_yMaxPlot = ymax;
  }
  return ok;
}

double ViewField::Efield(const double x, const double y, const double z,
                         const Parameter par) const {

  // Compute the field.
  double ex = 0., ey = 0., ez = 0., volt = 0.;
  int status = 0;
  Medium* medium = nullptr;
  if (!m_sensor) {
    m_component->ElectricField(x, y, z, ex, ey, ez, volt, medium, status);
  } else {
    m_sensor->ElectricField(x, y, z, ex, ey, ez, volt, medium, status);
  }
  if (m_useStatus && status != 0) return m_vBkg;
  switch (par) {
    case Parameter::Potential:
      return volt;
    case Parameter::Emag:
      return sqrt(ex * ex + ey * ey + ez * ez);
    case Parameter::Ex:
      return ex;
    case Parameter::Ey:
      return ey;
    case Parameter::Ez:
      return ez;
    default:
      break;
  }
  return volt;
}

double ViewField::Wfield(const double x, const double y, const double z,
                         const Parameter par,
                         const std::string& electrode, const double t) const {

  if (par == Parameter::Potential) {
    double v = 0.;
    if (m_sensor) {
      v = m_sensor->WeightingPotential(x, y, z, electrode);
      if (t > 0.) {
        v += m_sensor->DelayedWeightingPotential(x, y, z, t, electrode);
      }
    } else {
      v = m_component->WeightingPotential(x, y, z, electrode);
      if (t > 0.) {
        v += m_component->DelayedWeightingPotential(x, y, z, t, electrode);
      }
    }
    return std::max(0., v);
  }
  // TODO: delayed component.
  double ex = 0., ey = 0., ez = 0.;
  if (!m_sensor) {
    m_component->WeightingField(x, y, z, ex, ey, ez, electrode);
  } else {
    m_sensor->WeightingField(x, y, z, ex, ey, ez, electrode);
  }

  switch (par) {
    case Parameter::Emag:
      return sqrt(ex * ex + ey * ey + ez * ez);
    case Parameter::Ex:
      return ex;
    case Parameter::Ey:
      return ey;
    case Parameter::Ez:
      return ez;
    default:
      break;
  }
  return 0.;
} 

double ViewField::Bfield(const double x, const double y, const double z,
                         const Parameter par) const {

  // Compute the field.
  double bx = 0., by = 0., bz = 0.;
  int status = 0;
  if (!m_sensor) {
    m_component->MagneticField(x, y, z, bx, by, bz, status);
  } else {
    m_sensor->MagneticField(x, y, z, bx, by, bz, status);
  }
  if (m_useStatus && status != 0) return 0.;
  switch (par) {
    case Parameter::Bmag:
      return sqrt(bx * bx + by * by + bz * bz);
    case Parameter::Bx:
      return bx;
    case Parameter::By:
      return by;
    case Parameter::Bz:
      return bz;
    default:
      break;
  }
  return 0.;
}

void ViewField::PlotFieldLines(const std::vector<double>& x0,
                               const std::vector<double>& y0,
                               const std::vector<double>& z0,
                               const bool electron, const bool axis,
                               const short col) {
  
  if (x0.empty() || y0.empty() || z0.empty()) return;
  const size_t nLines = x0.size();
  if (y0.size() != nLines || z0.size() != nLines) {
    std::cerr << m_className << "::PlotLines:\n"
              << "    Mismatch in number of x/y/z coordinates.\n";
    return;
  }

  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::PlotFieldLines:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }

  auto pad = GetCanvas();
  pad->cd();
  pad->SetTitle("Field lines");
  // Determine the x-y range.
  const bool rangeSet = RangeSet(pad);
  if (axis || !rangeSet) {
    // Determine the plot limits.
    if (!SetPlotLimits()) {
      std::cerr << m_className << "::PlotFieldLines:\n"
                << "     Could not determine the plot limits.\n";
      return;
    }
  }
  if (axis) {
    auto frame = pad->DrawFrame(m_xMinPlot, m_yMinPlot,
                                m_xMaxPlot, m_yMaxPlot);
    frame->GetXaxis()->SetTitle(LabelX().c_str());
    frame->GetYaxis()->SetTitle(LabelY().c_str());
    pad->Update();
  } else if (!rangeSet) {
    SetRange(pad, m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);
  } 

  DriftLineRKF drift;
  Sensor sensor;
  if (m_sensor) {
    drift.SetSensor(m_sensor);
  } else {
    double xmin = 0., ymin = 0., zmin = 0.;
    double xmax = 0., ymax = 0., zmax = 0.;
    if (!m_component->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)) {
      if (m_userBox) {
        sensor.SetArea(m_xMinBox, m_yMinBox, m_zMinBox, 
                       m_xMaxBox, m_yMaxBox, m_zMaxBox);
      }
    }
    sensor.AddComponent(m_component);
    drift.SetSensor(&sensor);
  }
  const double lx = 0.01 * fabs(m_xMaxPlot - m_xMinPlot);
  const double ly = 0.01 * fabs(m_yMaxPlot - m_yMinPlot);
  drift.SetMaximumStepSize(std::min(lx, ly));
  for (size_t i = 0; i < nLines; ++i) {
    std::vector<std::array<float, 3> > xl;
    if (!drift.FieldLine(x0[i], y0[i], z0[i], xl, electron)) continue;
    DrawLine(xl, col, 1);
  }
  pad->Update();
}

bool ViewField::EqualFluxIntervals(
    const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1,
    std::vector<double>& xf, std::vector<double>& yf,
    std::vector<double>& zf,
    const unsigned int nPoints) const {

  if (nPoints < 2) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << "    Number of flux lines must be > 1.\n";
    return false;
  }

  // Set integration intervals.
  constexpr unsigned int nV = 5;
  // Compute the inplane vector normal to the track.
  const double xp = (y1 - y0) * m_plane[2] - (z1 - z0) * m_plane[1];
  const double yp = (z1 - z0) * m_plane[0] - (x1 - x0) * m_plane[2];
  const double zp = (x1 - x0) * m_plane[1] - (y1 - y0) * m_plane[0];
  // Compute the total flux, accepting positive and negative parts.
  double q = 0.;
  if (m_component) {
    q = m_component->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                       xp, yp, zp, 20 * nV, 0);
  } else {
    q = m_sensor->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                    xp, yp, zp, 20 * nV, 0);
  }
  const int isign = q > 0 ? +1 : -1;
  if (m_debug) {
    std::cout << m_className << "::EqualFluxIntervals:\n";
    std::printf("    Total flux: %15.e8\n", q);
  }
  // Compute the 1-sided flux in a number of steps.
  double fsum = 0.;
  unsigned int nOtherSign = 0;
  double s0 = -1.;
  double s1 = -1.;
  constexpr size_t nSteps = 1000;
  std::array<double, nSteps> sTab;
  std::array<double, nSteps> fTab;
  constexpr double ds = 1. / nSteps;
  const double dx = (x1 - x0) * ds;
  const double dy = (y1 - y0) * ds;
  const double dz = (z1 - z0) * ds;
  for (size_t i = 0; i < nSteps; ++i) {
    const double x = x0 + i * dx;
    const double y = y0 + i * dy;
    const double z = z0 + i * dz;
    if (m_component) {
      q = m_component->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                         xp, yp, zp, nV, isign);
    } else {
      q = m_sensor->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                      xp, yp, zp, nV, isign);
    }
    sTab[i] = (i + 1) * ds;
    if (q > 0) {
      fsum += q;
      if (s0 < -0.5) s0 = i * ds;
      s1 = (i + 1) * ds;
    }
    if (q < 0) ++nOtherSign;
    fTab[i] = fsum;
  }
  if (m_debug) {
    std::printf("    Used flux: %15.8e V. Start: %10.3f End: %10.3f\n",
                fsum, s0, s1);
  }
  // Make sure that the sum is positive.
  if (fsum <= 0) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << "    1-Sided flux integral is not > 0.\n";
    return false;
  } else if (s0 < -0.5 || s1 < -0.5 || s1 <= s0) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << "    No flux interval without sign change found.\n";
    return false;
  } else if (nOtherSign > 0) {
    std::cerr << m_className << "::EqualFluxIntervals:\n"
              << " The flux changes sign over the line.\n";
  }
  // Normalise the flux.
  const double scale = (nPoints - 1) / fsum;
  for (size_t i = 0; i < nSteps; ++i) fTab[i] *= scale;

  // Compute new cluster position.
  for (size_t i = 0; i < nPoints; ++i) {
    double s = std::min(s1, std::max(s0, Interpolate(sTab, fTab, i)));
    xf.push_back(x0 + s * (x1 - x0));
    yf.push_back(y0 + s * (y1 - y0));
    zf.push_back(z0 + s * (z1 - z0));
  }
  return true;
}

bool ViewField::FixedFluxIntervals(
    const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1,
    std::vector<double>& xf, std::vector<double>& yf,
    std::vector<double>& zf, const double interval) const {

  if (interval <= 0.) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    Flux interval must be > 0.\n";
    return false;
  }

  // Set integration intervals.
  constexpr unsigned int nV = 5;
  // Compute the inplane vector normal to the track.
  const double xp = (y1 - y0) * m_plane[2] - (z1 - z0) * m_plane[1];
  const double yp = (z1 - z0) * m_plane[0] - (x1 - x0) * m_plane[2];
  const double zp = (x1 - x0) * m_plane[1] - (y1 - y0) * m_plane[0];
  // Compute the total flux, accepting positive and negative parts.
  double q = 0.;
  if (m_component) {
    q = m_component->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                       xp, yp, zp, 20 * nV, 0);
  } else {
    q = m_sensor->IntegrateFluxLine(x0, y0, z0, x1, y1, z1, 
                                    xp, yp, zp, 20 * nV, 0);
  }
  const int isign = q > 0 ? +1 : -1;
  if (m_debug) {
    std::cout << m_className << "::FixedFluxIntervals:\n";
    std::printf("    Total flux: %15.8e V\n", q);
  }
  // Compute the 1-sided flux in a number of steps.
  double fsum = 0.;
  unsigned int nOtherSign = 0;
  double s0 = -1.;
  constexpr size_t nSteps = 1000;
  std::array<double, nSteps> sTab;
  std::array<double, nSteps> fTab;
  constexpr double ds = 1. / nSteps;
  const double dx = (x1 - x0) * ds;
  const double dy = (y1 - y0) * ds;
  const double dz = (z1 - z0) * ds;
  for (size_t i = 0; i < nSteps; ++i) {
    const double x = x0 + i * dx;
    const double y = y0 + i * dy;
    const double z = z0 + i * dz;
    if (m_component) {
      q = m_component->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                         xp, yp, zp, nV, isign);
    } else {
      q = m_sensor->IntegrateFluxLine(x, y, z, x + dx, y + dy, z + dz, 
                                      xp, yp, zp, nV, isign);
    }
    sTab[i] = (i + 1) * ds;
    if (q > 0) {
      fsum += q;
      if (s0 < -0.5) s0 = i * ds;
    }
    if (q < 0) ++nOtherSign;
    fTab[i] = fsum;
  }
  // Make sure that the sum is positive.
  if (m_debug) {
    std::printf("    Used flux: %15.8e V. Start offset: %10.3f\n", fsum, s0);
  }
  if (fsum <= 0) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    1-Sided flux integral is not > 0.\n";
    return false;
  } else if (s0 < -0.5) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    No flux interval without sign change found.\n";
    return false;
  } else if (nOtherSign > 0) {
    std::cerr << m_className << "::FixedFluxIntervals:\n"
              << "    Warning: The flux changes sign over the line.\n";
  }

  double f = 0.;
  while (f < fsum) {
    const double s = Interpolate(sTab, fTab, f);
    f += interval;
    xf.push_back(x0 + s * (x1 - x0));
    yf.push_back(y0 + s * (y1 - y0));
    zf.push_back(z0 + s * (z1 - z0));
  }
  return true;
}
}
