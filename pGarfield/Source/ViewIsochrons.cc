#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

#include <TAxis.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TH1F.h>

#include "Garfield/Sensor.hh"
#include "Garfield/Component.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewIsochrons.hh"

namespace {

double Interpolate(const std::vector<double>& y,
                   const std::vector<double>& x, const double xx) {

  const double tol = 1.e-6 * fabs(x.back() - x.front());
  if (xx < x.front()) return y.front();
  const auto it1 = std::upper_bound(x.cbegin(), x.cend(), xx);
  if (it1 == x.cend()) return y.back();
  const auto it0 = std::prev(it1);
  const double dx = (*it1 - *it0);
  if (dx < tol) return y[it0 - x.cbegin()];
  const double f = (xx - *it0) / dx;
  return y[it0 - x.cbegin()] * (1. - f) + f * y[it1 - x.cbegin()];
}

bool OnLine(const double x1, const double y1, const double x2, const double y2,
            const double u, const double v) {
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v)});
  epsx = std::max(1.e-10, epsx);
  epsy = std::max(1.e-10, epsy);

  if ((fabs(x1 - u) <= epsx && fabs(y1 - v) <= epsy) ||
      (fabs(x2 - u) <= epsx && fabs(y2 - v) <= epsy)) {
    // Point to be examined coincides with start or end.
    return true;
  } else if (fabs(x1 - x2) <= epsx && fabs(y1 - y2) <= epsy) {
    // The line (x1, y1) to (x2, y2) is in fact a point.
    return false;
  }
  double xc = 0., yc = 0.;
  if (fabs(u - x1) + fabs(v - y1) < fabs(u - x2) + fabs(v - y2)) {
    // (u, v) is nearer to (x1, y1).
    const double dx = (x2 - x1);
    const double dy = (y2 - y1);
    const double xl = ((u - x1) * dx + (v - y1) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x1;
      yc = y1;
    } else if (xl > 1.) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + xl * dx;
      yc = y1 + xl * dy;
    }
  } else {
    // (u, v) is nearer to (x2, y2).
    const double dx = (x1 - x2);
    const double dy = (y1 - y2);
    const double xl = ((u - x2) * dx + (v - y2) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x2;
      yc = y2;
    } else if (xl > 1.) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + xl * dx;
      yc = y2 + xl * dy;
    }
  }
  // See whether the point is on the line.
  if (fabs(u - xc) < epsx && fabs(v - yc) < epsy) {
    return true;
  }
  return false;
}

bool Crossing(const double x1, const double y1, const double x2,
              const double y2, const double u1, const double v1,
              const double u2, const double v2) {

  // Check for a point of one line located on the other line.
  if (OnLine(x1, y1, x2, y2, u1, v1) || OnLine(x1, y1, x2, y2, u2, v2) ||
      OnLine(u1, v1, u2, v2, x1, y1) || OnLine(u1, v1, u2, v2, x2, y2)) {
    return true;
  }
  // Matrix to compute the crossing point.
  std::array<std::array<double, 2>, 2> a;
  a[0][0] = y2 - y1;
  a[0][1] = v2 - v1;
  a[1][0] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u1), fabs(u2)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v1), fabs(v2)});
  epsx = std::max(epsx, 1.e-10);
  epsy = std::max(epsy, 1.e-10);
  if (fabs(det) < epsx * epsy) {
    // Parallel, non-touching.
    return false;
  }
  // Crossing, non-trivial lines: solve crossing equations.
  const double aux = a[1][1];
  a[1][1] = a[0][0] / det;
  a[0][0] = aux / det;
  a[1][0] = -a[1][0] / det;
  a[0][1] = -a[0][1] / det;
  // Compute crossing point.
  const double xc = a[0][0] * (x1 * y2 - x2 * y1) + a[1][0] * (u1 * v2 - u2 * v1);
  const double yc = a[0][1] * (x1 * y2 - x2 * y1) + a[1][1] * (u1 * v2 - u2 * v1);
  // See whether the crossing point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc)) {
    // Intersecting lines.
    return true;
  }
  // Crossing point not on both lines.
  return false;
}

}

namespace Garfield {

ViewIsochrons::ViewIsochrons() : ViewBase("ViewIsochrons") { }

void ViewIsochrons::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }

  m_sensor = s;
  m_component = nullptr;
}

void ViewIsochrons::SetComponent(Component* c) {
  if (!c) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }

  m_component = c;
  m_sensor = nullptr;
}

void ViewIsochrons::SetAspectRatioSwitch(const double ar) {

  if (ar < 0.) {
    std::cerr << m_className << "::SetAspectRatioSwitch: Value must be > 0.\n";
    return;
  }
  m_aspectRatio = ar;
}

void ViewIsochrons::SetLoopThreshold(const double thr) {

  if (thr < 0. || thr > 1.) {
    std::cerr << m_className << "::SetLoopThreshold:\n"
              << "    Value must be between 0 and 1.\n";
    return;
  }
  m_loopThreshold = thr;
}

void ViewIsochrons::SetConnectionThreshold(const double thr) {

  if (thr < 0. || thr > 1.) {
    std::cerr << m_className << "::SetConnectionThreshold:\n"
              << "    Value must be between 0 and 1.\n";
    return;
  }
  m_connectionThreshold = thr;
}

void ViewIsochrons::PlotIsochrons(const double tstep,
    const std::vector<std::array<double, 3> >& points, const bool rev, 
    const bool colour, const bool markers, const bool plotDriftLines) {

  if (!m_sensor && !m_component) {
    std::cerr << m_className << "::PlotIsochrons:\n"
              << "    Neither sensor nor component are defined.\n";
    return;
  }
  if (tstep <= 0.) {
    std::cerr << m_className << "::PlotIsochrons: Time step must be > 0.\n";
    return;
  }
  if (points.empty()) {
    std::cerr << m_className << "::PlotIsochrons:\n"
              << "    No starting points provided.\n";
    return;
  }
  if (!SetPlotLimits()) return;
  auto canvas = GetCanvas();
  canvas->cd();
  canvas->SetTitle("Isochrons");
  auto frame = canvas->DrawFrame(m_xMinPlot, m_yMinPlot, 
                                 m_xMaxPlot, m_yMaxPlot);
  frame->GetXaxis()->SetTitle(LabelX().c_str());
  frame->GetYaxis()->SetTitle(LabelY().c_str());
  canvas->Update();

  //-----------------------------------------------------------------------
  //   DRFEQT - The main routine (DRFEQT) accumulates equal drift time data
  //   DRFEQP   which is plotted as a set of contours in the entry DRFEQP.
  //-----------------------------------------------------------------------
  std::vector<std::vector<std::array<double, 3> > > driftLines;
  std::vector<std::array<double, 3> > startPoints;  
  std::vector<std::array<double, 3> > endPoints;  
  std::vector<int> statusCodes;
  // Accumulate drift lines.
  ComputeDriftLines(tstep, points, driftLines, startPoints, endPoints, 
                    statusCodes, rev);
  const unsigned int nDriftLines = driftLines.size();
  if (nDriftLines < 2) {
    std::cerr << m_className << "::PlotIsochrons: Too few drift lines.\n";
    return;
  }
  // Keep track of the largest number of contours.
  std::size_t nContours = 0;
  for (const auto& driftLine : driftLines) {
    nContours = std::max(nContours, driftLine.size());
  }
  if (nContours == 0) {
    std::cerr << m_className << "::PlotIsochrons: No contour lines.\n";
    return;
  }

  std::set<int> allStats;
  for (const auto stat : statusCodes) allStats.insert(stat);

  // DRFEQP
  if (m_debug) {
    std::cout << m_className << "::PlotIsochrons:\n"
              << "    Drawing " << nContours << " contours, "
              << nDriftLines << " drift lines.\n";
    std::printf("    Connection threshold:   %10.3f\n", m_connectionThreshold);
    std::printf("    Aspect ratio threshold: %10.3f\n", m_aspectRatio);
    std::printf("    Loop closing threshold: %10.3f\n", m_loopThreshold);
    if (m_sortContours) {
      std::cout << "    Sort contours.\n";
    } else {
      std::cout << "    Do not sort contours.\n";
    }
    if (m_checkCrossings) {
      std::cout << "    Check for crossings.\n";
    } else {
      std::cout << "    Do not check for crossings.\n";
    }
    if (markers) {
      std::cout << "    Mark isochron points.\n";
    } else {
      std::cout << "    Draw isochron lines.\n";
    }
  }

  // Loop over the equal time contours.
  TGraph graph;
  graph.SetLineColor(kGray + 2);
  graph.SetMarkerColor(kGray + 2);
  graph.SetLineStyle(m_lineStyle);
  graph.SetMarkerStyle(m_markerStyle);
  const double colRange = double(gStyle->GetNumberOfColors()) / nContours;
  for (unsigned int ic = 0; ic < nContours; ++ic) {
    if (colour) {
      const auto col = gStyle->GetColorPalette(int((ic + 0.99) * colRange));
      graph.SetLineColor(col);
      graph.SetMarkerColor(col);
    }
    for (int stat : allStats) {
      std::vector<std::pair<std::array<double, 4>, unsigned int> > contour;
      // Loop over the drift lines, picking up the points when OK.
      for (unsigned int k = 0; k < nDriftLines; ++k) {
        const auto& dl = driftLines[k]; 
        // Reject any undesirable combinations.
        if (statusCodes[k] != stat || ic >= dl.size()) continue;
        // Add the point to the contour line.
        std::array<double, 4> point = {dl[ic][0], dl[ic][1], dl[ic][2], 0.};
        contour.push_back(std::make_pair(point, k)); 
      }
      // Skip the plot of this contour if there are no points.
      if (contour.empty()) continue;
      bool circle = false;
      // If requested, sort the points on the contour line.
      if (m_sortContours && !markers) SortContour(contour, circle);
      // Plot this contour.
      if (markers) {
        // Simply mark the contours if this was requested.
        std::vector<double> xp;
        std::vector<double> yp;
        std::vector<double> zp;
        for (const auto& point : contour) {
          const double x = point.first[0];
          const double y = point.first[1];
          const double z = point.first[2];
          xp.push_back(m_proj[0][0] * x + m_proj[1][0] * y + z * m_plane[0]);
          yp.push_back(m_proj[0][1] * x + m_proj[1][1] * y + z * m_plane[1]);
          zp.push_back(m_proj[0][2] * x + m_proj[1][2] * y + z * m_plane[2]);
        }
        graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Psame");
        continue;
      }
      // Regular plotting.
      const double tolx = (m_xMaxPlot - m_xMinPlot) * m_connectionThreshold;
      const double toly = (m_yMaxPlot - m_yMinPlot) * m_connectionThreshold;
      // Flag to keep track if the segment is interrupted by a drift line
      // or if it is too long.
      bool gap = false;
      // Coordinates to be plotted.
      std::vector<double> xp;
      std::vector<double> yp;
      std::vector<double> zp;
      const unsigned int nP = contour.size();
      for (unsigned int i = 0; i < nP; ++i) { 
        gap = false;
        const auto x0 = contour[i].first[0];
        const auto y0 = contour[i].first[1];
        const auto z0 = contour[i].first[2];
        xp.push_back(m_proj[0][0] * x0 + m_proj[1][0] * y0 + z0 * m_plane[0]);
        yp.push_back(m_proj[0][1] * x0 + m_proj[1][1] * y0 + z0 * m_plane[1]);
        zp.push_back(m_proj[0][2] * x0 + m_proj[1][2] * y0 + z0 * m_plane[2]);
        if (i == nP - 1) break; 
        const auto x1 = contour[i + 1].first[0];
        const auto y1 = contour[i + 1].first[1];
        // Reject contour segments which are long compared with AREA.
        if (fabs(x1 - x0) > tolx || fabs(y1 - y0) > toly) gap = true;
        // Get the indices of the drift lines corresponding 
        // to these two points on the contour line.
        const auto i0 = contour[i].second;
        const auto i1 = contour[i + 1].second; 
        // Set the BREAK flag if it crosses some stored drift line segment.
        if (m_checkCrossings && !gap) {
          for (unsigned int k = 0; k < nDriftLines; ++k) {
            const auto& dl = driftLines[k];
            for (unsigned int jc = 0; jc < dl.size(); ++jc) {
              if ((i0 == k || i1 == k) && (jc == ic || jc + 1 == ic)) {
                continue;
              }
              const auto& p0 = dl[jc];
              const auto& p1 = jc == dl.size() - 1 ? endPoints[k] : dl[jc + 1];
              if (Crossing(p0[0], p0[1], p1[0], p1[1], x0, y0, x1, y1)) {
                gap = true;
                break;
              }
            }
            if (gap) break;
            if ((i0 == k || i1 == k) && ic == 0) continue;
            const auto& p0 = startPoints[k];
            if (Crossing(p0[0], p0[1], dl[0][0], dl[0][1], 
                         x0, y0, x1, y1)) {
              gap = true;
              break;
            }
          }
        }
        // If there has been a break, plot what we have already.
        if (gap) {
          if (xp.size() > 1) {
            // Plot line.
            graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Lsame");
          } else {
            // Plot markers.
            graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Psame");
          }
          xp.clear();
          yp.clear();
          zp.clear();
        }
      }
      // Plot the remainder.
      if (xp.size() > 1) {
        graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Lsame");
      } else if (!xp.empty()) {
        graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Psame");
      }
    }
  }

  gPad->Update();
  if (!plotDriftLines) return;

  graph.SetLineStyle(1);
  if (m_particle == Particle::Electron) {
    graph.SetLineColor(kOrange - 3);
  } else {
    graph.SetLineColor(kRed + 1);
  }
  for (unsigned int i = 0; i < nDriftLines; ++i) {
    std::vector<double> xp;
    std::vector<double> yp;
    const double x0 = startPoints[i][0];
    const double y0 = startPoints[i][1];
    const double z0 = startPoints[i][2];
    xp.push_back(m_proj[0][0] * x0 + m_proj[1][0] * y0 + z0 * m_plane[0]);
    yp.push_back(m_proj[0][1] * x0 + m_proj[1][1] * y0 + z0 * m_plane[1]);
    for (const auto& point : driftLines[i]) {
      const double x = point[0];
      const double y = point[1];
      const double z = point[2];
      xp.push_back(m_proj[0][0] * x + m_proj[1][0] * y + z * m_plane[0]);
      yp.push_back(m_proj[0][1] * x + m_proj[1][1] * y + z * m_plane[1]);
    }
    const double x1 = endPoints[i][0];
    const double y1 = endPoints[i][1];
    const double z1 = endPoints[i][2];
    xp.push_back(m_proj[0][0] * x1 + m_proj[1][0] * y1 + z1 * m_plane[0]);
    yp.push_back(m_proj[0][1] * x1 + m_proj[1][1] * y1 + z1 * m_plane[1]);
    graph.DrawGraph(xp.size(), xp.data(), yp.data(), "Lsame");
  }
}

void ViewIsochrons::ComputeDriftLines(const double tstep,
  const std::vector<std::array<double, 3> >& points,
  std::vector<std::vector<std::array<double, 3> > >& driftLines,
  std::vector<std::array<double, 3> >& startPoints,
  std::vector<std::array<double, 3> >& endPoints,
  std::vector<int>& statusCodes, const bool rev) {

  DriftLineRKF drift;
  Sensor sensor;
  if (m_sensor) {
    drift.SetSensor(m_sensor);
  } else {
    sensor.AddComponent(m_component);
    if (m_userBox) {
      sensor.SetArea(m_xMinBox, m_yMinBox, m_zMinBox,
                     m_xMaxBox, m_yMaxBox, m_zMaxBox);
    }
    drift.SetSensor(&sensor);
  }
  const double lx = 0.1 * fabs(m_xMaxPlot - m_xMinPlot);
  const double ly = 0.1 * fabs(m_yMaxPlot - m_yMinPlot);
  drift.SetMaximumStepSize(std::min(lx, ly));
  drift.EnableSignalCalculation(false);
  for (const auto& point : points) {
    if (m_particle == Particle::Electron) {
      if (m_positive) {
        drift.DriftPositron(point[0], point[1], point[2], 0.);
      } else {
        drift.DriftElectron(point[0], point[1], point[2], 0.);
      }
    } else {
      if (m_positive) {
        drift.DriftIon(point[0], point[1], point[2], 0.);
      } else {
        drift.DriftNegativeIon(point[0], point[1], point[2], 0.);
      }
    }
    const unsigned int nu = drift.GetNumberOfDriftLinePoints();
    // Check that the drift line has enough points.
    if (nu < 3) continue;
    int status = 0;
    double xf = 0., yf = 0., zf = 0., tf = 0.;
    drift.GetEndPoint(xf, yf, zf, tf, status);
    // Find the number of points to be stored.
    const unsigned int nSteps = static_cast<unsigned int>(tf / tstep);
    if (nSteps == 0) continue;
    std::vector<double> xu(nu, 0.);
    std::vector<double> yu(nu, 0.);
    std::vector<double> zu(nu, 0.);
    std::vector<double> tu(nu, 0.);
    for (unsigned int i = 0; i < nu; ++i) {
      drift.GetDriftLinePoint(i, xu[i], yu[i], zu[i], tu[i]);
    }
    if (rev) {
      for (auto& t : tu) t = tf - t;
      std::reverse(std::begin(xu), std::end(xu)); 
      std::reverse(std::begin(yu), std::end(yu)); 
      std::reverse(std::begin(zu), std::end(zu)); 
      std::reverse(std::begin(tu), std::end(tu)); 
    }
    std::vector<std::array<double, 3> > tab;
    // Interpolate at regular time intervals.
    for (unsigned int i = 0; i < nSteps; ++i) {
      const double t = (i + 1) * tstep;
      // tab.push_back(PLACO3(Interpolate(xu, tu, t),
      //                      Interpolate(yu, tu, t),
      //                      Interpolate(zu, tu, t)));
      std::array<double, 3> step = {Interpolate(xu, tu, t),
                                    Interpolate(yu, tu, t),
                                    Interpolate(zu, tu, t)};
      tab.push_back(step);
    }
    driftLines.push_back(std::move(tab));
    std::array<double, 3> start = {xu[0], yu[0], zu[0]};
    std::array<double, 3> end = {xu[nu - 1], yu[nu - 1], zu[nu - 1]};
    startPoints.push_back(std::move(start));
    endPoints.push_back(std::move(end));
    // Store the drift line return code.
    if (rev) {
      statusCodes.push_back(status);
    } else {
      statusCodes.push_back(0);
    }
  }
}

bool ViewIsochrons::SetPlotLimits() {

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

void ViewIsochrons::SortContour(
    std::vector<std::pair<std::array<double, 4>, unsigned int> >& contour,
    bool& circle) {

  if (contour.size() < 2) return;
  // First compute the centre of gravity.
  double xcog = 0.;
  double ycog = 0.;
  for (const auto& point : contour) {
    xcog += point.first[0];
    ycog += point.first[1];
  }
  const double scale = 1. / contour.size();
  xcog *= scale;
  ycog *= scale;
  // Compute angles wrt to the centre of gravity and principal axes.
  double sxx = 0.;
  double sxy = 0.;
  double syy = 0.;
  for (const auto& point : contour) {
    const double dx = point.first[0] - xcog;
    const double dy = point.first[1] - ycog;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  const double theta = 0.5 * atan2(2 * sxy, sxx - syy);
  const double ct = cos(theta);
  const double st = sin(theta);
  // Evaluate dispersions around the principal axes.
  sxx = 0.;
  syy = 0.;
  for (const auto& point : contour) {
    const double dx = point.first[0] - xcog;
    const double dy = point.first[1] - ycog;
    sxx += fabs(+ct * dx + st * dy);
    syy += fabs(-st * dx + ct * dy);
  }
  // Decide whether this is more linear or more circular.
  if (fabs(sxx) > m_aspectRatio * fabs(syy) || 
      fabs(syy) > m_aspectRatio * fabs(sxx)) {
    circle = false;
  } else {
    circle = true;
  }
  // Set a sorting coordinate accordingly.
  for (auto& point : contour) {
    const double dx = point.first[0] - xcog;
    const double dy = point.first[1] - ycog;
    point.first[3] = circle ? atan2(dy, dx) : ct * dx + st * dy;
  }
  // Sort the points.
  std::sort(contour.begin(), contour.end(), 
            [](const std::pair<std::array<double, 4>, int>& p1, 
               const std::pair<std::array<double, 4>, int>& p2) {
                 return p1.first[3] < p2.first[3]; }
           );
  if (!circle) return;
  // For circles, perhaps add the first point to the end of the list.
  // Compute breakpoint, total distance and maximum distance.
  double dsum = 0.;
  double dmax = -1.;
  unsigned int imax = 0;
  const unsigned int nPoints = contour.size();
  for (unsigned int j = 0; j < nPoints; ++j) {
    const auto& p1 = contour[j];
    const auto& p0 = j > 0 ? contour[j - 1] : contour.back();
    const double dx = p1.first[0] - p0.first[0];
    const double dy = p1.first[1] - p0.first[1];
    const double d = sqrt(dx * dx + dy * dy);
    if (j > 0) dsum += d;
    if (dmax < d) {
      dmax = d;
      imax = j;
    }
  }
  if (dmax < m_loopThreshold * dsum) {
    // If a true loop, close it.
    contour.push_back(contour[0]);
  } else {
    circle = false;
    if (imax > 0) {
      // Shift the points to make a line.
      std::rotate(contour.begin(), contour.begin() + imax, contour.end());
    }
  }
}

}
