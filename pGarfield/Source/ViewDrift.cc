#include <iostream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>

#include <TGraph.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TAxis.h>
#include <TAxis3D.h>
#include <TH1F.h>
#include <TView3D.h>
#include <TVirtualViewer3D.h>

#include "Garfield/Plotting.hh"
#include "Garfield/ViewDrift.hh"

namespace Garfield {

ViewDrift::ViewDrift() : ViewBase("ViewDrift") {
  m_driftLines.reserve(1000);
  m_tracks.reserve(100);
  m_exc.reserve(1000);
  m_ion.reserve(1000);
  m_att.reserve(1000);
}

void ViewDrift::Clear() {
  m_driftLines.clear();
  m_tracks.clear();

  m_exc.clear();
  m_ion.clear();
  m_att.clear();
}

void ViewDrift::SetClusterMarkerSize(const double size) {
  if (size > 0.) {
    m_markerSizeCluster = size;
  } else {
    std::cerr << m_className << "::SetClusterMarkerSize: Size must be > 0.\n";
  }
}

void ViewDrift::SetCollisionMarkerSize(const double size) {
  if (size > 0.) {
    m_markerSizeCollision = size;
  } else {
    std::cerr << m_className << "::SetCollisionMarkerSize: Size must be > 0.\n";
  }
}

void ViewDrift::GetDriftLine(const size_t i, 
                             std::vector<std::array<float, 3> >& driftLine,
                             bool& electron) const {
  driftLine.clear();
  if (i >= m_driftLines.size()) return;
  std::copy(m_driftLines[i].first.begin(), m_driftLines[i].first.end(),
            std::back_inserter(driftLine));
  if (m_driftLines[i].second == Particle::Electron) {
    electron = true;
  } else {
    electron = false;
  }
}

void ViewDrift::NewDriftLine(const Particle particle, const size_t np, 
    size_t& id, const float x0, const float y0, const float z0) {
  std::lock_guard<std::mutex> guard(m_mutex);
  // Create a new drift line and add it to the list.
  std::array<float, 3> p = {x0, y0, z0};
  constexpr size_t minSize = 1;
  std::vector<std::array<float, 3> > dl(std::max(minSize, np), p);
  m_driftLines.push_back(std::make_pair(std::move(dl), particle));
  // Return the index of this drift line.
  id = m_driftLines.size() - 1;
}

void ViewDrift::AddPhoton(const float x0, const float y0, const float z0, 
                          const float x1, const float y1, const float z1) {
  std::lock_guard<std::mutex> guard(m_mutex);
  std::array<float, 3> p0 = {x0, y0, z0};
  std::array<float, 3> p1 = {x1, y1, z1};
  m_photons.push_back({p0, p1});
}

void ViewDrift::NewChargedParticleTrack(const size_t np, size_t& id,
                                        const float x0, const float y0,
                                        const float z0) {
  std::lock_guard<std::mutex> guard(m_mutex);
  // Create a new track and add it to the list.
  constexpr size_t minSize = 1;
  std::vector<std::array<float, 3> > track(std::max(minSize, np));
  track[0] = {x0, y0, z0};
  m_tracks.push_back(std::move(track));
  // Return the index of this track.
  id = m_tracks.size() - 1;
}

void ViewDrift::SetDriftLinePoint(const size_t iL, const size_t iP,
                                  const float x, const float y,
                                  const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  if (iL >= m_driftLines.size() || iP >= m_driftLines[iL].first.size()) {
    std::cerr << m_className << "::SetDriftLinePoint: Index out of range.\n";
    return;
  }
  m_driftLines[iL].first[iP] = {x, y, z};
}

void ViewDrift::AddDriftLinePoint(const size_t iL, const float x,
                                  const float y, const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  if (iL >= m_driftLines.size()) {
    std::cerr << m_className << "::AddDriftLinePoint: Index out of range.\n";
    return;
  }
  std::array<float, 3> p = {x, y, z};
  m_driftLines[iL].first.push_back(std::move(p));
}

void ViewDrift::SetTrackPoint(const size_t iL, const size_t iP,
                              const float x, const float y, const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  if (iL >= m_tracks.size() || iP >= m_tracks[iL].size()) {
    std::cerr << m_className << "::SetTrackPoint: Index out of range.\n";
    return;
  }
  m_tracks[iL][iP] = {x, y, z};
}

void ViewDrift::AddTrackPoint(const size_t iL, const float x,
                              const float y, const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  if (iL >= m_tracks.size()) {
    std::cerr << m_className << "::AddTrackPoint: Index out of range.\n";
    return;
  }
  std::array<float, 3> p = {x, y, z};
  m_tracks[iL].push_back(std::move(p));
}

void ViewDrift::AddExcitation(const float x, const float y, const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  std::array<float, 3> p = {x, y, z};
  m_exc.push_back(std::move(p));
}

void ViewDrift::AddIonisation(const float x, const float y, const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  std::array<float, 3> p = {x, y, z};
  m_ion.push_back(std::move(p));
}

void ViewDrift::AddAttachment(const float x, const float y, const float z) {
  std::lock_guard<std::mutex> guard(m_mutex);
  std::array<float, 3> p = {x, y, z};
  m_att.push_back(std::move(p));
}

void ViewDrift::Plot(const bool twod, const bool axis, const bool snapshot) {
  if (twod) {
    Plot2d(axis, snapshot);
  } else {
    Plot3d(axis, false, snapshot);
  }
}

void ViewDrift::Plot2d(const bool axis, const bool snapshot) {
  auto pad = GetCanvas();
  pad->cd();
  pad->SetTitle("Drift lines");
  // Check if the canvas range has already been set.
  const bool rangeSet = RangeSet(pad);
  if (axis || !rangeSet) {
    // Determine the plot limits.
    if (!SetPlotLimits2d()) {
      std::cerr << m_className << "::Plot2d:\n"
                << "     Could not determine the plot limits.\n";
      return;
    }
  }
  if (axis) {
    auto frame = pad->DrawFrame(m_xMinPlot, m_yMinPlot, 
                                m_xMaxPlot, m_yMaxPlot);
    frame->GetXaxis()->SetTitle(LabelX().c_str());
    frame->GetYaxis()->SetTitle(LabelY().c_str());
  } else if (!rangeSet) {
    SetRange(pad, m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);
  } 
  if (snapshot) {
    std::vector<std::array<float, 3> > electrons;
    std::vector<std::array<float, 3> > holes;
    std::vector<std::array<float, 3> > ions;
    for (const auto& driftLine : m_driftLines) {
      if (driftLine.second == Particle::Electron) {
        electrons.push_back(driftLine.first.back());
      } else if (driftLine.second == Particle::Hole) {
        holes.push_back(driftLine.first.back());
      } else {
        ions.push_back(driftLine.first.back());
      }
    }
    DrawMarkers2d(electrons, m_colElectron, m_markerSizeCollision);
    DrawMarkers2d(holes, m_colHole, m_markerSizeCollision);
    DrawMarkers2d(ions, m_colIon, m_markerSizeCollision);
  } else {
    for (const auto& driftLine : m_driftLines) {
      const short lw = 1;
      if (driftLine.second == Particle::Electron) {
        DrawLine(driftLine.first, m_colElectron, lw);
      } else if (driftLine.second == Particle::Hole) {
        DrawLine(driftLine.first, m_colHole, lw);
      } else {
        DrawLine(driftLine.first, m_colIon, lw);
      }
    }
  }
  gPad->Update();

  for (const auto& track : m_tracks) {
    DrawLine(track, m_colTrack, 2);
    if (!m_drawClusters) continue;
    DrawMarkers2d(track, m_colTrack, m_markerSizeCluster);
  }

  TGraph gr;
  gr.SetLineColor(m_colPhoton);
  gr.SetLineStyle(2);
  for (const auto& photon : m_photons) {
    float xp0 = 0., yp0 = 0.;
    float xp1 = 0., yp1 = 0.;
    ToPlane(photon[0][0], photon[0][1], photon[0][2], xp0, yp0);
    ToPlane(photon[1][0], photon[1][1], photon[1][2], xp1, yp1);
    std::vector<float> xgr = {xp0, xp1};
    std::vector<float> ygr = {yp0, yp1};
    gr.DrawGraph(2, xgr.data(), ygr.data(), "Lsame"); 
  }

  if (!m_exc.empty()) {
    DrawMarkers2d(m_exc, m_colExcitation, m_markerSizeCollision);
  } 
  if (!m_ion.empty()) {
    DrawMarkers2d(m_ion, m_colIonisation, m_markerSizeCollision);
  }
  if (!m_att.empty()) {
    DrawMarkers2d(m_att, m_colAttachment, m_markerSizeCollision);
  }
 
  gPad->Update();
}

void ViewDrift::DrawMarkers2d(
    const std::vector<std::array<float, 3> >& points, const short col,
    const double size) {
  if (points.empty()) return;
  TGraph gr;
  gr.SetMarkerColor(col);
  gr.SetMarkerSize(size);
  gr.SetMarkerStyle(20);
  std::vector<float> xgr;
  std::vector<float> ygr;
  for (const auto& p : points) {
    if (!InBox(p)) continue;
    float xp = 0., yp = 0.;
    ToPlane(p[0], p[1], p[2], xp, yp);
    xgr.push_back(xp);
    ygr.push_back(yp); 
  }
  if (!xgr.empty()) {
    gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Psame");
  }
}

void ViewDrift::Plot3d(const bool axis, const bool ogl, 
                       const bool snapshot) {
  auto pad = GetCanvas();
  pad->cd();
  pad->SetTitle("Drift lines");
  if (!pad->GetView()) {
    if (!SetPlotLimits3d()) {
      std::cerr << m_className << "::Plot3d:\n"
                << "     Could not determine the plot limits.\n";
    }
    auto view = TView::CreateView(1, 0, 0);
    view->SetRange(m_xMinBox, m_yMinBox, m_zMinBox, 
                   m_xMaxBox, m_yMaxBox, m_zMaxBox);
    if (axis) view->ShowAxis();
    pad->SetView(view);
    if (ogl) pad->GetViewer3D("ogl");
  }

  if (snapshot) {
    std::vector<std::array<float, 3> > electrons;
    std::vector<std::array<float, 3> > holes;
    std::vector<std::array<float, 3> > ions;
    for (const auto& driftLine : m_driftLines) {
      if (driftLine.second == Particle::Electron) {
        electrons.push_back(driftLine.first.back());
      } else if (driftLine.second == Particle::Hole) {
        holes.push_back(driftLine.first.back());
      } else {
        ions.push_back(driftLine.first.back());
      }
    }
    DrawMarkers3d(electrons, m_colElectron, m_markerSizeCollision);
    DrawMarkers3d(holes, m_colHole, m_markerSizeCollision);
    DrawMarkers3d(ions, m_colIon, m_markerSizeCollision);
  } else {
    for (const auto& driftLine : m_driftLines) {
      std::vector<float> points;
      for (const auto& p : driftLine.first) {
        points.push_back(p[0]);
        points.push_back(p[1]);
        points.push_back(p[2]);
      }
      const int nP = driftLine.first.size();
      TPolyLine3D pl(nP, points.data());
      if (driftLine.second == Particle::Electron) {
        pl.SetLineColor(m_colElectron);
      } else if (driftLine.second == Particle::Hole) {
        pl.SetLineColor(m_colHole);
      } else {
        pl.SetLineColor(m_colIon);
      }
      pl.SetLineWidth(1);
      pl.DrawPolyLine(nP, points.data(), "same");
    }
  }
  for (const auto& track : m_tracks) {
    std::vector<float> points;
    for (const auto& p : track) {
      points.push_back(p[0]);
      points.push_back(p[1]);
      points.push_back(p[2]);
    }
    const int nP = track.size();
    TPolyLine3D pl(nP, points.data());
    pl.SetLineColor(m_colTrack);
    pl.SetLineWidth(1);
    pl.DrawPolyLine(nP, points.data(), "same");
    if (m_drawClusters) {
      DrawMarkers3d(track, m_colTrack, m_markerSizeCluster);
    }
  }

  if (!m_exc.empty()) {
    DrawMarkers3d(m_exc, m_colExcitation, m_markerSizeCollision);
  }
  if (!m_ion.empty()) {
    DrawMarkers3d(m_ion, m_colIonisation, m_markerSizeCollision);
  }
  if (!m_att.empty()) {
    DrawMarkers3d(m_att, m_colAttachment, m_markerSizeCollision);
  } 
  pad->Modified();
  pad->Update();

  if (axis) {
    TAxis3D* ax3d = TAxis3D::GetPadAxis();
    if (ax3d) {
      ax3d->SetLabelColor(kGray + 2);
      ax3d->SetAxisColor(kGray + 2);
      ax3d->SetXTitle("x");
      ax3d->SetYTitle("y");
      ax3d->SetZTitle("z");
    }
    pad->Update();
  }
}

void ViewDrift::DrawMarkers3d(
    const std::vector<std::array<float, 3> >& points, const short col,
    const double size) {

  const size_t nP = points.size();
  std::vector<float> xyz;
  for (size_t i = 0; i < nP; ++i) {
    xyz.push_back(points[i][0]);
    xyz.push_back(points[i][1]);
    xyz.push_back(points[i][2]);
  }
  TPolyMarker3D pm(nP, xyz.data());
  pm.SetMarkerColor(col);
  pm.SetMarkerSize(size);
  pm.DrawPolyMarker(nP, xyz.data(), 20, "same");
}

bool ViewDrift::SetPlotLimits2d() {

  if (m_userPlotLimits) return true;
  double xmin = 0., ymin = 0., xmax = 0., ymax = 0.;
  if (m_userBox) {
    if (PlotLimitsFromUserBox(xmin, ymin, xmax, ymax)) {
      m_xMinPlot = xmin;
      m_yMinPlot = ymin;
      m_xMaxPlot = xmax;
      m_yMaxPlot = ymax;
      return true;
    }
  }

  // Try to determine the limits from the drift lines themselves.
  std::array<double, 3> bbmin;
  std::array<double, 3> bbmax;
  bbmin.fill(std::numeric_limits<double>::max());
  bbmax.fill(-std::numeric_limits<double>::max());
  for (const auto& driftLine : m_driftLines) {
    for (const auto& p : driftLine.first) {
      for (unsigned int i = 0; i < 3; ++i) {
        bbmin[i] = std::min(bbmin[i], double(p[i])); 
        bbmax[i] = std::max(bbmax[i], double(p[i]));
      }
    }
  }
  for (const auto& track : m_tracks) {
    for (const auto& p : track) {
      for (unsigned int i = 0; i < 3; ++i) {
        bbmin[i] = std::min(bbmin[i], double(p[i])); 
        bbmax[i] = std::max(bbmax[i], double(p[i]));
      }
    }
  }
  return PlotLimits(bbmin, bbmax, 
                    m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);
}

bool ViewDrift::SetPlotLimits3d() {

  if (m_userBox) return true;
  if (m_driftLines.empty() && m_tracks.empty()) return false;
  // Try to determine the limits from the drift lines themselves.
  std::array<double, 3> bbmin;
  std::array<double, 3> bbmax;
  bbmin.fill(std::numeric_limits<double>::max());
  bbmax.fill(-std::numeric_limits<double>::max());
  for (const auto& driftLine : m_driftLines) {
    for (const auto& p : driftLine.first) {
      for (unsigned int i = 0; i < 3; ++i) {
        bbmin[i] = std::min(bbmin[i], double(p[i])); 
        bbmax[i] = std::max(bbmax[i], double(p[i]));
      }
    }
  }
  for (const auto& track : m_tracks) {
    for (const auto& p : track) {
      for (unsigned int i = 0; i < 3; ++i) {
        bbmin[i] = std::min(bbmin[i], double(p[i])); 
        bbmax[i] = std::max(bbmax[i], double(p[i]));
      }
    }
  }
  double r = 0.;
  for (size_t i = 0; i < 3; ++i) r += fabs(bbmax[i] - bbmin[i]);
  m_xMinBox = bbmin[0];
  m_yMinBox = bbmin[1];
  m_zMinBox = bbmin[2];
  m_xMaxBox = bbmax[0];
  m_yMaxBox = bbmax[1];
  m_zMaxBox = bbmax[2];
  if (fabs(m_xMaxBox - m_xMinBox) < 0.1 * r) {
    m_xMinBox -= 0.1 * r;
    m_xMaxBox += 0.1 * r;
  }
  if (fabs(m_yMaxBox - m_yMinBox) < 0.1 * r) {
    m_yMinBox -= 0.1 * r;
    m_yMaxBox += 0.1 * r;
  }
  if (fabs(m_zMaxBox - m_zMinBox) < 0.1 * r) {
    m_zMinBox -= 0.1 * r;
    m_zMaxBox += 0.1 * r;
  }
  return true;
}

}
