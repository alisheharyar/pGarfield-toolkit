#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <TH1F.h>
#include <TPolyLine.h>
#include <TGraph.h>
#include <TGeoSphere.h>
#include <TPolyLine3D.h>

#include "Garfield/ComponentCST.hh"
#include "Garfield/ComponentFieldMap.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"
#include "Garfield/TGeoTet.hh"
#include "Garfield/ViewFEMesh.hh"

namespace Garfield {

ViewFEMesh::ViewFEMesh() : ViewBase("ViewFEMesh") {}

ViewFEMesh::~ViewFEMesh() {
  Reset();
}

void ViewFEMesh::Reset() {
  for (auto it = m_volumes.begin(), end = m_volumes.end(); it != end; ++it) {
    if (*it) {
      TGeoShape* shape = (*it)->GetShape();
      if (shape) delete shape;
      delete *it;
    }
  }
  m_volumes.clear();
  for (auto it = m_media.begin(), end = m_media.end(); it != end; ++it) {
    if (*it) {
      TGeoMaterial* material = (*it)->GetMaterial();
      if (material) delete material;
      delete *it;
    }
  }
  m_media.clear();

  m_geoManager.reset(nullptr);
}
 
void ViewFEMesh::SetComponent(ComponentFieldMap* cmp) {
  if (!cmp) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }

  m_cmp = cmp;
}

// The plotting functionality here is ported from Garfield
//  with some inclusion of code from ViewCell.cc
bool ViewFEMesh::Plot(const bool twod) {
  if (!m_cmp) {
    std::cerr << m_className << "::Plot: Component is not defined.\n";
    return false;
  }

  double pmin = 0., pmax = 0.;
  if (!m_cmp->GetVoltageRange(pmin, pmax)) {
    std::cerr << m_className << "::Plot: Component is not ready.\n";
    return false;
  }

  auto pad = GetCanvas();
  pad->cd();

  if (!twod) {
    if (!m_cmp->m_is3d) {
      std::cerr << m_className << "::Plot:\n"
                << "    Cannot plot 2D mesh elements in 3D.\n";
      return false;
    }
    DrawElements3d();
    DrawDriftLines3d();
    return true;
  }

  if (!GetPlotLimits()) return false;

  if (!RangeSet(pad)) {
    SetRange(pad, m_xMinPlot, m_yMinPlot, m_xMaxPlot, m_yMaxPlot);
  }

  if (m_drawAxes) {
    if (!m_xaxis && !m_yaxis) {
      // Draw default axes.
      auto frame = pad->DrawFrame(m_xMinPlot, m_yMinPlot,
                                  m_xMaxPlot, m_yMaxPlot);
      if (m_xaxisTitle.empty()) {
        frame->GetXaxis()->SetTitle(LabelX().c_str());
      } else { 
        frame->GetXaxis()->SetTitle(m_xaxisTitle.c_str());
      }
      if (m_yaxisTitle.empty()) {
        frame->GetYaxis()->SetTitle(LabelY().c_str());
      } else {
        frame->GetYaxis()->SetTitle(m_yaxisTitle.c_str());
      }
    } else {
      // Draw custom axes.
      if (m_xaxis) m_xaxis->Draw();
      if (m_yaxis) m_yaxis->Draw();
    }
  }

  // Plot the mesh elements.
  auto cst = dynamic_cast<ComponentCST*>(m_cmp);
  if (cst) {
    std::cout << m_className << "::Plot: CST component. Calling DrawCST.\n";
    DrawCST(cst);
  } else {
    DrawElements2d();
  }
 
  DrawDriftLines2d();
 
  if (m_drawViewRegion && !m_viewRegionX.empty()) {
    TPolyLine poly;
    poly.SetLineColor(kSpring + 4);
    poly.SetLineWidth(3);
    std::vector<double> xv = m_viewRegionX;
    std::vector<double> yv = m_viewRegionY;
    // Close the polygon.
    xv.push_back(m_viewRegionX[0]);
    yv.push_back(m_viewRegionY[0]);
    poly.DrawPolyLine(xv.size(), xv.data(), yv.data(), "same");
  }
  gPad->Update();
  // Draw axes again so they are on top
  gPad->RedrawAxis("g");
  return true;
}

bool ViewFEMesh::GetPlotLimits() {

  if (m_userPlotLimits) {
    std::vector<double> xp = {m_xMinPlot, m_xMinPlot, m_xMaxPlot, m_xMaxPlot};
    std::vector<double> yp = {m_yMinPlot, m_yMaxPlot, m_yMaxPlot, m_yMinPlot};
    std::vector<double> xg(4, 0.);
    std::vector<double> yg(4, 0.);
    std::vector<double> zg(4, 0.);
    for (size_t i = 0; i < 4; ++i) {
      xg[i] = m_proj[0][0] * xp[i] + m_proj[1][0] * yp[i] + m_proj[2][0];
      yg[i] = m_proj[0][1] * xp[i] + m_proj[1][1] * yp[i] + m_proj[2][1];
      zg[i] = m_proj[0][2] * xp[i] + m_proj[1][2] * yp[i] + m_proj[2][2];
    }
    m_xMinBox = *std::min_element(xg.begin(), xg.end());
    m_xMaxBox = *std::max_element(xg.begin(), xg.end());
    m_yMinBox = *std::min_element(yg.begin(), yg.end());
    m_yMaxBox = *std::max_element(yg.begin(), yg.end());
    m_zMinBox = *std::min_element(zg.begin(), zg.end());
    m_zMaxBox = *std::max_element(zg.begin(), zg.end());
    m_viewRegionX = xp;
    m_viewRegionY = yp;
    return true;
  }
  if (!m_userBox) {
    // If not set by the user, get the bounding box of the component.
    if (!m_cmp) return false;
    if (!m_cmp->GetBoundingBox(m_xMinBox, m_yMinBox, m_zMinBox,
                               m_xMaxBox, m_yMaxBox, m_zMaxBox)) {
      std::cerr << m_className << "::GetPlotLimits:\n"
                << "    Bounding box of the component is not defined.\n"
                << "    Please set the limits explicitly (SetArea).\n";
      return false;
    }
    if (std::isinf(m_xMinBox) || std::isinf(m_xMaxBox) || 
        std::isinf(m_yMinBox) || std::isinf(m_yMaxBox) ||
        std::isinf(m_zMinBox) || std::isinf(m_zMaxBox)) {
      double x0 = 0., y0 = 0., z0 = 0.;
      double x1 = 0., y1 = 0., z1 = 0.;
      if (!m_cmp->GetElementaryCell(x0, y0, z0, x1, y1, z1)) {
        std::cerr << m_className << "::GetPlotLimits:\n"
                  << "    Cell boundaries are not defined.\n"
                  << "    Please set the limits explicitly (SetArea).\n";
      }
      if (std::isinf(m_xMinBox) || std::isinf(m_xMaxBox)) {
        m_xMinBox = x0;
        m_xMaxBox = x1;
      }
      if (std::isinf(m_yMinBox) || std::isinf(m_yMaxBox)) {
        m_yMinBox = y0;
        m_yMaxBox = y1;
      }
      if (std::isinf(m_zMinBox) || std::isinf(m_zMaxBox)) {
        m_zMinBox = z0;
        m_zMaxBox = z1;
      }
    }
  }
  // Determine the intersection of the viewing plane and the box.
  double xmin = 0., xmax = 0., ymin = 0., ymax = 0.;
  IntersectPlaneArea(xmin, ymin, xmax, ymax);
  if (m_viewRegionX.empty()) {
    std::cerr << m_className << "::GetPlotLimits: Empty view.\n"
              << "    Make sure the viewing plane (SetPlane)\n"
              << "    intersects with the bounding box.\n";
    return false;
  }
  m_xMinPlot = xmin;
  m_xMaxPlot = xmax;
  m_yMinPlot = ymin;
  m_yMaxPlot = ymax;
  return true;
}

void ViewFEMesh::SetPlane(const double fx, const double fy, const double fz,
                          const double x0, const double y0, const double z0) {
  if (fy * fy + fz * fz > 0) {
    SetPlane(fx, fy, fz, x0, y0, z0, 1, 0, 0);
  } else {
    SetPlane(fx, fy, fz, x0, y0, z0, 0, 1, 0);
  }
}

void ViewFEMesh::SetPlane(const double fx, const double fy, const double fz,
                          const double x0, const double y0, const double z0,
                          const double hx, const double hy, const double hz) {
  ViewBase::SetPlane(fx, fy, fz, x0, y0, z0, hx, hy, hz);
}

// Set the x-axis.
void ViewFEMesh::SetXaxis(TGaxis* ax) { m_xaxis = ax; }

// Set the y-axis.
void ViewFEMesh::SetYaxis(TGaxis* ay) { m_yaxis = ay; }

// Create default axes
void ViewFEMesh::CreateDefaultAxes() {
  // Create a new x and y axis.
  if (!GetPlotLimits()) {
    std::cerr << m_className << "::CreateDefaultAxes:\n"
              << "    Cannot determine the axis limits.\n";
    return;
  } 
  const double dx = std::abs(m_xMaxPlot - m_xMinPlot) * 0.1;
  const double dy = std::abs(m_yMaxPlot - m_yMinPlot) * 0.1;
  const double x0 = m_xMinPlot + dx;
  const double y0 = m_yMinPlot + dy;
  const double x1 = m_xMaxPlot - dx;
  const double y1 = m_yMaxPlot - dy;
  m_xaxis = new TGaxis(x0, y0, x1, y0, x0, x1, 2405, "x");
  m_yaxis = new TGaxis(x0, y0, x0, y1, y0, y1, 2405, "y");

  // Label sizes
  m_xaxis->SetLabelSize(0.025);
  m_yaxis->SetLabelSize(0.025);

  // Titles
  m_xaxis->SetTitleSize(0.03);
  m_xaxis->SetTitle(LabelX().c_str());
  m_yaxis->SetTitleSize(0.03);
  m_yaxis->SetTitle(LabelY().c_str());
}

// Use ROOT plotting functions to draw the mesh elements on the canvas.
// General methodology ported from Garfield
void ViewFEMesh::DrawElements2d() {
  // Get the map boundaries from the component.
  double mapxmax = m_cmp->m_mapmax[0];
  double mapxmin = m_cmp->m_mapmin[0];
  double mapymax = m_cmp->m_mapmax[1];
  double mapymin = m_cmp->m_mapmin[1];
  double mapzmax = m_cmp->m_mapmax[2];
  double mapzmin = m_cmp->m_mapmin[2];

  // Get the periodicities.
  double sx = mapxmax - mapxmin;
  double sy = mapymax - mapymin;
  double sz = mapzmax - mapzmin;
  const bool perX = m_cmp->m_periodic[0] || m_cmp->m_mirrorPeriodic[0];
  const bool perY = m_cmp->m_periodic[1] || m_cmp->m_mirrorPeriodic[1];
  const bool perZ = m_cmp->m_periodic[2] || m_cmp->m_mirrorPeriodic[2];

  // Get the plane information.
  const double fx = m_plane[0];
  const double fy = m_plane[1];
  const double fz = m_plane[2];
  const double dist = m_plane[3];

  // Construct single-column matrix for use as coordinate vector.
  TMatrixD xMat(3, 1);

  // Determine the number of periods present in the cell.
  const int nMinX = perX ? int(m_xMinBox / sx) - 1 : 0;
  const int nMaxX = perX ? int(m_xMaxBox / sx) + 1 : 0;
  const int nMinY = perY ? int(m_yMinBox / sy) - 1 : 0;
  const int nMaxY = perY ? int(m_yMaxBox / sy) + 1 : 0;
  const int nMinZ = perZ ? int(m_zMinBox / sz) - 1 : 0;
  const int nMaxZ = perZ ? int(m_zMaxBox / sz) + 1 : 0;

  bool cst = false;
  if (dynamic_cast<ComponentCST*>(m_cmp)) cst = true;

  // Loop over all elements.
  const auto nElements = m_cmp->GetNumberOfElements();
  for (size_t i = 0; i < nElements; ++i) {
    size_t mat = 0;
    bool driftmedium = false;
    std::vector<size_t> nodes;
    if (!m_cmp->GetElement(i, mat, driftmedium, nodes)) continue;
    // Do not plot the drift medium.
    if (driftmedium && !m_plotMeshBorders) continue;
    // Do not create polygons for disabled materials.
    if (m_disabledMaterial.count(mat) > 0 && m_disabledMaterial[mat]) {
      continue;
    }
    TGraph gr;
    const short col = m_colorMap.count(mat) != 0 ? m_colorMap[mat] : 1;
    gr.SetLineColor(col);
    if (m_colorMap_fill.count(mat) != 0) {
      gr.SetFillColor(m_colorMap_fill[mat]);
    } else {
      gr.SetFillColor(col);
    }
    gr.SetLineWidth(3);
    std::string opt = "";
    if (m_plotMeshBorders || !m_fillMesh) opt += "l";
    if (m_fillMesh) opt += "f";
    opt += "same";
    // Get the vertex coordinates in the basic cell.
    std::vector<double> vx0;
    std::vector<double> vy0;
    std::vector<double> vz0;
    const size_t nNodes = nodes.size();
    for (size_t j = 0; j < nNodes; ++j) {
      double xn = 0., yn = 0., zn = 0.;
      if (!m_cmp->GetNode(nodes[j], xn, yn, zn)) continue;
      vx0.push_back(xn);
      vy0.push_back(yn);
      vz0.push_back(zn);
    }
    if (vx0.size() != nNodes) {
      std::cerr << m_className << "::DrawElements2d:\n"
                << "    Error retrieving nodes of element " << i << ".\n";
      continue;
    }
    // Coordinates of vertices
    std::vector<double> vx(nNodes, 0.);
    std::vector<double> vy(nNodes, 0.);
    std::vector<double> vz(nNodes, 0.);
    // Loop over the periodicities in x.
    for (int nx = nMinX; nx <= nMaxX; nx++) {
      const double dx = sx * nx;
      // Determine the x-coordinates of the vertices.
      if (m_cmp->m_mirrorPeriodic[0] && nx != 2 * (nx / 2)) {
        for (size_t j = 0; j < nNodes; ++j) {
          vx[j] = mapxmin + (mapxmax - vx0[j]) + dx;
        }
      } else {
        for (size_t j = 0; j < nNodes; ++j) {
          vx[j] = vx0[j] + dx;
        }
      }

      // Loop over the periodicities in y.
      for (int ny = nMinY; ny <= nMaxY; ny++) {
        const double dy = sy * ny;
        // Determine the y-coordinates of the vertices.
        if (m_cmp->m_mirrorPeriodic[1] && ny != 2 * (ny / 2)) {
          for (size_t j = 0; j < nNodes; ++j) {
            vy[j] = mapymin + (mapymax - vy0[j]) + dy;
          }
        } else {
          for (size_t j = 0; j < nNodes; ++j) {
            vy[j] = vy0[j] + dy;
          }
        }

        // Loop over the periodicities in z.
        for (int nz = nMinZ; nz <= nMaxZ; nz++) {
          const double dz = sz * nz;
          // Determine the z-coordinates of the vertices.
          if (m_cmp->m_mirrorPeriodic[2] && nz != 2 * (nz / 2)) {
            for (size_t j = 0; j < nNodes; ++j) {
              vz[j] = mapzmin + (mapzmax - vz0[j]) + dz;
            }
          } else {
            for (size_t j = 0; j < nNodes; ++j) {
              vz[j] = vz0[j] + dz;
            }
          }

          // Store the x and y coordinates of the relevant mesh vertices.
          std::vector<double> vX;
          std::vector<double> vY;

          // Value used to determine whether a vertex is in the plane.
          const double pcf = std::max(
              {std::abs(vx[0]), std::abs(vy[0]), std::abs(vz[0]), 
               std::abs(fx), std::abs(fy), std::abs(fz), std::abs(dist)});
          const double tol = 1.e-4 * pcf;
          // First isolate the vertices that are in the viewing plane.
          std::vector<bool> in(nNodes, false);
          int cnt = 0;
          for (size_t j = 0; j < nNodes; ++j) {
            const double d = fx * vx[j] + fy * vy[j] + fz * vz[j] - dist;
            if (std::abs(d) < tol) {
              // Point is in the plane.
              in[j] = true;
              // Calculate the planar coordinates.
              double xp = 0., yp = 0.;
              ToPlane(vx[j], vy[j], vz[j], xp, yp);
              vX.push_back(xp);
              vY.push_back(yp);
            } else {
              if (d > 0.) {
                cnt += 1;
              } else { 
                cnt -= 1;
              }
            }
          }
          // Stop if all points are on the same side of the plane.
          if (std::abs(cnt) == (int)nNodes) continue;
          // Cut the sides that are not in the plane.
          if (cst) {
            const std::array<std::array<unsigned int, 3>, 8> neighbours = {{
              {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
              {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}
            }};
            for (size_t j = 0; j < nNodes; ++j) {
              for (unsigned int k : neighbours[j]) {
                if (in[j] || in[k]) continue;
                if (PlaneCut(vx[j], vy[j], vz[j], 
                             vx[k], vy[k], vz[k], xMat)) {
                  vX.push_back(xMat(0, 0));
                  vY.push_back(xMat(1, 0));
                }
              }
            } 
          } else { 
            // Tetrahedron.
            for (size_t j = 0; j < nNodes; ++j) {
              for (size_t k = j + 1; k < nNodes; ++k) {
                if (in[j] || in[k]) continue;
                if (PlaneCut(vx[j], vy[j], vz[j], 
                             vx[k], vy[k], vz[k], xMat)) {
                  vX.push_back(xMat(0, 0));
                  vY.push_back(xMat(1, 0));
                }
              }
            }
          }
          if (vX.size() < 3) continue;

          // Eliminate crossings of the polygon lines
          // (known as "butterflies" in Garfield).
          RemoveCrossings(vX, vY);

          // Create vectors to store the clipped polygon.
          std::vector<double> cX;
          std::vector<double> cY;

          // Clip the polygon to the view area.
          ClipToView(vX, vY, cX, cY);

          // If we have more than 2 vertices, add the polygon.
          if (cX.size() <= 2) continue;
          // Again eliminate crossings of the polygon lines.
          RemoveCrossings(cX, cY);
  
          // Draw the polygon.
          std::vector<float> xgr(cX.begin(), cX.end());
          std::vector<float> ygr(cY.begin(), cY.end());
          gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), opt.c_str());
        }  // end z-periodicity loop
      }    // end y-periodicity loop
    }      // end x-periodicity loop
  }        // end loop over elements

}

void ViewFEMesh::DrawElements3d() {

  // Get the map boundaries from the component.
  double mapxmax = m_cmp->m_mapmax[0];
  double mapxmin = m_cmp->m_mapmin[0];
  double mapymax = m_cmp->m_mapmax[1];
  double mapymin = m_cmp->m_mapmin[1];
  double mapzmax = m_cmp->m_mapmax[2];
  double mapzmin = m_cmp->m_mapmin[2];

  // Get the periodicities.
  double sx = mapxmax - mapxmin;
  double sy = mapymax - mapymin;
  double sz = mapzmax - mapzmin;
  const bool perX = m_cmp->m_periodic[0] || m_cmp->m_mirrorPeriodic[0];
  const bool perY = m_cmp->m_periodic[1] || m_cmp->m_mirrorPeriodic[1];
  const bool perZ = m_cmp->m_periodic[2] || m_cmp->m_mirrorPeriodic[2];

  // Set the plot limits.
  if (m_userBox) {
    if (std::isinf(m_xMinBox)) m_xMinBox = mapxmin;
    if (std::isinf(m_yMinBox)) m_yMinBox = mapymin;
    if (std::isinf(m_zMinBox)) m_zMinBox = mapzmin;
    if (std::isinf(m_xMaxBox)) m_xMaxBox = mapxmax;
    if (std::isinf(m_yMaxBox)) m_yMaxBox = mapymax;
    if (std::isinf(m_zMaxBox)) m_zMaxBox = mapzmax;
  } else {
    m_xMinBox = mapxmin;
    m_yMinBox = mapymin;
    m_zMinBox = mapzmin;
    m_xMaxBox = mapxmax;
    m_yMaxBox = mapymax;
    m_zMaxBox = mapzmax;
  } 
  // Determine the number of periods present in the cell.
  const int nMinX = perX ? int(m_xMinBox / sx) - 1 : 0;
  const int nMaxX = perX ? int(m_xMaxBox / sx) + 1 : 0;
  const int nMinY = perY ? int(m_yMinBox / sy) - 1 : 0;
  const int nMaxY = perY ? int(m_yMaxBox / sy) + 1 : 0;
  const int nMinZ = perZ ? int(m_zMinBox / sz) - 1 : 0;
  const int nMaxZ = perZ ? int(m_zMaxBox / sz) + 1 : 0;

  gGeoManager = nullptr;
  m_geoManager.reset(new TGeoManager("ViewFEMeshGeoManager", ""));
  TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0., 0., 0.);
  TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
  m_media.push_back(medVacuum);
  // Use silicon as "default" material.
  TGeoMaterial* matDefault = new TGeoMaterial("Default", 28.085, 14., 2.329);
  TGeoMedium* medDefault = new TGeoMedium("Default", 1, matDefault);
  const double hxw = std::max(std::abs(m_xMaxBox), std::abs(m_xMinBox));
  const double hyw = std::max(std::abs(m_yMaxBox), std::abs(m_yMinBox));
  const double hzw = std::max(std::abs(m_zMaxBox), std::abs(m_zMinBox));
  TGeoVolume* top = m_geoManager->MakeBox("Top", medVacuum, hxw, hyw, hzw);
  m_geoManager->SetTopVolume(top);
  m_volumes.push_back(top);

  // Loop over all elements.
  const auto nElements = m_cmp->GetNumberOfElements();
  for (size_t i = 0; i < nElements; ++i) {
    size_t mat = 0;
    bool driftmedium = false;
    std::vector<size_t> nodes;
    if (!m_cmp->GetElement(i, mat, driftmedium, nodes)) continue;
    // Do not plot the drift medium.
    if (driftmedium && !m_plotMeshBorders) continue;
    // Do not create polygons for disabled materials.
    if (m_disabledMaterial.count(mat) > 0 && m_disabledMaterial[mat]) {
      continue;
    }
    const short col = m_colorMap.count(mat) != 0 ? m_colorMap[mat] : 1;

    // Get the vertex coordinates in the basic cell.
    const size_t nNodes = nodes.size();
    if (nNodes != 4) continue;
    std::vector<std::array<double, 3> > v0;
    for (size_t j = 0; j < nNodes; ++j) {
      double xn = 0., yn = 0., zn = 0.;
      if (!m_cmp->GetNode(nodes[j], xn, yn, zn)) break;
      v0.push_back({xn, yn, zn});
    }
    if (v0.size() != nNodes) {
      std::cerr << m_className << "::DrawElements3d:\n"
                << "    Error retrieving nodes of element " << i << ".\n";
      continue;
    }
    // Coordinates of vertices
    std::array<std::array<double, 3>, 4> v;
    // Loop over the periodicities in x.
    for (int nx = nMinX; nx <= nMaxX; nx++) {
      const double dx = sx * nx;
      // Determine the x-coordinates of the vertices.
      if (m_cmp->m_mirrorPeriodic[0] && nx != 2 * (nx / 2)) {
        for (size_t j = 0; j < nNodes; ++j) {
          v[j][0] = mapxmin + (mapxmax - v0[j][0]) + dx;
        }
      } else {
        for (size_t j = 0; j < nNodes; ++j) {
          v[j][0] = v0[j][0] + dx;
        }
      }

      // Loop over the periodicities in y.
      for (int ny = nMinY; ny <= nMaxY; ny++) {
        const double dy = sy * ny;
        // Determine the y-coordinates of the vertices.
        if (m_cmp->m_mirrorPeriodic[1] && ny != 2 * (ny / 2)) {
          for (size_t j = 0; j < nNodes; ++j) {
            v[j][1] = mapymin + (mapymax - v0[j][1]) + dy;
          }
        } else {
          for (size_t j = 0; j < nNodes; ++j) {
            v[j][1] = v0[j][1] + dy;
          }
        }

        // Loop over the periodicities in z.
        for (int nz = nMinZ; nz <= nMaxZ; nz++) {
          const double dz = sz * nz;
          // Determine the z-coordinates of the vertices.
          if (m_cmp->m_mirrorPeriodic[2] && nz != 2 * (nz / 2)) {
            for (size_t j = 0; j < nNodes; ++j) {
              v[j][2] = mapzmin + (mapzmax - v0[j][2]) + dz;
            }
          } else {
            for (size_t j = 0; j < nNodes; ++j) {
              v[j][2] = v0[j][2] + dz;
            }
          }
          auto tet = new TGeoTet("Tet", v);
          std::string vname = "Tet" + std::to_string(m_volumes.size());
          TGeoVolume* vol = new TGeoVolume(vname.c_str(), tet, medDefault); 
          vol->SetLineColor(col);
          vol->SetTransparency(70.);
          m_volumes.push_back(vol);
          top->AddNodeOverlap(vol, 1, gGeoIdentity);
        }
      }
    }
  }
  m_geoManager->CloseGeometry();
  top->SetTransparency(100.);
  m_geoManager->SetTopVisible();
  m_geoManager->SetMaxVisNodes(m_volumes.size() + 1);
  m_geoManager->GetTopNode()->Draw("e");
}

void ViewFEMesh::DrawDriftLines2d() {
 
  if (!m_viewDrift) return;
  // Plot a 2D projection of the drift line.
  for (const auto& driftLine : m_viewDrift->m_driftLines) {
    TGraph gr;
    if (driftLine.second == Particle::Electron) {
      gr.SetLineColor(m_viewDrift->m_colElectron);
    } else if (driftLine.second == Particle::Hole) {
      gr.SetLineColor(m_viewDrift->m_colHole);
    } else {
      gr.SetLineColor(m_viewDrift->m_colIon);
    }
    std::vector<float> xgr;
    std::vector<float> ygr;
    // Loop over the points.
    for (const auto& point : driftLine.first) {
      // Project this point onto the plane.
      float xp = 0., yp = 0.;
      ToPlane(point[0], point[1], point[2], xp, yp);
      // Add this point if it is within the view.
      if (InView(xp, yp)) {
        xgr.push_back(xp);
        ygr.push_back(yp);
      }
    }
    if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "lsame");
    }
  }

}

void ViewFEMesh::DrawDriftLines3d() {

  if (!m_viewDrift) return;
  for (const auto& driftLine : m_viewDrift->m_driftLines) {
    std::vector<float> points;
    for (const auto& p : driftLine.first) {
      points.push_back(p[0]);
      points.push_back(p[1]);
      points.push_back(p[2]);
    }
    const int nP = driftLine.first.size();
    TPolyLine3D pl(nP, points.data());
    if (driftLine.second == Particle::Electron) {
      pl.SetLineColor(m_viewDrift->m_colElectron);
    } else if (driftLine.second == Particle::Hole) {
      pl.SetLineColor(m_viewDrift->m_colHole);
    } else {
      pl.SetLineColor(m_viewDrift->m_colIon);
    }
    pl.SetLineWidth(1);
    pl.DrawPolyLine(nP, points.data(), "same");
  }
}

void ViewFEMesh::DrawCST(ComponentCST* cst) {
  /*The method is based on ViewFEMesh::Draw, thus the first part is copied from
   * there.
   * At the moment only x-y, x-z, and y-z are available due to the simple
   * implementation.
   * The advantage of this method is that there is no element loop and thus it
   * is much
   * faster.
   */

  // Helper struct.
  struct PolygonInfo {
    double p1[2];
    double p2[2];
    double p3[2];
    double p4[2];
    int element;
  };

  // Get the map boundaries from the component
  double mapxmax = m_cmp->m_mapmax[0];
  double mapxmin = m_cmp->m_mapmin[0];
  double mapymax = m_cmp->m_mapmax[1];
  double mapymin = m_cmp->m_mapmin[1];
  double mapzmax = m_cmp->m_mapmax[2];
  double mapzmin = m_cmp->m_mapmin[2];

  // Get the periodicities.
  double sx = mapxmax - mapxmin;
  double sy = mapymax - mapymin;
  double sz = mapzmax - mapzmin;
  const bool perX = m_cmp->m_periodic[0] || m_cmp->m_mirrorPeriodic[0];
  const bool perY = m_cmp->m_periodic[1] || m_cmp->m_mirrorPeriodic[1];
  const bool perZ = m_cmp->m_periodic[2] || m_cmp->m_mirrorPeriodic[2];

  // Determine the number of periods present in the cell.
  const int nMinX = perX ? int(m_xMinBox / sx) - 1 : 0;
  const int nMaxX = perX ? int(m_xMaxBox / sx) + 1 : 0;
  const int nMinY = perY ? int(m_yMinBox / sy) - 1 : 0;
  const int nMaxY = perY ? int(m_yMaxBox / sy) + 1 : 0;
  const int nMinZ = perZ ? int(m_zMinBox / sz) - 1 : 0;
  const int nMaxZ = perZ ? int(m_zMaxBox / sz) + 1 : 0;

  std::vector<PolygonInfo> elements;
  int nMinU = 0, nMaxU = 0, nMinV = 0, nMaxV = 0;
  double mapumin = 0., mapumax = 0., mapvmin = 0., mapvmax = 0.;
  double su = 0., sv = 0.;
  bool mirroru = false, mirrorv = false;
  double uMin, vMin, uMax, vMax;
  unsigned int n_x, n_y, n_z;
  cst->GetNumberOfMeshLines(n_x, n_y, n_z);

  const double fx = m_plane[0];
  const double fy = m_plane[1];
  const double fz = m_plane[2];
  if (fx == 0 && fy == 0 && fz == 1) {
    // xy view
    std::cout << m_className << "::DrawCST: Creating x-y mesh view.\n";
    // Calculate the z position.
    unsigned int i, j, z;
    const double z0 = m_plane[3] * fz;
    if (!cst->Coordinate2Index(0, 0, z0, i, j, z)) {
      std::cerr << "    Could not determine the z-index of the plane.\n";
      return;
    }
    std::cout << "    The z-index of the plane is " << z << ".\n";
    nMinU = nMinX;
    nMaxU = nMaxX;
    nMinV = nMinY;
    nMaxV = nMaxY;
    uMin = m_xMinBox;
    uMax = m_xMaxBox;
    vMin = m_yMinBox;
    vMax = m_yMaxBox;

    mapumin = mapxmin;
    mapumax = mapxmax;
    mapvmin = mapymin;
    mapvmax = mapymax;
    su = sx;
    sv = sy;
    mirroru = perX;
    mirrorv = perY;
    for (unsigned int y = 0; y < (n_y - 1); y++) {
      for (unsigned int x = 0; x < (n_x - 1); x++) {
        auto elem = cst->Index2Element(x, y, z);
        double e_xmin, e_xmax, e_ymin, e_ymax, e_zmin, e_zmax;
        cst->GetElementBoundaries(elem, e_xmin, e_xmax, e_ymin, e_ymax,
                                        e_zmin, e_zmax);
        PolygonInfo tmp_info;
        tmp_info.element = elem;
        tmp_info.p1[0] = e_xmin;
        tmp_info.p2[0] = e_xmax;
        tmp_info.p3[0] = e_xmax;
        tmp_info.p4[0] = e_xmin;
        tmp_info.p1[1] = e_ymin;
        tmp_info.p2[1] = e_ymin;
        tmp_info.p3[1] = e_ymax;
        tmp_info.p4[1] = e_ymax;
        elements.push_back(std::move(tmp_info));
      }
    }
  } else if (fx == 0 && fy == -1 && fz == 0) {
    // xz-view
    std::cout << m_className << "::DrawCST: Creating x-z mesh view.\n";
    // Calculate the y position.
    unsigned int i = 0, j = 0, y = 0;
    const double y0 = m_plane[3] * fy;
    if (!cst->Coordinate2Index(0, y0, 0, i, y, j)) {
      std::cerr << "    Could not determine the y-index of the plane.\n";
      return;
    }
    std::cout << "    The y-index of the plane is " << y << ".\n";

    nMinU = nMinX;
    nMaxU = nMaxX;
    nMinV = nMinZ;
    nMaxV = nMaxZ;
    uMin = m_xMinBox;
    uMax = m_xMaxBox;
    vMin = m_zMinBox;
    vMax = m_zMaxBox;

    mapumin = mapxmin;
    mapumax = mapxmax;
    mapvmin = mapzmin;
    mapvmax = mapzmax;
    su = sx;
    sv = sz;
    mirroru = perX;
    mirrorv = perZ;
    for (unsigned int z = 0; z < (n_z - 1); z++) {
      for (unsigned int x = 0; x < (n_x - 1); x++) {
        auto elem = cst->Index2Element(x, y, z);
        double e_xmin, e_xmax, e_ymin, e_ymax, e_zmin, e_zmax;
        cst->GetElementBoundaries(elem, e_xmin, e_xmax, e_ymin, e_ymax,
                                        e_zmin, e_zmax);
        PolygonInfo tmp_info;
        tmp_info.element = elem;
        tmp_info.p1[0] = e_xmin;
        tmp_info.p2[0] = e_xmax;
        tmp_info.p3[0] = e_xmax;
        tmp_info.p4[0] = e_xmin;
        tmp_info.p1[1] = e_zmin;
        tmp_info.p2[1] = e_zmin;
        tmp_info.p3[1] = e_zmax;
        tmp_info.p4[1] = e_zmax;
        elements.push_back(std::move(tmp_info));
      }
    }
  } else if (fx == -1 && fy == 0 && fz == 0) {
    // yz-view
    std::cout << m_className << "::DrawCST: Creating z-y mesh view.\n";
    // Calculate the x position.
    unsigned int i, j, x;
    const double x0 = m_plane[3] * fx;
    if (!cst->Coordinate2Index(x0, 0, 0, x, i, j)) {
      std::cerr << "    Could not determine the x-index of the plane.\n";
      return;
    }
    std::cout << "    The x-index of the plane is " << x << ".\n";
    nMinU = nMinZ;
    nMaxU = nMaxZ;
    nMinV = nMinY;
    nMaxV = nMaxY;
    uMin = m_yMinBox;
    uMax = m_yMaxBox;
    vMin = m_zMinBox;
    vMax = m_zMaxBox;

    mapumin = mapzmin;
    mapumax = mapzmax;
    mapvmin = mapymin;
    mapvmax = mapymax;
    su = sz;
    sv = sy;
    mirroru = perZ;
    mirrorv = perY;
    for (unsigned int z = 0; z < (n_z - 1); z++) {
      for (unsigned int y = 0; y < (n_y - 1); y++) {
        auto elem = cst->Index2Element(x, y, z);
        double e_xmin, e_xmax, e_ymin, e_ymax, e_zmin, e_zmax;
        cst->GetElementBoundaries(elem, e_xmin, e_xmax, e_ymin, e_ymax,
                                        e_zmin, e_zmax);
        PolygonInfo tmp_info;
        tmp_info.element = elem;
        tmp_info.p1[0] = e_zmin;
        tmp_info.p2[0] = e_zmax;
        tmp_info.p3[0] = e_zmax;
        tmp_info.p4[0] = e_zmin;
        tmp_info.p1[1] = e_ymin;
        tmp_info.p2[1] = e_ymin;
        tmp_info.p3[1] = e_ymax;
        tmp_info.p4[1] = e_ymax;
        // Add the polygon to the mesh
        elements.push_back(std::move(tmp_info));
      }
    }
  } else {
    std::cerr << m_className << "::DrawCST:\n";
    std::cerr << "    The given plane name is not known.\n";
    std::cerr << "    Please choose one of the following: xy, xz, yz.\n";
    return;
  }

  std::cout << m_className << "::DrawCST:\n    " << elements.size()
            << " elements in the projection of the unit cell.\n";
  for (const auto& element : elements) {
    size_t mat = 0;
    bool driftmedium = false;
    std::vector<size_t> nodes;
    cst->GetElement(element.element, mat, driftmedium, nodes);
    // Do not plot the drift medium.
    if (driftmedium && !m_plotMeshBorders) continue;
    // Do not create polygons for disabled materials.
    if (m_disabledMaterial[mat]) continue;
    TGraph gr;
    const short col = m_colorMap.count(mat) > 0 ? m_colorMap[mat] : 1;
    gr.SetLineColor(col);
    if (m_colorMap_fill.count(mat) > 0) {
      gr.SetFillColor(m_colorMap_fill[mat]);
    } else {
      gr.SetFillColor(col);
    }
    if (m_plotMeshBorders)
      gr.SetLineWidth(3);
    else
      gr.SetLineWidth(1);

    for (int nu = nMinU; nu <= nMaxU; nu++) {
      for (int nv = nMinV; nv <= nMaxV; nv++) {
        // Add 4 points of the square
        float tmp_u[4], tmp_v[4];
        if (mirroru && nu != 2 * (nu / 2)) {
          // nu is odd
          tmp_u[0] = mapumin + (mapumax - element.p1[0]) + su * nu;
          tmp_u[1] = mapumin + (mapumax - element.p2[0]) + su * nu;
          tmp_u[2] = mapumin + (mapumax - element.p3[0]) + su * nu;
          tmp_u[3] = mapumin + (mapumax - element.p4[0]) + su * nu;
        } else {
          // nu is even
          tmp_u[0] = element.p1[0] + su * nu;
          tmp_u[1] = element.p2[0] + su * nu;
          tmp_u[2] = element.p3[0] + su * nu;
          tmp_u[3] = element.p4[0] + su * nu;
        }
        if (mirrorv && nv != 2 * (nv / 2)) {
          tmp_v[0] = mapvmin + (mapvmax - element.p1[1]) + sv * nv;
          tmp_v[1] = mapvmin + (mapvmax - element.p2[1]) + sv * nv;
          tmp_v[2] = mapvmin + (mapvmax - element.p3[1]) + sv * nv;
          tmp_v[3] = mapvmin + (mapvmax - element.p4[1]) + sv * nv;
        } else {
          tmp_v[0] = element.p1[1] + sv * nv;
          tmp_v[1] = element.p2[1] + sv * nv;
          tmp_v[2] = element.p3[1] + sv * nv;
          tmp_v[3] = element.p4[1] + sv * nv;
        }
        if (tmp_u[0] < uMin || tmp_u[1] > uMax || tmp_v[0] < vMin ||
            tmp_v[2] > vMax) {
          continue;
        }
        std::string opt = "";
        if (m_plotMeshBorders || !m_fillMesh) opt += "l";
        if (m_fillMesh) opt += "f";
        opt += "same";
        gr.DrawGraph(4, tmp_u, tmp_v, opt.c_str());
      }
    }
  }

}

// Removes duplicate points and line crossings by correctly ordering
//  the points in the provided vectors.
//
//  NOTE: This is a 2D version of the BUTFLD method in Garfield.  It
//   follows the same general algorithm.
//
// TODO: there is an algorithm which always sorts points correctly using cross
// product, see IntersectPlaneArea.
void ViewFEMesh::RemoveCrossings(std::vector<double>& x,
                                 std::vector<double>& y) {
  // Determine element dimensions
  double xmin = x[0], xmax = x[0];
  double ymin = y[0], ymax = y[0];
  for (int i = 1; i < (int)x.size(); i++) {
    if (x[i] < xmin) xmin = x[i];
    if (x[i] > xmax) xmax = x[i];
    if (y[i] < ymin) ymin = y[i];
    if (y[i] > ymax) ymax = y[i];
  }

  // First remove duplicate points
  double xtol = 1e-10 * std::abs(xmax - xmin);
  double ytol = 1e-10 * std::abs(ymax - ymin);
  for (int i = 0; i < (int)x.size(); i++) {
    for (int j = i + 1; j < (int)x.size(); j++) {
      if (std::abs(x[i] - x[j]) < xtol && std::abs(y[i] - y[j]) < ytol) {
        x.erase(x.begin() + j);
        y.erase(y.begin() + j);
        j--;
      }
    }
  }

  // No crossings with 3 points or less
  if (x.size() <= 3) return;

  // Save the polygon size so it is easily accessible
  int NN = x.size();

  // Keep track of the number of attempts
  int attempts = 0;

  // Exchange points until crossings are eliminated or we have attempted NN
  // times
  bool crossings = true;
  while (crossings && (attempts < NN)) {
    // Assume we are done after this attempt.
    crossings = false;

    for (int i = 1; i <= NN; i++) {
      for (int j = i + 2; j <= NN; j++) {
        // End the j-loop if we have surpassed N and wrapped around to i
        if ((j + 1) > NN && 1 + (j % NN) >= i) break;
        // Otherwise, detect crossings and attempt to eliminate them.
        // Determine if we have a crossing.
        double xc = 0., yc = 0.;
        if (!LinesCrossed(x[(i - 1) % NN], y[(i - 1) % NN], x[i % NN],
                          y[i % NN], x[(j - 1) % NN], y[(j - 1) % NN],
                          x[j % NN], y[j % NN], xc, yc)) {
          continue;
        }
        // Swap each point from i towards j with each corresponding point
        // from j towards i.
        for (int k = 1; k <= (j - i) / 2; k++) {
          double xs = x[(i + k - 1) % NN];
          double ys = y[(i + k - 1) % NN];
          x[(i + k - 1) % NN] = x[(j - k) % NN];
          y[(i + k - 1) % NN] = y[(j - k) % NN];
          x[(j - k) % NN] = xs;
          y[(j - k) % NN] = ys;

          // Force another attempt
          crossings = true;
        }
      }  // end loop over j
    }    // end loop over i

    // Increment the number of attempts
    attempts++;

  }  // end while(crossings)

  if (attempts > NN) {
    std::cerr << m_className << "::RemoveCrossings:\n    Warning: "
              << "Maximum attempts reached. Crossings not removed.\n";
  }
}

/// Return true if the specified point is in the view region.
bool ViewFEMesh::InView(const double x, const double y) const {
  // Test whether this vertex is inside the view.
  if (m_userPlotLimits) {
    return (x >= m_xMinPlot && x <= m_xMaxPlot && 
            y >= m_yMinPlot && y <= m_yMaxPlot);
  }
  bool edge = false;
  return IsInPolygon(x, y, m_viewRegionX, m_viewRegionY, edge);
}

//
// Determines whether the line connecting points (x1,y1) and (x2,y2)
// and the line connecting points (u1,v1) and (u2,v2) cross somewhere
// between the 4 points.  Sets the crossing point in (xc, yc).
//
// Ported from Garfield function CROSSD
//
bool ViewFEMesh::LinesCrossed(double x1, double y1, double x2, double y2,
                              double u1, double v1, double u2, double v2,
                              double& xc, double& yc) const {
  // Set the tolerances.
  double xtol = 1.0e-10 * std::max({std::abs(x1), std::abs(x2), std::abs(u1),
                                    std::abs(u2)});
  double ytol = 1.0e-10 * std::max({std::abs(y1), std::abs(y2), std::abs(v1),
                                    std::abs(v2)});
  if (xtol <= 0) xtol = 1.0e-10;
  if (ytol <= 0) ytol = 1.0e-10;

  // Compute the distances and determinant (dx,dy) x (du,dv).
  double dy = y2 - y1;
  double dv = v2 - v1;
  double dx = x1 - x2;
  double du = u1 - u2;
  double det = dy * du - dx * dv;

  // Check for crossing because one of the endpoints is on both lines.
  if (OnLine(x1, y1, x2, y2, u1, v1)) {
    xc = u1;
    yc = v1;
    return true;
  } else if (OnLine(x1, y1, x2, y2, u2, v2)) {
    xc = u2;
    yc = v2;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x1, y1)) {
    xc = x1;
    yc = y1;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x2, y2)) {
    xc = x2;
    yc = y2;
    return true;
  }
  // Check if the lines are parallel (zero determinant).
  if (std::abs(det) < xtol * ytol) return false;
  // No special case: compute point of intersection.

  // Solve crossing equations.
  xc = (du * (x1 * y2 - x2 * y1) - dx * (u1 * v2 - u2 * v1)) / det;
  yc = ((-1 * dv) * (x1 * y2 - x2 * y1) + dy * (u1 * v2 - u2 * v1)) / det;

  // Determine if this point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc))
    return true;

  // The lines do not cross if we have reached this point.
  return false;
}

// Determines whether the point (u,v) lies on the line connecting
// points (x1,y1) and (x2,y2).
//
// Ported from Garfield function ONLIND
//
bool ViewFEMesh::OnLine(double x1, double y1, double x2, double y2, double u,
                        double v) const {
  // Set the tolerances
  double xtol = 1.e-10 * std::max({std::abs(x1), std::abs(x2), std::abs(u)});
  double ytol = 1.e-10 * std::max({std::abs(y1), std::abs(y2), std::abs(v)});
  if (xtol <= 0) xtol = 1.0e-10;
  if (ytol <= 0) ytol = 1.0e-10;

  // To store the coordinates of the comparison point
  double xc = 0, yc = 0;

  // Check if (u,v) is the point (x1,y1) or (x2,y2)
  if ((std::abs(x1 - u) <= xtol && std::abs(y1 - v) <= ytol) ||
      (std::abs(x2 - u) <= xtol && std::abs(y2 - v) <= ytol)) {
    return true;
  }
  // Check if the line is actually a point
  if (std::abs(x1 - x2) <= xtol && std::abs(y1 - y2) <= ytol) {
    return false;
  }
  // Choose (x1,y1) as starting point if closer to (u,v)
  if (std::abs(u - x1) + std::abs(v - y1) <
      std::abs(u - x2) + std::abs(v - y2)) {
    // Compute the component of the line from (x1,y1) to (u,v)
    // along the line from (x1,y1) to (x2,y2)
    double dpar = ((u - x1) * (x2 - x1) + (v - y1) * (y2 - y1)) /
                  ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Determine the point on the line to which to compare (u,v)
    if (dpar < 0.0) {
      xc = x1;
      yc = y1;
    } else if (dpar > 1.0) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + dpar * (x2 - x1);
      yc = y1 + dpar * (y2 - y1);
    }
  } else { 
    // Choose (x2,y2) as starting point if closer to (u,v)
    // Compute the component of the line from (x2,y2) to (u,v)
    //  along the line from (x2,y2) to (x1,y1)
    double dpar = ((u - x2) * (x1 - x2) + (v - y2) * (y1 - y2)) /
                  ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Determine the point on the line to which to compare (u,v)
    if (dpar < 0.0) {
      xc = x2;
      yc = y2;
    } else if (dpar > 1.0) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + dpar * (x1 - x2);
      yc = y2 + dpar * (y1 - y2);
    }
  }

  // Compare the calculated point to (u,v)
  if (std::abs(u - xc) < xtol && std::abs(v - yc) < ytol) return true;

  return false;
}

// Ported from Garfield: determines the point of intersection, in planar
// coordinates, of a plane with the line connecting multiple points
// x1,y1,z1;x2,y2,z2: the world coordinates of the two points
// projMat;planeMat: the projection and plane matrices
// xMat: the resulting planar coordinates of the intersection point
bool ViewFEMesh::PlaneCut(double x1, double y1, double z1, double x2, double y2,
                          double z2, TMatrixD& xMat) {
  // Set up the matrix for cutting edges not in the plane
  TArrayD dataCut(9);
  TMatrixD cutMat(3, 3);
  dataCut[0] = m_proj[0][0];
  dataCut[1] = m_proj[1][0];
  dataCut[2] = x1 - x2;
  dataCut[3] = m_proj[0][1];
  dataCut[4] = m_proj[1][1];
  dataCut[5] = y1 - y2;
  dataCut[6] = m_proj[0][2];
  dataCut[7] = m_proj[1][2];
  dataCut[8] = z1 - z2;
  cutMat.SetMatrixArray(dataCut.GetArray());

  // Calculate the determinant of the cut matrix
  double cutDet =
      cutMat(0, 0) *
          (cutMat(1, 1) * cutMat(2, 2) - cutMat(1, 2) * cutMat(2, 1)) -
      cutMat(0, 1) *
          (cutMat(1, 0) * cutMat(2, 2) - cutMat(1, 2) * cutMat(2, 0)) +
      cutMat(0, 2) *
          (cutMat(1, 0) * cutMat(2, 1) - cutMat(1, 1) * cutMat(2, 0));

  // Do not proceed if the matrix is singular
  if (std::abs(cutDet) < 1e-20) return false;

  // Set up a coordinate vector (RHS of equation)
  TArrayD dataCoords(3);
  TMatrixD coordMat(3, 1);
  dataCoords[0] = x1 - m_proj[2][0];
  dataCoords[1] = y1 - m_proj[2][1];
  dataCoords[2] = z1 - m_proj[2][2];
  coordMat.SetMatrixArray(dataCoords.GetArray());

  // Invert the cut matrix and multiply to get the solution
  cutMat.SetTol(1e-20);
  cutMat.Invert();
  // Do not proceed if the matrix is singular
  if (!cutMat.IsValid()) return false;
  xMat = cutMat * coordMat;

  // Return success if the plane point is between the two vertices
  if (xMat(2, 0) < 0 || xMat(2, 0) > 1) return false;
  return true;
}

// Calculates view region and canvas dimensions based on projection plane
// and view area
bool ViewFEMesh::IntersectPlaneArea(double& xmin, double& ymin,
                                    double& xmax, double& ymax) {
  std::vector<TMatrixD> intersect_points;
  m_viewRegionX.clear();
  m_viewRegionY.clear();
  // Loop over box edges
  for (int i0 = 0; i0 < 2; ++i0) {
    for (int j0 = 0; j0 < 2; ++j0) {
      for (int k0 = 0; k0 < 2; ++k0) {
        for (int i1 = i0; i1 < 2; ++i1) {
          for (int j1 = j0; j1 < 2; ++j1) {
            for (int k1 = k0; k1 < 2; ++k1) {
              if (i1 - i0 + j1 - j0 + k1 - k0 != 1) continue;
              const double x0 = i0 ? m_xMinBox : m_xMaxBox;
              const double y0 = j0 ? m_yMinBox : m_yMaxBox;
              const double z0 = k0 ? m_zMinBox : m_zMaxBox;
              const double x1 = i1 ? m_xMinBox : m_xMaxBox;
              const double y1 = j1 ? m_yMinBox : m_yMaxBox;
              const double z1 = k1 ? m_zMinBox : m_zMaxBox;
              TMatrixD xMat(3, 1);
              if (!PlaneCut(x0, y0, z0, x1, y1, z1, xMat)) continue;
              if (m_debug) {
                std::cout << m_className << "::IntersectPlaneArea:\n"
                          << "    Intersection of plane at (" << xMat(0, 0) 
                          << ", " << xMat(1, 0) << ", " << xMat(2, 0) 
                          << ") with edge\n    (" 
                          << x0 << ", " << y0 << ", " << z0 << ")-(" 
                          << x1 << ", " << y1 << ", " << z1 << ")\n";
              }
              // Do not add same points (the case when plane contains an edge)
              bool skip = false;
              for (auto& p : intersect_points) {
                const double dx = xMat(0, 0) - p(0, 0);
                const double dy = xMat(1, 0) - p(1, 0);
                if (std::hypot(dx, dy) < 1e-10) {
                  skip = true;
                  break;
                }
              }
              if (!skip) intersect_points.push_back(xMat);
            }
          }
        }
      }
    }
  }
  if (intersect_points.size() < 3) {
    std::cerr << m_className << "::IntersectPlaneArea:\n"
              << "    WARNING: Empty intersection of view plane with area.\n";
    return false;
  }
  TMatrixD offset = intersect_points[0];
  xmin = xmax = intersect_points[0](0, 0);
  ymin = ymax = intersect_points[0](1, 0);
  // Remove crossings by sorting points rotation-wise.
  for (auto& p : intersect_points) p -= offset;
  std::sort(intersect_points.begin(), intersect_points.end(),
            [](const TMatrixD& a, const TMatrixD& b) -> bool {
              double cross_z = a(0, 0) * b(1, 0) - a(1, 0) * b(0, 0);
              return cross_z < 0;
            });
  for (auto& p : intersect_points) {
    p += offset;
    m_viewRegionX.push_back(p(0, 0));
    m_viewRegionY.push_back(p(1, 0));
    xmin = std::min(p(0, 0), xmin);
    ymin = std::min(p(1, 0), ymin);
    xmax = std::max(p(0, 0), xmax);
    ymax = std::max(p(1, 0), ymax);
  }
  return true;
}

// Ported from Garfield (function INTERD):
// Returns true if the point (x,y) is inside of the specified polygon.
// x: the x-coordinate
// y: the y-coordinate
// px: the x-vertices of the polygon
// py: the y-vertices of the polygon
// edge: a variable set to true if the point is located on the polygon edge
bool ViewFEMesh::IsInPolygon(double x, double y, 
                             const std::vector<double>& px,
                             const std::vector<double>& py, bool& edge) const {
  // Get the number and coordinates of the polygon vertices.
  const size_t pN = px.size();

  // Handle the special case of less than 2 vertices.
  if (pN < 2) return false;
  // Handle the special case of exactly 2 vertices (a line).
  if (pN == 2) return OnLine(px[0], py[0], px[1], py[1], x, y);

  // Set the minimum and maximum coordinates of all polygon vertices.
  double px_min = px[0], py_min = py[0];
  double px_max = px[0], py_max = py[0];
  for (size_t i = 0; i < pN; i++) {
    px_min = std::min(px_min, px[i]);
    py_min = std::min(py_min, py[i]);
    px_max = std::max(px_max, px[i]);
    py_max = std::max(py_max, py[i]);
  }

  // Set the tolerances
  double xtol = 1.0e-10 * std::max(std::abs(px_min), std::abs(px_max));
  double ytol = 1.0e-10 * std::max(std::abs(py_min), std::abs(py_max));
  if (xtol <= 0) xtol = 1.0e-10;
  if (ytol <= 0) ytol = 1.0e-10;

  // If we have essentially one x value, check to see if y is in range.
  if (std::abs(px_max - px_min) < xtol) {
    edge = (y > (py_min - ytol) && y < (py_max + ytol) &&
            std::abs(px_max + px_min - 2 * x) < xtol);
    return false;
  }
  // If we have essentially one y value, check to see if x is in range.
  if (std::abs(py_max - py_min) < ytol) {
    edge = (x > (px_min - xtol) && x < (px_max + xtol) &&
            std::abs(py_max + py_min - 2 * y) < ytol);
    return false;
  }

  // Set "infinity" points.
  double xinf = px_min - std::abs(px_max - px_min);
  double yinf = py_min - std::abs(py_max - py_min);

  // Loop until successful or maximum iterations (100) reached.
  int niter = 0;
  bool done = false;
  int ncross = 0;

  while (!done && niter < 100) {
    // Assume we will finish on this loop.
    done = true;

    // Loop over all edges, counting the number of edges crossed by a line
    // extending from (x, y) to (xinf, yinf).
    ncross = 0;
    for (size_t i = 0; (done && i < pN); i++) {
      // Determine whether the point lies on the edge.
      if (OnLine(px[i % pN], py[i % pN], px[(i + 1) % pN], py[(i + 1) % pN], x,
                 y)) {
        edge = true;
        return false;
      }

      // Determine whether this edge is crossed; if so increment the counter.
      double xc = 0., yc = 0.;
      if (LinesCrossed(x, y, xinf, yinf, px[i % pN], py[i % pN],
                       px[(i + 1) % pN], py[(i + 1) % pN], xc, yc))
        ncross++;

      // Ensure this vertex is not crossed by the line from (x,y)
      //  to (xinf,yinf); if so recompute (xinf,yinf) and start over.
      if (OnLine(x, y, xinf, yinf, px[i], py[i])) {
        // Recompute (xinf,yinf).
        xinf = px_min - RndmUniform() * std::abs(px_max - xinf);
        yinf = py_min - RndmUniform() * std::abs(py_max - yinf);

        // Start over.
        done = false;
        niter++;
      }
    }
  }

  // If we failed to finish iterating, return false.
  if (niter >= 100) {
    std::cerr << m_className << "::IsInPolygon: Unable to determine whether ("
              << x << ", " << y << ") is inside a polygon. Returning false.\n";
    return false;
  }

  // Point is inside for an odd, nonzero number of crossings.
  return (ncross != 2 * (ncross / 2));
}

// Ported from Garfield (method GRCONV):
// Clip the specified polygon to the view region; return the clipped polygon.
// px: the x-vertices of the polygon
// py: the y-vertices of the polygon
// cx: to contain the x-vertices of the clipped polygon
// cy: to contain the y-vertices of the clipped polygon
void ViewFEMesh::ClipToView(std::vector<double>& px, std::vector<double>& py,
                            std::vector<double>& cx, std::vector<double>& cy) {
  // Get the number and coordinates of the polygon vertices.
  int pN = (int)px.size();

  // Clear the vectors to contain the final polygon.
  cx.clear();
  cy.clear();

  // Set up the view vertices.
  const auto& vx = m_viewRegionX;
  const auto& vy = m_viewRegionY;
  const int vN = m_viewRegionX.size();

  // Do nothing if we have less than 2 points.
  if (pN < 2) return;

  // Loop over the polygon vertices.
  for (int i = 0; i < pN; i++) {
    // Flag for skipping check for edge intersection.
    bool skip = false;

    // Loop over the view vertices.
    for (int j = 0; j < vN; j++) {
      // Determine whether this vertex lies on a view edge:
      // if so add the vertex to the final polygon.
      if (OnLine(vx[j % vN], vy[j % vN], vx[(j + 1) % vN], vy[(j + 1) % vN],
                 px[i], py[i])) {
        // Add the vertex.
        cx.push_back(px[i]);
        cy.push_back(py[i]);

        // Skip edge intersection check in this case.
        skip = true;
      }

      // Determine whether a corner of the view area lies on this edge:
      // if so add the corner to the final polygon.
      if (OnLine(px[i % pN], py[i % pN], px[(i + 1) % pN], py[(i + 1) % pN],
                 vx[j], vy[j])) {
        // Add the vertex.
        cx.push_back(vx[j]);
        cy.push_back(vy[j]);

        // Skip edge intersection check in this case.
        skip = true;
      }
    }

    // If we have not skipped the edge intersection check, look for an
    // intersection between this edge and the view edges.
    if (skip) continue;
    // Loop over the view vertices.
    for (int j = 0; j < vN; j++) {
      // Check for a crossing with this edge;
      // if one exists, add the crossing point.
      double xc = 0., yc = 0.;
      if (LinesCrossed(vx[j % vN], vy[j % vN], vx[(j + 1) % vN],
                       vy[(j + 1) % vN], px[i % pN], py[i % pN],
                       px[(i + 1) % pN], py[(i + 1) % pN], xc, yc)) {
        // Add a vertex.
        cx.push_back(xc);
        cy.push_back(yc);
      }
    }
  }

  // Find all view field vertices inside the polygon.
  for (int j = 0; j < vN; j++) {
    // Test whether this vertex is inside the polygon.
    // If so, add it to the final polygon.
    bool edge = false;
    if (IsInPolygon(vx[j], vy[j], px, py, edge)) {
      // Add the view vertex.
      cx.push_back(vx[j]);
      cy.push_back(vy[j]);
    }
  }

  // Find all polygon vertices inside the box.
  for (int i = 0; i < pN; i++) {
    // Test whether this vertex is inside the view.
    // If so, add it to the final polygon.
    bool edge = false;
    if (IsInPolygon(px[i], py[i], vx, vy, edge)) {
      // Add the polygon vertex.
      cx.push_back(px[i]);
      cy.push_back(py[i]);
    }
  }
}
}  // namespace Garfield
