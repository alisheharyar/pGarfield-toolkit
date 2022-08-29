#include <iostream>
#include <cstdio>
#include <cmath>
#include <limits>

#include <TROOT.h>
#include <TGraph.h>

#include "Garfield/Sensor.hh"
#include "Garfield/Component.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/ViewBase.hh"

namespace {

bool Invert(std::array<std::array<double, 3>, 3>& a) {

  // Compute cofactors.
  const double c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  const double c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
  const double c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
  const double c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
  const double c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
  const double c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
  const double c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
  const double c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
  const double c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  const double t1 = fabs(a[0][0]);
  const double t2 = fabs(a[1][0]);
  const double t3 = fabs(a[2][0]);
  double det = 0.;
  double pivot = 0.;
  if (t2 < t1 && t3 < t1) {
    pivot = a[0][0];
    det = c22 * c33 - c23 * c32;
  } else if (t1 < t2 && t3 < t2) {
    pivot = a[1][0];
    det = c13 * c32 - c12 * c33;
  } else {
    pivot = a[2][0];
    det = c23 * c12 - c22 * c13;
  }
  if (det == 0.) return false;
  const double s = pivot / det;
  a[0][0] = s * c11;
  a[0][1] = s * c21;
  a[0][2] = s * c31;
  a[1][0] = s * c12;
  a[1][1] = s * c22;
  a[1][2] = s * c32;
  a[2][0] = s * c13;
  a[2][1] = s * c23;
  a[2][2] = s * c33;
  return true;
}

std::string Fmt(const double x) {
  char buf[100];
  sprintf(buf, "%g", x);
  return std::string(buf);
}

}

namespace Garfield {

ViewBase::ViewBase(const std::string& name) :
    m_className(name) { 

  plottingEngine.SetDefaultStyle();
}

TPad* ViewBase::GetCanvas() {
  if (!m_pad) {
    std::string name = FindUnusedCanvasName("c" + m_className);
    if (!m_canvas) m_canvas.reset(new TCanvas(name.c_str(), ""));
    m_pad = m_canvas.get();
  }
  return m_pad;
}

bool ViewBase::RangeSet(TVirtualPad* pad) {
  if (!pad) return false;
  if (pad->GetListOfPrimitives()->GetSize() == 0 && 
      pad->GetX1() == 0 && pad->GetX2() == 1 && 
      pad->GetY1() == 0 && pad->GetY2() == 1) {
    return false;
  }
  return true;
}

void ViewBase::SetRange(TVirtualPad* pad, const double x0, const double y0,
                        const double x1, const double y1) {
  if (!pad) return;
  const double bm = pad->GetBottomMargin();
  const double lm = pad->GetLeftMargin();
  const double rm = pad->GetRightMargin();
  const double tm = pad->GetTopMargin();
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  pad->Range(x0 - dx * (lm / (1. - rm - lm)),
             y0 - dy * (bm / (1. - tm - lm)),
             x1 + dx * (rm / (1. - rm - lm)),
             y1 + dy * (tm / (1. - tm - lm)));

}

void ViewBase::SetArea(const double xmin, const double ymin, 
                       const double xmax, const double ymax) {
  // Check range, assign if non-null.
  if (xmin == xmax || ymin == ymax) {
    std::cerr << m_className << "::SetArea: Null area is not permitted.\n"
              << "      " << xmin << " < x < " << xmax << "\n"
              << "      " << ymin << " < y < " << ymax << "\n";
    return;
  }
  m_xMinPlot = std::min(xmin, xmax);
  m_yMinPlot = std::min(ymin, ymax);
  m_xMaxPlot = std::max(xmin, xmax);
  m_yMaxPlot = std::max(ymin, ymax);
  m_userPlotLimits = true;
}

void ViewBase::SetArea(const double xmin, const double ymin, const double zmin,
                       const double xmax, const double ymax,
                       const double zmax) {
  // Check range, assign if non-null
  if (xmin == xmax || ymin == ymax || zmin == zmax) {
    std::cerr << m_className << "::SetArea: Null area range not permitted.\n";
    return;
  }
  m_xMinBox = std::min(xmin, xmax);
  m_yMinBox = std::min(ymin, ymax);
  m_zMinBox = std::min(zmin, zmax);
  m_xMaxBox = std::max(xmin, xmax);
  m_yMaxBox = std::max(ymin, ymax);
  m_zMaxBox = std::max(zmin, zmax);
  m_userBox = true;
}

void ViewBase::SetPlane(const double fx, const double fy, const double fz,
                        const double x0, const double y0, const double z0) {
  // Calculate two in-plane vectors for the normal vector
  const double fnorm = sqrt(fx * fx + fy * fy + fz * fz);
  if (fnorm > 0 && fx * fx + fz * fz > 0) {
    const double fxz = sqrt(fx * fx + fz * fz);
    m_proj[0][0] = fz / fxz;
    m_proj[0][1] = 0;
    m_proj[0][2] = -fx / fxz;
    m_proj[1][0] = -fx * fy / (fxz * fnorm);
    m_proj[1][1] = (fx * fx + fz * fz) / (fxz * fnorm);
    m_proj[1][2] = -fy * fz / (fxz * fnorm);
    m_proj[2][0] = x0;
    m_proj[2][1] = y0;
    m_proj[2][2] = z0;
  } else if (fnorm > 0 && fy * fy + fz * fz > 0) {
    const double fyz = sqrt(fy * fy + fz * fz);
    m_proj[0][0] = (fy * fy + fz * fz) / (fyz * fnorm);
    m_proj[0][1] = -fx * fz / (fyz * fnorm);
    m_proj[0][2] = -fy * fz / (fyz * fnorm);
    m_proj[1][0] = 0;
    m_proj[1][1] = fz / fyz;
    m_proj[1][2] = -fy / fyz;
    m_proj[2][0] = x0;
    m_proj[2][1] = y0;
    m_proj[2][2] = z0;
  } else {
    std::cout << m_className << "::SetPlane:\n"
              << "    Normal vector has zero norm. No new projection set.\n";
  }

  // Store the plane description
  m_plane[0] = fx;
  m_plane[1] = fy;
  m_plane[2] = fz;
  m_plane[3] = fx * x0 + fy * y0 + fz * z0;

  UpdateProjectionMatrix();
}

void ViewBase::SetPlane(const double fx, const double fy, const double fz,
                        const double x0, const double y0, const double z0,
                        const double hx, const double hy, const double hz) {

  const double fnorm = sqrt(fx * fx + fy * fy + fz * fz);
  if (fnorm < Small) {
    std::cout << m_className << "::SetPlane:\n"
              << "    Normal vector has zero norm. No new projection set.\n";
    return;
  }
  // Normalise the vector.
  const double wx = fx / fnorm;
  const double wy = fy / fnorm;
  const double wz = fz / fnorm;
  // Store the plane description.
  m_plane[0] = wx;
  m_plane[1] = wy;
  m_plane[2] = wz;
  m_plane[3] = wx * x0 + wy * y0 + wz * z0;

  double d = hx * wx + hy * wy + hz * wz;
  double ux = hx - d * wx;
  double uy = hy - d * wy;
  double uz = hz - d * wz;
  double unorm = std::sqrt(ux * ux + uy * uy + uz * uz);
  if (unorm < 1.e-10) {  
    // Wrong in-plane x hint (close to norm).
    if (fy * fy + fz * fz > 0) {
      // Taking global x as in-plane x hint.
      ux = 1;
      uy = 0;
      uz = 0;  
    } else {
      // Taking global y as in-plane x hint.
      ux = 0;
      uy = 1;
      uz = 0;  
    }
    d = ux * wx + uy * wy + uz * wz;
    ux -= d * wx;
    uy -= d * wy;
    uz -= d * wz;
    unorm = std::sqrt(ux * ux + uy * uy + uz * uz);
  }
  ux /= unorm;
  uy /= unorm;
  uz /= unorm;

  m_prmat[0][0] = ux;
  m_prmat[1][0] = uy;
  m_prmat[2][0] = uz;
  // In-plane y = cross product [z,x]
  m_prmat[0][1] = wy * uz - wz * uy;
  m_prmat[1][1] = wz * ux - wx * uz;
  m_prmat[2][1] = wx * uy - wy * ux;
  m_prmat[0][2] = wx;
  m_prmat[1][2] = wy;
  m_prmat[2][2] = wz;

  for (unsigned int i = 0; i < 3; ++i) {
    m_proj[0][i] = m_prmat[i][0];
    m_proj[1][i] = m_prmat[i][1];
  }
  m_proj[2][0] = x0;
  m_proj[2][1] = y0;
  m_proj[2][2] = z0;
  if (!Invert(m_prmat)) {
    std::cerr << m_className << "::SetPlane:\n"
              << "    Inversion failed; reset to default.\n";
    SetPlaneXY();
  }
  if (m_debug) {
    std::cout << m_className << "::SetPlane:\n    PRMAT:\n";
    for (size_t i = 0; i < 3; ++i) {
      std::printf("  %10.5f  %10.5f  %10.5f\n", 
                  m_prmat[i][0], m_prmat[i][1], m_prmat[i][2]);
    }
    std::cout << "    PROJ:\n";
    for (size_t i = 0; i < 3; ++i) {
      std::printf("  %10.5f  %10.5f  %10.5f\n", 
                  m_proj[i][0], m_proj[i][1], m_proj[i][2]);
    }
    std::cout << "   PLANE:\n";
    std::printf("  %10.5f  %10.5f  %10.5f  %10.5f\n", 
                m_plane[0], m_plane[1], m_plane[2], m_plane[3]);
  }       
}

void ViewBase::Rotate(const double theta) {
  // Rotate the axes
  double auxu[3], auxv[3];
  const double ctheta = cos(theta);
  const double stheta = sin(theta);
  for (int i = 0; i < 3; ++i) {
    auxu[i] = ctheta * m_proj[0][i] - stheta * m_proj[1][i];
    auxv[i] = stheta * m_proj[0][i] + ctheta * m_proj[1][i];
  }
  for (int i = 0; i < 3; ++i) {
    m_proj[0][i] = auxu[i];
    m_proj[1][i] = auxv[i];
  }

  UpdateProjectionMatrix();
}

void ViewBase::SetPlaneXY() {
  m_proj = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 0}}};
  m_plane = {0, 0, 1, 0};
  m_prmat = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
}

void ViewBase::SetPlaneXZ() {
  m_proj = {{{1, 0, 0}, {0, 0, 1}, {0, 0, 0}}};
  m_plane = {0, 1, 0, 0};
  UpdateProjectionMatrix();
}

void ViewBase::SetPlaneYZ() {
  m_proj = {{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};
  m_plane = {1, 0, 0, 0};
  UpdateProjectionMatrix();
}

void ViewBase::SetPlaneZX() {
  m_proj = {{{0, 0, 1}, {1, 0, 0}, {0, 0, 0}}};
  m_plane = {0, 1, 0, 0};
  UpdateProjectionMatrix();
}

void ViewBase::SetPlaneZY() {
  m_proj = {{{0, 0, 1}, {0, 1, 0}, {0, 0, 0}}};
  m_plane = {1, 0, 0, 0};
  UpdateProjectionMatrix();
}

std::string ViewBase::FindUnusedFunctionName(const std::string& s) {
  int idx = 0;
  std::string fname = s + "_0";
  while (gROOT->GetListOfFunctions()->FindObject(fname.c_str())) {
    ++idx;
    fname = s + "_" + std::to_string(idx);
  }
  return fname;
}

std::string ViewBase::FindUnusedHistogramName(const std::string& s) {
  int idx = 0;
  std::string hname = s + "_0";
  while (gDirectory->GetList()->FindObject(hname.c_str())) {
    ++idx;
    hname = s + "_" + std::to_string(idx);
  }
  return hname;
}

std::string ViewBase::FindUnusedCanvasName(const std::string& s) {
  int idx = 0;
  std::string hname = s + "_0";
  while (gROOT->GetListOfCanvases()->FindObject(hname.c_str())) {
    ++idx;
    hname = s + "_" + std::to_string(idx);
  }
  return hname;
}

void ViewBase::UpdateProjectionMatrix() {

  m_prmat[0][0] = m_proj[0][0];
  m_prmat[1][0] = m_proj[0][1];
  m_prmat[2][0] = m_proj[0][2];
  m_prmat[0][1] = m_proj[1][0];
  m_prmat[1][1] = m_proj[1][1];
  m_prmat[2][1] = m_proj[1][2];
  const double vnorm = sqrt(m_plane[0] * m_plane[0] +
                            m_plane[1] * m_plane[1] +
                            m_plane[2] * m_plane[2]);
  if (vnorm <= 0.) {
    std::cerr << m_className << "::UpdateProjectionMatrix:\n"
              << "    Zero norm vector; reset to default.\n";
    SetPlaneXY();
    return;
  }
  m_prmat[0][2] = m_plane[0] / vnorm;
  m_prmat[1][2] = m_plane[1] / vnorm;
  m_prmat[2][2] = m_plane[2] / vnorm;
  if (!Invert(m_prmat)) {
    std::cerr << m_className << "::UpdateProjectionMatrix:\n"
              << "    Inversion failed; reset to default.\n";
    SetPlaneXY();
  }
}

void ViewBase::Clip(const std::array<float, 3>& x0, 
                    const std::array<float, 3>& x1,
                    std::array<float, 3>& xc) const {

  xc.fill(0.);
  const bool in0 = InBox(x0);
  const bool in1 = InBox(x1);
  if (in0 == in1) return;
  xc = in0 ? x1 : x0;
  const std::array<float, 3> dx = {x1[0] - x0[0], x1[1] - x0[1], 
                                   x1[2] - x0[2]};
  std::array<double, 3> bmin = {m_xMinBox, m_yMinBox, m_zMinBox};
  std::array<double, 3> bmax = {m_xMaxBox, m_yMaxBox, m_zMaxBox};
  for (size_t i = 0; i < 3; ++i) {
    if (dx[i] == 0. || (xc[i] >= bmin[i] && xc[i] <= bmax[i])) continue;
    const double b = xc[i] < bmin[i] ? bmin[i] : bmax[i];
    const double s = (b - xc[i]) / dx[i];
    xc[i] = b;
    for (size_t j = 0; j < 3; ++j) {
      if (j != i) xc[j] += dx[j] * s;
    }
  }
}

void ViewBase::DrawLine(const std::vector<std::array<float, 3> >& xl,
                        const short col, const short lw) {

  const size_t nP = xl.size();
  if (nP < 2) return;

  TGraph gr;
  gr.SetLineColor(col);
  gr.SetLineWidth(lw);
 
  std::vector<float> xgr;
  std::vector<float> ygr;
  auto x0 = xl[0];
  bool in0 = InBox(x0);
  if (in0) {
    float xp = 0., yp = 0.;
    ToPlane(x0[0], x0[1], x0[2], xp, yp);
    xgr.push_back(xp);
    ygr.push_back(yp);
  }
  for (unsigned int j = 1; j < nP; ++j) {
    auto x1 = xl[j];
    bool in1 = InBox(x1);
    if (in1 != in0) {
      float xp = 0., yp = 0.;
      std::array<float, 3> xc;
      Clip(x0, x1, xc);
      ToPlane(xc[0], xc[1], xc[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp);
    } 
    if (in1) {
      float xp = 0., yp = 0.;
      ToPlane(x1[0], x1[1], x1[2], xp, yp);
      xgr.push_back(xp);
      ygr.push_back(yp);
    } else if (!xgr.empty()) {
      gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Lsame");
      xgr.clear();
      ygr.clear();
    }
    x0 = x1;
    in0 = in1;
  }
  if (!xgr.empty()) {
    gr.DrawGraph(xgr.size(), xgr.data(), ygr.data(), "Lsame");
  }
}

std::string ViewBase::LabelX() {

  std::string xLabel = "";
  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[0][0] - 1) < tol) {
    xLabel = "#it{x}";
  } else if (fabs(m_proj[0][0] + 1) < tol) {
    xLabel = "#minus#it{x}";
  } else if (fabs(m_proj[0][0]) > tol) {
    xLabel = Fmt(m_proj[0][0]) + " #it{x}";
  }

  // y portion
  if (!xLabel.empty()) {
    if (m_proj[0][1] < -tol) {
      xLabel += " #minus ";
    } else if (m_proj[0][1] > tol) {
      xLabel += " #plus ";
    }
    if (fabs(m_proj[0][1] - 1) < tol || fabs(m_proj[0][1] + 1) < tol) {
      xLabel += "#it{y}";
    } else if (fabs(m_proj[0][1]) > tol) {
      xLabel += Fmt(fabs(m_proj[0][1])) + " #it{y}";
    }
  } else {
    if (fabs(m_proj[0][1] - 1) < tol) {
      xLabel = "#it{y}";
    } else if (fabs(m_proj[0][1] + 1) < tol) {
      xLabel = "#minus#it{y}";
    } else if (fabs(m_proj[0][1]) > tol) {
      xLabel = Fmt(m_proj[0][1]) + " #it{y}";
    }
  }

  // z portion
  if (!xLabel.empty()) {
    if (m_proj[0][2] < -tol) {
      xLabel += " #minus ";
    } else if (m_proj[0][2] > tol) {
      xLabel += " #plus ";
    }
    if (fabs(m_proj[0][2] - 1) < tol || fabs(m_proj[0][2] + 1) < tol) {
      xLabel += "#it{z}";
    } else if (fabs(m_proj[0][2]) > tol) {
      xLabel += Fmt(fabs(m_proj[0][2])) + " #it{z}";
    }
  } else {
    if (fabs(m_proj[0][2] - 1) < tol) {
      xLabel = "#it{z}";
    } else if (fabs(m_proj[0][2] + 1) < tol) {
      xLabel = "#minus#it{z}";
    } else if (fabs(m_proj[0][2]) > tol) {
      xLabel = Fmt(m_proj[0][2]) + " #it{z}";
    }
  }

  // Unit
  xLabel += " [cm]";
  return xLabel;

}

std::string ViewBase::LabelY() {

  std::string yLabel = "";
  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_proj[1][0] - 1) < tol) {
    yLabel = "#it{x}";
  } else if (fabs(m_proj[1][0] + 1) < tol) {
    yLabel = "#minus#it{x}";
  } else if (fabs(m_proj[1][0]) > tol) {
    yLabel = Fmt(m_proj[1][0]) + " #it{x}";
  }

  // y portion
  if (!yLabel.empty()) {
    if (m_proj[1][1] < -tol) {
      yLabel += " #minus ";
    } else if (m_proj[1][1] > tol) {
      yLabel += " #plus ";
    }
    if (fabs(m_proj[1][1] - 1) < tol || fabs(m_proj[1][1] + 1) < tol) {
      yLabel += "#it{y}";
    } else if (fabs(m_proj[1][1]) > tol) {
      yLabel += Fmt(fabs(m_proj[1][1])) + " #it{y}";
    }
  } else {
    if (fabs(m_proj[1][1] - 1) < tol) {
      yLabel = "#it{y}";
    } else if (fabs(m_proj[1][1] + 1) < tol) {
      yLabel = "#minus#it{y}";
    } else if (fabs(m_proj[1][1]) > tol) {
      yLabel = Fmt(m_proj[1][1]) + " #it{y}";
    }
  }

  // z portion
  if (!yLabel.empty()) {
    if (m_proj[1][2] < -tol) {
      yLabel += " #minus ";
    } else if (m_proj[1][2] > tol) {
      yLabel += " #plus ";
    }
    if (fabs(m_proj[1][2] - 1) < tol || fabs(m_proj[1][2] + 1) < tol) {
      yLabel += "#it{z}";
    } else if (fabs(m_proj[1][2]) > tol) {
      yLabel += Fmt(fabs(m_proj[1][2])) + " #it{z}";
    }
  } else {
    if (fabs(m_proj[1][2] - 1) < tol) {
      yLabel = "#it{z}";
    } else if (fabs(m_proj[1][2] + 1) < tol) {
      yLabel = "#minus#it{z}";
    } else if (fabs(m_proj[1][2]) > tol) {
      yLabel = Fmt(m_proj[1][2]) + " #it{z}";
    }
  }

  // Unit
  yLabel += " [cm]";
  return yLabel;
}

std::string ViewBase::PlaneDescription() {

  std::string description;

  constexpr double tol = 1.e-4;
  // x portion
  if (fabs(m_plane[0] - 1) < tol) {
    description = "x";
  } else if (fabs(m_plane[0] + 1) < tol) {
    description = "-x";
  } else if (fabs(m_plane[0]) > tol) {
    description = Fmt(m_plane[0]) + " x";
  }

  // y portion
  if (!description.empty()) {
    if (m_plane[1] < -tol) {
      description += " - ";
    } else if (m_plane[1] > tol) {
      description += " + ";
    }
    if (fabs(m_plane[1] - 1) < tol || fabs(m_plane[1] + 1) < tol) {
      description += "y";
    } else if (fabs(m_plane[1]) > tol) {
      description += Fmt(fabs(m_plane[1])) + " y";
    }
  } else {
    if (fabs(m_plane[1] - 1) < tol) {
      description = "y";
    } else if (fabs(m_plane[1] + 1) < tol) {
      description = "-y";
    } else if (fabs(m_plane[1]) > tol) {
      description = Fmt(m_plane[1]) + " y";
    }
  }

  // z portion
  if (!description.empty()) {
    if (m_plane[2] < -tol) {
      description += " - ";
    } else if (m_plane[2] > tol) {
      description += " + ";
    }
    if (fabs(m_plane[2] - 1) < tol || fabs(m_plane[2] + 1) < tol) {
      description += "z";
    } else if (fabs(m_plane[2]) > tol) {
      description += Fmt(fabs(m_plane[2])) + " z";
    }
  } else {
    if (fabs(m_plane[2] - 1) < tol) {
      description = "z";
    } else if (fabs(m_plane[2] + 1) < tol) {
      description = "-z";
    } else if (fabs(m_plane[2]) > tol) {
      description = Fmt(m_plane[2]) + " z";
    }
  }

  // Constant
  description += " = " + Fmt(m_plane[3]);
  return description;
}

bool ViewBase::PlotLimits(Sensor* sensor, 
                          double& xmin, double& ymin, 
                          double& xmax, double& ymax) const {

  if (!sensor) return false;
  // Try to get the area/bounding box from the sensor/component.
  std::array<double, 3> bbmin;
  std::array<double, 3> bbmax;
  if (!sensor->GetArea(bbmin[0], bbmin[1], bbmin[2], 
                       bbmax[0], bbmax[1], bbmax[2])) {
    std::cerr << m_className << "::PlotLimits:\n"
              << "    Sensor area is not defined.\n"
              << "    Please set the plot limits explicitly (SetArea).\n";
    return false;
  }
  return PlotLimits(bbmin, bbmax, xmin, ymin, xmax, ymax);
}

bool ViewBase::PlotLimits(Component* cmp, 
                          double& xmin, double& ymin, 
                          double& xmax, double& ymax) const {

  if (!cmp) return false;
  // Try to get the area/bounding box from the sensor/component.
  std::array<double, 3> bbmin;
  std::array<double, 3> bbmax;
  if (!cmp->GetBoundingBox(bbmin[0], bbmin[1], bbmin[2], 
                           bbmax[0], bbmax[1], bbmax[2])) {
    std::cerr << m_className << "::PlotLimits:\n"
              << "    Bounding box of the component is not defined.\n"
              << "    Please set the plot limits explicitly (SetArea).\n";
    return false;
  }
  if (std::isinf(bbmin[0]) || std::isinf(bbmax[0]) ||
      std::isinf(bbmin[1]) || std::isinf(bbmax[1]) ||
      std::isinf(bbmin[2]) || std::isinf(bbmax[2])) {
    std::array<double, 3> cellmin = {0., 0., 0.};
    std::array<double, 3> cellmax = {0., 0., 0.};
    if (!cmp->GetElementaryCell(cellmin[0], cellmin[1], cellmin[2],
                                cellmax[0], cellmax[1], cellmax[2])) {
      std::cerr << m_className << "::PlotLimits:\n"
                << "    Cell boundaries are not defined.\n"
                << "    Please set the plot limits explicitly (SetArea).\n";
      return false;
    }
    for (size_t i = 0; i < 3; ++i) {
      if (std::isinf(bbmin[i]) || std::isinf(bbmax[i])) {
        bbmin[i] = cellmin[i];
        bbmax[i] = cellmax[i];
      }
    }
  }
  return PlotLimits(bbmin, bbmax, xmin, ymin, xmax, ymax);
}

bool ViewBase::PlotLimitsFromUserBox(double& xmin, double& ymin,
                                     double& xmax, double& ymax) const {

  std::array<double, 3> bbmin = {m_xMinBox, m_yMinBox, m_zMinBox};
  std::array<double, 3> bbmax = {m_xMaxBox, m_yMaxBox, m_zMaxBox};
  return PlotLimits(bbmin, bbmax, xmin, ymin, xmax, ymax);
}

bool ViewBase::PlotLimits(std::array<double, 3>& bbmin,
                          std::array<double, 3>& bbmax,
                          double& xmin, double& ymin,
                          double& xmax, double& ymax) const {
  constexpr double tol = 1.e-4;
  double umin[2] = {-std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max()};
  double umax[2] = {std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max()};
  for (unsigned int i = 0; i < 3; ++i) {
    bbmin[i] -= m_proj[2][i];
    bbmax[i] -= m_proj[2][i];
    for (unsigned int j = 0; j < 2; ++j) {
      if (fabs(m_proj[j][i]) < tol) continue;
      const double t1 = bbmin[i] / m_proj[j][i];
      const double t2 = bbmax[i] / m_proj[j][i];
      const double tmin = std::min(t1, t2);
      const double tmax = std::max(t1, t2);
      if (tmin > umin[j] && tmin < umax[j]) umin[j] = tmin;
      if (tmax < umax[j] && tmax > umin[j]) umax[j] = tmax;
    }
  }
  xmin = umin[0];
  xmax = umax[0];
  ymin = umin[1];
  ymax = umax[1];
  return true;
}

}
