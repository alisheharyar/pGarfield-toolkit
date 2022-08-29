#include "Garfield/ComponentParallelPlate.hh"

#include <TF1.h>
#include <TF2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "Garfield/GarfieldConstants.hh"

namespace Garfield {

ComponentParallelPlate::ComponentParallelPlate() : Component("ParallelPlate") {}

void ComponentParallelPlate::Setup(const int N, std::vector<double> eps,
                                   std::vector<double> d, const double V,
                                   std::vector<int> sigmaIndex) {

  // Here I switch conventions with the z-axis the direction of drift.
  std::vector<double> placeHolder(N + 1, 0);

  const int Nholder1 = eps.size();
  const int Nholder2 = d.size();
  if (N != Nholder1 || N != Nholder2) {
    std::cout << m_className
              << "::Inconsistency between the number of layers, permittivities "
                 "and thicknesses given.\n";
    return;
  } else if (N < 2) {

    std::cout << m_className
              << "::Setup:: Number of layers must be larger then 1.\n";
    return;
  }

  if (m_debug)
    std::cout << m_className << "::Setup:: Loading parameters.\n";
  m_epsHolder = eps;
  m_eps = placeHolder;

  m_dHolder = d;
  m_d = placeHolder;
  m_N = N + 1;
  m_V = V;

  m_sigmaIndex = sigmaIndex;

  std::vector<double> m_zHolder(N + 1);
  m_zHolder[0] = 0;
  for (int i = 1; i <= N; i++) {
    m_zHolder[i] = m_zHolder[i - 1] + m_dHolder[i - 1];

    if (m_debug)
      std::cout << m_className << "Setup:: layer " << i
                << ":: z = " << m_zHolder[i]
                << ", epsr = " << m_epsHolder[i - 1] << ".\n";
  }
  m_z = m_zHolder;

  if (m_debug)
    std::cout << m_className << "Setup:: Constructing matrices.\n";
  constructGeometryMatrices(m_N);

  if (m_debug)
    std::cout << m_className
              << "Setup:: Computing weighting potential functions.\n";
  setHIntegrand();
  setwpStripIntegrand();
  setwpPixelIntegrand();

  std::cout << m_className << "Setup:: Geometry with N = " << N
            << " layers set.\n";
}

bool ComponentParallelPlate::GetBoundingBox(double &x0, double &y0, double &z0,
                                            double &x1, double &y1,
                                            double &z1) {
  // If a geometry is present, try to get the bounding box from there.

  // Here I switch conventions back with the y-axis the direction of drift.
  if (m_geometry) {
    if (m_geometry->GetBoundingBox(x0, y0, z0, x1, y1, z1))
      return true;
  }
  z0 = -std::numeric_limits<double>::infinity();
  x0 = -std::numeric_limits<double>::infinity();
  z1 = +std::numeric_limits<double>::infinity();
  x1 = +std::numeric_limits<double>::infinity();

  y0 = 0.;
  y1 = m_z.back();
  return true;
}

double ComponentParallelPlate::IntegratePromptPotential(const Electrode &el,
                                                        const double x,
                                                        const double y,
                                                        const double z) {
  switch (el.ind) {
  case structureelectrode::Plane: {
    return wpPlane(z);
    break;
  }
  case structureelectrode::Pixel: {
    m_wpPixelIntegral.SetParameters(x, y, el.xpos, el.ypos, el.lx, el.ly,
                                    z); //(x,y,x0,y0,lx,ly,z)
    int im;
    double epsm;
    getLayer(z, im, epsm);
    double upLim = m_upperBoundIntegration;
    if (z == 0 || m_upperBoundIntegration / z > 200) {
      upLim = 200;
    } else {
      upLim *= 1 / z;
    }
    return m_wpPixelIntegral.Integral(0, upLim, 0, upLim, 1.e-12);
    break;
  }
  case structureelectrode::Strip: {
    m_wpStripIntegral.SetParameters(x, el.xpos, el.lx, z); //(x,x0,lx,z)
    int im;
    double epsm;
    getLayer(z, im, epsm);
    double upLim = m_upperBoundIntegration;
    if (z == 0 || m_upperBoundIntegration / z > 200) {
      upLim = 200;
    } else {
      upLim *= 1 / z;
    }
    return m_wpStripIntegral.Integral(0, upLim, 1.e-12);
    break;
  }
  default: {
    std::cerr << m_className << "::IntegratePromptPotential:\n"
              << "    Unknown electrode type.\n";
    return 0.;
  }
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double &ex,
                                           double &ey, double &ez, Medium *&m,
                                           int &status) {
  // Here I switch conventions back with the y-axis the direction of drift.

  ex = ey = ez = 0;

  int im = -1;
  double epsM = -1;
  if (!getLayer(y, im, epsM)) {
    if (m_debug)
      std::cout << m_className << "::ElectricField: Not inside geometry.\n";
    status = -6;
    return;
  }

  ey = constEFieldLayer(im);

  m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;

  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (epsM == 1) {
    status = 0;
  } else {
    status = -5;
  }
}

void ComponentParallelPlate::ElectricField(const double x, const double y,
                                           const double z, double &ex,
                                           double &ey, double &ez, double &v,
                                           Medium *&m, int &status) {
  // Here I switch conventions back with the y-axis the direction of drift.

  ex = ey = ez = v = 0;

  int im = -1;
  double epsM = -1;
  if (!getLayer(y, im, epsM)) {
    if (m_debug)
      std::cout << m_className << "::ElectricField: Not inside geometry.\n";
    status = -6;
    return;
  }

  ey = constEFieldLayer(im);

  v = -m_V - (y - m_z[im - 1]) * constEFieldLayer(im);
  for (int i = 1; i <= im - 1; i++) {
    v -= (m_z[i] - m_z[i - 1]) * constEFieldLayer(i);
  }

  m = m_geometry ? m_geometry->GetMedium(x, y, z) : m_medium;

  if (!m) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField: No medium at (" << x << ", "
                << y << ", " << z << ").\n";
    }
    status = -6;
    return;
  }

  if (epsM == 1) {
    status = 0;
  } else {
    status = -5;
  }
}

bool ComponentParallelPlate::GetVoltageRange(double &vmin, double &vmax) {
  if (m_V == 0)
    return false;

  if (m_V < 0) {
    vmin = m_V;
    vmax = 0;
  } else {
    vmin = 0;
    vmax = m_V;
  }
  return true;
}

double ComponentParallelPlate::WeightingPotential(const double x,
                                                  const double y,
                                                  const double z,
                                                  const std::string &label) {

  // Here I switch conventions back with the y-axis the direction of drift.

  double ret = 0.;

  for (auto &electrode : m_readout_p) {
    if (electrode.label == label) {
      double yin = y;
      if (!electrode.formAnode)
        yin = m_z.back() - y;
      if (!electrode.m_usegrid) {
        ret += IntegratePromptPotential(electrode, z, x, yin);
      } else {
        ret += FindWeightingPotentialInGrid(electrode, z, x, yin);
      }
    }
  }
  return ret;
}

void ComponentParallelPlate::Reset() {
  m_readout.clear();
  m_readout_p.clear();

  m_cMatrix.clear();
  m_vMatrix.clear();
  m_gMatrix.clear();
  m_wMatrix.clear();

  m_sigmaIndex.clear();
  m_eps.clear();
  m_d.clear();
  m_z.clear();

  m_N = 0;
  m_V = 0;

  m_medium = nullptr;
}

void ComponentParallelPlate::UpdatePeriodicity() {
  if (m_debug) {
    std::cerr << m_className << "::UpdatePeriodicity:\n"
              << "    Periodicities are not supported.\n";
  }
}

void ComponentParallelPlate::AddPixel(double x, double z, double lx_input,
                                      double lz_input, const std::string &label,
                                      bool anode) {

  // Here I switch conventions back with the y-axis the direction of drift.

  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it != m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPixel:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode pixel;
  pixel.label = label;
  pixel.ind = structureelectrode::Pixel;
  pixel.xpos = z;
  pixel.ypos = x;
  pixel.lx = lz_input;
  pixel.ly = lx_input;

  pixel.formAnode = anode;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(pixel));
  std::cout << m_className << "::AddPixel: Added pixel electrode.\n";
}

void ComponentParallelPlate::AddStrip(double z, double lz_input,
                                      const std::string &label, bool anode) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it != m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddStrip:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode strip;
  strip.label = label;
  strip.ind = structureelectrode::Strip;
  strip.xpos = z;
  strip.lx = lz_input;

  strip.formAnode = anode;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(strip));

  std::cout << m_className << "::AddStrip: Added strip electrode.\n";
}

void ComponentParallelPlate::AddPlane(const std::string &label, bool anode) {
  const auto it = std::find(m_readout.cbegin(), m_readout.cend(), label);
  if (it != m_readout.end() && m_readout.size() > 0) {
    std::cerr << m_className << "::AddPlane:\n"
              << "Note that the label " << label << " is already in use.\n";
  }
  Electrode plate;
  plate.label = label;
  plate.ind = structureelectrode::Plane;

  plate.formAnode = anode;

  m_readout.push_back(label);
  m_readout_p.push_back(std::move(plate));

  std::cout << m_className << "::AddPlane: Added plane electrode.\n";
}

Medium *ComponentParallelPlate::GetMedium(const double x, const double y,
                                          const double z) {
  if (m_geometry) {
    return m_geometry->GetMedium(x, y, z);
  } else if (m_medium) {
    return m_medium;
  }
  return nullptr;
}

bool ComponentParallelPlate::Nsigma(
    int N, std::vector<std::vector<int>> &sigmaMatrix) {
  int nCol = N - 1;
  int nRow = pow(2, N - 1);
  // array to store binary number
  std::vector<int> binaryNum(nCol, 0);

  for (int i = 0; i < nRow; i++) {
    if (decToBinary(i, binaryNum)) {
      sigmaMatrix.push_back(binaryNum);
      std::reverse(sigmaMatrix[i].begin(), sigmaMatrix[i].end());
      std::for_each(sigmaMatrix[i].begin(), sigmaMatrix[i].end(),
                    [](int &n) { n = 1 - 2 * n; });
    }
  }
  return true;
}

bool ComponentParallelPlate::Ntheta(
    int N, std::vector<std::vector<int>> &thetaMatrix,
    std::vector<std::vector<int>> &sigmaMatrix) {
  int nCol = N - 1;
  int nRow = pow(2, N - 1);

  std::vector<int> thetaRow(nCol, 1);
  std::vector<int> thetaRowReset(nCol, 1);

  for (int i = 0; i < nRow; i++) {
    for (int j = 0; j < nCol; j++) {
      for (int l = j; l < nCol; l++)
        thetaRow[j] *= sigmaMatrix[i][l];
    }
    thetaMatrix.push_back(thetaRow);
    thetaRow = thetaRowReset;
  }
  return true;
}

void ComponentParallelPlate::constructGeometryMatrices(const int N) {

  int nRow = N;

  std::vector<std::vector<int>> sigmaMatrix;
  std::vector<std::vector<int>> thetaMatrix;

  for (int n = 1; n <= nRow; n++) {
    // building sigma and theta matrices
    Nsigma(n, sigmaMatrix);
    Ntheta(n, thetaMatrix, sigmaMatrix);
    // store solution for row n-1
    m_sigmaMatrix.push_back(sigmaMatrix);
    m_thetaMatrix.push_back(thetaMatrix);
    // reset
    sigmaMatrix.clear();
    thetaMatrix.clear();
  }
}

void ComponentParallelPlate::constructGeometryFunction(const int N) {

  int nRow = N;
  int nCol = pow(2, N - 1);
  // reset
  m_cMatrix.clear();
  m_vMatrix.clear();
  m_gMatrix.clear();
  m_wMatrix.clear();

  std::vector<double> cHold(nCol, 1);
  std::vector<double> vHold(nCol, 0);
  std::vector<double> gHold(nCol, 1);
  std::vector<double> wHold(nCol, 0);

  for (int n = 1; n <= nRow; n++) {
    int ix1 = 0;
    int ix2 = 0;

    for (int i = 0; i < nCol; i++) {
      // cyclic permutation over the rows of sigma
      if (ix1 == pow(2, n - 1))
        ix1 = 0;
      if (ix2 == pow(2, N - n))
        ix2 = 0;
      // normalization
      cHold[i] *= 1 / pow(2, n - 1);
      gHold[i] *= 1 / pow(2, N - n);
      // summation for c and v
      if (n > 1) {
        for (int j = 0; j < n - 1; j++) {
          cHold[i] *= (m_eps[j] + m_sigmaMatrix[n - 1][ix1][j] * m_eps[j + 1]) /
                      m_eps[j + 1];
          vHold[i] += (m_thetaMatrix[n - 1][ix1][j] - 1) * m_d[j];
        }
      }
      // summation for g and w
      for (int j = 0; j < N - n; j++) {
        gHold[i] *= (m_eps[N - j - 1] +
                     m_sigmaMatrix[N - n][ix2][j] * m_eps[N - j - 2]) /
                    m_eps[N - j - 2];
        wHold[i] += (m_thetaMatrix[N - n][ix2][j] - 1) * m_d[N - 1 - j];
      }
      ix1++;
      ix2++;
    }

    // store solution for row n
    m_cMatrix.push_back(cHold);
    m_vMatrix.push_back(vHold);
    m_gMatrix.push_back(gHold);
    m_wMatrix.push_back(wHold);

    // reset
    std::fill(cHold.begin(), cHold.end(), 1);
    std::fill(vHold.begin(), vHold.end(), 0);
    std::fill(gHold.begin(), gHold.end(), 1);
    std::fill(wHold.begin(), wHold.end(), 0);
  }
}

void ComponentParallelPlate::setHIntegrand() {
  auto hFunction = [=](double *k, double * /*p*/) {
    double kk = k[0];
    double z = k[1];

    double hNorm = 0;
    double h = 0;

    int im = -1;
    double epsM = -1;
    if (!getLayer(z, im, epsM))
      return 0.;
    LayerUpdate(z, im, epsM);

    for (int i = 0; i < pow(2, m_N - im - 1); i++) {
      h += m_gMatrix[im][i] * sinh(kk * (m_wMatrix[im][i] + m_z.back() - z));
    }
    for (int i = 0; i < pow(2, m_N - 1); i++) {
      hNorm += m_cMatrix[m_N - 1][i] *
               sinh(kk * (m_vMatrix[m_N - 1][i] + m_z.back()));
    }
    return h * m_eps[0] / (m_eps[m_N - 1] * hNorm);
  };
  TF2 *hF = new TF2("hFunction", hFunction, 0, m_upperBoundIntegration, 0,
                    m_z.back(), 0);

  hF->Copy(m_hIntegrand);

  delete hF;
}

void ComponentParallelPlate::setwpPixelIntegrand() {
  auto intFunction = [=](double *k, double *p) {
    double kx = k[0];
    double ky = k[1];

    double K = sqrt(kx * kx + ky * ky);

    double x = p[0];
    double y = p[1];
    double x0 = p[2];
    double y0 = p[3];
    double wx = p[4];
    double wy = p[5];
    double z = p[6];

    double sol = cos(kx * (x - x0)) * sin(kx * wx / 2) * cos(ky * (y - y0)) *
                 sin(ky * wy / 2) * m_hIntegrand.Eval(K, z) / (kx * ky);

    return 4 * sol / (Pi * Pi);
  };

  TF2 *wpPixelIntegrand =
      new TF2("wpPixelIntegrand", intFunction, 0, m_upperBoundIntegration, 0,
              m_upperBoundIntegration, 7);
  wpPixelIntegrand->SetNpx(
      10000); // increasing number of points the function is evaluated on
  wpPixelIntegrand->SetNpy(10000);
  wpPixelIntegrand->Copy(m_wpPixelIntegral);

  delete wpPixelIntegrand;
}

void ComponentParallelPlate::setwpStripIntegrand() {
  auto intFunction = [=](double *k, double *p) {
    double kk = k[0];
    double x = p[0];
    double x0 = p[1];
    double wx = p[2];
    double z = p[3];
    double sol =
        cos(kk * (x - x0)) * sin(kk * wx / 2) * m_hIntegrand.Eval(kk, z) / kk;
    return 2 * sol / Pi;
  };
  TF1 *wpStripIntegrand =
      new TF1("wpStripIntegrand", intFunction, 0, m_upperBoundIntegration, 4);
  wpStripIntegrand->SetNpx(
      1000); // increasing number of points the function is evaluated on
  wpStripIntegrand->Copy(m_wpStripIntegral);

  delete wpStripIntegrand;
}

bool ComponentParallelPlate::decToBinary(int n, std::vector<int> &binaryNum) {
  int L = binaryNum.size();
  // counter for binary array
  int i = 0;
  while (n > 0) {
    if (i + 1 > L) {
      std::cerr << m_className
                << "::decToBinary: Size of binary exceeds amount of colomb.\n";
      return false; // Triggered if binary expression is larger then n.
    }
    // storing remainder in binary array
    binaryNum[i] = n % 2;
    n = n / 2;
    i++;
  }
  return true; // Succesfully
}

void ComponentParallelPlate::SetWeightingPotentialGrid(
    const double xmin, const double xmax, const double xsteps,
    const double ymin, const double ymax, const double ysteps,
    const double zmin, const double zmax, const double zsteps,
    const std::string &label) {

  for (auto &electrode : m_readout_p) {
    if (electrode.label == label) {
      if (electrode.m_usegrid) {
        std::cerr << m_className
                  << "::SetWeightingPotentialGrid: Overwriting grid.\n";
      }

      if (electrode.grid.SetMesh(xsteps, ysteps, zsteps, xmin, xmax, ymin, ymax,
                                 zmin, zmax)) {
        std::cerr << m_className << "::SetWeightingPotentialGrid: Mesh set for "
                  << label << ".\n";
      }

      electrode.grid.SaveWeightingField(this, label, label + "map", "xyz");

      LoadWeightingPotentialGrid(label);
    }
  }
}

void ComponentParallelPlate::SetWeightingPotentialGrids(
    const double xmin, const double xmax, const double xsteps,
    const double ymin, const double ymax, const double ysteps,
    const double zmin, const double zmax, const double zsteps) {

  for (auto &electrode : m_readout_p) {
    SetWeightingPotentialGrid(xmin, xmax, xsteps, ymin, ymax, ysteps, zmin,
                              zmax, zsteps, electrode.label);
  }
}

double ComponentParallelPlate::FindWeightingPotentialInGrid(Electrode &el,
                                                            const double x,
                                                            const double y,
                                                            const double z) {
  return el.grid.WeightingPotential(y, z, x, el.label);
}

} // namespace Garfield
