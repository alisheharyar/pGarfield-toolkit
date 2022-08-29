/* given a matrix 'a[1...m][1...n]', this routine computes its singular
   value decomposition, [A] = [U][W](transpose)[V]. The matrix [U] replaces
   'a' on output. The diagonal matrix of singular values [W] is output as a
   vector 'w[1...n]'. The matrix [V] (not its transpose) is output as
   'v[1...n][1...n]'. 'm' must be greater than or equal to 'n'; if it is
   smaller, then 'a' should be fitted upto square with zero rows */

#include <math.h>
#include "NR.h"

#ifdef __cplusplus
namespace neBEM {
#endif

static double at, bt, ct;
static double maxarg1, maxarg2;

#define PYTHAG(a, b)                              \
  ((at = fabs(a)) > (bt = fabs(b))                \
       ? (ct = bt / at, at * sqrt(1.0 + ct * ct)) \
       : (bt ? (ct = at / bt, bt * sqrt(1.0 + ct * ct)) : 0.0))
#define MAX(a, b) \
  (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define SIGNS(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX_ITS 1000

void svdcmp(double **a, int m, int n, double *w, double **v) {
  int flag, i, its, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  double *rv1;

  if (m < n) nrerror("SVDCMP: A must be augmented with extra zeros.");

  rv1 = dvector(1, n);

  for (i = 1; i <= n; i++) /* householder reduction to bidiagonal form */
  {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i <= m) {
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : scale)
#endif
      for (k = i; k <= m; k++) scale += fabs(a[k][i]);
      if (scale) {
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : s)
#endif
        for (k = i; k <= m; k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }

        f = a[i][i];
        g = -SIGNS(sqrt(s), f);
        h = f * g - s;
        a[i][i] = f - g;

        if (i != n) {
          for (j = l; j <= n; j++) {
            s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : s)
#endif
            for (k = i; k <= m; k++) s += a[k][i] * a[k][j];
            f = s / h;
#ifdef _OPENMP
#pragma omp parallel for private(k)  // OMPCheck
#endif
            for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
          }
        }
#ifdef _OPENMP
#pragma omp parallel for private(k)
#endif
        for (k = i; k <= m; k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if (i <= m && i != n) {
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : scale)
#endif
      for (k = l; k <= n; k++) scale += fabs(a[i][k]);
      if (scale) {
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : s)
#endif
        for (k = l; k <= n; k++) {
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = -SIGNS(sqrt(s), f);
        h = f * g - s;
        a[i][l] = f - g;
#ifdef _OPENMP
#pragma omp parallel for private(k)
#endif
        for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
        if (i != m) {
          for (j = l; j <= m; j++) {
            s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : s)
#endif
            for (k = l; k <= n; k++) s += a[j][k] * a[i][k];
#ifdef _OPENMP
#pragma omp parallel for private(k)  // OMPCheck
#endif
            for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
          }
        }
#ifdef _OPENMP
#pragma omp parallel for private(k)
#endif
        for (k = l; k <= n; k++) a[i][k] *= scale;
      }
    }
    anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

  for (i = n; i >= 1; i--) /* acumulation of right-hand transformations */
  {
    if (i < n) {
      if (g) {
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
        for (j = l; j <= n; j++) /* double division to avoid possible */
          v[j][i] = (a[i][j] / a[i][l]) / g; /* underflow */
        for (j = l; j <= n; j++) {
          s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : s)
#endif
          for (k = l; k <= n; k++) s += a[i][k] * v[k][j];
#ifdef _OPENMP
#pragma omp parallel for private(k)  // OMPCheck
#endif
          for (k = l; k <= n; k++) v[k][j] += s * v[k][i];
        }
      }
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  for (i = n; i >= 1; i--) /* accumulation of the left-hand transformations */
  {
    l = i + 1;
    g = w[i];
    if (i < n) {
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for (j = l; j <= n; j++) a[i][j] = 0.0;
    }
    if (g) {
      g = 1.0 / g;
      if (i != n) {
        for (j = l; j <= n; j++) {
          s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : s)
#endif
          for (k = l; k <= m; k++) s += a[k][i] * a[k][j];
          f = (s / a[i][i]) * g;
#ifdef _OPENMP
#pragma omp parallel for private(k)  // OMPCheck
#endif
          for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
        }
      }
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for (j = i; j <= m; j++) a[j][i] *= g;
    } else {
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for (j = i; j <= m; j++) a[j][i] = 0.0;
    }

    ++a[i][i];
  }

  for (k = n; k >= 1; k--) /* diagonalization of the bidiagonal form */
    for (its = 1; its <= MAX_ITS; its++) /* loop over alllowed values */
    {
      flag = 1;
      for (l = k; l >= 1; l--) /* test for splitting */
      {
        nm = l - 1; /* note that rv1 is always zero */
        if ((double)(fabs(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }
        if ((double)(fabs(w[nm]) + anorm) == anorm) break;
      }
      if (flag) {
        c = 0.0; /* cancellation of rv1[l], if l > 1 */
        s = 1.0;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          if ((double)(fabs(f) + anorm) == anorm) break;
          g = w[i];
          h = PYTHAG(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g * h;
          s = (-f * h);
#ifdef _OPENMP
#pragma omp parallel for private(j, y, z)  // OMPCheck
#endif
          for (j = 1; j <= m; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z * s;
            a[j][i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) /* convergence */
      {
        if (z < 0.0) /* singular value is made nonnegative */
        {
          w[k] = -z;
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
          for (j = 1; j <= n; j++) v[j][k] = (-v[j][k]);
        }
        break;
      }

      if (its == MAX_ITS)
        nrerror("no convergence in 1000 SVDCMP iterations.\n");

      x = w[l]; /* shift from bottom 2-by-2 minor */
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGNS(g, f))) - h)) / x;

      c = s = 1.0; /* next QR transformation */
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g *= c;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
#ifdef _OPENMP
#pragma omp parallel for private(jj, x, z)  // OMPCheck
#endif
        for (jj = 1; jj <= n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }

        z = PYTHAG(f, h);
        w[j] = z; /* rotation can be arbitrary if z = 0 */
        if (z) {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
#ifdef _OPENMP
#pragma omp parallel for private(jj, y, z)  // OMPCheck
#endif
        for (jj = 1; jj <= m; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }

  free_dvector(rv1, 1, n);
}

#ifdef __cplusplus
} // namespace
#endif
