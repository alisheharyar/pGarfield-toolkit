/* get the inverse of a matrix using LU decomposition */
#include <math.h>
#include <stdio.h>

#include "NR.h"

#define TINY 1.0e-20

#ifdef __cplusplus
namespace neBEM {
#endif

void ludcmp(double **a, int n, int *index, double *d) {
  int i, j, k;
  int imax = 1;
  double big, sum, dum;
  double *vv; /* vv stores the implicit scaling of each row */

  vv = dvector(1, n); /* No looping yet. */
  *d = 1.0;           /* loop over rows to get implicit scaling info. */
#ifdef _OPENMP
#pragma omp parallel for private(j, big)
#endif
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if (fabs(a[i][j]) > big) big = fabs(a[i][j]);
    if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
    vv[i] = 1.0 / big;
  }

  for (j = 1; j <= n; j++) /* outer loop of Crout's method */
  {
    for (i = 1; i < j; i++) {
      sum = a[i][j];
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(- : sum)
#endif
      for (k = 1; k < i; k++)  // OMPChk: parallelization may not help much
        sum = sum - a[i][k] * a[k][j];
      a[i][j] = sum;
    }

    big = 0.0; /* search for largest pivotal element. */
    for (i = j; i <= n; i++) {
      sum = a[i][j];
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(- : sum)
#endif
      for (k = 1; k < j; k++)  // OMPChk: parallelization may not help much
        sum = sum - a[i][k] * a[k][j];
      a[i][j] = sum;

      dum = vv[i] * fabs(sum);
      if (dum >= big) { /* is the figure of merit better */
        big = dum;      /* than the best? */
        imax = i;
      }
    }

    if (j != imax) {
#ifdef _OPENMP
#pragma omp parallel for private(k, dum)
#endif
      for (k = 1; k <= n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);       /* change the parity of d. */
      vv[imax] = vv[j]; /* interchange the scale factor. */
    }
    index[j] = imax;
    if (a[j][j] == 0.0) /* for some applications on singular */
      a[j][j] = TINY;   /* matrices, it is desirable to */
                        /* substitute TINY for zero. */

    if (j != n) /* finally, divide by the pivot element. */
    {
      dum = 1.0 / a[j][j];
#ifdef _OPENMP
#pragma omp parallel for private(i)  // OMPCheck: may not help much
#endif
      for (i = j + 1; i <= n; i++) a[i][j] = a[i][j] * dum;
    }
  }

  free_dvector(vv, 1, n);
}

void lubksb(double **a, int n, int *index, double *b) {
  int i, ii = 0, ip, j;
  double sum;

  for (i = 1; i <= n; i++) /* forward substitution. */
  {
    ip = index[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii) {
#ifdef _OPENMP
#pragma omp parallel for private(j) reduction(- : sum)
#endif
      for (j = ii; j <= i; j++) sum = sum - a[i][j] * b[j];
    } else if (sum)
      ii = i;
    b[i] = sum;
  }

  for (i = n; i >= 1; i--) /* back substitution. */
  {
    sum = b[i];
#ifdef _OPENMP
#pragma omp parallel for private(j) reduction(- : sum)
#endif
    for (j = i + 1; j <= n; j++) sum = sum - a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

/*
void main(void)
{
double	**a, **y, d, *col;
int	n, i, j, *index;
char	filename[50];
FILE	*fp;

printf("Give filename: ");
scanf("%s", filename);
fp = fopen(filename, "r");
fscanf(fp, "%d", &n);
a = dmatrix(1, n, 1, n);
y = dmatrix(1, n, 1, n);
col = dvector(1, n);
index = ivector(1, n);
for(i = 1; i <= n; ++i)
        for(j = 1; j <= n; ++j)
                fscanf(fp, "%le", &a[i][j]);
fclose(fp);
fp = fopen("luc1.out", "w");
for(i = 1; i <= n; ++i)
        {
        for(j = 1; j <= n; ++j)
                fprintf(fp, "%f ", a[i][j]);
        fprintf(fp, "\n");
        }
fclose(fp);

ludcmp(a, n, index, &d);	Decompose the matrix just once.

fp = fopen("luc2.out", "w");
for(i = 1; i <= n; ++i)
        {
        for(j = 1; j <= n; ++j)
                fprintf(fp, "%f ", a[i][j]);
        fprintf(fp, "\n");
        }
fclose(fp);

for(j = 1; j <= n; j++)		Find inverse by columns.
        {
        for(i = 1; i <= n; i++)
                col[i] = 0.0;
        col[j] = 1.0;

        lubksb(a, n, index, col);
        for(i = 1; i <= n; i++)
                y[i][j] = col[i];
        }

fp = fopen("luc.out", "w");
for(i = 1; i <= n; ++i)
        {
        for(j = 1; j <= n; ++j)
                fprintf(fp, "%f ", y[i][j]);
        fprintf(fp, "\n");
        }
fclose(fp);
}
*/

#ifdef __cplusplus
}  // namespace
#endif
