#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "NR.h"

#define NR_END 1
#define FREE_ARG char *

#ifdef __cplusplus
namespace neBEM {
#endif

void nrerror(const char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr, " run_time error......\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, ".....now exiting to system.....\n");
  exit(1);
}

float *fvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl...nh] */
{
  float *v;

  v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
  if (!v) nrerror("allocation failure in vector()");
  return v - nl + NR_END;
}

int *ivector(long nl, long nh)
/* allocate a int vector with subscript range v[nl...nh] */
{
  int *v;

  v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v - nl + NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl...nh] */
{
  unsigned char *v;

  v = (unsigned char *)malloc(
      (size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector()");
  return v - nl + NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl...nh] */
{
  unsigned long *v;

  v = (unsigned long *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
  if (!v) nrerror("allocation failure in lvector()");
  return v - nl + NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl...nh] */
{
  double *v;

  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v - nl + NR_END;
}

float **fmatrix(long nrl, long nrh, long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl...nrh][ncl...nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointer to rows */
  return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointer to row */
  m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
  if (!m) nrerror("allocation faliore 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate row and set pointer to them */
  m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl]) nrerror("allocation faliore 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointer to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an int with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointer to row */
  m = (int **)malloc((size_t)((nrow + NR_END) * sizeof(int *)));
  if (!m) nrerror("allocation faliore 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate row and set pointer to them */
  m[nrl] = (int *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
  if (!m[nrl]) nrerror("allocation falior 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointer to rows */
  return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long /*oldch*/,
                  long newrl, long newcl)
/* point a submatrix [newrl...][newcl...] to a [oldrl...][oldcl...] */
{
  long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
  float **m;

  /* allocate array of pointers to rows */
  m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
  if (!m) nrerror(" allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /*set pointers to rows */
  for (i = oldrl, j = newrl; i <= oldrh; i++, j++) m[j] = a[i] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
   declared in the standard C manner as a[nrow][ncol], where nrh-nrl+1 and
   ncol=nch-ncl+1. the routine should be called with the address &a[0][0]
   as the first argument. */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
  if (!m) nrerror(" allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /*set pointers to rows */
  m[nrl] = a - ncl;
  for (i = 1, j = nrl + 1; i < nrow; i++, j++) m[j] = m[j - 1] + ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][mdl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t = (float ***)malloc((size_t)((nrow + NR_END) * sizeof(float **)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (float **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float *)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] =
      (float *)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
  }
  /* return pointer to array of pointers to rows */
  return t;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)

/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][mdl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t = (double ***)malloc((size_t)((nrow + NR_END) * sizeof(double **)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] =
      (double **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double *)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = (double *)malloc(
      (size_t)((nrow * ncol * ndep + NR_END) * sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
  }
  /* return pointer to array of pointers to rows */
  return t;
}

// from f4tensor.txt
double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                    long nwl, long nwh)
// allocate a double 4tensor with range
// t[nrl..nrh][ncl..nch][ndl..ndh][nwl..nwh]
{
  long i, j, k, nrow = nrh - nrl + 1, ncol = nch - ncl + 1,
                ndep = ndh - ndl + 1, nwid = nwh - nwl + 1;
  double ****t;

  // allocate pointers to pointers to pointers to rows
  t = (double ****)malloc((size_t)((nrow + NR_END) * sizeof(double ***)));
  if (!t) nrerror("allocation failure 1 in d4tensor()");
  t += NR_END;
  t -= nrl;

  // allocate pointers to pointers to rows and set pointers to them
  t[nrl] =
      (double ***)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double **)));
  if (!t[nrl]) nrerror("allocation failure 2 in d4tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  // allocate pointers to rows and set pointers to them
  t[nrl][ncl] = (double **)malloc(
      (size_t)((nrow * ncol * ndep + NR_END) * sizeof(double *)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d4tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  // allocate rows and set pointers to them
  t[nrl][ncl][ndl] = (double *)malloc(
      (size_t)((nrow * ncol * ndep * nwid + NR_END) * sizeof(double)));
  if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in d4tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nwl;

  for (i = nrl; i <= nrh; i++) {
    if (i > nrl) {
      t[i] = t[i - 1] + ncol;
      t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
      t[i][ncl][ndl] = t[i - 1][ncl][ndl] + ncol * ndep * nwid;
    }
    for (j = ncl; j <= nch; j++) {
      if (j > ncl) {
        t[i][j] = t[i][j - 1] + ndep;
        t[i][j][ndl] = t[i][j - 1][ndl] + ndep * nwid;
      }

      for (k = ndl; k <= ndh; k++) {
        if (k > ndl) t[i][j][k] = t[i][j][k - 1] + nwid;
      }
    }
  }

  // return pointer to pointer to array of pointers to rows
  return t;
}

/* from f4tensor2.c - inability to start from non-zero index
double ****d4tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh, int
nsl, int nsh)
// for example:  d4tensor(0,NX-1,0,NY-1,0,NZ-1,0,NS)
// allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
{
        int i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,nsext=nsh-nsl+1;
        double ****tt;

        // allocate pointers to pointers to rows
        tt=(double ****) malloc((size_t)((nrow)*sizeof(double***)));
        if (!tt) nrerror("allocation failure 1 in d4tensor()");
        tt -= nrl;

        // allocate pointers to rows and set pointers to them
        tt[nrl]=(double ***) malloc((size_t)((nrow*ncol)*sizeof(double**)));
        if (!tt[nrl]) nrerror("allocation failure 2 in d4tensor()");
        tt[nrl] -= ncl;

        // allocate rows and set pointers to them
        tt[nrl][ncl]=(double **)
malloc((size_t)((nrow*ncol*ndep)*sizeof(double*))); if (!tt[nrl][ncl])
nrerror("allocation failure 3 in d4tensor()"); tt[nrl][ncl] -= ndl;

        // allocate rows and set pointers to them
        tt[nrl][ncl][ndl]=(double *)
malloc((size_t)((nrow*ncol*ndep*nsext)*sizeof(double))); if (!tt[nrl][ncl][ndl])
nrerror("allocation failure 4 in d4tensor()"); tt[nrl][ncl][ndl] -= nsl;


        for(i=nrl+1;i<=nrh;i++) tt[i]=tt[i-1]+ncol;
        //////////////////////////////////////////////////////////////
        for(j=ncl+1;j<=nch;j++) tt[nrl][j]=tt[nrl][j-1]+ndep;
        for(i=nrl+1;i<=nrh;i++)
        {
                tt[i][ncl]=tt[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++)
                tt[i][j]=tt[i][j-1]+ndep;
        }
           ////////////////////////////////////////////////////////////////
        for(k=ndl+1;k<=ndh;k++)  tt[nrl][ncl][k]=tt[nrl][ncl][k-1]+nsext;   //
(0;0;0~ndh).

                for(i=nrl+1;i<=nrh;i++) //  (1~nrh;0;0~ndh).
                {
              tt[i][ncl][ndl]=tt[i-1][ncl][ndl]+ncol*ndep*nsext;
                        for(k=ndl+1;k<=ndh;k++)
                          tt[i][ncl][k]=tt[i][ncl][k-1]+nsext;
                }
                for(j=ncl+1;j<=nch;j++) //   (0;1~nch;0~ndh).
                {
                  tt[nrl][j][ndl]=tt[nrl][j-1][ndl]+ndep*nsext;
                        for(k=ndl+1;k<=ndh;k++)
                          tt[nrl][j][k]=tt[nrl][j][k-1]+nsext;
                }

        for(i=nrl+1;i<=nrh;i++)                                               //
(1~nrh;1~nch;0~ndh).
        {
      tt[i][ncl][ndl]=tt[i-1][ncl][ndl]+ncol*ndep*nsext;
          for(j=ncl+1;j<=nch;j++)
          {
                  tt[i][j][ndl]=tt[i][j-1][ndl]+ndep*nsext;
                for(k=ndl+1;k<=ndh;k++)
                  tt[i][j][k]=tt[i][j][k-1]+nsext;
          }
        }
        // return pointer to array of pointers to rows
        return tt;
}


// original version
double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long
ndh,long nd4l,long nd4h)

// allocate a double d4tensor with range
t[nrl..nrh][ncl..nch][ndl..ndh][nd4l..nd4h]
{
   long i, j, k, nrow=nrh-nrl+1, ncol=nch-ncl+1,
ndep=ndh-ndl+1,nd4th=nd4h-nd4l+1; double ****t;

   // allocate pointers to pointers to pointers to rows
   t=(double ****) malloc((size_t)((nrow+NR_END)*sizeof(double***)));
   if(!t) nrerror("allocation failure 1 in d4tensor()");
   t += NR_END;
   t -= nrl;

   // allocate pointers to pointers to rows and set pointers to them
   t[nrl]=(double ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
   if(!t[nrl]) nrerror("allocation failure 2 in d4tensor()");
   t[nrl] += NR_END;
   t[nrl] -= ncl;

   // allocate pointers to rows and set pointers to them
   t[nrl][ncl]=(double **)
malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double*))); if(!t[nrl][ncl])
nrerror("allocation failure 3 in d4tensor()"); t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;

   // allocate  rows and set pointers to them
   t[nrl][ncl][ndl]=(double *)
malloc((size_t)((nrow*ncol*ndep*nd4th+NR_END)*sizeof(double)));
   if(!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in d4tensor()");
   t[nrl][ncl][ndl] += NR_END;
   t[nrl][ncl][ndl] -= nd4l;

   for(k=ndl+1; k<=ndh; k++) t[nrl][ncl][k]=t[nrl][ncl][k-1]+nd4th;
       for(j=ncl+1; j<=nch; j++)
        {
         t[nrl][j]=t[nrl][j-1]+ndep;
         t[nrl][j][ndl]=t[nrl][j-1][ndl]+ndep*nd4th;
         for(k=ndl+1; k<=ndh; k++)t[nrl][j][k]=t[nrl][j][k-1]+nd4th;
        }
     for(i=nrl+1; i<=nrh; i++)
      {
       t[i]=t[i-1]+ncol;
       t[i][ncl]=t[i-1][ncl]+ncol*ndep;
       for(j=ncl+1; j<=nch; j++)
        {
         t[i][j]=t[i][j-1]+ndep;
         t[i][j][ndl]=t[i][j-1][ndl]+ndep*nd4th;
         for(k=ndl+1; k<=ndh; k++)t[i][j][k]=t[i][j][k-1]+nd4th;
        }

      }
   // return pointer to array of pointers to rows
   return t;
}
*/

/* free a float vector  allocated with vector() */
void free_fvector(float *v, long nl, long /*nh*/) {
  free((FREE_ARG)(v + nl - NR_END));
}

/* free an int vector  allocated with vector() */
void free_ivector(int *v, long nl, long /*nh*/) {
  free((FREE_ARG)(v + nl - NR_END));
}

/* free an unsigned char vector  allocated with vector() */
void free_cvector(unsigned char *v, long nl, long /*nh*/) {
  free((FREE_ARG)(v + nl - NR_END));
}

/* free an unsigned long vector  allocated with vector() */
void free_lvector(unsigned long *v, long nl, long /*nh*/) {
  free((FREE_ARG)(v + nl - NR_END));
}

/* free a double vector  allocated with vector() */
void free_dvector(double *v, long nl, long /*nh*/) {
  free((FREE_ARG)(v + nl - NR_END));
}

/* free a float matrix allocated by matrix() */
void free_fmatrix(float **m, long nrl, long /*nrh*/, long ncl, long /*nch*/) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long /*nrh*/, long ncl, long /*nch*/) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

/* free an int matrix allocated by imatrix() */
void free_imatrix(int **m, long nrl, long /*nrh*/, long ncl, long /*nch*/) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

/* free a submatrix allocated by submatrix() */
void free_submatrix(float **b, long nrl, long /*nrh*/, long /*ncl*/,
                    long /*nch*/) {
  free((FREE_ARG)(b + nrl - NR_END));
}

/* free a matrix allocated by convert_matrix() */
void free_convert_matrix(float **b, long nrl, long /*nrh*/, long /*ncl*/,
                         long /*nch*/) {
  free((FREE_ARG)(b + nrl - NR_END));
}

/* free a f3tensor allocated by f3tensor() */
void free_f3tensor(float ***t, long nrl, long /*nrh*/, long ncl, long /*nch*/,
                   long ndl, long /*ndh*/) {
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

/* free a d3tensor allocated by d3tensor() */
void free_d3tensor(double ***t, long nrl, long /*nrh*/, long ncl, long /*nch*/,
                   long ndl, long /*ndh*/) {
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

// from f4tensor.txt
/* free a double d4tensor allocated by d4tensor() */
void free_d4tensor(double ****t, long nrl, long /*nrh*/, long ncl, long /*nch*/,
                   long ndl, long /*ndh*/, long nwl, long /*nwh*/) {
  free((FREE_ARG)(t[nrl][ncl][ndl] + nwl - NR_END));
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

/* from f4tensor2.c - inability to start from non-zero index
void free_d4tensor(double ****tt, int nrl, int nrh, int ncl, int nch,
        int ndl, int ndh, int nsl, int nsh)
// free a double d4tensor allocated by d4tensor()
{
        free((FREE_ARG) (tt[nrl][ncl][ndl]+nsl));
        free((FREE_ARG) (tt[nrl][ncl]+ndl));
        free((FREE_ARG) (tt[nrl]+ncl));
        free((FREE_ARG) (tt+nrl));
}

// original version
void free_d4tensor(double ****t,long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh, long nd4l, long nd4h)
// free a d4tensor allocated by d4tensor()
{
   free((FREE_ARG) (t[nrl][ncl][ndl]+nd4l-NR_END));
   free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
   free((FREE_ARG) (t[nrl]+ncl-NR_END));
   free((FREE_ARG) (t+nrl-NR_END));
}
*/

#ifdef __cplusplus
} // namespace
#endif
