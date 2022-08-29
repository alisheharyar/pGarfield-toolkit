#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#ifdef DEFINE_NRGLOBAL
#define NRGLOBAL
#else
#define NRGLOBAL extern
#endif

#ifdef __cplusplus
namespace neBEM {
#endif

/*
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) ==0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1)>(dmaxarg2) ?\
                   (dmaxarg1):(dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1)>(dminarg2) ?\
                   (dminarg1):(dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1)>(maxarg2) ?\
                   (maxarg1):(maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1)>(minarg2) ?\
                   (minarg1):(minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b)(lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1)>(lmaxarg2) ?\
                  (lmaxarg1):(lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b)(lminarg1=(a),lminarg2=(b),(lminarg1)>(lminarg2) ?\
(lminarg1):(lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1)>(imaxarg2) ?\
                   (imaxarg1):(imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1)>(iminarg2) ?\
                   (iminarg1):(iminarg2))

#define SIGN(a,b) ((b)>0.0 ? fabs(a) : -fabs(a))
*/

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

NRGLOBAL void nrerror(const char error_text[]);

NRGLOBAL float *fvector(long nl, long nh);

NRGLOBAL int *ivector(long nl, long nh);

NRGLOBAL unsigned char *cvector(long nl, long nh);

NRGLOBAL double *dvector(long nl, long nh);

NRGLOBAL float **fmatrix(long nrl, long nrh, long ncl, long nch);

NRGLOBAL double **dmatrix(long nrl, long nrh, long ncl, long nch);

NRGLOBAL int **imatrix(long nrl, long nrh, long ncl, long nch);

NRGLOBAL float **submatrix(float **a, long oldrl, long oldrh, long oldcl,
                           long oldch, long newrl, long newcl);

NRGLOBAL float **convert_matrix(float *a, long nrl, long nrh, long ncl,
                                long nch);

NRGLOBAL int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl,
                         long ndh);

NRGLOBAL float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl,
                           long ndh);

NRGLOBAL double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl,
                            long ndh);

NRGLOBAL double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl,
                             long ndh, long nd4l, long nd4h);

NRGLOBAL void free_fvector(float *v, long nl, long nh);

NRGLOBAL void free_ivector(int *v, long nl, long nh);

NRGLOBAL void free_cvector(unsigned char *v, long nl, long nh);

NRGLOBAL void free_lvector(unsigned long *v, long nl, long nh);

NRGLOBAL void free_dvector(double *v, long nl, long nh);

NRGLOBAL void free_fmatrix(float **m, long nrl, long nrh, long ncl, long nch);

NRGLOBAL void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

NRGLOBAL void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

NRGLOBAL void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);

NRGLOBAL void free_convert_matrix(float **b, long nrl, long nrh, long ncl,
                                  long nch);

NRGLOBAL void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
                            long ndl, long ndh);

NRGLOBAL void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
                            long ndl, long ndh);

NRGLOBAL void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
                            long ndl, long ndh);

NRGLOBAL void free_d4tensor(double ****t, long nrl, long nrh, long ncl,
                            long nch, long ndl, long ndh, long nd4l, long nd4h);

#else /*ANSI*/
/* traditional - K&R */

NRGLOBAL void nrerror();
NRGLOBAL float *vector();

#endif /* ANSI */

#ifdef __cplusplus
}  // namespace
#endif

#endif /* _NR_UTILS_H_ */

#ifndef _NR_H_
#define _NR_H_

#ifdef __cplusplus
namespace neBEM {
#endif

NRGLOBAL void gaussj(double **a, int n, double *b, int m);

NRGLOBAL void ludcmp(double **a, int matsize, int *b, double *c);

NRGLOBAL void lubksb(double **a, int matsize, int *b, double *c);

NRGLOBAL void svdcmp(double **a, int matrow, int matcol, double *w, double **v);

NRGLOBAL void svbksb(double **a, double *w, double **v, int matrow, int matcol,
                     double *V, double *ChDen);

#ifdef __cplusplus
}  // namespace
#endif

#endif /* NR */
