#ifndef BLAS_H
#define BLAS_H

#include <cassert>
#include <cstdlib>
#include <cmath>


inline unsigned GE(unsigned i, unsigned j, unsigned n)
{
  return i+j*n;
}

const double D_ZERO =  0.0;
const double D_ONE  =  1.0;
const double D_MONE = -1.0;
const float S_ZERO  =  0.0;
const float S_ONE   =  1.0;
const float S_MONE  = -1.0;

const double D_PREC = 1e-16;

const unsigned N_ONE = 1;
const int N_MONE = -1;
const char JOB_STR[] = "NTOSVULCRA";

extern "C"
{
  /******************************************************************/
  //double precision real
  /******************************************************************/
  unsigned idamax_(const unsigned*, const double*, const unsigned*);
  void dcopy_(const unsigned*, const double*, const unsigned*,
              double*, const unsigned*);
  void daxpy_(const unsigned*, const double*, const double*,
              const unsigned*, double*, const unsigned*);
  void dscal_(const unsigned*, const double*, const double*,
              const unsigned*);
  double ddot_(const unsigned*, const double*, const unsigned*,
               const double*, const unsigned*);
  double dnrm2_(const unsigned*, const double* const, const unsigned*);

  void dgtsv_(const unsigned*, const unsigned*, const double*,
              const double*, const double*, const double*, const unsigned*,
              const int*);
  void dgemm_(const char*, const char*, const unsigned*, const unsigned*,
              const unsigned*, const double*, const double*, const unsigned*,
              const double*, const unsigned*, const double*, double*,
              const unsigned*);
  void dger_(const unsigned*, const unsigned*, const double*, const double*,
             const unsigned*, double*, const unsigned*, const double*,
             const unsigned*);
  void dgemv_(const char*, const unsigned*, const unsigned*, const double*,
              const double*, const unsigned*, const double*, const unsigned*,
              const double*, double*, const unsigned*);
  void dorgqr_(const unsigned*, const unsigned*, const unsigned*,
               double*, const unsigned*, double*, double*,
               const unsigned*, int*);
  void dormqr_(const char*, const char*, const unsigned*, const unsigned*,
               const unsigned*, double*, const unsigned*, double*,
               double*, const unsigned*, double*, const unsigned*,
               int*);
  void dsyev_(const char*, const char*, const unsigned*, double*,
              const unsigned*, double*, double*, const unsigned*, int*);
  void dgeqrf_(const unsigned*, const unsigned*, double*, const unsigned*,
               double*, double*, const unsigned*, int*);
  void dgeqp3_(const unsigned*, const unsigned*, const double*,
               const unsigned*, const unsigned*, const double*, const double*,
               const unsigned*, int*);
  void dgeqpf_(const unsigned*, const unsigned*, const double*,
               const unsigned*, const unsigned*, const double*, const double*,
               const unsigned*, int*);
  void dgesvd_(const char*, const char*, const unsigned*, const unsigned*,
               double*, const unsigned*, double*, double*, const unsigned*,
               double*, const unsigned*, double*, const unsigned*, int*);
  void dgetrf_(const unsigned*, const unsigned*, double*, const unsigned*,
               unsigned*, int*);
  void dgetrs_(const char*, const unsigned*, const unsigned*, double*,
               const unsigned*, const unsigned*, double*, const unsigned*,
               int*);
  void dgetri_(const unsigned*, double*, const unsigned*, unsigned*, double*,
               const unsigned*, int*);
  void dspmv_(const char*, const unsigned*, const double*,
              const double*, const double*, const unsigned*,
              const double*, double*, const unsigned*);
  void dsptrf_(const char*, const unsigned*, const double*, int*, int*);
  void dsptri_(const char*, const unsigned*, const double*, int*, double*, int*);
  void dpotrf_(const char*, const unsigned*, const double*, const unsigned*,
               int*);
  void dpotri_(const char*, const unsigned*, const double*, const unsigned*,
               int*);
  void dpptrf_(const char*, const unsigned*, const double*, int*);
  void dpptri_(const char*, const unsigned*, const double*, int*);
  void dtptrs_(const char*, const char*, const char*, const unsigned*,

               const unsigned*, double*, const unsigned*, double*,
               double*, const unsigned*, double*, const unsigned*,
               int*);
  void dtpsv_(const char*, const char*, const char*, const unsigned*,
              double*, double*, const unsigned*);
  void dtrtrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, const double*, const unsigned*,
               double*, const unsigned*, int*);
  void dtrsm_(const char*, const char*, const char*, const char*,
              const unsigned*, const unsigned*, const double*, const double*,
              const unsigned*, double*, const unsigned*);
  void dtpmv_(const char*, const char*, const char*, const unsigned*,
              const double*, double*, const unsigned*);
  void dlacpy_(const char*, const unsigned*, const unsigned*, const double*,
               const unsigned*, double*, const unsigned*);
  void dlaset_(const char*, const unsigned*, const unsigned*, const double*,
               const double*, double*, const unsigned*);
  void dtrmm_(const char *, const char *, const char *, const char *,
              unsigned *, unsigned *, const double *, double *, unsigned *,
              double *, unsigned *);
  void dswap_(const unsigned*, double*, const unsigned*, double*,
	      const unsigned*);

  /******************************************************************/
  //single precision real
  /******************************************************************/
  unsigned isamax_(const unsigned*, const float*, const unsigned*);
  void scopy_(const unsigned*, const float*, const unsigned*,
              float*, const unsigned*);
  void saxpy_(const unsigned*, const float*, const float*,
              const unsigned*, float*, const unsigned*);
  void sscal_(const unsigned*, const float*, const float*,
              const unsigned*);
  float sdot_(const unsigned*, const float*, const unsigned*,
              const float*, const unsigned*);
  float snrm2_(const unsigned*, const float* const, const unsigned*);

  void sgtsv_(const unsigned*, const unsigned*, const float*,
              const float*, const float*, const float*, const unsigned*,
              const int*);
  void sgemm_(const char*, const char*, const unsigned*, const unsigned*,
              const unsigned*, const float*, const float*, const unsigned*,
              const float*, const unsigned*, const float*, float*,
              const unsigned*);
  void sger_(const unsigned*, const unsigned*, const float*, const float*,
             const unsigned*, float*, const unsigned*, const float*,
             const unsigned*);
  void sgemv_(const char*, const unsigned*, const unsigned*, const float*,
              const float*, const unsigned*, const float*, const unsigned*,
              const float*, float*, const unsigned*);
  void sorgqr_(const unsigned*, const unsigned*, const unsigned*,
               float*, const unsigned*, float*, float*,
               const unsigned*, int*);
  void sormqr_(const char*, const char*, const unsigned*, const unsigned*,
               const unsigned*, float*, const unsigned*, float*,
               float*, const unsigned*, float*, const unsigned*,
               int*);
  void ssyev_(const char*, const char*, const unsigned*, float*,
              const unsigned*, float*, float*, const unsigned*, int*);
  void stpsv_(const char*, const char*, const char*, const unsigned*,
              float*, float*, const unsigned*);
  void sgeqrf_(const unsigned*, const unsigned*, float*, const unsigned*,
               float*, float*, const unsigned*, int*);
  void sgesvd_(const char*, const char*, const unsigned*, const unsigned*,
               float*, const unsigned*, float*, float*, const unsigned*,
               float*, const unsigned*, float*, const unsigned*, int*);
  void sgetrf_(const unsigned*, const unsigned*, float*, const unsigned*,
               unsigned*, int*);
  void sgetri_(const unsigned*, float*, const unsigned*, unsigned*, float*,
               const unsigned*, int*);
  void sspmv_(const char*, const unsigned*, const float*,
              const float*, const float*, const unsigned*,
              const float*, float*, const unsigned*);
  void strsm_(const char*, const char*, const char*, const char*,
              const unsigned*, const unsigned*, const float*, const float*,
              const unsigned*, float*, const unsigned*);
  void ssptrf_(const char*, const unsigned*, const float*, int*, int*);
  void ssptri_(const char*, const unsigned*, const float*, int*, float*, int*);
  void spotrf_(const char*, const unsigned*, const float*, const unsigned*,
               int*);
  void spotri_(const char*, const unsigned*, const float*, const unsigned*,
               int*);
  void spptrf_(const char*, const unsigned*, const float*, int*);
  void spptri_(const char*, const unsigned*, const float*, int*);
  void stptrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, float*, float*, const unsigned*, int*);
  void strtrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, const float*, const unsigned*,
               float*, const unsigned*, int*);
  void stpmv_(const char*, const char*, const char*, const unsigned*,
              const float*, float*, const unsigned*);
  void sswap_(const unsigned*, float*, const unsigned*, float*,
	      const unsigned*);
}


namespace blas
{
  inline void swap(const unsigned n, double* x, const unsigned incx,
		   double* y, const unsigned incy )
  {
    dswap_(&n, x, &incx, y, &incy);
  }

  inline void swap(const unsigned n, float* x, const unsigned incx,
		   float* y, const unsigned incy )
  {
    sswap_(&n, x, &incx, y, &incy);
  }

  inline void laset(const unsigned m, const unsigned n, const double a,
		    const double b, double* A, unsigned ldA)
  {
    dlaset_(JOB_STR, &m, &n, &a, &b, A, &ldA);
  }
  inline void lasetu(const unsigned m, const unsigned n, const double a,
		     const double b, double* A, unsigned ldA)
  {
    dlaset_(JOB_STR+5, &m, &n, &a, &b, A, &ldA);
  }
  inline void lasetl(const unsigned m, const unsigned n, const double a,
		     const double b, double* A, unsigned ldA)
  {
    dlaset_(JOB_STR+6, &m, &n, &a, &b, A, &ldA);
  }

  inline unsigned maxi(const unsigned n, double* const v)
  {
    return idamax_(&n, v, &N_ONE);
  }
  inline unsigned maxi(const unsigned n, float* const v)
  {
    return isamax_(&n, v, &N_ONE);
  }

  inline void load(const unsigned n, double e, double* const v)
  {
    for (unsigned i=0; i<n; ++i) v[i] = e;
  }
  inline void load(const unsigned n, float e, float* const v)
  {
    for (unsigned i=0; i<n; ++i) v[i] = e;
  }

  inline void setzero(const unsigned n, unsigned* const v)
  {
    for (unsigned i=0; i<n; ++i) v[i] = 0;
  }
  inline void setzero(const unsigned n, double* const v)
  {
    load(n, D_ZERO, v);
  }
  inline void setzero(const unsigned n, float* const v)
  {
    load(n, S_ZERO, v);
  }
  inline double nrm2(const unsigned n, const double* const v)
  {
    return dnrm2_(&n, v, &N_ONE);
  }
  inline float nrm2(const unsigned n, const float* const v)
  {
    return snrm2_(&n, v, &N_ONE);
  }

  inline void copy(const unsigned n, const double*const orig, double* dest)
  {
    dcopy_(&n, orig, &N_ONE, dest, &N_ONE);
  }
  inline void copy(const unsigned n, const float*const orig, float* dest)
  {
    scopy_(&n, orig, &N_ONE, dest, &N_ONE);
  }

  inline void lacpy(const unsigned m, const unsigned n, double* A,
		    const unsigned ldA, double* B, const unsigned ldB)
  {
    dlacpy_(JOB_STR, &m, &n, A, &ldA, B, &ldB);
  }
  inline void lacpyu(const unsigned m, const unsigned n, double* A,
		     const unsigned ldA, double* B, const unsigned ldB)
  {
    dlacpy_(JOB_STR+5, &m, &n, A, &ldA, B, &ldB);
  }


  inline void copy(const unsigned n, double* orig, const unsigned inco,
		   double* dest, const unsigned incd)
  {
    dcopy_(&n, orig, &inco, dest, &incd);
  }
  inline void copy(const unsigned n, float* orig, const unsigned inco,
		   float* dest, const unsigned incd)
  {
    scopy_(&n, orig, &inco, dest, &incd);
  }

  // Conv.Copy double2float
  inline void copy(const unsigned n, double* orig, float* dest)
  {
    for (unsigned i=0; i<n; i++) dest[i] = (float) orig[i];
  }

  // Conv.Copy float2double
  inline void copy(const unsigned n, float* orig, double* dest)
  {
    for (unsigned i=0; i<n; i++) dest[i] = (double) orig[i];
  }

  // Scalar product conj(x)*y
  inline double scpr(const unsigned n, const double* const v1,
		     const double* const v2)
  {
    return ddot_(&n, v1, &N_ONE, v2, &N_ONE);
  }
  inline float scpr(const unsigned n, const float* const v1,
		    const float* const v2)
  {
    return sdot_(&n, v1, &N_ONE, v2, &N_ONE);
  }

  inline double sqrsum(const unsigned n, double* const v)
  {
    return ddot_(&n, v, &N_ONE, v, &N_ONE);
  }

  inline float sqrsum(const unsigned n, float* const v)
  {
    return sdot_(&n, v, &N_ONE, v, &N_ONE);
  }

  inline void add(const unsigned n, double* const x, double* const y)
  {
    daxpy_(&n, &D_ONE, x, &N_ONE, y, &N_ONE);
  }
  inline void add(const unsigned n, float* const x, float* const y)
  {
    saxpy_(&n, &S_ONE, x, &N_ONE, y, &N_ONE);
  }

  inline void axpy(const unsigned n, const double d, const double* const x,
		   double* const y)
  {
    daxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
  }
  inline void axpy(const unsigned n, const float d, const float* const x,
		   float* const y)
  {
    saxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
  }

  inline void scal(const unsigned n, const double d, double* const x, unsigned incx)
  { dscal_(&n, &d, x, &incx); }
  inline void scal(const unsigned n, const float d, float* const x, unsigned incx)
  { sscal_(&n, &d, x, &incx); }

  template<class T>
    inline void scal(const unsigned n, const T d, T* const x)
    { scal(n, d, x, N_ONE); }

  template<class T> inline void normalize(unsigned n, T* x)
    {
      T s = 1.0/blas::nrm2(n, x);
      blas::scal(n, s, x);
    }

  template<class T> inline void mkOrth(unsigned n, const T* v1, T* v2)
    {
      T s = -blas::scpr(n, v1, v2);
      blas::axpy(n, s, v1, v2);
    }

  template<class T>
    inline void scal(const unsigned n, const T d, T* const x, T* const y)
    {
      for (unsigned i=0; i<n; ++i) y[i] = d * x[i];
    }




  // y = d Ax
  inline void gemv(const unsigned m, const unsigned n, double d,
		   const double* A, double *x, double *y)
  {
    dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE);
  }
  inline void gemv(const unsigned m, const unsigned n, float d, const float* A,
		   float *x, float *y)
  {
    sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE);
  }

  // y += d Ax
  inline void gemva(const unsigned m, const unsigned n, double d, const double* A,
		    const double *x, double *y)
  {
    dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);
  }
  inline void gemva(const unsigned m, const unsigned n, float d, const float* A,
		    const float *x, float *y)
  {
    sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);
  }

  // y = d A^H x
  inline void gemhv(const unsigned m, const unsigned n, double d, const double* A,
		    const double *x, double *y)
  {
    dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE);
  }
  inline void gemhv(const unsigned m, const unsigned n, float d, const float* A,
		    const float *x, float *y)
  {
    sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE);
  }

  // y += d A^H x
  inline void gemhva(const unsigned m, const unsigned n, double d,
		     const double* A, unsigned ldA, const double *x, unsigned incx,
		     double *y, unsigned incy)
  {
    dgemv_(JOB_STR+1, &m, &n, &d, A, &ldA, x, &incx, &D_ONE, y, &incy);
  }
  inline void gemhva(const unsigned m, const unsigned n, float d,
		     const float* A, unsigned ldA, const float *x, unsigned incx,
		     float *y, unsigned incy)
  {
    sgemv_(JOB_STR+1, &m, &n, &d, A, &ldA, x, &incx, &S_ONE, y, &incy);
  }

  inline void gemhva(const unsigned m, const unsigned n, double d, const double* A,
		     const double *x, double *y)
  { gemhva(m, n, d, A, m, x, N_ONE, y, N_ONE); }
  inline void gemhva(const unsigned m, const unsigned n, float d, const float* A,
		     const float *x, float *y)
  { gemhva(m, n, d, A, m, x, N_ONE, y, N_ONE); }

  // y += d A x (A symm. dense in packed format)
  inline void gemsva(const unsigned n, double d, double* A,
		     double *x, double *y)
  {
    dspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &D_ONE, y, &N_ONE);
  }
  inline void gemsva(const unsigned n, float d, float* A,
		     float *x, float *y)
  {
    sspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &S_ONE, y, &N_ONE);
  }

  // sovles Ax=B, A is a triangluar Matrix
  inline void gtsv(const unsigned* n, const double* DiagLower,
		   const double* Diag, const double* DiagUpper,
		   const double* B, const int* INFO)
  {
    dgtsv_(n, &N_ONE, DiagLower, Diag, DiagUpper, B, n, INFO);
  }
  inline void gtsv(const unsigned* n, const float* DiagLower,
		   const float* Diag, const float* DiagUpper,
		   const float* B, const int* INFO)
  {
    sgtsv_(n, &N_ONE, DiagLower, Diag, DiagUpper, B, n, INFO);
  }

  // C = d A B, A is m x p, B is p x n
  inline void gemm(const unsigned m, const unsigned p, const unsigned n,
		   const double d, const double* const A, const unsigned ldA,
		   const double* const B, const unsigned ldB,
		   double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &D_ZERO, C, &ldC);
  }
  inline void gemm(const unsigned m, const unsigned p, const unsigned n,
		   const float d, const float* const A, const unsigned ldA,
		   const float* const B, const unsigned ldB,
		   float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &S_ZERO, C, &ldC);
  }

  // C += d A B, A is m x p, B is p x n
  inline void gemma(const unsigned m, const unsigned p, const unsigned n,
		    const double d, const double* const A, const unsigned ldA,
		    const double* const B, const unsigned ldB,
		    double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &D_ONE, C, &ldC);
  }
  inline void gemma(const unsigned m, const unsigned p, const unsigned n,
		    const float d, const float* const A, const unsigned ldA,
		    const float* const B, const unsigned ldB,
		    float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &S_ONE, C, &ldC);
  }

  // C = d A^H B, A is m x p, B is m x n
  inline void gemhm(const unsigned m, const unsigned p, const unsigned n,
		    const double d, const double* A, const unsigned ldA,
		    const double *B, const unsigned ldB,
		    double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
	   &D_ZERO, C, &ldC);
  }
  inline void gemhm(const unsigned m, const unsigned p, const unsigned n,
		    const float d, const float* const A, const unsigned ldA,
		    const float* const B, const unsigned ldB,
		    float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
	   &S_ZERO, C, &ldC);
  }

  // C += d A^H B, A is m x p, B is m x n
  inline void gemhma(unsigned m, unsigned p, unsigned n, double d,
		     const double* const A, const unsigned ldA, const double* const B,
		     const unsigned ldB, double* C, unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
	   &D_ONE, C, &ldC);
  }
  inline void gemhma(unsigned m, unsigned p, unsigned n, float d,
		     const float* const A, const unsigned ldA, const float* const B,
		     const unsigned ldB, float* C, unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
	   &S_ONE, C, &ldC);
  }

  // C = d A B^H, A is m x p, B is n x p
  inline void gemmh(const unsigned m, const unsigned p, const unsigned n,
		    const double d, const double* const A, const unsigned ldA,
		    const double* const B, const unsigned ldB,
		    double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &D_ZERO, C, &ldC);
  }
  inline void gemmh(const unsigned m, const unsigned p, const unsigned n,
		    const float d, const float* const A, const unsigned ldA,
		    const float* const B, const unsigned ldB,
		    float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &S_ZERO, C, &ldC);
  }

  // C += d A B^H, A is m x p, B is n x p
  inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
		     const double d, const double* const A, const unsigned ldA,
		     const double* const B, const unsigned ldB,
		     double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &D_ONE, C, &ldC);
  }
  inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
		     const float d, const float* const A, const unsigned ldA,
		     const float* const B, const unsigned ldB,
		     float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &S_ONE, C, &ldC);
  }
  inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
		     const double* const A, const unsigned ldA,
		     const double* const B, const unsigned ldB,
		     double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &D_ONE, A, &ldA, B, &ldB,
	   &D_ONE, C, &ldC);
  }
  inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
		     const float* const A, const unsigned ldA,
		     const float* const B, const unsigned ldB,
		     float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &S_ONE, A, &ldA, B, &ldB,
	   &S_ONE, C, &ldC);
  }

  // C = d A^H B^H, A is p x m, B is n x p
  inline void gemhmh(const unsigned m, const unsigned p, const unsigned n,
		     const double d, const double* const A, const unsigned ldA,
		     const double* const B, const unsigned ldB,
		     double* C, const unsigned ldC)
  {
    dgemm_(JOB_STR+1, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &D_ZERO, C, &ldC);
  }
  inline void gemhmh(const unsigned m, const unsigned p, const unsigned n,
		     const float d, const float* const A, const unsigned ldA,
		     const float* const B, const unsigned ldB,
		     float* C, const unsigned ldC)
  {
    sgemm_(JOB_STR+1, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
	   &S_ZERO, C, &ldC);
  }

  //C += d*AB, A is mxm (packed upper half is stored), B is mxn and regular matrix
  inline void sygemma(const unsigned m, const unsigned n,
		      const double* const A, const double* const B,
		      const double d, double* const C)
  {
    for(unsigned i=0;i<m;i++){
      for(unsigned j=0;j<n;j++){
	for(unsigned k=i;k<m;k++){
	  if(i==k){
	    C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
	  }else{
	    C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
	    C[j*m+k] += d*A[i+k*(k+1)/2]*B[i+j*m];
	  }
	}
      }
    }
  }
  inline void sygemma(const unsigned m, const unsigned n,
		      const float* const A, const float* const B,
		      const float d, float* const C)
  {
    for(unsigned i=0;i<m;i++){
      for(unsigned j=0;j<n;j++){
	for(unsigned k=i;k<m;k++){
	  if(i==k){
	    C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
	  }else{
	    C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
	    C[j*m+k] += d*A[i+k*(k+1)/2]*B[i+j*m];
	  }
	}
      }
    }
  }

  //C += d*AB, A is mxn and regular matrix, B is nxn (packed upper half is stored)
  inline void gesymma(const unsigned m, const unsigned n,
		      const double* const A, const double* const B,
		      const double d, double* const C)
  {
    for(unsigned i=0;i<m;i++){
      for(unsigned j=0;j<n;j++){
	for(unsigned k=j;k<n;k++){
	  if(j==k)
	    C[j*m+i] += d*A[i+k*m]*B[k+j*(j+1)/2];
	  else{
	    C[j*m+i] += d*A[i+k*m]*B[j+k*(k+1)/2];
	    C[k*m+i] += d*A[i+j*m]*B[j+k*(k+1)/2];
	  }
	}
      }
    }
  }
  inline void gesymma(const unsigned m, const unsigned n,
		      const float* const A, const float* const B,
		      const float d, float* const C)
  {
    for(unsigned i=0;i<m;i++){
      for(unsigned j=0;j<n;j++){
	for(unsigned k=j;k<n;k++){
	  if(j==k)
	    C[j*m+i] += d*A[i+k*m]*B[k+j*(j+1)/2];
	  else{
	    C[j*m+i] += d*A[i+k*m]*B[j+k*(k+1)/2];
	    C[k*m+i] += d*A[i+j*m]*B[j+k*(k+1)/2];
	  }
	}
      }
    }
  }

  // C += d A^H A, C is a symm. matrix (packed upper half is stored), A is mxn
  inline void symhm(const unsigned m, const unsigned n, const double* const A,
		    const double d, double* C)
  {
    for (unsigned j=0; j<n; ++j) {
      for (unsigned i=0; i<=j; ++i) {
	double sum = 0.0;
	for (unsigned k=0; k<m; ++k) sum += A[k+i*m] * A[k+j*m];
	C[i+j*(j+1)/2] += d * sum;
      }
    }
  }
  inline void symhm(const unsigned m, const unsigned n, const float* const A,
		    const float d, float* C)
  {
    for (unsigned j=0; j<n; ++j) {
      for (unsigned i=0; i<=j; ++i) {
	float sum = 0.0;
	for (unsigned k=0; k<m; ++k) sum += A[k+i*m] * A[k+j*m];
	C[i+j*(j+1)/2] += d * sum;
      }
    }
  }

  // C += d A A^H, C is a symm. matrix (packed upper half is stored), A is mxn
  inline void symmh(unsigned m, unsigned n, double* A, double d, double* C)
  {
    for (unsigned k=0; k<n; ++k) {
      for (unsigned j=0; j<=n; ++j) {
	double e = d * A[j+k*m];
	for (unsigned i=0; i<j; ++i) C[i+j*(j+1)/2] += e * A[i+k*m];
      }
    }
  }
  inline void symmh(unsigned m, unsigned n, float* A, float d, float* C)
  {
    for (unsigned k=0; k<n; ++k) {
      for (unsigned j=0; j<n; ++j) {
	float e = d * A[j+k*m];
	for (unsigned i=0; i<=j; ++i) C[i+j*(j+1)/2] += e * A[i+k*m];
      }
    }
  }

  // Singular Value Decomposition
  inline int gesvdS(unsigned m, unsigned n, double* A, double* S,
		    double* U, unsigned ldU, double* VT, unsigned ldVT,
		    unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+3, JOB_STR+3, &m, &n, A, &m, S, U, &ldU, VT, &ldVT,
	    wk, &nwk, &INF);
    return INF;
  }
  inline int gesvd(unsigned m, unsigned n, double* A, double* S,
		   double* VT, unsigned ldVT, unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
	    wk, &nwk, &INF);
    return INF;
  }
  inline int gesvd(unsigned m, unsigned n, float* A, float* S,
		   float* VT, unsigned ldVT, unsigned nwk, float* wk)
  {
    int INF;
    sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
	    wk, &nwk, &INF);
    return INF;
  }

  inline int gesvd(unsigned m, unsigned n, double* A, double* S,
		   double* U, unsigned ldU, double* VT, unsigned ldVT,
		   unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+9, JOB_STR+9, &m, &n, A, &m, S, U, &ldU, VT, &ldVT,
	    wk, &nwk, &INF);
    return INF;
  }

  // compute singular values
  inline int svals(unsigned m, unsigned n, double* A, double* S,
		   unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n, wk, &nwk, &INF);
    return INF;
  }
  inline int svals(unsigned m, unsigned n, float* A, float* S,
		   unsigned nwk, float* wk)
  {
    int INF;
    sgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n, wk, &nwk, &INF);
    return INF;
  }

  // triangular factorisation
  inline int getrf(const unsigned n, double* A, unsigned* ipiv)
  {
    int INF;
    dgetrf_(&n, &n, A, &n, ipiv, &INF);
    return INF;
  }
  inline int getrf(const unsigned n, float* A, unsigned* ipiv)
  {
    int INF;
    sgetrf_(&n, &n, A, &n, ipiv, &INF);
    return INF;
  }

  // upper triangular packed MV
  inline void utrpv(const unsigned n, double* A, double* x)
  {
    dtpmv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, x, &N_ONE);
  }
  inline void utrpv(const unsigned n, float* A, float* x)
  {
    stpmv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, x, &N_ONE);
  }

  // lower triangular packed transpose MV
  inline void ltrphv(const unsigned n, double* A, double* x)
  {
    dtpmv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, x, &N_ONE);
  }
  inline void ltrphv(const unsigned n, float* A, float* x)
  {
    stpmv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, x, &N_ONE);
  }

  // QR factorisation
  inline int geqrf(const unsigned m, const unsigned n, double* A,
		   double* tau, unsigned nwk, double* wk)
  {
    int INF;
    dgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int geqrf(const unsigned m, const unsigned n, float* A,
		   float* tau, unsigned nwk, float* wk)
  {
    int INF;
    sgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int geqrf(const unsigned m, const unsigned n, double* A,
		   const unsigned ldA, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    dgeqrf_(&m, &n, A, &ldA, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int geqrf(const unsigned m, const unsigned n, float* A,
		   const unsigned ldA, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    sgeqrf_(&m, &n, A, &ldA, tau, wk, &nwk, &INF);
    return INF;
  }

  // Multiply a general Matrix with the Q-Matrix (QR factorization), Q C
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
		   double* A, double* tau, double* C,
		   unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
		   float* A, float* tau, float* C,
		   unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
		   double* A, const unsigned ldA, double* tau, double* C,
		   const unsigned ldC, unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
	    &INF);
    return INF;
  }
  inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
		   float* A, const unsigned ldA, float* tau, float* C,
		   const unsigned ldC, unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
	    &INF);
    return INF;
  }

  // Q^H C
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
		    double* A, const unsigned ldA, double* tau, double* C,
		    const unsigned ldC, unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
	    &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
		    float* A, const unsigned ldA, float* tau, float* C,
		    const unsigned ldC, unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
	    &INF);
    return INF;
  }

  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
		    double* A, const unsigned ldA, double* tau, double* C,
		    unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &m, wk, &nwk,
	    &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
		    float* A, const unsigned ldA, float* tau, float* C,
		    unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &m, wk, &nwk,
	    &INF);
    return INF;
  }

  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
		    double* A, double* tau, double* C,
		    unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
		    float* A, float* tau, float* C,
		    unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  inline int morqr(const unsigned m, const unsigned n, const unsigned p,
		   double* A, const unsigned ldA, double* tau, double* C,
		   const unsigned ldC, unsigned nwk, double* wk)
  {
    int INF;
    dormqr_(JOB_STR+8, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
	    &INF);
    return INF;
  }
  inline int morqr(const unsigned m, const unsigned n, const unsigned p,
		   float* A, const unsigned ldA, float* tau, float* C,
		   const unsigned ldC, unsigned nwk, float* wk)
  {
    int INF;
    sormqr_(JOB_STR+8, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
	    &INF);
    return INF;
  }

  inline void ger(unsigned M, unsigned N, double d, double* X, unsigned INCX,
		  double* y, unsigned INCY, double* A, unsigned LDA)
  {
    dger_(&M, &N, &d, X, &INCX, y, &INCY, A, &LDA);
  }

  inline void ger(unsigned M, unsigned N, float d, float* X, unsigned INCX,
		  float* y, unsigned INCY, float* A, unsigned LDA)
  {
    sger_(&M, &N, &d, X, &INCX, y, &INCY, A, &LDA);
  }

  // return Q-Matrix (QR factorization) in A
  inline int orgqr(const unsigned m, const unsigned n, double* A, double* tau,
		   unsigned nwk, double* wk)
  {
    int INF;
    dorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int orgqr(const unsigned m, const unsigned n, float* A, float* tau,
		   unsigned nwk, float* wk)
  {
    int INF;
    sorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }

  // B=A^H, A in R^(m,n)
  inline void transpose(unsigned m, unsigned n, double* A, double* B)
  {
    for (unsigned i=0; i<m; ++i) dcopy_(&n, A+i, &m, B+i*n, &N_ONE);
  }
  inline void transpose(unsigned m, unsigned n, float* A, float* B)
  {
    for (unsigned i=0; i<m; ++i) scopy_(&n, A+i, &m, B+i*n, &N_ONE);
  }

  // product of an upper triangular matrix U and a matrix A, A:=U A
  inline void utrgemm(unsigned m, unsigned n, double* U, unsigned ldU,
		      double* A, unsigned ldA)
  {
    dtrmm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR, &m, &n, &D_ONE, U, &ldU,
	   A, &ldA);
  }


  // A += d U U^H, where U is upper triangular in packed storage
  // only the upper triangular part of A is computed
  inline void putrmmh(double d, unsigned n, double* U, double* A)
  {
    for (unsigned j=0; j<n; ++j) {
      for (unsigned l=j; l<n; ++l) {
	unsigned idl = l*(l+1)/2;
	double e = d * U[j+idl];
	for (unsigned i=0; i<=j; ++i) A[i+j*(j+1)/2] += e * U[i+idl];
      }
    }
  }

  inline void putrmmh(float d, unsigned n, float* U, float* A)
  {
    for (unsigned j=0; j<n; ++j) {
      for (unsigned l=j; l<n; ++l) {
	unsigned idl = l*(l+1)/2;
	float e = d * U[j+idl];
	for (unsigned i=0; i<=j; ++i) A[i+j*(j+1)/2] += e * U[i+idl];
      }
    }
  }

  inline void fill0_ltr(unsigned n, double* A)
  {
    for (unsigned j=0; j<n; ++j) {
      *A++ = (double) j;
      for (unsigned i=j+1; i<n; ++i) *A++ = D_ZERO;
    }
  }
  inline void fill0_ltr(unsigned n, float* A)
  {
    for (unsigned j=0; j<n; ++j) {
      *A++ = (float) j;
      for (unsigned i=j+1; i<n; ++i) *A++ = S_ZERO;
    }
  }

  // fill nxn matrix A with identity
  inline void fillId(unsigned n, double *A)
  {
    for (unsigned j=0; j<n; ++j) {
      unsigned i = 0;
      for (; i<j; ++i) *A++ = D_ZERO;
      *A++ = D_ONE;
      ++i;
      for (; i<n; ++i) *A++ = D_ZERO;
    }
  }
  inline void fillId(unsigned n, float *A)
  {
    for (unsigned j=0; j<n; ++j) {
      unsigned i = 0;
      for (; i<j; ++i) *A++ = S_ZERO;
      *A++ = S_ONE;
      ++i;
      for (; i<n; ++i) *A++ = S_ZERO;
    }
  }

  // fill nxn upper triang. matrix A with identity (packed storage)
  inline void fillId_utr(unsigned n, double *A)
  {
    for (unsigned j=0; j<n; ++j) {
      for (unsigned i=0; i<j; ++i) *A++ = D_ZERO;
      *A++ = D_ONE;
    }
  }
  inline void fillId_utr(unsigned n, float *A)
  {
    for (unsigned j=0; j<n; ++j) {
      for (unsigned i=0; i<j; ++i) *A++ = S_ZERO;
      *A++ = S_ONE;
    }
  }

  // fill nxn normalized lower triang. matrix A with identity (packed storage)
  inline void fillId_ltr(unsigned n, double *A)
  {
    for (unsigned i=0; i<n; ++i) {
      for (unsigned j=0; j<i; ++j) *A++ = D_ZERO;
      // for pivoting, a ltr is assumed to have ones on the diagonal
      *A++ = (double) i;
    }
  }
  inline void fillId_ltr(unsigned n, float *A)
  {
    for (unsigned i=0; i<n; ++i) {
      for (unsigned j=0; j<i; ++j) *A++ = S_ZERO;
      // for pivoting, a ltr is assumed to have ones on the diagonal
      *A++ = (float) i;
    }
  }

}


namespace lapack
{

  // general inversion
  inline void geinv(const unsigned n, double* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    double* const work = new double[lwork];
    assert(work!=NULL);

    int INFO = blas::getrf(n, A, ipiv);
    assert(INFO==0);

    dgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv(const unsigned n, float* const A)
  {
    unsigned* const ipiv = new unsigned[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    float* const work = new float[lwork];
    assert(work!=NULL);

    int INFO = blas::getrf(n, A, ipiv);
    assert(INFO==0);

    sgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }

  inline void geinv_sym(const unsigned n, double* const A)
  {
    int* const ipiv = new int[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    double* const work = new double[lwork];
    assert(work!=NULL);

    int INFO;
    dsptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
    assert(INFO==0);

    dsptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }
  inline void geinv_sym(const unsigned n, float* const A)
  {
    int* const ipiv = new int[n];
    assert(ipiv!=NULL);

    const unsigned lwork = 4*n;
    float* const work = new float[lwork];
    assert(work!=NULL);

    int INFO;
    ssptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
    assert(INFO==0);

    ssptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
    assert(INFO==0);

    delete [] work;
    delete [] ipiv;
  }

  // packed triangular factorisation of positive definite matrix
  inline int pptrf(const unsigned n, double* A)
  {
    int inf;
    dpptrf_(JOB_STR+5, &n, A, &inf);
    return inf;
  }
  inline int pptrf(const unsigned n, float* A)
  {
    int inf;
    spptrf_(JOB_STR+5, &n, A, &inf);
    return inf;
  }

  // packed triangular factorization of hermitian matrix
  inline int sptrf(const unsigned n, double* A, int* ipiv)
  {
    int inf;
    dsptrf_(JOB_STR+5, &n, A, ipiv, &inf);
    return inf;
  }
  inline int sptrf(const unsigned n, float* A, int* ipiv)
  {
    int inf;
    ssptrf_(JOB_STR+5, &n, A, ipiv, &inf);
    return inf;
  }

  // lower triangular solve
  inline void ltrs(const unsigned n, double* A,
		   const unsigned p, double* B, const unsigned ldB)
  {
    //  dtptrs_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      dtpsv_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
  }
  inline void ltrs(const unsigned n, float* A,
		   const unsigned p, float* B, const unsigned ldB)
  {
    // stptrs_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      stpsv_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
  }

  // lower triangular transpose solve
  inline void ltrhs(const unsigned n, double* A,
		    const unsigned p, double* B, const unsigned ldB)
  {
    //  dtptrs_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      dtpsv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
  }
  inline void ltrhs(const unsigned n, float* A,
		    const unsigned p, float* B, const unsigned ldB)
  {
    //  stptrs_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      stpsv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
  }

  // unit upper triangular solve (with L and R stored in one matrix)
  // XR=B, R is pxp, B is nxp
  inline void utrcs(const unsigned p, const double* LR, const unsigned ldLR,
		    const unsigned n, double* X, const unsigned ldX)
  {
    dtrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
	   LR, &ldLR, X, &ldX);
  }

  inline void utrcs(const unsigned p, const float* LR, const unsigned ldLR,
		    const unsigned n, float* X, const unsigned ldX)
  {
    strsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
	   LR, &ldLR, X, &ldX);
  }

  // unit upper triangular solve (with L and R stored in one matrix)
  // RX=B, R is nxn, B is nxp
  inline void utlcs(const unsigned n, float* LR, const unsigned ldLR,
		    const unsigned p, float* X, const unsigned ldX)
  {
    strsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
	   LR, &ldLR, X, &ldX);
  }
  inline void utlcs(const unsigned n, double* LR, const unsigned ldLR,
		    const unsigned p, double* X, const unsigned ldX)
  {
    dtrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
	   LR, &ldLR, X, &ldX);
  }

  // unit lower triangular solve (with L and R stored in one matrix)
  // XL=B, L is pxp, B is nxp
  inline void ltrcs(const unsigned p, float* LR, const unsigned ldLR,
		    const unsigned n, float* X, const unsigned ldX)
  {
    strsm_(JOB_STR+8, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
	   LR, &ldLR, X, &ldX);
  }
  inline void ltrcs(const unsigned p, double* LR, const unsigned ldLR,
		    const unsigned n, double* X, const unsigned ldX)
  {
    dtrsm_(JOB_STR+8, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
	   LR, &ldLR, X, &ldX);
  }


  // unit lower triangular transposed solve (with L and R stored in one matrix)
  // XL^T=B, L is pxp, B is nxp
  inline void ltrhcs(const unsigned p, double* LR, const unsigned ldLR,
		     const unsigned n, double* X, const unsigned ldX)
  {
    dtrsm_(JOB_STR+8, JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, &p, &D_ONE,
	   LR, &ldLR, X, &ldX);
  }

  // unit lower triangular solve (with L and R stored in one matrix)
  // LX=B, L is nxn, B is nxp
  inline void ltlcs(const unsigned n, float* LR, const unsigned ldLR,
		    const unsigned p, float* X, const unsigned ldX)
  {
    strsm_(JOB_STR+6, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
	   LR, &ldLR, X, &ldX);
  }
  inline void ltlcs(const unsigned n, double* LR, const unsigned ldLR,
		    const unsigned p, double* X, const unsigned ldX)
  {
    dtrsm_(JOB_STR+6, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
	   LR, &ldLR, X, &ldX);
  }

  // upper triangular solve
  inline void utrs(const unsigned n, double* A,
		   const unsigned p, double* B, const unsigned ldB)
  {
    //  dtptrs_(JOB_STR+5, JOB_STR, JOB_STR, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      dtpsv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, B+i*ldB, &N_ONE);
  }
  inline void utrs(const unsigned n, float* A,
		   const unsigned p, float* B, const unsigned ldB)
  {
    //stptrs_(JOB_STR+5, JOB_STR, JOB_STR, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      stpsv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, B+i*ldB, &N_ONE);
  }

  // upper triangluar transpose solve
  inline void utrhs(const unsigned n, double* A,
		    const unsigned p, double* B, const unsigned ldB)
  {
    //dtptrs_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      dtpsv_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, A, B+i*ldB, &N_ONE);
  }
  inline void utrhs(const unsigned n, float* A,
		    const unsigned p, float* B, const unsigned ldB)
  {
    //stptrs_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, &p, A, B, &ldB, &inf);
    for (unsigned i=0; i<p; ++i)
      stpsv_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, A, B+i*ldB, &N_ONE);
  }
}

#endif
