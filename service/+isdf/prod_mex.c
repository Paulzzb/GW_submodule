#include "mex.h"
#include "blas.h"
#include <stddef.h>

static int all_real(const mxArray *a, const mxArray *b, const mxArray *c, const mxArray *d) {
  return !mxIsComplex(a) && !mxIsComplex(b) && !mxIsComplex(c) && !mxIsComplex(d);
}

static mxComplexDouble *to_complex_buffer(const mxArray *a, mwSize n, int *must_free) {
  mwSize i;
  if (mxIsComplex(a)) {
    *must_free = 0;
    return mxGetComplexDoubles(a);
  }
  *must_free = 1;
  {
    const double *r = mxGetDoubles(a);
    mxComplexDouble *buf = (mxComplexDouble *)mxCalloc(n, sizeof(mxComplexDouble));
    for (i = 0; i < n; ++i) {
      buf[i].real = r[i];
      buf[i].imag = 0.0;
    }
    return buf;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mxArray *Psi, *psi, *Phi, *phi;
  mwSize m, n, k1, k2, i;
  ptrdiff_t m_, n_, k1_, k2_, lda, ldb, ldc;
  char trans_n = 'N', trans_c = 'C';
  double alpha = 1.0;
  double beta = 0.0;
  double alpha_c[2] = {1.0, 0.0};
  double beta_c[2] = {0.0, 0.0};

  if (nrhs != 4) {
    mexErrMsgIdAndTxt("isdf:prod_mex:nrhs", "Expected 4 inputs: Psi, psi, Phi, phi.");
  }
  if (nlhs > 1) {
    mexErrMsgIdAndTxt("isdf:prod_mex:nlhs", "Expected one output.");
  }

  Psi = prhs[0];
  psi = prhs[1];
  Phi = prhs[2];
  phi = prhs[3];

  if (!mxIsDouble(Psi) || !mxIsDouble(psi) || !mxIsDouble(Phi) || !mxIsDouble(phi)) {
    mexErrMsgIdAndTxt("isdf:prod_mex:type", "prod_mex supports double real/complex inputs only.");
  }

  m = mxGetM(Psi);
  k1 = mxGetN(Psi);
  n = mxGetM(psi);
  if (mxGetN(psi) != k1) {
    mexErrMsgIdAndTxt("isdf:prod_mex:size", "size(Psi,2) must equal size(psi,2).");
  }
  if (mxGetM(Phi) != m) {
    mexErrMsgIdAndTxt("isdf:prod_mex:size", "size(Psi,1) must equal size(Phi,1).");
  }
  if (mxGetM(phi) != n) {
    mexErrMsgIdAndTxt("isdf:prod_mex:size", "size(psi,1) must equal size(phi,1).");
  }
  k2 = mxGetN(Phi);
  if (mxGetN(phi) != k2) {
    mexErrMsgIdAndTxt("isdf:prod_mex:size", "size(Phi,2) must equal size(phi,2).");
  }

  m_ = (ptrdiff_t)m;
  n_ = (ptrdiff_t)n;
  k1_ = (ptrdiff_t)k1;
  k2_ = (ptrdiff_t)k2;
  lda = (ptrdiff_t)m;
  ldb = (ptrdiff_t)n;
  ldc = (ptrdiff_t)m;

  if (all_real(Psi, psi, Phi, phi)) {
    double *M1 = (double *)mxCalloc(m * n, sizeof(double));
    double *M2 = (double *)mxCalloc(m * n, sizeof(double));
    plhs[0] = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);
    {
      double *out = mxGetDoubles(plhs[0]);
      dgemm(&trans_n, &trans_c, &m_, &n_, &k1_, &alpha,
            mxGetDoubles(Psi), &lda, mxGetDoubles(psi), &ldb, &beta, M1, &ldc);
      dgemm(&trans_n, &trans_c, &m_, &n_, &k2_, &alpha,
            mxGetDoubles(Phi), &lda, mxGetDoubles(phi), &ldb, &beta, M2, &ldc);
      for (i = 0; i < m * n; ++i) {
        out[i] = M1[i] * M2[i];
      }
    }
    mxFree(M1);
    mxFree(M2);
    return;
  }

  {
    mxComplexDouble *M1, *M2;
    mxComplexDouble *PsiC, *psiC, *PhiC, *phiC;
    int free_PsiC = 0, free_psiC = 0, free_PhiC = 0, free_phiC = 0;

    PsiC = to_complex_buffer(Psi, m * k1, &free_PsiC);
    psiC = to_complex_buffer(psi, n * k1, &free_psiC);
    PhiC = to_complex_buffer(Phi, m * k2, &free_PhiC);
    phiC = to_complex_buffer(phi, n * k2, &free_phiC);

    plhs[0] = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxCOMPLEX);
    M1 = mxGetComplexDoubles(plhs[0]);
    M2 = (mxComplexDouble *)mxCalloc(m * n, sizeof(mxComplexDouble));

    zgemm(&trans_n, &trans_c, &m_, &n_, &k1_, alpha_c, (double *)PsiC, &lda, (double *)psiC, &ldb, beta_c, (double *)M1, &ldc);
    zgemm(&trans_n, &trans_c, &m_, &n_, &k2_, alpha_c, (double *)PhiC, &lda, (double *)phiC, &ldb, beta_c, (double *)M2, &ldc);

    for (i = 0; i < m * n; ++i) {
      double ar = M1[i].real;
      double ai = -M1[i].imag;
      double br = M2[i].real;
      double bi = M2[i].imag;
      M1[i].real = ar * br - ai * bi;
      M1[i].imag = ar * bi + ai * br;
    }

    mxFree(M2);
    if (free_PsiC) mxFree(PsiC);
    if (free_psiC) mxFree(psiC);
    if (free_PhiC) mxFree(PhiC);
    if (free_phiC) mxFree(phiC);
  }
}
