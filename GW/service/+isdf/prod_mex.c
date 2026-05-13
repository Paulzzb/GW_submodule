#include "mex.h"
#include "blas.h"
#include <stddef.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mxArray *Psi, *psi, *Phi, *phi;
  mwSize m, n, k1, k2, i;
  ptrdiff_t m_, n_, k1_, k2_, lda, ldb, ldc;
  char trans_n = 'N', trans_c = 'C';
  float alpha_c[2] = {1.0f, 0.0f};
  float beta_c[2] = {0.0f, 0.0f};
  mxComplexSingle *M1, *M2;

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

  if (!mxIsSingle(Psi) || !mxIsSingle(psi) || !mxIsSingle(Phi) || !mxIsSingle(phi)) {
    mexErrMsgIdAndTxt("isdf:prod_mex:type", "prod_mex supports single only.");
  }
  if (!mxIsComplex(Psi) || !mxIsComplex(psi) || !mxIsComplex(Phi) || !mxIsComplex(phi)) {
    mexErrMsgIdAndTxt("isdf:prod_mex:type", "prod_mex supports single complex only.");
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

  plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxCOMPLEX);
  M1 = mxGetComplexSingles(plhs[0]);
  M2 = (mxComplexSingle *)mxCalloc(m * n, sizeof(mxComplexSingle));

  /* M1_tmp = Psi * psi^H */
  cgemm(&trans_n, &trans_c, &m_, &n_, &k1_, alpha_c, (float *)mxGetComplexSingles(Psi), &lda, (float *)mxGetComplexSingles(psi), &ldb, beta_c, (float *)M1, &ldc);
  /* M2 = Phi * phi^H */
  cgemm(&trans_n, &trans_c, &m_, &n_, &k2_, alpha_c, (float *)mxGetComplexSingles(Phi), &lda, (float *)mxGetComplexSingles(phi), &ldb, beta_c, (float *)M2, &ldc);

  /* out = conj(M1_tmp) .* M2 */
  for (i = 0; i < m * n; ++i) {
    float ar = M1[i].real;
    float ai = -M1[i].imag;
    float br = M2[i].real;
    float bi = M2[i].imag;
    M1[i].real = ar * br - ai * bi;
    M1[i].imag = ar * bi + ai * br;
  }

  mxFree(M2);
}
