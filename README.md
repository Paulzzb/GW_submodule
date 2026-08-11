# A MATLAB toolbox for low-rank planewave GW calculation

A MATLAB toolbox for **low-rank** planewave **GW** calculations (COHSEX and full-frequency paths), with ISDF-based compression.

It originated as a module inside [KSSOLV](https://doi.org/10.1021/acs.jpca.1c03762) and is now a **standalone** codebase.

---

## Status (beta)

This is a **beta / test** release. Many parts are still under active development and may change without notice. Interfaces and defaults are not fully stabilized yet — sorry about that.

---

## Quick start

From the repository root in MATLAB:

```matlab
QPstartup
```

Then see the small demos under [`examples/`](examples/README.md), or the user guides:

- English: [`doc_user_EN/README.md`](doc_user_EN/README.md)
- 中文: [`doc_user_ZH/README.md`](doc_user_ZH/README.md)

---

## Paper reproduction (Benchmarks)

To reproduce systems from the paper, see [`benchmarks/`](benchmarks/README.md).

That tree ships **Quantum ESPRESSO inputs only** (no `*.save` / wavefunctions). Run QE yourself, then point this code at the resulting save directory. For a tiny end-to-end demo on a laptop, prefer `examples/` instead.

---

## Citation

If you use this code, please cite:

- Z. Zhou, H. Ma, W. Wu, W. Gao, J. Yang, M. Shao, and W. Hu,
  “A fast low-rank inversion algorithm of dielectric matrix in GW approximation,”
  [arXiv:2403.12340](https://arxiv.org/abs/2403.12340) (2024).

---

## Contact

Questions and bug reports: **zbzhou21@m.fudan.edu.cn**

---

## Third-party elliptic routines

Cauchy ISDF uses Jacobi / complete elliptic integrals, which are implemented in
**Tobin A. Driscoll’s Schwarz–Christoffel Toolbox** (`ellipkkp`, `ellipjc`; BSD-3-Clause).
These are **not** original to this project. See
[`service/+isdf/+Cauchy/NOTICE`](service/+isdf/+Cauchy/NOTICE).

Where only the real complete integral \(K(k)\) is needed, MATLAB’s built-in
`ellipke(k^2)` is numerically equivalent; this tree uses the vendored
`ellipkkp` for consistency with `ellipjc`.

---

## License

BSD 3-Clause. See [`LICENSE`](LICENSE). Contributors: see [`AUTHORS`](AUTHORS).
Third-party notices: [`service/+isdf/+Cauchy/NOTICE`](service/+isdf/+Cauchy/NOTICE).
