# benchmarks/ — paper systems (QE inputs only)

No `*.save` / wavefunctions are shipped. Run Quantum ESPRESSO yourself, then point
this GW code at the resulting save directory.

## Inventory

| Dir | System | Atoms | Cell (bohr) | UPF in-tree | Notes |
|-----|--------|------:|-------------|-------------|-------|
| `Si8/` | Si conventional | 8 | 10.3344 | via `pseudo_dir='../Si64'` | Gamma; nscf `nbnd=2048` |
| `Si64/` | Si 2×2×2 of Si8 | 64 | 20.6688 | `Si_ONCV_PBE-1.0.upf` | Gamma; nscf `nbnd=16384`; `ecutwfc=50` |
| `STO3/` | SrTiO₃ primitive | 5 | 7.455 | via `../STO3_8` | prefix `STO2`; nscf `nbnd=1800`; `ecutwfc=40` |
| `STO3_8/` | SrTiO₃ 2×2×2 | 40 | 14.91 | Sr/Ti/O ONCV PBE | prefix `STO2`; nscf `nbnd=24000`; `ecutwfc=80` |

Each case typically has: `scf.in`, `nscf.in`, `pp_in` (pw2bgw for `vxc.dat` / optional WFN).

## Suggested QE workflow (per case)

```bash
cd benchmarks/<case>
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out
pw2bgw.x < pp_in > pp_out
# then either:
#   - use QE *.save directly with this GW code (groundstate_type='qe'), or
#   - run pw2bgw.x < pp_in if you need BGW-style WFN + vxc.dat
```

## Still missing / known gaps

- Per-case README with QE version, walltime hints, and expected `outdir` layout for GW
- MATLAB / GW `test` namelist stubs under each case (how to set `groundstate_dir`)
- `Si8/` has no local UPF (depends on `../Si64`)
- Folder name `STO3` vs QE `prefix='STO2'` is historical; do not rename prefix without updating `pp_in`
- `Si8` uses `ecutwfc=40`, `Si64` uses `50` — intentional only if paper says so; otherwise align
- Huge `nbnd` values (especially Si64 / STO3_8) need large memory; document machine requirements
