# Si_gamma QE save (not shipped in git)

Gamma-only Si needs a **thin** groundstate for the public suite (small
`wfc*.hdf5`, well under GitHub’s 100 MB limit). Until that dataset is ready,
copy your local QE `*.save` contents here:

- `data-file-schema.xml`
- `wfc*.hdf5`
- `charge-density.hdf5` (if used)
- `vxc.dat`
- pseudopotentials as needed

`example/run_all` will SKIP this case when the groundstate is missing.
`Si_k/qe.save` remains the small in-repo multi-k example.
