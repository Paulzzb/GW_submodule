test_nok — minimal layout like test_Si for qp_cohsex / gw_cohsex_multi_k smoke test
==================================================================================

Run (MATLAB), from this folder:
  test_input          % optional: (re)build ./SAVE from ./test via input_driver
  test_gw_cohsex_multi_k

The bundled Si.save/ here may be incomplete for QE load (e.g. missing wfc*.hdf5).
If input_driver fails or GW stages error out, copy from a working profile:

  xcopy /E /I ..\test_Si\SAVE .\SAVE
  xcopy /E /I ..\test_Si\Si.save .\Si.save

(Linux/mac: cp -r ../test_Si/SAVE . && cp -r ../test_Si/Si.save .)

Then run test_gw_cohsex_multi_k again.
