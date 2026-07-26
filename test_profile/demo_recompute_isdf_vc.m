function id_new = demo_recompute_isdf_vc(case_dir, output_dir, inputfile)
%DEMO_RECOMPUTE_ISDF_VC  Wrapper: demo_recompute_isdf(..., 'vc', ...).
  id_new = demo_recompute_isdf(case_dir, 'vc', output_dir, inputfile);
end
