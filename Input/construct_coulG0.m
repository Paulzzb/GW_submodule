function coulG0 = construct_coulG0(data, config);
% Depending on config, construct the v(q+G) when |q+G| ~ 0
trun_method = config.CUTOFFS.coulomb_truncation_method;
trun_factor = config.CUTOFFS.coulomb_truncation_parameter;

fourpi = 4*pi;
eightpi = 8*pi;
switch trun_method 
  case 0
    warning('GW:construct_vcoul0:trunc_method', 'None truncation not implemented yet');
    warning('GW:construct_vcoul0:trunc_method', 'Set default value to 0.0');
    coulG0 = 0.0;
  case 2
    coulG0 = fourpi * trunc_factor^2;
  otherwise
    msg = sprintf('Unsupported truncation type: %d', trun_method);
    GWerror(msg);
end

  
end
