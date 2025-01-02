function opt = setISDF(opt, mol, options_in)

  if ~isfield(options_in, 'optionsISDF')
    optionsISDF = [];
  else
  %   opt.ISDFCauchy.isdfoptions.seed = 0;
  %   opt.ISDFCauchy.exxmethod = 'kmeans'; 
	%   opt.ISDFCauchy.isdfoptions.weight = 'add'; 
	%   opt.ISDFCauchy.isdfoptions.sys = mol;
	%   opt.ISDFCauchy.isdfoptions.init = 'wrs';
  % %  opt.ISDFCauchy.froErr = 1e-8; optionsISDFCauchy.MaxIter = 10;
  % else
    optionsISDF = options_in.optionsISDF
  end
	
  optionsISDF.isISDF = true;
 
  if ~isfield(optionsISDF, 'exxmethod')
    optionsISDF.exxmethod = 'qrcp';
  end
  if ~isfield(optionsISDF, 'isdfoptions')
    optionsISDF.isdfoptions = [];
  end
  if ~isfield(optionsISDF.isdfoptions, 'seed')
    optionsISDF.isdfoptions.seed = 0;
  end
  if ~isfield(optionsISDF.isdfoptions, 'init')
    optionsISDF.isdfoptions.init = 'wrs';
  end
  if ~isfield(optionsISDF.isdfoptions, 'weight')
    optionsISDF.isdfoptions.weight = 'add';
  end
  if ~isfield(optionsISDF.isdfoptions, 'sys')
    optionsISDF.isdfoptions.sys = mol;
  end
  opt.ISDFCauchy = optionsISDF;
  % set informations of ISDF	
  
	if ~isfield(options_in, 'vcrank_ratio')
  	opt.ISDFCauchy.vcrank_ratio = 8;
  else
  	opt.ISDFCauchy.vcrank_ratio = options_in.vcrank_ratio;
  end
  
	if ~isfield(options_in, 'vsrank_ratio')
  	opt.ISDFCauchy.vsrank_ratio = 8;
  else
  	opt.ISDFCauchy.vsrank_ratio = options_in.vsrank_ratio;
  end
  
	if ~isfield(options_in, 'ssrank_ratio')
  	opt.ISDFCauchy.ssrank_ratio = 8;
  else
  	opt.ISDFCauchy.ssrank_ratio = options_in.ssrank_ratio;
  end
  
	% set information of Cauchy
	if ~isfield(options_in, 'isCauchy')
	  opt.ISDFCauchy.isCauchy = true;
	else
    opt.ISDFCauchy.isCauchy = options_in.isCauchy;
  end
  
	if opt.ISDFCauchy.isCauchy
    if isfield(options_in, 'optionsCauchy')
      optionsCauchy = options.optionsCauchy;
    else
      optionsCauchy = [];
    if ~isfield(optionsCauchy, 'froErr')
      optionsCauchy.froErr = 1e-6; 
    end
    if ~isfield(optionsCauchy, 'MaxIter')
			optionsCauchy.MaxIter = 10;
    end
    opt.ISDFCauchy.optionsCauchy = optionsCauchy;
  end  

end
