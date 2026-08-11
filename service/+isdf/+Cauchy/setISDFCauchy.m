function opt = setISDFCauchy(data, config)

def = isdf.Cauchy.default_ISDFCauchy();

opt = def;

% Read data from config and data
opt.isISDF = config.ISDF.isisdf;
opt.is_helper = config.ISDF.is_helper;
opt.exxmethod = config.ISDF.exxmethod;
opt.vcrank_ratio = config.ISDF.isdf_ratio_type1;
opt.vsrank_ratio = config.ISDF.isdf_ratio_type2;
opt.ssrank_ratio = config.ISDF.isdf_ratio_type3;

% Switch lives on &ISDF iscauchy (defaults froErr / MaxIter).
copt = isdf.cauchy_opts(config);
opt.isCauchy = copt.isCauchy;
opt.froErr = copt.froErr;
opt.MaxIter = copt.MaxIter;

opt.isdfoptions.weight = config.ISDF.weight;
opt.isdfoptions.seed = config.ISDF.seed;
opt.isdfoptions.init = config.ISDF.init;
opt.isdfoptions.sys = data.sys;

end % EOF
