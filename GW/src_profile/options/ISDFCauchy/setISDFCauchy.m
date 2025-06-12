function opt = setISDFCauchy(data, config)

def = default_ISDFCauchy();

opt = def;

% Read data from config and data
opt.isISDF = config.ISDF.isisdf;
opt.exxmethod = config.ISDF.exxmethod;
opt.vcrank_ratio = config.ISDF.isdf_ratio_type1;
opt.vsrank_ratio = config.ISDF.isdf_ratio_type2;
opt.ssrank_ratio = config.ISDF.isdf_ratio_type3;
opt.isCauchy = 0;
% opt.froErr = config.
% opt.MaxIter = config.

opt.isdfoptions.weight = config.ISDF.weight;
opt.isdfoptions.seed = config.ISDF.seed;
opt.isdfoptions.init = config.ISDF.init;
opt.isdfoptions.sys = data.sys;
% opt.isdfoptions.rank = data.sys;

end % EOF