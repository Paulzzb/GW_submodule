function obj = fromConfig(config, sys)
% fromConfig - Construct GWOptions from parsed config struct

    obj = GWOptions();

    if isfield(config, 'CONTROL')
        obj.Constant = optionsConstant(config.CONTROL, sys);
    end

    if isfield(config, 'SYSTEM')
        obj = setGWCal(obj, config.SYSTEM);
    end

    if isfield(config, 'ISDF') && isfield(config.ISDF, 'isISDF') && config.ISDF.isISDF
        obj = setISDF(obj, config.ISDF, sys);
    else
        obj.ISDFCauchy.isISDF = false;
    end

    obj.Groundstate = struct();
    keys = {'isGW', 'isBSE', 'frequency_dependence', ...
            'coulomb_truncation', 'nv', 'nc', 'amin', 'input', 'inputfile'};

    for i = 1:length(keys)
        k = keys{i};
        if isfield(config.SYSTEM, k)
            obj.Groundstate.(k) = config.SYSTEM.(k);
        end
    end
end
