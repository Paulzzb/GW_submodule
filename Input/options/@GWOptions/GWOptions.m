classdef GWOptions
    properties (SetAccess = public)
        ISDFCauchy
        Constant
        GWCal
        Groundstate
    end

    methods
        function obj = GWOptions(data, config)
            % Empty constructor with default structure
            obj.ISDFCauchy = setISDFCauchy(data, config);
            obj.Constant = setConstant(data, config);
            obj.GWCal = setGWCal(data, config);
            obj.Groundstate = struct();
        end
    end
end
