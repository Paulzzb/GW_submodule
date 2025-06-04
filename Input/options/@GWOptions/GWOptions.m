classdef GWOptions
    properties (SetAccess = public)
        ISDFCauchy
        Constant
        GWCal
        Groundstate
    end

    methods
        function obj = GWOptions()
            % Empty constructor with default structure
            obj.ISDFCauchy = struct();
            obj.Constant = struct();
            obj.GWCal = struct();
            obj.Groundstate = struct();
        end
    end
end
