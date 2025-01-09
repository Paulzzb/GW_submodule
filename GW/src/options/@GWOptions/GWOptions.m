classdef GWOptions
  properties (SetAccess = protected)

  end

  properties (SetAccess = private)
  
  end

  properties (SetAccess = public)
    ISDFCauchy;
    Constant;
    GWCal;
    Groundstate
  end

%   properties (Constant)
%     nv; nv_ener; nv_oper; nc; nc_ener; nc_oper; n_oper; n_ener; ...
%     ng; nr; ne; vol; Fouriercell; ...
%     TOL_ZERO; TOL_SMALL; INF; bdot; ...
%     ry2ev; ...
%     qk; 
%   
%   end
  
  methods
    function options_out = GWOptions(options_in, sys)
      % [options_in, sys] = init(mol, options_in);
      options_out.Constant = optionsConstant(options_in, sys);
      options_out = setGWCal(options_out, options_in);
      if (options_in.isISDF == true)
        options_out = setISDF(options_out, options_in, sys);
      else
        options_out.ISDFCauchy.isISDF = false;
      end
        
      options_out.Groundstate = struct();
      
      options_out.Groundstate.isGW = options_in.isGW;
      options_out.Groundstate.isBSE = options_in.isBSE;
      options_out.Groundstate.frequency_dependence ...
      = options_in.frequency_dependence;
      options_out.Groundstate.coulomb_truncation ...
      = options_in.coulomb_truncation; 
      options_out.Groundstate.nv = options_in.nv;
      options_out.Groundstate.nc = options_in.nc;
      options_out.Groundstate.amin = options_in.amin;
      
      options_out.Groundstate.input = options_in.input;
      if isfield(options_in, 'inputfile')
        options_out.Groundstate.inputfile = options_in.inputfile;
      end
      if isfield(options_in, 'options_kssolv')
        options_out.Groundstate.options_kssolv = options_in.options_kssolv;
      end

    return
    end % function GWOptions 

  end

end
