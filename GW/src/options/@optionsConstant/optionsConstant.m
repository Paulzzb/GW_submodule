classdef optionsConstant
  properties (SetAccess = protected)

    nv; nv_ener; nv_oper; nc; nc_ener; nc_oper; n_oper; n_ener; ...
    ng; nr; ne; vol; ...
    TOL_ZERO; TOL_SMALL; INF; ...
    ry2ev;  
    ZERO; CZERO; 
  
    % For gpp 
    limitone; limittwo; sexcut;
  end
 
  methods
    function Constants = optionsConstant(options_in, sys)
      Constants.nv = options_in.nv;
      Constants.nc = options_in.nc;
      Constants.nv_oper = options_in.nv_oper;
      Constants.nc_oper = options_in.nc_oper;
      Constants.nv_ener = options_in.nv_ener;
      Constants.nc_ener = options_in.nc_ener;
      Constants.n_oper = options_in.n_oper;
      Constants.n_ener = options_in.n_ener;

      Constants.ng = sys.ng;
      Constants.nr = sys.nr;
      Constants.ne = sys.ne;
      Constants.vol = sys.vol;

      Constants.TOL_ZERO = 1e-12;
      Constants.TOL_SMALL = 1e-6;
      Constants.INF = 1e+12;
      Constants.ry2ev = 13.60569253;
      Constants.CZERO = complex(0.0, 0.0);
      Constants.ZERO = 0.0;

      if options_in.frequency_dependence == 1
        if isfield(options_in, 'gpp_brodening')
          limittwo = gpp_brodening.^2;
        else
          limittwo = 0.25;
        end
        if isfield(options_in, 'gpp_sexcutoff')
          sexcut = options_in.gpp_sexcutoff;
        else
          sexcut = 4.0;
        end
        Constants.limitone = 0.25 * 1e+6;
        Constants.limittwo = limittwo;
        Constants.sexcut = sexcut;
      end
      if options_in.frequency_dependence == 2
        ;
      end
    end
  end

end 
