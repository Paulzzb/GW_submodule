classdef QPenergy

  properties
    bandindices
    freq_dep
    freq_dep_method
    ev
    Ex
    Esex_x
    Ecoh
    Eres
    Eint
    Sig
    Vxc
    Eqp
    
    TOL_DEGENERACY = 1e-6 * 13.60569253;
    fout = 'GWenergy.out';
  end

  methods
    function obj = QPenergy(GWinfo, GWoptions)
      ry2ev = GWoptions.Constant.ry2ev;
      nv = GWoptions.Constant.nv;
      nv_ener = GWoptions.Constant.nv_ener;
      nc_ener = GWoptions.Constant.nc_ener;
      
      obj.freq_dep = GWoptions.GWCal.freq_dep;
      obj.freq_dep_method = GWoptions.GWCal.freq_dep_method;
      obj.fout = GWoptions.GWCal.fileName;

      obj.bandindices = nv-nv_ener+1:nv+nc_ener;
      obj.ev = GWinfo.ev(obj.bandindices) * ry2ev;
      obj.Vxc = GWinfo.Vxc(obj.bandindices) * ry2ev;
      % check degeneracy
    end % function GWenergy

    function check_degeneracy(obj, ev)
      % Check degeneracy according to ev
      n = size(ev, 1);
      if (n > 1)
        e1 = obj.ev(1:end-1, :);
        e2 = obj.ev(2:end, :);
      end
    end % function degeneracy


  end % methods
end % classdef