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
    fout;

    kpt         (:, 3)        = [0,0,0]
    kptweights  (:, 1)        = 1

  end

  methods
    function obj = QPenergy(GWinfo, config)
      defcon = default_constant();
      deffile = filename_map();
      ry2ev = defcon.ry2ev;
      
      obj.freq_dep = config.FREQUENCY.frequency_dependence;
      obj.freq_dep_method = config.FREQUENCY.frequency_dependence_method;
      obj.fout = deffile.QPenergy;

      obj.bandindices = config.SYSTEM.energy_band_index_min:config.SYSTEM.energy_band_index_max;
      obj.ev = GWinfo.ev(obj.bandindices) * ry2ev;
      obj.Vxc = GWinfo.Vxc(obj.bandindices) * ry2ev;
    end % function 
  end % methods
end % classdef