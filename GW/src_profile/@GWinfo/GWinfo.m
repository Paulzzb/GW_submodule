classdef GWinfo 
  properties (SetAccess = public)
    coulG;
    coulG_list;
    coulG0;
    supercell;
    bdot;
    vol
    ne 
    gvec
    gvec_list
    gvecrho
    Vxc
    rho
    idxnz
    ev
    psig
    psir
    Ggrid4psig 
    occupation

    % nkpts
    % kpts
    % nspin
    % nspinor

  end

  properties (SetAccess = protected)

  end
  
  methods
    function GWinfo = GWinfo()
      % if nargin == 0
      %   return;
      % end
      GWinfo = GWinfo();
      % if nargin < 2
      %   error('Number of inputs less than 2');
      % end
      % GWinfo = gwsetup(GWinfo, mol, options);
    end
  end % method

end % classdef
