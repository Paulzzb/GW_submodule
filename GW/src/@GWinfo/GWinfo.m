classdef GWinfo 
  properties (SetAccess = public)
    coulG;
    coulG0;
    supercell;
    bdot;
    qk;
    vol
    ne 
    gvec
    gvec2
    gvecrho
    Vxc
    rho
    idxnz
    ev
    Z
    psir
    % delete below later
    ntot 
    nv 
    aqs
  end

  properties (SetAccess = protected)

  end
  
  methods
    function GWinfo = GWinfo(mol, options)
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
