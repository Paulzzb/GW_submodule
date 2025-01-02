classdef ksinfo 
  properties (SetAccess = public)
		F ;
    mol ;
    coulG;
    coulG0
    bdot
    qk 
    ntot 
    vol
    ne 
    nv 
    gvec
    gvec2
		gvecrho
		Vxc
    rho
		idxnz
		ev
		Z
		aqs
	end

	properties (SetAccess = protected)

	end
	
	methods
	  function ksinfo = ksinfo(mol, options)
      if nargin == 0
				return;
			end
			ksinfo = ksinfo();
			if nargin < 2
				error('Number of inputs less than 2');
			end
      ksinfo = gwsetup(ksinfo, mol, options);
		  % ksinfo = setksinfo(ksinfo, options);
	  end
  end % method

end % classdef
