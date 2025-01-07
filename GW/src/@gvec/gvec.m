classdef gvec
  properties (SetAccess = public)
	  components
		index_vec
		ng
		nfftgridpts
		fftgrid
		idxnz
	end

	properties (SetAccess = protected)
    from
	end
	
	methods
	  function gvec_out = gvec(mol, varargin)
      if nargin == 0
				return;
			end
			nvar = length(varargin);
     % gvec = gvec();
      ecut = [];
			dir = [];
			method = 'kssolv';
     % if 
		  for it = 1:2:nvar
		    attr_name = varargin{it};
		    value = varargin{it+1};
        if mod(nvar, 2) == 1
	      	error('Wrong input for gvec.set');
	      end
				if strcmp(attr_name, 'dir')
					dir = value;
					continue;
				end
				if strcmp(attr_name, 'ecut')
					ecut = value;
				  continue;
				end
				if strcmp(attr_name, 'method')
          method = value;
				end
		    gvec_out.(attr_name) = value;
      end	
      
			if isempty(ecut)
				ecut = mol.ecut;
			end

      gvec_out.from = method;

			switch lower(method)
				case 'kssolv'
			    grid = Ggrid(mol, ecut);
	        ng = grid.ng;
	        coulG = [grid.gkx, grid.gky, grid.gkz] * mol.supercell / (2*pi);
	        coulG = round(coulG);
	        gvec_out.components = coulG;
	        gvec_out.index_vec = zeros(mol.n1 * mol.n2 * mol.n3, 1);
	        gvec_out.index_vec(grid.idxnz) = (1:ng)';
	        gvec_out.ng = ng;
	        gvec_out.nfftgridpts = prod([mol.n1, mol.n2, mol.n3]);
	        gvec_out.fftgrid = [mol.n1, mol.n2, mol.n3];
					gvec_out.idxnz = grid.idxnz;
          % gvec_out = finalize_kssolv(mol, ecut, gvec_out);
				case 'bgw'
				  % gvec_out = finalize_bgw(mol, dir, gvec_out);
				otherwise
				  error('Method for gvec is not supported!')
			end
    end % function gvec
	end % method

end % classdef
