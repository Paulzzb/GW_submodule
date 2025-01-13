function gvec_out = set(mol, varargin)
  
  nvar = length(varargin);
  if mod(nvar, 2) == 1
    error('Wrong input for gvec.set');
  end
  
  for it = 1:2:nvar
    attr_name = varargin{it};
    value = varargin{it+1};
    
    gvec_out.(attr_name) = value;
  end

end
