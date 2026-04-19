function varargout = isdf_schur_update(mode, varargin)
%ISDF_SCHUR_UPDATE  Persistent state for Schur / Gram refresh.
%   CCH and invCCH are both Nisdf-by-Nisdf (centroid Gram and its inverse).
%
%   isdf_schur_update('init', Nisdf)
%   isdf_schur_update('init', Nisdf, fftgrid_fine)
%   isdf_schur_update('init', Nisdf, fftgrid_fine, Nremain)
%   [CCH, invCCH, Nremain, fftgrid_fine] = isdf_schur_update('get')
%   isdf_schur_update('clear')

  persistent Nisdf Nadd CCH invL_CCH L_CCH Psi_on_grid Phi_on_grid ...
             Nisdfmax

  if nargin < 1
    error('isdf_schur_update:mode', 'First argument ''mode'' is required.');
  end

  switch lower(mode)
    case 'update'
      % varargin{1}: Nadd
      % varargin{2}: indices for new isdf points, in fftgrid_fine
      % 
      Nadd = round(double(varargin{1}));
      if (Nadd - varargin{1}) > 1e-6
        error('isdf_schur_update:update', 'Nadd is not an integer.');
      end
      Psi_on_new_grid = double(varargin{2}); % size: Nadd, nb1, nspin
      Phi_on_new_grid = double(varargin{3}); % size: Nadd, nb2, nspin
      % Update process
      C2C2H = isdf_prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);
      % C1C2H = isdf_prod(Psi_on_grid(1:Nisdf, :), Psi_on_new_grid, Phi_on_grid(1:Nisdf, :), Phi_on_new_grid); 
      C2C1H = isdf_prod(Psi_on_new_grid, Psi_on_grid(1:Nisdf, :), Phi_on_new_grid, Phi_on_grid(1:Nisdf, :)); 
      L_CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf) = C2C1H / L_CCH(1:Nisdf, 1:Nisdf)';
      S = C2C2H - L_CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf)*L_CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf)';
      L_CCH(Nisdf+1:Nisdf+Nadd, Nisdf+1:Nisdf+Nadd) = chol(S, "lower");
      % Seems like I don't need CCH ...
      CCH(Nisdf+1:Nisdf+Nadd, Nisdf+1:Nisdf+Nadd) = C2C2H;
      CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf) ...
      = isdf_prod(Psi_on_grid(1:Nisdf, :), Psi_on_new_grid, Phi_on_grid(1:Nisdf, :), Phi_on_new_grid); 
      CCH(1:Nisdf, Nisdf+1:Nisdf+Nadd) = CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf)';
      invL_CCH(Nisdf+1:Nisdf+Nadd, Nisdf+1:Nisdf+Nadd) = inv(L_CCH(Nisdf+1:Nisdf+Nadd, Nisdf+1:Nisdf+Nadd));
      invL_CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf) = ...
      - invL_CCH(Nisdf+1:Nisdf+Nadd, Nisdf+1:Nisdf+Nadd) * L_CCH(Nisdf+1:Nisdf+Nadd, 1:Nisdf) * invL_CCH(1:Nisdf, 1:Nisdf);
      ...
      Psi_on_grid(Nisdf+1:Nisdf+Nadd, :) = Psi_on_new_grid;
      Phi_on_grid(Nisdf+1:Nisdf+Nadd, :) = Phi_on_new_grid;
      Nisdf = Nisdf + Nadd;

    case 'init'
      % varargin{1}: Nisdf
      if nargin < 2 || isempty(varargin{1})
        error('isdf_schur_update:init', '''init'' requires scalar or nonempty Nisdf.');
      end
      Nisdf = double(varargin{1}(1));
      Nisdfmax = Nisdf + 200; % TODO: make it a parameter
      %
      nb1 = size(varargin{2}, 2);
      nb2 = size(varargin{3}, 2);
      Psi_on_grid = zeros(Nisdfmax, nb1);
      Phi_on_grid = zeros(Nisdfmax, nb2);
      Psi_on_grid(1:Nisdf, :) = double(varargin{2}); % size: Nisdf, nb1
      Phi_on_grid(1:Nisdf, :) = double(varargin{3}); % size: Nisdf, nb2
      
      CCH = zeros(Nisdfmax, Nisdfmax);
      L_CCH = zeros(Nisdfmax, Nisdfmax);
      invL_CCH = zeros(Nisdfmax, Nisdfmax);
      CCH(1:Nisdf, 1:Nisdf) = isdf_prod(Psi_on_grid(1:Nisdf, :), Psi_on_grid(1:Nisdf, :), ...
                      Phi_on_grid(1:Nisdf, :), Phi_on_grid(1:Nisdf, :));
      if condest(CCH(1:Nisdf, 1:Nisdf)) < 1e+12
        L_CCH(1:Nisdf, 1:Nisdf) = chol(CCH(1:Nisdf, 1:Nisdf), "lower");
        invL_CCH(1:Nisdf, 1:Nisdf) = inv(L_CCH(1:Nisdf, 1:Nisdf));
      else
        warning('isdf_schur_update:init', 'CCH is ill-conditioned, using pseudoinverse.');
        L_CCH(1:Nisdf, 1:Nisdf) = chol(CCH(1:Nisdf, 1:Nisdf), "lower");
        invL_CCH(1:Nisdf, 1:Nisdf) = inv(L_CCH(1:Nisdf, 1:Nisdf));
      end
      Nadd = Nisdf;
    case 'get'
      varargout{1} = Nisdf;
      varargout{2} = Nadd; 
      varargout{3} = invL_CCH;
      if nargout >= 2
        varargout{4} = L_CCH;
      end
      if nargout >= 5
        varargout{5} = Psi_on_grid;
      end
      if nargout >= 6
        varargout{6} = Phi_on_grid;
      end
      if nargout >= 7
        error('isdf_schur_update:get', 'too much output arguments.');
      end
    case 'clear'
      clear Nisdf Nadd CCH invL_CCH L_CCH Psi_on_grid Phi_on_grid Nisdfmax

    otherwise
      error('isdf_schur_update:mode', 'Unknown mode ''%s''.', mode);
  end
end
