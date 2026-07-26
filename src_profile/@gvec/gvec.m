classdef gvec
  % Class contains information about the grid vectors in reciprocal space.
  % 
  % Normally, the user doesn't need to know the above this class,
  % since they should not manipulate the grid vectors directly.
  % 
  % Theorically, the grid vectors in reciprocal space is n1 x n2 x n3,
  % the same as the grid in real space. After cutting off the grids with
  %  a given 'ecut', it remains 'ng' grids.
  % The above process generates a map:
  %    1:ng -> 1:n1*n2*n3 -> (1:n1) x (1:n2) x (1:n3)
  %    ─(multiplies 2*pi*inv(supercell)')─> reciprocal grids
  % 
  % @gvec class contains the following properties:
  %     ng        : number of grid vectors
  %     fftgrid   : exact [n1, n2, n3]
  %     nfftgridpts: n1*n2*n3
  %     idxnz     : the mapping of 1:ng -> 1:n1*n2*n3
  %     components: the mapping of 1:ng -> (1:n1) x (1:n2) x (1:n3)
  % 
  % To generate a @gvec
  %         gvec_out = gvec(gvecinput),
  % Please contains the following information in gvecinput:
  %     n1(int, essential), n2(int, essential), n3(int, essential),
  %     ecut(double, esential), 
  %     supercell(3*3, double, essential),
  %     Ggrid_forming_method(string, 'kssolv' or 'qe'),
  % Or an alternative way is to provide ng * 3 double array, containing
  % the grid vectors in (1:n1) x (1:n2) x (1:n3)

  % Zhengbang Zhou, Fudan University, 01-13-2025

  % properties (SetAccess = protected)
  properties (SetAccess = public)
    components
    ng
    nfftgridpts
    fftgrid
    idxnz
    
    % with respecting to gmapping :
    % gvec.component(ii, :) - Go_list(jj, :) = fullfftgrid(G0_mappling(ii, jj))
    nGo          (1,1) {mustBeInteger}% number of Go vector
    Go_list      (:,:)                % [nGo x 3] G0 vectors
    G0_mapping   (:,:) {mustBeInteger}% [nGo x 1]

    % g_rot   (:, :) {mustBeInteger}       % [ng x nrot]
    g_rot                             % cell [nrot, 1]
  end

  
  methods
    function gvec_out = gvec(gvecinput)
      if nargin == 0
        return;
      end

      requiredFields = {'n1', 'n2', 'n3', 'supercell'};
      for i = 1:numel(requiredFields)
        fieldName = requiredFields{i};
        if ~isfield(gvecinput, fieldName)
          error('%s is not specified', fieldName);
        end
      end

      n1 = gvecinput.n1;
      n2 = gvecinput.n2;
      n3 = gvecinput.n3;
      gvec_out.fftgrid = [n1, n2, n3];
      gvec_out.nfftgridpts = n1 * n2 * n3;

      
      if ~isfield(gvecinput, 'qpoint') 
        gvecinput.qpoint = [0, 0, 0]; % Then it is suitable for original non-kpoints code.
      end

      if isfield(gvecinput, 'ecut') 
        flagecut= 1;
      else
        flagecut = 0;
      end

      if isfield(gvecinput, 'Ggrid') % Provide Ggrid directly
        flagGgrid = 1;
      else
        flagGgrid = 0;
      end

      if (~flagecut && ~flagGgrid)
        error('Please provide either ecut or Ggrid directly');
      end

      if flagGgrid
        gvec_out.components = gvecinput.Ggrid;
        ng = size(gvecinput.Ggrid, 1);
        gvec_out.ng = ng;
        idxnz = findidxnz(gvec_out.components, [n1, n2, n3]); 
        gvec_out.idxnz = idxnz;
      else
        % Generate it based on ecut and Ggrid_forming_method
        gvec_out = mill(gvec_out, gvecinput);
        ng = gvec_out.ng;
      end      

    end % function gvec
  end % method

  % methods (Static, Access = private)
  %   function idxnz = findidxnz(grid, n1n2n3)
  %   % Find the mapping from 3-D grid indices to 1-D indices,
  %   % Use 'our mapping rules'
  %   n1 = n1n2n3(1);
  %   n2 = n1n2n3(2);
  %   n3 = n1n2n3(3);
    
  %   idxnz = 1 + (grid(:, 1) + (grid(:, 1) < 0) * n1) ...
  %           + (grid(2) + (grid(2) < 0) * n2) * n1 ... 
  %           + (grid(3) + (grid(3) < 0) * n3) * n1 * n2;
    
  %   end % EOF
  % end

end % classdef
