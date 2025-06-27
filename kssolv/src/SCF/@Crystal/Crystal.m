classdef Crystal < Molecule
    % CRYSTAL KSSOLV class for crystal
    %    cry = CRYSTAL(str1,field1,str2,field2,...) returns a crystal
    %    class of the given fields with respect to the name strings.
    %
    %    The crystal class contains the following fields.
    %        Field     Explaination
    %      ----------------------------------------------------------
    %        kpts        The array of k-points
    %        nkpts       The number of k-points
    %        wks         The weights for k-points
    %        scfkpts     Used in hybrid functional band calculation
    %   
    %    See also Atom.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = protected)
        kpts
        nkpts
        wks
        nkxyz = 0
        scfkpts % used in hybrid functional band calculation
		nqs
    end
    properties (SetAccess = public)
		modifyidxnz = false
		qeidxnz
    end
    methods
    	function cry = Crystal(varargin)
        	if nargin == 0
                return;
            end
            cry = set(cry,varargin{:});
            cry = finalize(cry);
        end
    end
end
