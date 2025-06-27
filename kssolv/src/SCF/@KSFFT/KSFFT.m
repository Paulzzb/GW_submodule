classdef KSFFT
	% KSFFT KSSOLV class for Fast Fourier transform
	%    F = KSFFT() returns an empty KSFFT
	%
	%    F = KSFFT(mol) returns KSFFT with respect to the molecule.
	%
	%    F = KSFFT(mol,ecut) returns KSFFT with respect to the molecule
	%    and a specified energy cutoff.
	%
	%    KSFFT constructs an inverse Fourier transform
	%    that converts Fourier coefficients of a wavefunction 
	%    (saved as a vector) to a wavefunction in real space
	%    saved as a vector
	%
	%    Once constructed, F can be used as a matrix
	%
	%    y = F*x gives the inverse Fourier transform of x
	%    x = F'*y gives the Fourier transform of y
	%
	%    F contains the volume factor
	%
	%    The KSFFT class contains the following fields.
	%        Field       Explaination
	%      ----------------------------------------------------------
	%        n1,n2,n3    The mesh size for FFT
	%        vol         The volumn of the supercell
	%        idxnz       Indices for non-zero entries in the n1,n2,n3
    %        forward     Whether the FFT is forward, i.e., from real space
    %                    to G space
    %        inverse     Whether the FFT is inverse, i.e., from G space
    %                    to real space
	%
	%    See also KSFFT.

	%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
	%                          Stanford University and Lawrence Berkeley
	%                          National Laboratory
	%  This file is distributed under the terms of the MIT License.

	properties (SetAccess = protected)
		n1
    	n2
    	n3
    	vol
    	idxnz
    	forward
    	inverse
	end
	methods
		function F = KSFFT(varargin)
			switch (nargin)
    			case 0
        			F.n1 = [];
       				F.n2 = [];
        			F.n3 = [];
                    F.vol = [];
        			F.idxnz = [];
        			F.forward = 1;
        			F.inverse = 0;
    			case 1
        			mol = varargin{1};
        			if ( isa(mol,'Molecule') )
            			% the input arguement is a Molecule object
            			n1 = get(mol,'n1');
            			n2 = get(mol,'n2');
            			n3 = get(mol,'n3');
            			F.n1 = n1;
            			F.n2 = n2;
            			F.n3 = n3;
            			grid = Ggrid(mol);
            			F.vol = det(get(mol,'supercell'));
            			F.idxnz = get(grid,'idxnz');
						if (isprop(mol,'modifyidxnz')&&(mol.modifyidxnz))
							F.idxnz=mol.qeidxnz;
						end
                        F.forward = 1;
                    	F.inverse = 0;
        			else
            			error('The input must be a Molecule object');
        			end
    			case 2
        			mol = varargin{1};
					ecut= varargin{2};
					if abs(ecut-mol.ecut)>1e-16 && abs(ecut-mol.ecut2)>1e-16
						fprintf('##############################################################\n')
						fprintf('Warning! The ecut may be wrong!\n')
						fprintf('input ecut:%f, mol.ecut:%f\n',ecut,mol.ecut);
					end
        			if ( isa(mol,'Molecule') )
            			% the input arguement is a Molecule object
            			n1 = get(mol,'n1');
            			n2 = get(mol,'n2');
            			n3 = get(mol,'n3');
            			F.n1 = n1;
            			F.n2 = n2;
            			F.n3 = n3;
            			grid = Ggrid(mol,ecut);
            			F.vol = det(get(mol,'supercell'));
            			F.idxnz = get(grid,'idxnz');
                        F.forward = 1;
                    	F.inverse = 0;
        			else
            			error('The input must be a Molecule object');
        			end
    			otherwise
        			error('Must have one arguement');
			end
		end
	end
end
