classdef Molecule
    % MOLECULE KSSOLV class for molecule
    %    mol = MOLECULE(str1,field1,str2,field2,...) returns a molecule
    %    class of the given fields with respect to the name strings.
    %
    %    The molecule class contains the following fields.
    %        Field     Explaination
    %      ----------------------------------------------------------
    %        name        Molecule name
    %        supercell   The cell of the molecule
    %        atoms       The atom list for distinct atoms
    %        alist       List of atom indices in atoms
    %        xyzlist     The list of x,y,z location of the atoms
    %        ecut        Energy cut
    %        ecut2       It is usually 4 times Energy cut (TODO:not clear)
    %        n1,n2,n3    Number of discretization points in each dimension
    %        vext        External potential on the molecule
    %        nspin       Number of density components
    %        tot_mag     Total magnetic moments
    %        smear       Broadening method for occupations
    %        temperature The temperature of the system
    %        natoms      The number of atoms for each distinct type of atom
    %        vol         The volumn of the system
    %        lspinorb    Whether including spin-orbit coupling or not
    %        initmag     Initial magnetic moments of atoms
    %        noncolin    Whether spin-noncollinear calculation or not
    %        domag       Whether open switch for magnetization calculation or not
    %        lsda        Whether spin-polarized calculation or not
    %        alat        The crystal contant
    %        nel         The number of electrons in the molecule
    %
    %    See also Atom.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    properties (SetAccess = protected)
        name
        supercell
        atoms
        alist
        ecut
        ecut2
        n1
        n2
        n3
        vext
        nspin
        tot_mag
        smear
        temperature
        natoms
        vol
        lspinorb
        initmag
        noncolin
        domag
        lsda
        alat 
    end
    properties (SetAccess = public)
        nel
        nelup
        neldw
	nbnd
        abnd % whether set nbnd or not
        funct
        ppvar
        xyzlist
	efermi
        xyzforce
	info
	zeta
        lsign
        ux
    end
    methods
        function mol = Molecule(varargin)
            if nargin == 0
                return;
            end
            mol = set(mol,varargin{:});
            mol = finalize(mol);
        end
    end
end
