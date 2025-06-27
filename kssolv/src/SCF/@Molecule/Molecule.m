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
    %        vol         The volumn of the system
    %        atoms       The atom list for distinct atoms
    %        alist       List of atom indices in atoms
    %        natoms      The number of atoms for each distinct type of atom
    %        xyzlist     The list of x,y,z location of the atoms
    %        ecut        Energy cut
    %        ecut2       It is usually 4 times Energy cut (TODO:not clear)
    %        n1,n2,n3    Number of discretization points in each dimension
    %        vext        External potential on the molecule
    %        smear       Broadening method for occupations
    %        temperature The temperature of the system
    %        alat        The crystal contant
    %        nspin       Number of density components
    %        tot_mag     Total magnetic moments
    %        lspinorb    Whether including spin-orbit coupling or not
    %        initmag     Initial magnetic moments of atoms
    %        noncolin    Whether spin-noncollinear calculation or not
    %        domag       Whether open switch for magnetization calculation or not
    %        lsda        Whether spin-polarized calculation or not
    %        nel         The number of electrons in the molecule
    %        nbnd        The total number of bands
    %        extranbnd   The number of empty bands
    %        
    %    See also Atom.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = protected)
        name
        supercell
        vol
        atoms
        alist
        natoms
        ecut
        ecut2
        n1
        n2
        n3
        vext
        smear
        temperature
        alat
        nspin
        tot_mag
        lspinorb
        initmag
        noncolin
        domag
        lsda
    end
    properties (SetAccess = public)
        nel
        nelup
        neldw
        nbnd
    	extranbnd
        funct
        ppvar
        xyzlist
        efermi
        xyzforce
        info
    	zeta
        lsign
        ux
        omega
        exx_lr
        exx_sr
        dft_c
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
