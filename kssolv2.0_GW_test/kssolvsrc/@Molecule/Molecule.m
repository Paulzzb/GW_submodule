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
    %        xyzlist     The list of x,y,z positions of the atoms
    %        ecutwfc     Energy cut for wavefunction
    %        ecutrho     Energy cut for density and potential (4*ecutwfc)
    %        gridwfc     Grid corresponds to ecutwfc
    %        gridrho     Grid corresponds to ecutrho
    %        n1,n2,n3    Number of discretization points in each dimension
    %        vext        External potential on the molecule
    %        nspin       Number of spins
    %        temperature The temperature of the system
    %        natoms      The number of atoms for each distinct type of atom
    %        vol         The volumn of the system
    %        nel         The number of electrons in the molecule
    %
    %    See also Atom.
    
    properties (SetAccess = protected)
        name
        supercell
        atoms
        alist
        xyzlist
        ecutwfc
        ecutrho
        gridwfc
        gridrho
        n1
        n2
        n3
        vext
        nspin
        nspinor
        temperature
        natoms
        nbnd
        vol
        nel
    end
    properties (SetAccess = public)
        funct
        ppvar
        xyzforce
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
