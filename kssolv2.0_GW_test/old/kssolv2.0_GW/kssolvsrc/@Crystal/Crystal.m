classdef Crystal < Molecule
    % CRYSTAL KSSOLV class for crystal
    %    cry = CRYSTAL(str1,field1,str2,field2,...) returns a crystal
    %    class of the given fields with respect to the name strings.
    %
    %    The crystal class contains the following fields.
    %        Field     Explaination
    %      ----------------------------------------------------------
    %        name        Crystal name
    %        supercell   The cell of the crystal
    %        atoms       The atom list for distinct atons
    %        alist       List of atom indices in atoms
    %        xyzlist     The list of x,y,z location of the atoms
    %        ecutwfc     Energy cut for wavefunction
    %        ecutrho     Energy cut for density and potential (4*ecutwfc)
    %        gridwfc     Grid corresponds to ecutwfc
    %        gridrho     Grid corresponds to ecutrho
    %        n1,n2,n3    Number of discretization points in each dimension
    %        vext        External potential on the crystal
    %        nspin       Number of spins
    %        kpts        The array of k-points
    %        wks         The weights for k-points
    %        nkpts       The number of k-points
    %        temperature The temperature of the system
    %        natoms      The number of atoms for each distinct type of atom
    %        vol         The volumn of the system
    %        nel         The number of electrons in the crystal
    %
    %    See also Atom.
    
    properties (SetAccess = protected)
        kpts
        wks
        nkpts = 0
        bvec
        bdot
    end
    properties (SetAccess = public)
        idxnz
        vxc
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
