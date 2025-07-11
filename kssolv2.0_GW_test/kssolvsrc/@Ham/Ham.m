classdef Ham
    % HAM KSSOLV class for Hamiltonion
    %    H = HAM() returns an empty Hamiltonion.
    %
    %    X = HAM(mol) returns Hamiltonion with respect to the molecule.
    %
    %    H = HAM(mol,rho) returns Hamiltonion with respect to the molecule
    %    and rho as the initial density.
    %
    %    The Ham class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        n1,n2,n3    Number of discretization points in each dimension
    %        gkin        Kinetic energy in Fourier space for wave function,
    %                    i.e., gkin = gk^2/2/me
    %        idxnz       Indices for non-zero entries in the n1,n2,n3
    %        vtot        Total local potential
    %        vion        Local potential from the pseudo potential
    %        vext        External potential
    %        vnp         Density-dependent potential (Hartree + xc)
    %        vnlmat      Nonlocal potential from the pseudo potential,
    %                    which is known as the beta in the KB format
    %        vnlsign     Nonlocal potential form the pseudo potential,
    %                    which is known as the middle matrix in KB format
    %        rho         Density function in real space
    %        ishybrid    Hybrid xc functionals (Now ony HSE06 is available)
    %        vexx        Exact exchange potential
    %        eband       Eigenvalue of Hamiltonian
    %
    %    See also Atom, Molecule, Wavefun.
    
    %properties ( SetAccess = ?BlochHam )
    properties ( SetAccess = public ) % make it public for now
        n1
        n2
        n3
        gkin
        idxnz
        idxnz2
        vion
        vnlmat
        vnlsign
        ishybrid
    end
    properties (SetAccess = public)
        rho
        vtot
        vnp
        vext
        vexx
        eband
    end
    methods
        function H = Ham(mol,rho)
            if nargin == 0
                H.n1 = 0;
                H.n2 = 0;
                H.n3 = 0;
                return;
            end
            
            H.n1 = mol.n1;
            H.n2 = mol.n2;
            H.n3 = mol.n3;
            grid  = mol.gridwfc;
            H.gkin  = grid.gkk/(2*meDef());
            H.idxnz = grid.idxnz;
            
            if isempty(mol.ppvar)
                mol.ppvar  = PpVariable(mol);
                ne = 0;
                for it = 1:length(mol.atoms)
                    ne = ne + mol.ppvar.venums(it)*mol.natoms(it);
                end
                if ne ~= mol.nel
                    warning(['The number of valence electrons in Pp ' ...
                        'file is different fron Periodic Table']);
                    mol = set(mol,'nel',ne);
                end
            end
            
            % compute the local ionic potential
            if nargin == 1
                [H.vion,rho] = getvloc(mol);
            else
                [H.vion,~] = getvloc(mol);
            end
            
            % renormalize charge density
            rho = corrcharge(mol,rho);
            
            % compute nonlocal ionic potential
            [H.vnlmat,H.vnlsign] = getvnl(mol);
            
            % Calculate Hartree and exchange-correlation potential
            % from the known density rho
            [vhart,vxc,~,rho]=getVhxc(mol,rho);
            
            % Save a copy of the Hartree and exchange-correlation
            % separately for DCM updating purposes
            H.vnp = vhart+vxc;
            
            % there may be other external potential
            H.vext = mol.vext;
            
            % vtot does not include non-local ionici potential
            H.vtot = getVtot(H.vion, H.vext, vhart, vxc);
            H.rho  = rho;
            H.ishybrid = false;
        end
    end
end
