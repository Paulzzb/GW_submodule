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
    %        nspin       Number of density components
    %        gkin        Kinetic energy in Fourier space for wave function,
    %                    i.e., gkin = gk^2/2/me
    %        idxnz       Indices for non-zero entries in the n1,n2,n3 with ecut
    %        idxnz2      Indices for non-zero entries in the n1,n2,n3 with ecut2
    %        vion        Local potential from the pseudo potential
    %        vnlmat      Nonlocal potential from the pseudo potential,
    %                    which is known as the beta in the KB format
    %        vnlsign     Nonlocal potential form the pseudo potential,
    %                    which is known as the middle matrix in KB format
    %        ishybrid    Hybrid xc functionals (Now ony HSE06 is available)
    %        rho         Density function in real space
    %        vtot        Total local potential
    %        vnp         Density-dependent potential (Hartree + xc)
    %        vext        External potential
    %        vexx        Exact exchange potential
    %        eband       Eigenvalue of Hamiltonian
    %        ik          Index of k-point
    %        wks         Weight of k-point
    %        lspinorb    Whether add spin-orbit coupling effect
    %        dv          Difference between potential of this and previous iteration
    %
    %    See also Atom, Molecule, Wavefun.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    
    properties ( SetAccess = public ) % make it public for now
        n1
        n2
        n3
        nspin
        gkin
        idxnz
        idxnz2
        vion
        vnlmat
        vnlsign
        ishybrid
        rho
        vtot
        vnp
        vext
        vexx
        eband
	    ik
	    wks
        zeta % used in spin-noncollinear H*X
	    F
	    vol
	    exxmethod
        lspinorb
        dv % used to correct force
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
            H.nspin = mol.nspin;
            grid  = Ggrid(mol);
            H.gkin  = grid.gkk/(2*meDef());
            H.idxnz = grid.idxnz;
            H.ishybrid = mol.exx_sr ~= 0 || mol.exx_lr ~= 0;
            H.lspinorb = mol.lspinorb;
            
            if isempty(mol.ppvar)
                mol.ppvar  = PpVariable(mol);
                ne = 0;
                for it = 1:length(mol.atoms)
                    ne = ne + mol.ppvar.venums(it)*mol.natoms(it);
                end
                if ne ~= mol.nel
                    warning(['The number of valence electrons in Pp ' ...
                        'file is different from Periodic Table']);
                    mol = set(mol,'nel',ne);
                end
            end
            
            % compute the local ionic potential
            if nargin == 1
                [H.vion,rho] = getvloc(mol);
            else
                [H.vion,~] = getvloc(mol);
            end

	    	H.ik = -1;
            
            % renormalize charge density
            rho = corrcharge(mol,rho);
            
            % compute nonlocal ionic potential
            [H.vnlmat,H.vnlsign] = getvnl(mol);
            
            % Calculate Hartree and exchange-correlation potential
            % from the known density rho
            [vhart,vxc,~,rho]=getVhxc(mol,rho);
            
            % Save a copy of the Hartree and exchange-correlation
            % separately for DCM updating purposes
            if mol.nspin == 1||(mol.nspin == 4 && ~mol.domag)
                H.vnp = vhart + vxc;
            elseif mol.nspin == 2
                H.vnp = cell(2,1);
                H.vnp{1} = vhart + vxc{1};
                H.vnp{2} = vhart + vxc{2};
            elseif mol.nspin == 4 && mol.domag
                H.vnp = cell(4,1);
                H.vnp{1} = vhart + vxc{1};
                for i = 2:4
                    H.vnp{i} = vxc{i};
                end
            end
            
            % there may be other external potential
            H.vext = mol.vext;
            
            % vtot does not include non-local ionici potential
            H.vtot = getVtot(mol, H.vion, H.vext, vhart, vxc);
            H.rho  = rho;
        end
    end
end
