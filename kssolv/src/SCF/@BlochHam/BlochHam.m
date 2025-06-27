classdef BlochHam
    % BLOCHHAM KSSOLV class for Bloch Hamiltonion
    %    H = BLOCHHAM() returns an empty Bloch Hamiltonion.
    %
    %    X = BLOCHHAM(cry) returns Bloch Hamiltonion with respect to the
    %    crystal.
    %
    %    H = BLOCHHAM(cry,rho) returns Bloch Hamiltonion with respect to
    %    the crystal and rho as the initial rho.
    %
    %    The BlochHam class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        n1,n2,n3    Number of discretization points in each dimension
    %        nkpts       Number of k-points
    %        wks         Weights of k-points
    %        gkincell    Cell of kinetic energy in Fourier space for each
    %                    k-points
    %        idxnz       Indices for non-zero entries in the n1,n2,n3
    %        vion        Local potential from the pseudo potential
    %        vext        External potential
    %        vnlmatcell  Cell of nonlocal potential from the pseudo
    %                    potential for each k-points, which is known as the
    %                    beta in the KB format
    %        vnlsigncell Cell of nonlocal potential form the pseudo
    %                    potential for each k-points, which is known as the
    %                    middle matrix in KB format
    %        lspinorb    Whether including spin-orbit coupling or not
    %        rho         Density function in real space
    %        vtot        Total local potential
    %        vnp         Density-dependent potential (Hartree + xc)
    %        vexx        Exact exchange potential
    %        nspin       Number of density components
    %        eband       Eigenvalues of Hamiltonian
    %           dv       Difference between potential of this and previous iteration
    %        ishybrid    Whether hybrid functional calculations or not
    %
    %    See also Atom, Molecule, Wavefun, Ham, BlochWavefun.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    properties (SetAccess = protected)
        n1
        n2
        n3
        nkpts
        wks
        gkincell
        idxnz
        vion
        vext
        vnlmatcell
        vnlsigncell
        lspinorb
    end
    properties (SetAccess = public)
        rho
        vtot
        vnp
        vexx
        nspin
        eband
        dv 
        ishybrid
    end
    methods
        function BH = BlochHam(cry,rho)
            if nargin == 0
                BH.n1 = 0;
                BH.n2 = 0;
                BH.n3 = 0;
                BH.nkpts = 0;
                return;
            end
            
            BH.n1 = cry.n1;
            BH.n2 = cry.n2;
            BH.n3 = cry.n3;
            BH.nspin = cry.nspin;
            BH.nkpts = cry.nkpts;
            BH.wks = cry.wks;
            BH.ishybrid = cry.exx_sr ~= 0 || cry.exx_lr ~= 0;
            BH.lspinorb = cry.lspinorb;
            
            BH.gkincell = cell(BH.nkpts,1);
            grid  = Ggrid(cry);
            if cry.modifyidxnz
                idxnz = cry.qeidxnz;
                gkx  =  grid.kkxyz(idxnz,1);
                gky  =  grid.kkxyz(idxnz,2);
                gkz  =  grid.kkxyz(idxnz,3);
                for ik = 1:BH.nkpts
                    xyz = [gkx+cry.kpts(ik,1) ...
                        gky+cry.kpts(ik,2) ...
                        gkz+cry.kpts(ik,3)];
                    BH.gkincell{ik} = sum(xyz.*xyz,2)/(2*meDef());
                end
                BH.idxnz = idxnz;
            else
                for ik = 1:BH.nkpts
                    xyz = [grid.gkx+cry.kpts(ik,1) ...
                        grid.gky+cry.kpts(ik,2) ...
                        grid.gkz+cry.kpts(ik,3)];
                    BH.gkincell{ik} = sum(xyz.*xyz,2)/(2*meDef());
                end
                BH.idxnz = grid.idxnz;
            end
            
            % compute the local ionic potential
            if nargin == 1
                [BH.vion,rho] = getvloc(cry);
            else
                [BH.vion,~] = getvloc(cry);
            end
            
            % compute nonlocal ionic potential
            [BH.vnlmatcell,BH.vnlsigncell] = getvnlcell(cry);
            
            % Calculate Hartree and exchange-correlation potential
            % from the known density rho
            [vhart,vxc,~,rho]=getVhxc(cry,rho);
            
            % Save a copy of the Hartree and exchange-correlation
            % separately for DCM updating purposes
            if cry.nspin == 1
                BH.vnp = vhart+vxc;
            elseif cry.nspin == 2
                BH.vnp = cell(2,1);
                BH.vnp{1} = vhart+vxc{1};
                BH.vnp{2} = vhart+vxc{2};
            elseif cry.nspin == 4
                BH.vnp = cell(cry.nspin,1);
                BH.vnp{1} = vhart + vxc{1};
                for i = 2:cry.nspin
                    BH.vnp{i} = vxc{i};
                end
            end
            
            % there may be other external potential
            BH.vext = cry.vext;
            
            % vtot does not include non-local ionici potential
            BH.vtot = getVtot(cry, BH.vion, BH.vext, vhart, vxc);
            BH.rho  = rho;
            
        end
    end
end
