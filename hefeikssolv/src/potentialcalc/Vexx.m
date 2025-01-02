function VexxPsi = Vexx(Y, X, mol, options)
% Exact Exchange Operator
if isa(mol,'Crystal')
    VexxPsi = Vexxk(Y, X, mol, options);
else %Molecule
    VexxPsi = Vexx0(Y, X, mol, options);
end
%fprintf('Vexx time: %f\n',t);
VexxPsi = -(0.25*mol.vol) * VexxPsi;
end

function VexxPsi = Vexxk(Psi, BPhi, mol, options)
assert(isprop(Psi,'ik'),'k index of Psi is empty!');
ik	= Psi.ik;

switch options.exxmethod
    case 'normal'
        fprintf('------------Calculation of Vexx*X in standard method------------\n');
        n1=mol.n1;n2=mol.n2;n3=mol.n3;
        n123  = n1*n2*n3;
        F = KSFFT(mol);
        F2 = KSFFT(mol,options.exxcut);
        exxgkk = options.exxgkk;
        % LSDA case: select realwavefunctions labeled with same spin
        if mol.lsda
            BPhi = BPhi{Psi.ispin};
        end
        wks = BPhi.wks;

        if ~mol.noncolin
	    Psi = F'*Psi.psi;
	    VexxPsi = zeros(size(Psi));
            nk_phi = size(BPhi.psi,3);
        else
            npw = length(Psi.idxnz);
            Psi_up = F'*Psi.psi(1:npw,:);
            Psi_dw = F'*Psi.psi(npw+1:end,:);
            VexxPsi_up = zeros(size(Psi_up));
            VexxPsi_dw = zeros(size(Psi_dw));
            nk_phi = size(BPhi.psi{1},3);
        end
      
	for ik1=1:nk_phi
	    %facb equivalent to facb in exx.f90 of QE
	    facb = squeeze(exxgkk(:,ik,ik1));
    	    nocc_max = BPhi.nocc_max(ik1);
	    occ = BPhi.occ(1:nocc_max,ik1);
            if ~mol.noncolin
	        phi = BPhi.psi(:,1:nocc_max,ik1);
	        for i = 1 : size(Psi,2)
		    prod_pair = conj(phi).*Psi(:,i);
                    kernel = Poisson_solver(mol,prod_pair,facb);
		    VexxPsi(:,i) = VexxPsi(:,i) + wks(ik1)*((phi.*kernel)*occ);
	        end
            else
                phi_up=BPhi.psi{1}(:,1:nocc_max,ik1);
                phi_dw=BPhi.psi{2}(:,1:nocc_max,ik1);
                for i = 1 : size(Psi.psi,2)
                    prod_pair = conj(phi_up).*Psi_up(:,i) + ...
                                conj(phi_dw).*Psi_dw(:,i);
                    kernel = Poisson_solver(mol,prod_pair,facb);
                    VexxPsi_up(:,i) = VexxPsi_up(:,i) + wks(ik1)*((phi_up.*kernel)*occ);
                    VexxPsi_dw(:,i) = VexxPsi_dw(:,i) + wks(ik1)*((phi_dw.*kernel)*occ);
                end
	    end
        end
         
        if ~mol.noncolin
             VexxPsi = F*VexxPsi;
        else
             VexxPsi = [F*VexxPsi_up;F*VexxPsi_dw];
        end                   
    case {'ace','qrcp','kmeans','nscf','fake-scf'}
	if ~mol.lsda
            zeta = mol.zeta{ik};
            VexxPsi = zeta*(zeta'*Psi);
        else
            zeta = mol.zeta{Psi.ispin}{ik};
            VexxPsi = zeta*(zeta'*Psi);
        end
    otherwise
	error('Do not know how to calculate Vexx.')
    end
end

function VexxPsi = Vexx0(Psi, Phi, mol, options)
switch options.exxmethod
    case 'normal'
        fprintf('------------Calculation of Vexx*X in standard method------------\n');
        n1=mol.n1;n2=mol.n2;n3=mol.n3;
        n123  = n1*n2*n3;
        F = KSFFT(mol);
        F2 = KSFFT(mol,options.exxcut);
        exxgkk = options.exxgkk;

        if mol.lsda
            Phi = Phi{Psi.ispin};
        end
        
        if ~mol.noncolin
            Psi = F'*Psi.psi;
            VexxPsi = zeros(size(Psi));
        else
            npw = length(Psi.idxnz);
            Psi_up = F'*Psi.psi(1:npw,:);
            Psi_dw = F'*Psi.psi(npw+1:end,:);
            VexxPsi_up = zeros(size(Psi_up));
            VexxPsi_dw = zeros(size(Psi_dw));
        end
        
        if ~mol.noncolin
            for i = 1:size(Psi,2)
                prod_pair = conj(Phi.psi).*Psi(:,i);
                kernel = Poisson_solver(mol,prod_pair,exxgkk);
                VexxPsi(:,i) = (Phi.psi.*kernel)*Phi.occ;
            end
            VexxPsi = F*VexxPsi;
        else  
            phi_up = Phi.psi{1};
            phi_dw = Phi.psi{2};
            for i = 1 : size(Psi_up,2)
                prod_pair = conj(phi_up).*Psi_up(:,i) +...
                            conj(phi_dw).*Psi_dw(:,i); 
                kernel = Poisson_solver(mol,prod_pair,exxgkk);
                VexxPsi_up(:,i) = (phi_up.*kernel)*Phi.occ;
                VexxPsi_dw(:,i) = (phi_dw.*kernel)*Phi.occ;
            end
            VexxPsi = [F*VexxPsi_up;F*VexxPsi_dw];
        end
    case {'ace','qrcp','kmeans'}
        if mol.nspin == 1||mol.nspin == 4
            zeta = mol.zeta;
            VexxPsi = zeta*(zeta'*Psi);
        elseif mol.nspin == 2
            zeta = mol.zeta{Psi.ispin};
            VexxPsi = zeta*(zeta'*Psi);
        end
    otherwise
        error('Do not know how to calculate Vexx.')
end
end
