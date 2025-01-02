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
        if ~mol.noncolin
            % LSDA case: select the same-spin realwavefunction
            if mol.lsda
                BPhi = BPhi{Psi.spin};
            end
            F	= KSFFT(mol);
            Psi = F'*Psi.psi;
            F2	= KSFFT(mol,options.exxcut);
            wks = BPhi.wks;
            nk_phi = size(BPhi.psi,3);
            exxgkk	= options.exxgkk;
            VexxPsi = zeros(size(Psi));
            for ik1=1:nk_phi
                %facb equivalent to facb in exx.f90 of QE
                facb = squeeze(exxgkk(:,ik,ik1));
                nocc_max = BPhi.nocc_max(ik1);
                occ=BPhi.occ(1:nocc_max,ik1);
                phi=BPhi.psi(:,1:nocc_max,ik1);
                for i = 1 : size(Psi,2)
                    A = phi.*(F2'*(facb.*(F2*(conj(phi).*Psi(:,i)))));
                    VexxPsi(:,i) = VexxPsi(:,i)+wks(ik1)*sum(A.*occ',2);
                end
            end
            VexxPsi = F*VexxPsi;
        else
            F	= KSFFT(mol);
            npw = length(Psi.idxnz);
            Psi_up = F'*Psi.psi(1:npw,:);
            Psi_dw = F'*Psi.psi(npw+1:end,:);
            F2	= KSFFT(mol,options.exxcut);
            wks = BPhi.wks;
            nk_phi = size(BPhi.psi{1},3);
            exxgkk	= options.exxgkk;
            VexxPsi_up = zeros(size(Psi_up));
            VexxPsi_dw = zeros(size(Psi_dw));
            for ik1=1:nk_phi
                %facb equivalent to facb in exx.f90 of QE
                facb = squeeze(exxgkk(:,ik,ik1));
                nocc_max = BPhi.nocc_max(ik1);
                occ=BPhi.occ(1:nocc_max,ik1);
                phi_up=BPhi.psi{1}(:,1:nocc_max,ik1);
                phi_dw=BPhi.psi{2}(:,1:nocc_max,ik1);
                for i = 1 : size(Psi_up,2)
                    P_r = F2'*(facb.*(F2*(conj(phi_up).*Psi_up(:,i))))...
                        +F2'*(facb.*(F2*(conj(phi_dw).*Psi_dw(:,i))));
                    A = phi_up.*P_r;
                    B = phi_dw.*P_r;
                    VexxPsi_up(:,i) = VexxPsi_up(:,i)+wks(ik1)*sum(A.*occ',2);
                    VexxPsi_dw(:,i) = VexxPsi_dw(:,i)+wks(ik1)*sum(B.*occ',2);
                end
            end
            VexxPsi_up = F*VexxPsi_up;
            VexxPsi_dw = F*VexxPsi_dw;
            VexxPsi    = [VexxPsi_up;VexxPsi_dw];
        end  
        case {'ace'}
        if mol.nspin == 2
            zeta = mol.zeta{ik+(Psi.spin-1)*mol.nkpts};
            VexxPsi = zeta*(zeta'*Psi);
        else
            zeta = mol.zeta{ik};
            VexxPsi = zeta*(zeta'*Psi);
        end
	case {'qrcp','kmeans','nscf','fake-scf'}
        if mol.nspin == 2
            zeta = mol.zeta{Psi.spin}(:,:,ik);
            VexxPsi = zeta*(zeta'*Psi);
        else
            zeta = mol.zeta(:,:,ik);
            VexxPsi = zeta*(zeta'*Psi);
        end
	otherwise
		error('Do not know how to calculate Vexx.')
end
end

function VexxPsi = Vexx0(Psi, Phi, mol, options)
switch options.exxmethod
	case 'normal'
        if ~mol.noncolin
             fprintf('------------Calculation of Vexx*X in standard method------------\n');
            if mol.lsda
                Phi = Phi{Psi.spin};
            end
            F   = KSFFT(mol);
            Psi = F'*Psi.psi;
            F2  = KSFFT(mol,options.exxcut);
            exxgkk = options.exxgkk;
            VexxPsi = zeros(size(Psi));
            n1=mol.n1;n2=mol.n2;n3=mol.n3;
            exxgkk3D=zeros(n1,n2,n3);
            exxgkk3D(get(F2,'idxnz'))=exxgkk;

            %for i = 1 : size(Psi,2)
             %   A = Phi.psi.*(F2'*(exxgkk.*(F2*(conj(Phi.psi).*Psi(:,i)))));
              %  VexxPsi(:,i) = sum(A.*Phi.occ',2);
            %end
            %VexxPsi = F*VexxPsi;
            for i = 1:size(Psi,2)
                D=zeros(size(Phi.psi));
                for j = 1 : size(Phi.psi,2)
                    B=reshape(conj(Phi.psi(:,j)).*Psi(:,i),n1,n2,n3);
                    B=fft3(B);
                    C=exxgkk3D.*B;
                    C=ifft3(C);
                    D(:,j)=C(:);
                end
                A=Phi.psi.*D;
                VexxPsi(:,i) = A*Phi.occ;
            end
            VexxPsi = F*VexxPsi;
        else
             fprintf('------------Calculation of Vexx*X in standard method------------\n');
             F   = KSFFT(mol);
             npw = length(Psi.idxnz);
             n123 = mol.n1*mol.n2*mol.n3;
             Psi_up = F'*Psi.psi(1:npw,:); 
             Psi_dw = F'*Psi.psi(npw+1:end,:);
             Phi_up = Phi.psi{1};
             Phi_dw = Phi.psi{2};
             F2  = KSFFT(mol,options.exxcut);
             exxgkk = options.exxgkk;
             VexxPsir = zeros(n123*2,ncols(Psi));
             VexxPsi  = zeros(npw*2,ncols(Psi));
             %for i = 1 : ncols(Psi)
              %   A = F2'*(exxgkk.*(F2*(conj(Phi_up).*Psi_up(:,i)+conj(Phi_dw).*Psi_dw(:,i))));
               %  VexxPsir(1:n123,i) = sum((Phi_up.*A).*Phi.occ',2);
                % VexxPsir(n123+1:end,i) = sum((Phi_dw.*A).*Phi.occ',2);
            %end
            %clear A;
            n1=mol.n1;n2=mol.n2;n3=mol.n3;
            exxgkk3D=zeros(n1,n2,n3);
            exxgkk3D(get(F2,'idxnz'))=exxgkk;

    for i = 1 : size(Psi_up,2)
        D=zeros(size(Phi_up));
        for j = 1 : size(Phi_up,2)
            B=reshape(conj(Phi_up(:,j)).*Psi_up(:,i)+conj(Phi_dw(:,j)).*Psi_dw(:,i),n1,n2,n3);
            B=fft3(B);
            C=exxgkk3D.*B;
            C=ifft3(C);
            D(:,j)=C(:);
        end
        A_up=Phi_up.*D;
        A_dw=Phi_dw.*D;
        VexxPsir(1:n123,i) = A_up*Phi.occ;
        VexxPsir(n123+1:end,i) = A_dw*Phi.occ;
    end

    clear A_up;
    clear A_dw;

            VexxPsi(1:npw,:) = F*VexxPsir(1:n123,:);
            VexxPsi(npw+1:end,:) = F*VexxPsir(n123+1:end,:);
        end
	case {'ace','qrcp','kmeans'}
        if mol.nspin == 1||mol.nspin == 4
            zeta = mol.zeta;
            VexxPsi = zeta*(zeta'*Psi);
        elseif mol.nspin == 2
            zeta = mol.zeta{Psi.spin};
            VexxPsi = zeta*(zeta'*Psi);
        end
	otherwise
		error('Do not know how to calculate Vexx.')
end
end
