function [zeta_mu,ind_mu] = samp(Phi, Psi, options)

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

if m1 ~= m2,
	error('Wrong inputs: row dimensions of Phi and Psi do not match\n');
end

opt=options.isdfoptions;
rk=opt.rank;
rng(opt.seed,'twister');
switch lower(options.exxmethod)
    case 'default'
        %Obtain the rank by QRCP in the first run, then change to K-means
        %method.
        
	case 'qrcp'
		%Do dimension reduction on Phi and Psi via random Guassian projection
		rsamp = 1.5*rk;
		r1 = min( ceil(sqrt((n1/n2)*rsamp)), n1 );
		r2 = min( ceil(sqrt((n2/n1)*rsamp)), n2 );
		G1 = randn(n1,r1);
		G2 = randn(n2,r2);
		PhiG = Phi*G1;
		PsiG = Psi*G2;
		BG = prod_states(PhiG, PsiG);
		[~, R, e] = qr(BG',0);
		clear BG;
		ind_mu = e(1:rk);
		zeta_mu(e,:) = [eye(rk,rk), R(1:rk,1:rk)\R(1:rk,rk+1:end)]';
	case 'kmeans'
		%Compute electron density 
		switch lower(opt.weight)
			case 'prod'
				weight = sum((abs(Phi)).^2, 2).*sum((abs(Psi)).^2, 2);
			case 'add'
				weight = sum((abs(Phi)).^2, 2) + sum((abs(Psi)).^2, 2);
			case 'power'
				weight = (sum((abs(Phi)).^2, 2).*sum((abs(Psi)).^2, 2)).^(opt.power/2);
			case 'hf'
				Z = prod_states(Phi, Psi);
				weight = HF_weight(Z,opt.mol,rk);
			otherwise
				error('Unknown weight kind %s.',opt.weight);
		end
		weight_square=weight.^2;
		
		% Do kmeans
		ind_mu = k_means(rk, weight_square, options);
otherwise
	error('Unknown sample method %s.',opt.samp);
end

function weight = HF_weight(Z,mol,rk)
	grid = Ggrid(mol);
	gkk = grid.gkk;
	F = KSFFT(mol);
	T = randn(size(Z,2),round(1.5*rk));
	ZT = Z*T;
	Zk=F*ZT;
	FZ=zeros(size(Zk));
	for i=1:size(gkk,1)
		kk=gkk(i);
		if kk~=0
			FZ(i,:)=Zk(i,:)*(4*pi/kk);
		end
	end
	G=F'*FZ;
	weight=zeros(size(Z,1),1);
	for i=1:size(Z,1)
		weight(i)=G(i,:)*ZT(i,:)';
	end
return
return
