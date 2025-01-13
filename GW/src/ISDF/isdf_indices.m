function [ind_mu] = isdf_indices(Phi, Psi, options)
% COMPARE to isdf.m, this code generate INDICES of ISDF only.

% The function isdf() performs the
% Interpolative Separable Density Fitting
% decomposition:
%   phi_i(r)*psi_j(r) = sum_{mu} zeta_mu(r)*phi_i(r_mu)*psi_j(r_mu).
%
% Inputs:
%   Phi: A set of wavefunctions (size: M*N1).
%   Psi: Another set of wavefunctions (size: M*N2).
%   option: A structure.
%  fields in option:
%     rk: rank of product states [ Phi_i*Psi_j ], 1 <= rk <= min(M, N1*N2).
%    samp: sample method, 'qrcp' or 'kmeans'
%    seed: random seed controls random number generator
%    weight: weight type used in kmeans, 'multiply' or 'add'
%    init: initialization method in kmeans
%    points: cluster points in kmeans
% Outputs:
%   ind_mu: indices of interpolation points

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

if m1 ~= m2
    error('Wrong inputs: row dimensions of Phi and Psi do not match\n');
end

opt=options.isdfoptions;
rng(opt.seed,'twister');

switch lower(options.exxmethod)
    case {'default'}
        %Obtain the rank by QRCP in the first run, then change to K-means
        %method.
        %Guess an initial rank
        rk = min(max(n1,n2)*16,n1*n2);
        for iter=1:4
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
            %This threshold controls the accuracy of ISDF decomposition.
            threshold=max(abs(diag(R)))*options.isdf_threshold;
            n=sum(abs(diag(R))>threshold);
            if n>=rk
                %Initial guess of rk is too small
                rk = rk*2;
                if iter==4
                    fprintf('Warning: ISDF rank is not large enough, error may be large!\n');
                end
            else
                rk = n;
                fprintf('ISDF rank is set to %i.\n',rk);
                break
            end
        end
        ind_mu = e(1:rk);
        options.exxmethod = 'kmeans';
        options.isdfoptions.rank = rk;
    case {'qrcp'}
        %Do dimension reduction on Phi and Psi via random Guassian projection
        rk=opt.rank;
        rsamp = 1.2*rk;
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
    case 'kmeans'
        %Compute electron density
        rk=opt.rank;
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
        tic
        ind_mu = k_means(rk, weight_square, options);
        toc
        gcp
        tic
        ind_mu = k_means_hpc(rk, weight_square, options);
        toc
        % Construct the interpolation basis.
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

