function [zeta_mu, ind_mu, options] = isdf(Phi, Psi, ind_mu_in, options, ispin)
% The function isdf() performs the
% Interpolative Separable Density Fitting
% decomposition:
%   phi_i(r)*psi_j(r) = sum_{mu} zeta_mu(r)*phi_i(r_mu)*psi_j(r_mu).
%
% Inputs:
%   Phi: A set of wavefunctions (size: M*N1).
%   Psi: Another set of wavefunctions (size: M*N2).
%   ind_mu_in: Indices of selected interpolation points in previous ISDF process
%   option: A structure.
%	fields in option:
%   	rk: rank of product states [ Phi_i*Psi_j ], 1 <= rk <= min(M, N1*N2).
%		samp: sample method, 'qrcp' or 'kmeans'
%		seed: random seed controls random number generator
%		weight: weight type used in kmeans, 'multiply' or 'add'
%		init: initialization method in kmeans
%		points: cluster points in kmeans
% Outputs:
%   zeta_mu: interpolation vectors (basis) for density fitting
%   ind_mu: indices of interpolation points

if isempty(ispin)
    ispin = 0;
else

if ~iscell(Phi)
    [m1, n1] = size(Phi);
    [m2, n2] = size(Psi);
else
    Phi{1} = conj(Phi{1});
    Phi{2} = conj(Phi{2});
    [m1, n1] = size(Phi{1});
    [m2, n2] = size(Psi{1});
end

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
            if ~iscell(Phi)
                PhiG = Phi*G1;
                PsiG = Psi*G2;
                BG = prod_states(PhiG, PsiG);
            else
                PhiGup = Phi{1}*G1;
                PhiGdw = Phi{2}*G1;
                PsiGup = Psi{1}*G2;
                PsiGdw = Psi{2}*G2;
                BG = prod_states(PhiGup,PsiGup) + prod_states(PhiGdw,PsiGdw);
            end
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
        %zeta_mu(e,:) = [eye(rk,rk), R(1:rk,1:rk)\R(1:rk,rk+1:end)]';
        % the accuracy of this method to determine interpolative vectors is bad
        % although it's faster. So use the standard LS method.
        zeta_mu = isdf_kernel(Phi, Psi, ind_mu, options);
        options.exxmethod = 'kmeans';
        if ispin == 0
            options.isdfoptions.rank = rk;            
        else
            if isempty(options.isdfoptions.rank)
                options.isdfoptions.rank = zeros(1,2);
            end
            options.isdfoptions.rank(ispin) = rk;
        end
    case {'qrcp'}
        if ~opt.fixip
            fprintf('Sample for interpolation points for ISDF\n');
            %Do dimension reduction on Phi and Psi via random Guassian projection
            if ispin == 0
                rk = opt.rank;
            else
                rk = opt.rank(ispin);
            end

            rsamp = 1.5*rk;
            r1 = min( ceil(sqrt((n1/n2)*rsamp)), n1 );
            r2 = min( ceil(sqrt((n2/n1)*rsamp)), n2 );
            G1 = randn(n1,r1);
            G2 = randn(n2,r2);
            if ~iscell(Phi)
                PhiG = Phi*G1;
                PsiG = Psi*G2;
                BG = prod_states(PhiG, PsiG);
            else
                PhiGup = Phi{1}*G1;
                PhiGdw = Phi{2}*G1;
                PsiGup = Psi{1}*G2;
                PsiGdw = Psi{2}*G2;                
                BG = prod_states(PhiGup,PsiGup) + prod_states(PhiGdw,PsiGdw);
            end
            [~, ~, e] = qr(BG',0);
            clear BG;
            ind_mu = e(1:rk);
        else
            fprintf('Keep the interpolation points for ISDF unchanged\n');
            if ispin == 0
                ind_mu = ind_mu_in;
            else
                ind_mu = ind_mu_in{ispin};
            end
        end
        fprintf('The number of interpolation points is %d\n',length(ind_mu));
        % zeta_mu(e,:) = [eye(rk,rk), R(1:rk,1:rk)\R(1:rk,rk+1:end)]';
        % the accuracy of this method to determine interpolative vectors is bad
        % although it's faster. So use the standard LS method.
        zeta_mu = isdf_kernel(Phi, Psi, ind_mu, options);            
    case 'kmeans'
        if ~opt.fixip
            fprintf('Sample for interpolation points for ISDF\n');
            % Compute weight function for grid points
            if ispin == 0
                rk = opt.rank;
            else
                rk = opt.rank(ispin);
            end

            switch lower(opt.weight)
                case 'prod'
                    if ~iscell(Phi)
                        weight = sum((abs(Phi)).^2, 2).*sum((abs(Psi)).^2, 2);
                    else
                        weight = sum((abs(Phi{1})).^2, 2).*sum((abs(Psi{1})).^2, 2)...
                               + sum((abs(Phi{2})).^2, 2).*sum((abs(Psi{2})).^2, 2);
                    end
                case 'add'
                    if ~iscell(Phi)
                        weight = sum((abs(Phi)).^2, 2) + sum((abs(Psi)).^2, 2);
                    else
                        weight = sum((abs(Phi{1})).^2, 2) + sum((abs(Psi{1})).^2, 2)...
                               + sum((abs(Phi{2})).^2, 2) + sum((abs(Psi{2})).^2, 2);
                    end
                case 'power'
                    if ~iscell(Phi)
                        weight = (sum((abs(Phi)).^2, 2).*sum((abs(Psi)).^2, 2)).^(opt.power/2);
                    else
                    	weight = (sum((abs(Phi{1})).^2, 2).*sum((abs(Psi{1})).^2, 2)).^(opt.power/2)...
                           	+ (sum((abs(Phi{2})).^2, 2).*sum((abs(Psi{2})).^2, 2)).^(opt.power/2);
                	end
            	case 'hf'
                	if ~iscell(Phi)
                    	Z = prod_states(Phi, Psi);
                	else
                    	Z = prod_states(Phi{1}, Psi{1}) + prod_states(Phi{2}, Psi{2});
                	end
                	weight = HF_weight(Z,opt.mol,rk);
            	case 'mix'
                	if ~iscell(Phi)
                    	weight = (sum((abs(Phi)).^2, 2).*sum((abs(Psi)).^2, 2)).^(1/2);
                	else
                    	weight = (sum((abs(Phi{1})).^2, 2).*sum((abs(Psi{1})).^2, 2)).^(1/2)...
                           	+ (sum((abs(Phi{2})).^2, 2).*sum((abs(Psi{2})).^2, 2)).^(1/2);
                	end
                	weight = weight.^2;
            	otherwise
                	error('Unknown weight kind %s.',opt.weight);
        	end
            % Do kmeans
            ind_mu = k_means(rk, weight, options);
        else
            fprintf('Keep the interpolation points for ISDF unchanged\n');
            if ispin == 0
                ind_mu = ind_mu_in;
            else
                ind_mu = ind_mu_in{ispin};
            end
    	end
        fprintf('The number of interpolation points is %d\n',length(ind_mu));
        % Construct the interpolation basis.
        zeta_mu = isdf_kernel(Phi, Psi, ind_mu, options);
    otherwise
        error('Unknown sample method %s.',opt.samp);
    end
end
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
end

