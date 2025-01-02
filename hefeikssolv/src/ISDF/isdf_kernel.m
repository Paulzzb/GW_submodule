function zeta_mu = isdf_kernel(Phi, Psi, ind_mu, options)
% Least-squares fitting for interpolation vectors zeta_mu
if ~iscell(Phi)
    [m1, ~] = size(Phi);
    [m2, ~] = size(Psi);
else
    [m1, ~] = size(Phi{1});
    [m2, ~] = size(Psi{2});
end

if m1 ~= m2
    error('Wrong inputs: row dimensions of Phi and Psi do not match!\n');
end
m = m1;
if length(ind_mu) > m
    error('Wrong inputs: ind_mu is too long!\n');
end

rk=options.isdfoptions.rank;
if ~iscell(Phi)
    phi=Phi(ind_mu,:);
    psi=Psi(ind_mu,:);
    C1=(Phi*phi').*(Psi*psi');
    C2=(phi*phi').*(psi*psi');
    flag=sum(abs(C2-C2'),'all');
    if flag > 1e-16
        fprintf('Error introduced in C2\n');
        C2=(C2+C2')/2;
    end
    zeta_mu=C1/C2;
else
    % Fitting procedures of twice-sampling ISDF for decomposing 
    % the sum of orbital pairs
    phi = cell(2,1);
    psi = cell(2,1);
    for is = 1:2
        phi{is} = Phi{is}(ind_mu,:);
        psi{is} = Psi{is}(ind_mu,:);
    end
    C1 = zeros(m,rk);
    C2 = zeros(rk,rk);
    for i = 1:2
        for j = 1:2
            C1 = C1 + (Phi{i}*phi{j}').*(Psi{i}*psi{j}');
        end
    end
    for i = 1:2
        C2 = C2 + (phi{i}*phi{i}').*(psi{i}*psi{i}');
    end
    C2_updw = (phi{1}*phi{2}').*(psi{1}*psi{2}');
    C2 = C2 + C2_updw + C2_updw';
    flag=sum(abs(C2-C2'),'all');
    if flag > 1e-16
        fprintf('Error introduced in C2\n');
        C2=(C2+C2')/2;
    end
    zeta_mu = C1/C2;
end
