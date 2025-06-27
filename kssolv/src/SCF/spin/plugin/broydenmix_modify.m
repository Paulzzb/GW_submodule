function [vnew,df,dv,cdf] = broydenmix_modify(mol, vin, vout, beta, df, dv, cdf, iter, ...
 maxmix, brank)
% modified broyden mixing method.
% rho is stored in nspin cell
% brank is the number of the old density information which is stored
% in df,dv,cdf and transmitted to next mixing.
vin = transform_rho(mol,vin,'real2rec');
vout = transform_rho(mol,vout,'real2rec');
if mol.nspin == 2
    vin = transform_rho(mol,vin,'updw2sum');
    vout = transform_rho(mol,vout,'updw2sum');
end

if mol.nspin ~= 1
    for i = 1:mol.nspin
        vin{i} = vin{i}*2;
        vout{i} = vout{i}*2;
    end
else
     vin = vin*2;
     vout = vout*2;
end

if (brank > maxmix)
    fprintf('brank = %d, maxmix = %d\n', brank, maxmix);
    error('brank too big\n');
end

vout = rho_mix(mol,-1,vout,vin);
if isempty(df)
    df = cell(brank,1);
end
if isempty(dv)
    dv = cell(brank,1);
end
if isempty(cdf)
    cdf = cell(2,1);
end

iter_used = min(iter-1,brank);
ipos = iter-1-floor(abs(iter-2)/brank)*brank;
if iter > 1
    df{ipos} = rho_mix(mol,-1,cdf{1},vout);
    dv{ipos} = rho_mix(mol,-1,cdf{2},vin);
end
cdf{1} = vout;
cdf{2} = vin;
if iter_used > 0
    betamix = zeros(iter_used,iter_used);
    %calbeta = tic;
    for i = 1:iter_used
        for j = i:iter_used
            betamix(i,j) = rho_ddot(mol,df{j},df{i});
            betamix(j,i) = betamix(i,j);
        end
    end
    %fprintf(" Betamix calculating time = %20.3e\n",toc(calbeta));
    betamix = inv(betamix);
    work = zeros(iter_used,1);
    %calwork = tic;
    for i = 1:iter_used
        work(i) = rho_ddot(mol,df{i},vout);
    end
    %fprintf(" Work calculating time = %20.3e\n",toc(calwork));
    %calvout = tic;
    for i = 1:iter_used
        gamma0 = betamix(:,i)'*work;
        vin = rho_mix(mol,-gamma0,vin,dv{i});
        vout = rho_mix(mol,-gamma0,vout,df{i});
    end
    %fprintf(" Vout calculating time = %20.3e\n",toc(calvout));
end
vnew = rho_mix(mol,beta,vin,vout); 
vnew = transform_rho(mol,vnew,'rec2real');
if mol.nspin == 2
    vnew = transform_rho(mol,vnew,'sum2updw');
end

if mol.nspin ~= 1
    for i = 1:mol.nspin
        vnew{i} = vnew{i}/2;
    end
else
    vnew = vnew/2;
end
end

    





            
    

    



        
