function [Q, piv, k] = scdm_rand_transform(U,tol,r,rho,ortho)


maxiter = 5;
NN = size(U,1);
count = ceil(4*r*log(r));

rhoabs = abs(rho);
rhoabs = rhoabs/sum(rhoabs);
rhosum = cumsum(rhoabs);

[~, I] = histc(rand(1,count),[0 rhosum']);
I = unique(I);

if length(I) < r
    error('too few columns were sampled, this should not happen');
end
count = length(I);
for k = 1:maxiter
    
    if count == NN
        [Q, R, idx] = qr(U',0);
    else
        [Q, R, idx] = qr(U(I,:)',0);
    end
    idx = idx(1:r);
    piv = I(idx);
    sval = svd(U(piv,:));
    
    if sval(end)/sval(1) > tol || count == NN
        break;
    else
        count = count + r;
        if count > NN
            count = NN;
        end
        
        [~, new] = histc(rand(1,r),[0 rhosum']);
        I = unique([I new]);
    end
    
end

if ~ortho
    Q = U(piv,:)';
end


