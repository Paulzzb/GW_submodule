function expint = EXPINT(n, x)

    % EXPINT(1, x) is equal to -ei(-x)
    
    maxit=200;
    eps=1E-12;
    big=1.797693127679088E296;
    euler = 0.577215664901532860606512;
    
    if sum(~( n >= 0 & x >= 0 & (x > 0 | n > 1) ))
        fprintf('error for inputed n and x...');
    end

    expint = zeros(size(x));
    
    if n == 0
        expint = exp(-x)./x;
         return;
    end
  
    nm1 = n - 1;
    
    idxs = x == 0;
    expint(idxs) = 1/nm1;
    
    idxl = x > 1;
    b = x(idxl)+n;
    c = big;
    d = 1./b;
    h = d;
    for i = 1:maxit
        a = -i*(nm1+i);
        b = b+2;
        d = 1./(a*d+b);
        c = b+a./c;
        del = c.*d;
        h = h.*del;
        if abs(del-1) <= eps
            break;
        end
    end
    if i > maxit
        error('beyond max iteration');
    end
    expint(idxl) = h.*exp(-x(idxl));
    
    idxm = ~idxl & ~idxs;
    if nm1 ~= 0
        expint(idxm) = 1/nm1;
    else
        expint(idxm) = -log(x(idxm))-euler;
    end
    fact = 1;
    for i = 1:maxit
        fact = -fact/i.*x(idxm);
        if i ~= nm1
            del = -fact/(i-nm1);
        else
            iarsum = 0;
            for k = 1:nm1
                iarsum = iarsum + 1/k;
            end
            del = fact.*(-log(x(idxm))-euler+iarsum);
        end
        expint(idxm) = expint(idxm) + del;
        if abs(del) < abs(expint(idxm))*eps
            break;
        end
    end
    if i > maxit
        error('beyond max iteration');
    end
    
end