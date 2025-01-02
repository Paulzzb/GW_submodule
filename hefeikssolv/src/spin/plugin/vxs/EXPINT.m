function expint = EXPINT(n, x)
    maxit=200;
    eps=1E-12;
    big=1.797693127679088E296;
    euler = 0.577215664901532860606512;
    
    if (~((n >= 0)&&(x >= 0.0)&&((x > 0.0)||(n > 1))))
        fprintf('error for inputed n and x...');
    end

    if n == 0
        expint = exp(-x)/x;
         return;
    end
  
    nm1 = n-1;
    if x == 0.0
        expint = 1.0/nm1;
    elseif x > 1.0
         b = x+n;
         c = big;
         d = 1.0/b;
         h = d;
         for i = 1:maxit
             a = -i*(nm1+i);
             b = b+2.0;
             d = 1.0/(a*d+b);
             c = b+a/c;
             del = c*d;
             h = h*del;
             if abs(del-1.0) <= eps
                 break;
             end
         end
         if i > maxit
             fprintf('beyond max iteration');
         end

         expint = h*exp(-x);
    else
         if nm1 ~= 0
            expint = 1.0/nm1;
         else
            expint = -log(x)-euler;
         end
         fact = 1.0;
         for i = 1:maxit
            fact = -fact*x/i;
            if i ~= nm1
               del = -fact/(i-nm1);
            else
                iarsum = 0.0;
                for k = 1:nm1
                    iarsum = iarsum + 1.0/k;
                end 
                del = fact*(-log(x)-euler+iarsum);
            end
            expint = expint+del;
            if abs(del) < abs(expint)*eps
                break;
            end
         end
        
         if i > maxit
              fprintf('beyond max iteration');
         end
    end
end