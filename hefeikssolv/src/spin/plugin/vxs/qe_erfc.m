function erfc = qe_erfc(x)   
    p2 =[ 3.004592610201616E2 4.519189537118719E2...
             3.393208167343437E2 1.529892850469404E2...
             4.316222722205674E1 7.211758250883094...
             5.641955174789740E-1 -1.368648573827167E-7];
    q2 = [3.004592609569833E2 7.909509253278980E2...
              9.313540948506096E2 6.389802644656312E2...
              2.775854447439876E2 7.700015293522947E1...
              1.278272731962942E1 1.000000000000000];
    p3 = [-2.996107077035422E-3 -4.947309106232507E-2...
             -2.269565935396869E-1 -2.786613086096478E-1...
             -2.231924597341847E-2];
    q3 = [1.062092305284679E-2 1.913089261078298E-1...
              1.051675107067932 1.987332018171353...
              1.000000000000000];

    pim1 = 0.56418958354775629; 
    ax = abs(x);  
    
    if ax > 26.0
       erfc = 0.0;  
    elseif ax > 4.0  
       x2 = x^2;  
       xm2 = (1.0 / ax)^2;  
       erfc = (1.0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) ...
            + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) ...
            ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * ...
            (q3 (4) + xm2 * q3 (5) ) ) ) ) );
    elseif ax > 0.47 
       x2 = x^2;  
       erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) ...
            + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) ...
            + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * ...
            (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * ...
            (q2 (7) + ax * q2 (8) ) ) ) ) ) ) );
    else  
       erfc = 1.0 - qe_erf (ax);  
    end
    
    if x < 0.0 
        erfc = 2.0 - erfc;  
    end 
end 