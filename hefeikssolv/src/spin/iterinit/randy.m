function [num,info] = randy(info,varargin)
%  generate a number between 0 and 1 stochasticly.
m = 714025;
ia = 1366;
ic = 150889;
ntab = 97;
rm = 1.0/m;

if nargin==2
    info.idum = min(abs(irand),ic);
    info.first = true;
end

if info.first
    info.idum = 0;
    info.ir = zeros(ntab,1);
    info.first = false;
    info.idum = mod(ic-info.idum,m);
    
    for j=1:ntab
        info.idum = mod(ia*info.idum+ic,m);
        info.ir(j) = info.idum;
    end
    info.idum = mod(ia*info.idum+ic,m);
    info.iy = info.idum;
end

j = floor(1+(ntab*info.iy)/m);
if j>ntab&&j<1
    error('randy,j out of range');
end
info.iy = info.ir(j);
num = info.iy*rm;
info.idum = mod(ia*info.idum+ic,m);
info.ir(j) = info.idum;
end
    
    
