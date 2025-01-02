function dos = getdos(ev, nx, sigma)

% Initialization 
[m,n]   = size(ev);
ne      = m;
emin    = min(ev);
emax    = max(ev);
xgrid   = (emax - emin)/(nx - 1);
dos     = zeros(nx,2);
  
for ix = 1 : nx
  dos(ix,1) = emin + (ix - 1)*xgrid;
end

for ie = 1 : ne
  for ix = 1 : nx
    x = emin + (ix - 1)*xgrid - ev(ie);
    %if(Gaussian)
    dos(ix,2) = dos(ix,2) + 1/(sigma*sqrt(2*pi))*exp(-x^2/(2*sigma^2));
    %if(Lorentzian) 
    %dos(ix,2) = dos(ix,2) + sigma/(pi*(x^2+sigma^2));
  end
end

