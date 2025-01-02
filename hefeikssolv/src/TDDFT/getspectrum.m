function spectrum = getspectrum(mol, ksinfo, nvbands, ncbands, emax, nw, sigma)

% Initialization 
nv      = ksinfo.nv;
Z       = ksinfo.Z;
ev      = ksinfo.ev;
F       = ksinfo.F;
vol     = ksinfo.vol;
ntot    = ksinfo.ntot;
rho     = ksinfo.rho;

n1      = mol.n1;
n2      = mol.n2;
n3      = mol.n3;
[ng,nr] = size(F);
wgrid   = (emax - 0.0)/(nw - 1);

[x,y,z] = ndgrid(0:n1-1,0:n2-1,0:n3-1);
xyz = [x(:)/n1 y(:)/n2 z(:)/n3]*mol.supercell;
x = reshape(xyz(:,1),n1,n2,n3);
y = reshape(xyz(:,2),n1,n2,n3);
z = reshape(xyz(:,3),n1,n2,n3);

for iv = 1 : nv+ncbands
  Z(:,iv) = Z(:,iv) / (norm(Z(:,iv)));
end

psivr = zeros(nr, nvbands);
psicr = zeros(nr, ncbands);
spectrum = zeros(nw,2);

for iw = 1:nw
  spectrum(iw,1) = (iw - 1)*wgrid;
end

for iv = 1:nvbands
  psivr(:,iv) = conj(F'*Z(:,nv-nvbands+iv));
  for ic = 1:ncbands
    de = ev(nv+ic) - ev(nv-nvbands+iv);
    if(de < emax)
      psicr(:,ic) = F'*Z(:,nv+ic);
      psixr = 0.0;
      psiyr = 0.0;
      psizr = 0.0;
      for ir = 1 : nr
        psixr = psixr +  psivr(ir,iv)*x(ir)*psicr(ir,ic); 
        psiyr = psiyr +  psivr(ir,iv)*y(ir)*psicr(ir,ic); 
        psizr = psizr +  psivr(ir,iv)*z(ir)*psicr(ir,ic); 
      end
      psivcr = (abs(psixr))^2 + (abs(psiyr))^2 + (abs(psizr))^2;
      for iw = 1:nw
        w = (iw - 1)*wgrid - de;
         %if(Gaussian)
         spectrum(iw,2) = spectrum(iw,2) + 1/(sigma*sqrt(2*pi))*exp(-w^2/(2*sigma^2))*psivcr;
         %if(Lorentzian) 
         %spectrum(iw,2) = spectrum(iw,2) + sigma/(pi*(w^2+sigma^2));
      end
    end
  end
end

factor = 1.0;

spectrum = factor*spectrum;
