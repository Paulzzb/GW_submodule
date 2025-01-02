function [uc,v1c_up,v1c_dw,v2c] = VGCCorrelation_pbc_spin(rho,zeta,grho,iflag)
%PBE correlation for LSDA 
  ga = 0.031091;
  be = [0.06672455060314922, 0.046];
  third=1/3; pi34=0.6203504908994;
  % pi34=(3/4pi)^(1/3)
  xkf=1.919158292677513; xks=1.128379167095513;
  % xkf=(9 pi/4)^(1/3) , xks=sqrt(4/pi)

  rs = pi34 ./ rho.^third;
  
  [ec, vc_up, vc_dn] = VCorrelation_pw_spin(rs,zeta);
  kf = xkf ./ rs;
  ks = xks * sqrt(kf);
  
  fz = 0.5*( (1+zeta).^(2/3) + (1-zeta).^(2/3) );
  fz2 = fz .* fz;
  fz3 = fz2 .* fz;

  dfz = ( (1+zeta).^(-1/3) - (1 - zeta).^(-1/3) ) / 3;
  
  t  = sqrt(grho) ./ (2 * fz .* ks .* rho);
  expe = exp( - ec ./ (fz3 * ga) );
  af   = be(iflag) / ga * (1 ./ (expe-1) );
  bfup = expe .* (vc_up - ec) ./ fz3;
  bfdw = expe .* (vc_dn - ec) ./ fz3;
  
  y  = af .* t .* t;
  xy = (1 + y) ./ (1 + y + y .* y);
  qy = y .* y .* (2 + y) ./ (1 + y + y .* y).^2;
  s1 = 1 + be(iflag) / ga .* t .* t .* xy;
  
  h0 = fz3 .* ga .* log(s1);
  
  dh0up = be(iflag) * t .* t .* fz3 ./ s1 .* ( -7/3 * xy - qy .* ...
          (af .* bfup / be(iflag)-7/3) );
  
  dh0dw = be(iflag) * t .* t .* fz3 ./ s1 .* ( -7/3 * xy - qy .* ...
          (af .* bfdw / be(iflag)-7/3) );
  
  dh0zup =   (3 .* h0 ./ fz - be(iflag) * t .* t .* fz2 ./ s1 .*  ...
             (2 .* xy - qy .* (3 * af .* expe .* ec ./ fz3 ./ ...
             be(iflag)+2) ) ) .* dfz .* (1 - zeta);
  
  dh0zdw = - (3 .* h0 ./ fz - be(iflag) * t .* t .* fz2 ./ s1 .*  ...
             (2 .* xy - qy .* (3 * af .* expe .* ec ./ fz3 ./ ...
             be(iflag)+2) ) ) .* dfz .* (1 + zeta);
  
  ddh0 = be(iflag) .* fz ./ (2 * ks .* ks .* rho) .* (xy - qy) ./ s1;
  
   uc     = rho .* h0;
  
  v1c_up = h0 + dh0up + dh0zup;
  v1c_dw = h0 + dh0dw + dh0zdw;

  v2c    = ddh0; 
end

