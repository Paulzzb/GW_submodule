function [vc,uc] = VCorrelation_pz_spin(rs,zeta)
%Calculate VCorrelation for linear spin system
%by PZ functional  
p43 = 4.0/3.0;
third = 1.0/3.0;
vc = cell(2,1);

[vcu,ucu] = pz(rs,1);%unpolarized part
[vcp,ucp] = pz_polarized(rs);%polarized part

fz = ((1.0+zeta).^p43 + (1.0-zeta).^p43 - 2.0)/(2.0^p43-2.0);
dfz = p43*((1.0+zeta).^third-(1.0-zeta).^third)/(2.0^p43-2.0);
  
uc = ucu + fz .* (ucp - ucu);
vc{1} = vcu + fz .* (vcp - vcu) + (ucp - ucu) .* dfz .* ( 1.0 - zeta);
vc{2} = vcu + fz .* (vcp - vcu) + (ucp - ucu) .* dfz .* (-1.0 - zeta);

end

