function [ux, lsign] = compute_ux(mag)
% if the initial magnetizations of all atoms are in the same direction,
% lsign = 1 and ux is the unit axis vector which is used for the GGA
% functinal calculation.
% else, lsign = 0.
ux = 0;
lsign = false;
na = length(mag.amag);
m_loc = zeros(3,na);
m_loc(1,:) = (mag.amag).*sin(mag.theta).*cos(mag.phi);
m_loc(2,:) = (mag.amag).*sin(mag.theta).*sin(mag.phi);
m_loc(3,:) = (mag.amag).*cos(mag.theta);

idx = (abs(mag.amag) > 1e-6);
m_loc = m_loc(:,idx);
if sum(idx)
    lsign = true;
    ux = m_loc(:,1);
end
for i = 2:sum(idx)
    lsign = lsign&&is_parallel(ux,m_loc(:,i));
end
if lsign
    ux = ux/sqrt(ux'*ux);
    fprintf("Fixed quantization axis for GGA.\n");
end
end


