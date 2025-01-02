function [vx,ux] = VExchange_sla_spin(rho,zeta)
%Calculate the VExchange by slater with alpha = 2/3
%rho is the total denisty, zeta is spin polarizability  
f = - 1.10783814957303361;
%f = -9/8*(3/pi)^(1/3)
alpha = 2.0/3.0;
third = 1.0/3.0;
p43 = 4.0/3.0;
vx = cell(2,1);

rho13 = ((1.0+zeta).*rho).^third;
uxup = f*alpha*rho13;
vx{1} = p43*f*alpha*rho13;
  
rho13 = ((1.0-zeta).*rho).^third;
uxdw = f*alpha*rho13;
vx{2} = p43*f*alpha*rho13;
  
ux = 0.5*((1.0+zeta).*uxup+(1.0-zeta).*uxdw);
  
end

