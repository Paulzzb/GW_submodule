function gw_fullfreq_cd(GWinfo, options)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling gwCalculation.
% 
% Parameters
%   Input:
%     GWinfo: class @GWinfo, contains ground-state information.
%     options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%     GWoutput (directly saved) : contains the calculated self-energies.
% 
% Structure
%   The code contains the following steps 
% 
%   if use_ISDF
%     call gw_fullfreq_cd_res_ISDF
%     call gw_fullfreq_cd_int_ISDF
%   else
%     call gw_fullfreq_cd_res
%     call gw_fullfreq_cd_int
%   end
%   Calculate the exchange part of self energies.
%   Generate GWoutput

if (options.ISDFCauchy.isISDF)
  Eres = gw_fullfreq_cd_res_ISDF(GWinfo, options);
  Eint = gw_fullfreq_cd_int_ISDF(GWinfo, options);
else
  Eres = gw_fullfreq_cd_res(GWinfo, options);
  Eint = gw_fullfreq_cd_int(GWinfo, options);
end

Ex = gw_x(GWinfo, options);
Vxc = GWinfo.Vxc * options.Constant.ry2ev;
ev  = GWinfo.ev(1:length(Ex))  * options.Constant.ry2ev;
Sigma = ev - Vxc + Ex + Eres + Eint
GWenergy.Eres = Eres;
GWenergy.Eint = Eint;
GWenergy.Ex = Ex;
GWenergy.Vxc = Vxc;
GWenergy.ev = ev;
save(options.GWCal.fileName, 'GWenergy')

end
