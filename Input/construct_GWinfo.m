function GWinfor = struct2GWinfo(data, config)
% Convert from data to @GWinfor
%         data contains fields: 'rhor', 'Vxc', 'ev', 'psig', 'sys', 'occupation',
%                               'reciprocal_grid_info'.

% Step 1: Extract data
rhor = data.rhor;
Vxc = data.Vxc;
ev = data.ev;
psig = data.psig;
sys = data.sys;
occupation = data.occupation;
reciprocal_grid_info = data.reciprocal_grid_info;

% Step 2: Convert to @GWinfo
% @GWinfo contains fields: coulG, coulG0, supercell, bdot, qk, ntot(number of total band calculated),
% vol, ne, nv, gvec, gvec2(double Ecut), gvecrho(for density), Vxc, rho(density in reciprocal space),
% ev, Z
% abandoned: aqs, idxnz
ha2ry = 2.0;
ry2ev = 13.60569253;
GWinfor = GWinfo();

% Compute Coulomb interaction
disp("construct_Dcoul.m is still ongoing")
return
% coulG = construct_vcoul(data, config);
% coulG0 = construct_coulG0(data, config);

end