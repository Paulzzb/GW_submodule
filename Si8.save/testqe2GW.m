cd ..\GW\
gw_startup();
cd ..\hefeikssolv\
KSSOLV_startup();
cd ..\Si8.save\
options_in = [];
options_in.inputfile = './';
options_in.nv_ener = 4;
options_in.nc_ener = 4;
options_in.nc = 16;
[sys, options__] = read_qe_gw_bgw('./');
[GWinput, options] = qe2GW(options_in);

gwCalculation(GWinput, options);