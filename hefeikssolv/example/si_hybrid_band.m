%Calculate PBE & HSE energy and HSE band structure
addpath('~/matlab/dev');
KSSOLV_startup
name='si8';
nk=2;
%nk=[4,3,3];
args.nk=nk;
args.ecut=3;
%args.ng=40;
args.funct='PBE';
%args.temperature=300;
[sys,options]=feval([name,'_setup'],args);
%sys.nbnd = 30;
options.exxmethod='kmeans';
options.aceconv=true;

[sys, H, X, info] = scf(sys,options);

isdfoptions.rank=sys.nel*16;
isdfoptions.weight='power';
isdfoptions.power=1;
isdfoptions.init='random';
isdfoptions.seed=0;
isdfoptions.sys=sys;
options.isdfoptions=isdfoptions;

options.X0=X;
options.rho0=H.rho;
sys.funct='HSE06';
[sys, H, X, info] = scf(sys,options);
rho=H.rho;
save('kms','sys','X','rho','options');

kpts = {0,0,0,21,'G';
        -0.5,0,0,11,'X';
        %0,0.5,0,11,'Y'
    };
options.maxphiiter=4;
[sys,ebands,BH,BX]=eband(sys,options,kpts,rho,X);
save_bands(ebands,sys.efermi,kpts,sys);
