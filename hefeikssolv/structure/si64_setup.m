function [sys,options]=si8_setup(args)
	a1 = Atom('Si');
	atomlist = [ ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1 ...
	];
	C = 20.5250*eye(3);
	coefs = [
	     0.000000000         0.000000000         0.000000000
	     0.000000000         0.250000000         0.250000000
	     0.250000000         0.000000000         0.250000000
	     0.250000000         0.250000000         0.000000000
	     0.375000000         0.125000000         0.375000000
	     0.125000000         0.125000000         0.125000000
	     0.125000000         0.375000000         0.375000000
	     0.375000000         0.375000000         0.125000000
	     0.500000000         0.000000000         0.000000000
	     0.500000000         0.250000000         0.250000000
	     0.750000000         0.000000000         0.250000000
	     0.750000000         0.250000000         0.000000000
	     0.875000000         0.125000000         0.375000000
	     0.625000000         0.125000000         0.125000000
	     0.625000000         0.375000000         0.375000000
	     0.875000000         0.375000000         0.125000000
	     0.000000000         0.500000000         0.000000000
	     0.000000000         0.750000000         0.250000000
	     0.250000000         0.500000000         0.250000000
	     0.250000000         0.750000000         0.000000000
	     0.375000000         0.625000000         0.375000000
	     0.125000000         0.625000000         0.125000000
	     0.125000000         0.875000000         0.375000000
	     0.375000000         0.875000000         0.125000000
	     0.500000000         0.500000000         0.000000000
	     0.500000000         0.750000000         0.250000000
	     0.750000000         0.500000000         0.250000000
	     0.750000000         0.750000000         0.000000000
	     0.875000000         0.625000000         0.375000000
	     0.625000000         0.625000000         0.125000000
	     0.625000000         0.875000000         0.375000000
	     0.875000000         0.875000000         0.125000000
	     0.000000000         0.000000000         0.500000000
	     0.000000000         0.250000000         0.750000000
	     0.250000000         0.000000000         0.750000000
	     0.250000000         0.250000000         0.500000000
	     0.375000000         0.125000000         0.875000000
	     0.125000000         0.125000000         0.625000000
	     0.125000000         0.375000000         0.875000000
	     0.375000000         0.375000000         0.625000000
	     0.500000000         0.000000000         0.500000000
	     0.500000000         0.250000000         0.750000000
	     0.750000000         0.000000000         0.750000000
	     0.750000000         0.250000000         0.500000000
	     0.875000000         0.125000000         0.875000000
	     0.625000000         0.125000000         0.625000000
	     0.625000000         0.375000000         0.875000000
	     0.875000000         0.375000000         0.625000000
	     0.000000000         0.500000000         0.500000000
	     0.000000000         0.750000000         0.750000000
	     0.250000000         0.500000000         0.750000000
	     0.250000000         0.750000000         0.500000000
	     0.375000000         0.625000000         0.875000000
	     0.125000000         0.625000000         0.625000000
	     0.125000000         0.875000000         0.875000000
	     0.375000000         0.875000000         0.625000000
	     0.500000000         0.500000000         0.500000000
	     0.500000000         0.750000000         0.750000000
	     0.750000000         0.500000000         0.750000000
	     0.750000000         0.750000000         0.500000000
	     0.875000000         0.625000000         0.875000000
	     0.625000000         0.625000000         0.625000000
	     0.625000000         0.875000000         0.875000000
	     0.875000000         0.875000000         0.625000000
	];
	xyzlist = coefs*C';
	name='si64';
	nk=0;
	ng=[];
	ecut=12.5;
    if isfield(args,'nk')
        nk=args.nk;
    end
    if isfield(args,'ng')
        ng=args.ng;
    end
    if isfield(args,'ecut')
        ecut=args.ecut;
    end
    if nk==0
        sys = Molecule('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
            'ecut',ecut, 'name',name,'funct','HSE06');
    else
        if numel(nk)==1 & nk<0
            [kpts,wks]=getk(nk);
            sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
                'ecut',ecut, 'name',name, 'funct','HSE06','kpts',kpts, 'wks',wks);
        else
            if numel(nk)==1
                nk=[nk,nk,nk];
            end
            sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
                'ecut',ecut, 'name',name, 'autokpts', nk,'funct','HSE06');
        end
    end
	options = setksopt();
	options.scftol = 1e-7;
	options.phitol = 1e-7;
	options.maxscfiter=25;
	options.maxphiiter=20;
	%options.betamix=0.7;
	%options.mixdim=8;
end
