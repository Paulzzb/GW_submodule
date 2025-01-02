function [sys,options]=si_setup(args)
	a1 = Atom('Si');
	ncell=args.ncell;
	N=prod(ncell);
	atomlist = zeros(1,N*8,'Atom');
	for i = 1:N*8
		atomlist(i)=a1;
	end
	C = 10.216*diag(ncell);
	coefs = [
	    0.000000000         0.000000000         0.000000000
	    0.000000000         0.500000000         0.500000000
	    0.500000000         0.000000000         0.500000000
	    0.500000000         0.500000000         0.000000000
	    0.750000000         0.250000000         0.750000000
	    0.250000000         0.250000000         0.250000000
	    0.250000000         0.750000000         0.750000000
	    0.750000000         0.750000000         0.250000000
	];
	for i = 1:3
		pos = coefs;
		v = [zeros(1,i-1),1,zeros(1,3-i)];
        u = [ones(1,i-1),ncell(i),ones(1,3-i)];
		for j=1:ncell(i)-1
			pos = [pos;coefs+v*j];
		end
		coefs = pos./u;
    end
	xyzlist = coefs*C';
	name=['si',num2str(N*8)];
	nk=0;
	ng=[];
	ecut=12.5;
	funct='PBE';
	if isfield(args,'nk')
		nk=args.nk;
	end
	if isfield(args,'ng')
		ng=args.ng;
	end
	if isfield(args,'ecut')
		ecut=args.ecut;
	end
	if isfield(args,'funct')
		funct=args.funct;
	end
	if nk==0
		sys = Molecule('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
			'ecut',ecut, 'name',name,'funct',funct);
	else
		if numel(nk)==1 & nk<0
			[kpts,wks]=getk(nk);
	    	sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    	    	'ecut',ecut, 'name',name, 'funct',funct,'kpts',kpts, 'wks',wks);
		else
			if numel(nk)==1
				nk=[nk,nk,nk];
			end
	    	sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    	    	'ecut',ecut, 'name',name, 'autokpts', nk,'funct',funct);
		end
	end
	options = setksopt();
	options.scftol = 1e-7;
	options.phitol = 1e-7;
	options.maxscfiter=25;
	options.maxphiiter=20;
end
