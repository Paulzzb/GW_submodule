function [sys,options]=benzene_setup(args)
	a1 = Atom('H');
	a2 = Atom('C');
	atomlist = [a1 a1 a1 a1 a1 a1 a2 a2 a2 a2 a2 a2];
	C = 22.676722*eye(3);
	coefs = [
		 0.292420000		 0.500000000		 0.500000000
		 0.396239996		 0.320270002		 0.500000000
		 0.603760004		 0.320270002		 0.500000000
		 0.707579970		 0.500000000		 0.500000000
		 0.603760004		 0.679729998		 0.500000000
		 0.396239996		 0.679729998		 0.500000000
	
		 0.383379996		 0.500000000		 0.500000000
		 0.441709995		 0.399040014		 0.500000000
		 0.558290005		 0.399040014		 0.500000000
		 0.616620004		 0.500000000		 0.500000000
		 0.558290005		 0.600960016		 0.500000000
		 0.441709995		 0.600960016		 0.500000000
	];
	xyzlist = coefs*C';
	name='benzene';
	nk=0;
	ng=[];
	ecut=15;
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
		if numel(nk)==1
			nk=[nk,nk,nk];
		end
		sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
			'ecut',ecut, 'name',name, 'autokpts', nk,'funct',funct);
	end
	options = setksopt();
    options.scftol = 1e-7;
    options.phitol = 1e-7;
    options.maxscfiter=25;
    options.maxphiiter=20;
end

%qe HSE result
%total energy -74.77989878 Ry
%Fock energy   -4.81897599 Ry
