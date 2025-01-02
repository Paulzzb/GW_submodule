function [sys,options]=structure_setup(args)
	if ~(isfield(args,'ng1') & isfield(args,'ng2') & isfield(args,'ng3'))
		if isfield(args,'ng')
        	args.ng1=args.ng;
        	args.ng2=args.ng;
        	args.ng3=args.ng;
		else
        	args.ng1=[];
        	args.ng2=[];
        	args.ng3=[];
    end
	if ~isfield(args,'nk')
        args.nk=0;
    end
	if ~isfield(args,'ecut')
        args.ecut=10;
    end
	if ~isfield(args,'funct')
        args.funct='PBE';
    end
	if args.nk==0
		sys = Molecule('supercell',args.C,'n1',args.ng1,'n2',args.ng2,'n3',args.ng3, 'atomlist',args.atomlist, 'xyzlist' ,args.xyzlist, 'ecut',args.ecut, 'name',args.name,'funct',args.funct);
	else
		if numel(args.nk)==1
			args.nk=[args.nk,args.nk,args.nk];
		end
		sys = Crystal('supercell',args.C,'n1',args.ng1,'n2',args.ng2,'n3',args.ng3, 'atomlist',args.atomlist, 'xyzlist' ,args.xyzlist, 'ecut',args.ecut, 'name',args.name, 'autokpts', args.nk,'funct',args.funct);
	end
	options = setksopt();
    options.scftol = 1e-7;
    options.phitol = 1e-7;
    options.maxscfiter=25;
    options.maxphiiter=20;
end
