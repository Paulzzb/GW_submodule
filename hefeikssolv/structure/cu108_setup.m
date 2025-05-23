function [sys,options]=cu32_setup(args)
	%kssolvpptype('pz-hgh');
    a1 = Atom('Cu');
	atomlist = [ ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1, a1, a1, ...
	a1, a1, a1, a1, a1, a1, a1, a1...
	];
	C = 20.4923778*eye(3);
	coefs = [
	     0.000000000         0.000000000         0.000000000
	     0.000000000         0.166669995         0.166669995
	     0.166669995         0.000000000         0.166669995
	     0.166669995         0.166669995         0.000000000
	     0.333330005         0.000000000         0.000000000
	     0.333330005         0.166669995         0.166669995
	     0.500000000         0.000000000         0.166669995
	     0.500000000         0.166669995         0.000000000
	     0.666670024         0.000000000         0.000000000
	     0.666670024         0.166669995         0.166669995
	     0.833329976         0.000000000         0.166669995
	     0.833329976         0.166669995         0.000000000
	     0.000000000         0.333330005         0.000000000
	     0.000000000         0.500000000         0.166669995
	     0.166669995         0.333330005         0.166669995
	     0.166669995         0.500000000         0.000000000
	     0.333330005         0.333330005         0.000000000
	     0.333330005         0.500000000         0.166669995
	     0.500000000         0.333330005         0.166669995
	     0.500000000         0.500000000         0.000000000
	     0.666670024         0.333330005         0.000000000
	     0.666670024         0.500000000         0.166669995
	     0.833329976         0.333330005         0.166669995
	     0.833329976         0.500000000         0.000000000
	     0.000000000         0.666670024         0.000000000
	     0.000000000         0.833329976         0.166669995
	     0.166669995         0.666670024         0.166669995
	     0.166669995         0.833329976         0.000000000
	     0.333330005         0.666670024         0.000000000
	     0.333330005         0.833329976         0.166669995
	     0.500000000         0.666670024         0.166669995
	     0.500000000         0.833329976         0.000000000
	     0.666670024         0.666670024         0.000000000
	     0.666670024         0.833329976         0.166669995
	     0.833329976         0.666670024         0.166669995
	     0.833329976         0.833329976         0.000000000
	     0.000000000         0.000000000         0.333330005
	     0.000000000         0.166669995         0.500000000
	     0.166669995         0.000000000         0.500000000
	     0.166669995         0.166669995         0.333330005
	     0.333330005         0.000000000         0.333330005
	     0.333330005         0.166669995         0.500000000
	     0.500000000         0.000000000         0.500000000
	     0.500000000         0.166669995         0.333330005
	     0.666670024         0.000000000         0.333330005
	     0.666670024         0.166669995         0.500000000
	     0.833329976         0.000000000         0.500000000
	     0.833329976         0.166669995         0.333330005
	     0.000000000         0.333330005         0.333330005
	     0.000000000         0.500000000         0.500000000
	     0.166669995         0.333330005         0.500000000
	     0.166669995         0.500000000         0.333330005
	     0.333330005         0.333330005         0.333330005
	     0.333330005         0.500000000         0.500000000
	     0.500000000         0.333330005         0.500000000
	     0.500000000         0.500000000         0.333330005
	     0.666670024         0.333330005         0.333330005
	     0.666670024         0.500000000         0.500000000
	     0.833329976         0.333330005         0.500000000
	     0.833329976         0.500000000         0.333330005
	     0.000000000         0.666670024         0.333330005
	     0.000000000         0.833329976         0.500000000
	     0.166669995         0.666670024         0.500000000
	     0.166669995         0.833329976         0.333330005
	     0.333330005         0.666670024         0.333330005
	     0.333330005         0.833329976         0.500000000
	     0.500000000         0.666670024         0.500000000
	     0.500000000         0.833329976         0.333330005
	     0.666670024         0.666670024         0.333330005
	     0.666670024         0.833329976         0.500000000
	     0.833329976         0.666670024         0.500000000
	     0.833329976         0.833329976         0.333330005
	     0.000000000         0.000000000         0.666670024
	     0.000000000         0.166669995         0.833329976
	     0.166669995         0.000000000         0.833329976
	     0.166669995         0.166669995         0.666670024
	     0.333330005         0.000000000         0.666670024
	     0.333330005         0.166669995         0.833329976
	     0.500000000         0.000000000         0.833329976
	     0.500000000         0.166669995         0.666670024
	     0.666670024         0.000000000         0.666670024
	     0.666670024         0.166669995         0.833329976
	     0.833329976         0.000000000         0.833329976
	     0.833329976         0.166669995         0.666670024
	     0.000000000         0.333330005         0.666670024
	     0.000000000         0.500000000         0.833329976
	     0.166669995         0.333330005         0.833329976
	     0.166669995         0.500000000         0.666670024
	     0.333330005         0.333330005         0.666670024
	     0.333330005         0.500000000         0.833329976
	     0.500000000         0.333330005         0.833329976
	     0.500000000         0.500000000         0.666670024
	     0.666670024         0.333330005         0.666670024
	     0.666670024         0.500000000         0.833329976
	     0.833329976         0.333330005         0.833329976
	     0.833329976         0.500000000         0.666670024
	     0.000000000         0.666670024         0.666670024
	     0.000000000         0.833329976         0.833329976
	     0.166669995         0.666670024         0.833329976
	     0.166669995         0.833329976         0.666670024
	     0.333330005         0.666670024         0.666670024
	     0.333330005         0.833329976         0.833329976
	     0.500000000         0.666670024         0.833329976
	     0.500000000         0.833329976         0.666670024
	     0.666670024         0.666670024         0.666670024
	     0.666670024         0.833329976         0.833329976
	     0.833329976         0.666670024         0.833329976
	     0.833329976         0.833329976         0.666670024
	];
    xyzlist = coefs*C';
    name='cu108';
    nk=0;
    ng=[];
    ecut=10;
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
    options.maxscfiter=40;
    options.maxphiiter=20;
end
