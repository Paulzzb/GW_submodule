function [sys,options]=cu4_setup(args)
    a1 = Atom('Cu');
    atomlist = [ a1, a1, a1, a1 ];
    C = 6.8308*eye(3);
    coefs = [
         0.000000000         0.000000000         0.000000000
         0.000000000         0.500000000         0.500000000
         0.500000000         0.000000000         0.500000000
         0.500000000         0.500000000         0.000000000
    ];
    xyzlist = coefs*C';
    name='cu4';
    nk=0;
    ng=[];
    ecut=10;
    funct='PBE';
	temperature=0;
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
    if isfield(args,'temperature')
        temperature=args.temperature;
    end
    if nk==0
        sys = Molecule('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
            'ecut',ecut, 'name',name,'funct',funct,'temperature',temperature);
    else
        if numel(nk)==1 & nk<0
            [kpts,wks]=getk(nk);
            sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
                'ecut',ecut, 'name',name, 'funct',funct,'kpts',kpts, 'wks',wks,'temperature',temperature);
        else
            if numel(nk)==1
                nk=[nk,nk,nk];
            end
            sys = Crystal('supercell',C,'n1',ng,'n2',ng,'n3',ng, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
                'ecut',ecut, 'name',name, 'autokpts', nk,'funct',funct,'temperature',temperature);
        end
    end
    options = setksopt();
    options.scftol = 1e-7;
    options.phitol = 1e-7;
    options.maxscfiter=40;
    options.maxphiiter=20;
end
