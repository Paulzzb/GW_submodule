function [sys,options]=sih4_20_setup(args)
	a1 = Atom('Si');
	a2 = Atom('H');
	atomlist = [a1 a2 a2 a2 a2];
	C = 20*eye(3);
	coefs = [
		 0.0      0.0       0.0
		 0.0805   0.0805    0.0805
		-0.0805  -0.0805    0.0805
		 0.0805  -0.0805   -0.0805
		-0.0805   0.0805   -0.0805
	];
	xyzlist = coefs*C';

	args.C=C;
	args.atomlist=atomlist;
	args.xyzlist=xyzlist;
	args.name='sih4';

	[sys,options]=structure_setup(args);
end
