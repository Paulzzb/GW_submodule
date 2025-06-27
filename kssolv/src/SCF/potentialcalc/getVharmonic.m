function [vHarmonic] = getVharmonic(mol,w)
% GETVHARMONIC Computes a Harmonic external potential.
%    [vHarmonic] = getvHarmonic(mol,w) computes a Harmonic external 
%    potential by V = (1/2) w**2 * (r-r0)**2 where r0 is the center 
%    of the cell.
%
%    mol       --- a Molecule object
%    w         --- osciallator strength

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

C = mol.supercell;
r0 = 0.5*(C(:,1) + C(:,2) + C(:,3));

n1  = mol.n1;
n2  = mol.n2;
n3  = mol.n3;
vol = mol.vol;

% Harmonic potential in real space.
vHarmonic = [];
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            r = (i1-1)*C(:,1)/n1 + (i2-1)*C(:,2)/n2 + (i3-1)*C(:,3)/n3;
            dr = r-r0;
            dr2 = norm(dr);
	    vHarmonic(i1,i2,i3) = 0.5*w*w*dr2*dr2;
        end
    end
end            

end
