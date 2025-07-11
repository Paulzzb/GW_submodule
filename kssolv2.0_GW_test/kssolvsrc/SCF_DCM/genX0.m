function X0 = genX0(Obj,nXcols)
% GENX0 generates the initial wave functions.
%    X0 = GENX0(Obj,ncols) generates the initial wave functions for
%    either a Molecule/Crystal or Ham/BlochHam object with nXcols columns
%    per k-point;
%
%   See also Molecule, Crystal, Wavefun, BlochWavefun.

if nargin == 1
    if isa( Obj, 'Crystal' )
        nXcols = Obj.nel/2*Obj.nspin*ones(Obj.nkpts,1);
    else
        nXcols = Obj.nel/2*Obj.nspin;
    end
end

if isa( Obj, 'Ham') || isa(Obj,'BlochHam')
   if nargin < 2
      error('Need the number of vectors, if the first argument is an Ham object');
   end

   if isa(Obj, 'BlochHam')  
      nXcols = nXcols*ones(Obj.nkpts,1);
   end
end

nkpts = numel(nXcols);

n1   = Obj.n1;
n2   = Obj.n2;
n3   = Obj.n3;

if  isa(Obj,'Molecule') || isa(Obj,'Crystal')
    idxnz = Obj.gridwfc.idxnz;
else
    idxnz = Obj.idxnz;
end;

Qcell = cell(nkpts,1);

for ik = 1:nkpts
    psir = randn(n1,n2,n3,nXcols(ik));
    psif = reshape(fft3(psir),n1*n2*n3,[]);
    psif = psif(idxnz,:);
    [Qcell{ik},~]=qr(psif,0);
end

if isa( Obj, 'Crystal' ) || isa( Obj, 'BlochHam')
    X0 = BlochWavefun(Qcell,n1,n2,n3,idxnz,Obj.wks);
else
    X0 = Wavefun(Qcell{1},n1,n2,n3,idxnz);
end

end
