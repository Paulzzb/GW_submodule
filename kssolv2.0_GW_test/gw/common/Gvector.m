classdef Gvector
% convert the order of wave vectors so that they correspond to those in QE
% (ekin order)
    properties (SetAccess = protected)
        index_vec
        mill
        nr
        ng
        n1
        n2
        n3
        ecut
        
        
        
    end
    methods
        function gvector = Gvector(sigrid,sys)
            
            
            if nargin == 0
                return;
            end
            if isa(sigrid,'Ggrid') & isa(sys,'Crystal')
                
                sigkx=sigrid.gkx;
                sigky=sigrid.gky;
                sigkz=sigrid.gkz;
                mill=[sigkx sigky sigkz]/(2 * pi * inv(sys.supercell).');
                siG= round(mill);
                FFTgrid=[sys.n1 sys.n2 sys.n3];
                Gkk=sum(siG.*siG,2);
                n1=size(siG,1);
                for i=1:n1
                    GK(i,1)=siG(i,3)+sys.n3*(siG(i,2)+sys.n2*siG(i,1));
                end
                G=[siG Gkk GK];
                A=sortrows(G,[4 5]);
                A=A(:,1:3);
                index_vec=zeros(size(sigrid.kkxyz,1),1);
                if  (2*max(abs(A(:,1)))> sys.n1) | (2*max(abs(A(:,2)))> sys.n2) ...
                        | (2*max(abs(A(:,3)))> sys.n3 )
                    warning("gvectors must be in the interval [-FFTgrid/2, FFTgrid/2]")
                end
                for ig=1:sigrid.ng%%%Compute index_vec indices relating G-vectors in reduced coordinates to positions in the FFT grid，index_vec(FFT排序address)会给
                    iadd=round(((A(ig,1)+ floor(FFTgrid(1)/2)) *FFTgrid(2)+A(ig,2)+ floor(FFTgrid(2)/2)) *FFTgrid(3)+...
                        A(ig,3) + floor(FFTgrid(3)/2) +1);
                    index_vec(iadd,1)=ig;
                end
            end
            
            gvector.index_vec=index_vec;
            gvector.mill=A;
            gvector.ecut=sigrid.ecut;
            gvector.ng=sigrid.ng;
            gvector.nr=size(sigrid.kkxyz,1);
            gvector.n1=sys.n1;
            gvector.n2=sys.n2;
            gvector.n3=sys.n3;
        end
        
    end
end

