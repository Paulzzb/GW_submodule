function [kpts2,idxnz]=kgrid(sys)
    %Generate Fourier convolution kgrids for exxgkk
    nks=sys.nkxyz;
    nkx=nks(1);
    nky=nks(2);
    nkz=nks(3);

    [I,J,K] = ndgrid((0:nkx-1)+((0:nkx-1) >= nkx/2)*(nkx-1), ...
        (0:nky-1)+((0:nky-1) >= nky/2)*(nky-1), ...
        (0:nkz-1)+((0:nkz-1) >= nkz/2)*(nkz-1));
    iks=reshape(cat(4,I,J,K),[],3);
    idxnz=sum(iks.*[1,2*nks(1)-1,(2*nks(1)-1)*(2*nks(2)-1)],2)+1;

    nkx=2*nkx-1;
    nky=2*nky-1;
    nkz=2*nkz-1;
    [I,J,K] = ndgrid((0:nkx-1)-((0:nkx-1) >= nkx/2)*nkx, ...
        (0:nky-1)-((0:nky-1) >= nky/2)*nky, ...
        (0:nkz-1)-((0:nkz-1) >= nkz/2)*nkz);
    kpts2=2*pi*reshape(cat(4,I,J,K),[],3)./nks/(sys.supercell.');
return
