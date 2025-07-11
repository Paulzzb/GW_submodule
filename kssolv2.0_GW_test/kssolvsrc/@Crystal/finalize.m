function cry = finalize(cry)
% CRYSTAL/FINALIZE Finalize function for crystal class
%    cry = FINALIZE(cry) returns a crystal class of with finalized fields.
%
%    See also Crystal.

cry = finalize@Molecule(cry);

if isempty(cry.kpts)
    cry.kpts = [0,0,0];
else
    % NOTE: the transpose is very important.
    %Creci = 2*pi*inv(cry.supercell)';
    %cry.kpts = cry.kpts*Creci;
end

if numel(cry.wks) ~= cry.nkpts
    cry.wks = ones(cry.nkpts,1)/cry.nkpts;
end

%Value of reciprocal space
cry.bvec = 2*pi*inv(cry.supercell)';
cry.bdot = cry.bvec * cry.bvec.';

end
