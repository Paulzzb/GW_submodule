function CoulombG(mill, supercell, amin, truncation)

% Calculate the coulomb potential
ng = size(mill, 1);
coulG = zeros(ng, 1);

Creci = 2 * pi * inv(supercell)';
kkxyz = mill * Creci;
gkk = sum(kkxyz.^2, 2);

switch truncation
	case 2 % spherical_truncation
    for j = 1:ng
      if ( abs(gkk(j)) ~= 0 )
        coulG(j) = 8.0*pi/(gkk(j));
        coulG(j) = coulG(j,4)*(1-cos(sqrt(gkk(j))*amin));
      else
        %coulG(j) = 4.0*pi*amin^2/2;
        coulG(j) = 0.0;
      end
    end 
  case 6 %cell_slab_truncation
    zc = supercell(3,3) / 2;
    for j = 1:ng
      if ( abs(gkk(j)) ~= 0 )
          coulG(j,4) = 8.0*pi/(gkk(j));
          gkkxy = sqrt(grid.gkx(j)^2+grid.gky(j)^2);
    %      grid.gkz
          coulG(j,4) = coulG(j,4)*(1-exp(-gkkxy*zc)*cos(grid.gkz(j)*zc));
      else
          %coulG(j) = 4.0*pi*amin^2/2;
          coulG(j,4) = 0.0;
      end;
    end;
  otherwise
	  fprintf('options.coulomb_truncation = %d is not supported.\n', ...
		         options.coulomb_truncation);
    error();
end;



end % EOF