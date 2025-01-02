function [X,ev,options] = updateX(mol, H, X, prec, options)
%
% usage: [X,ev] = updateX(mol, H, X, prec, options);
%
% purupse: update the wavefunctions by computing the invariant
%          subspace associated with the lowest eigenvalues of H
%
global version;
%global time;
eigmethod = options.eigmethod;
verbose   = options.verbose;
if (any(strcmp(options.verbose,{'off';'OFF';'Off'})))
   verbose = 0;
else
   verbose = 1;
end
ncol = ncols(X);

switch lower(eigmethod) 
   case {'lobpcg'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      if version == 'CPU'
%         if time == 'on'
            t1 = clock;
	    [X, ev, lvec, rvec] = lobpcg_cpu(H, X, prec, cgtol, maxcgiter, verbose);
%	 else
%	    t1 = clock;
%            [X, ev, lvec, rvec] = lobpcg_cpu(H, X, prec, cgtol, maxcgiter,verbose);
%            timeforLobpcgCPU = etime(clock,t1)
%         end
      elseif version == 'GPU'
	 t1 = clock;
	 [X, ev, lvec, rvec] = lobpcg_gpu(H, X, prec, cgtol, maxcgiter,verbose);
         timeforLobpcgGPU = etime(clock,t1)
      end
      
   case {'lobpcg2'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      [X, ev, lvec, rvec] = lobpcg2(H, X, prec, cgtol, maxcgiter,verbose);
   case {'davidson'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      if version == 'CPU'
	 t2 = clock;
	 [X, ev] = davidson_cpu(H, X, prec, cgtol, maxcgiter,verbose);
	 timeforDavidsonCPU = etime(clock,t2)
      elseif version == 'GPU'
         t2 = clock;
         [X, ev] = davidson_gpu(H, X, prec, cgtol, maxcgiter,verbose);
         timeforDavidsonGPU = etime(clock,t2)
      end
   case {'davidson2'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      if version == 'CPU'
         t3 = clock;
         [X, ev] = davidson2_cpu(H, X, prec, cgtol, maxcgiter,verbose);	 
	 timeforDavidson2CPU = etime(clock,t3)
      elseif version == 'GPU'
	 t3 = clock;
	 [X, ev] = davidson2_gpu(H, X, prec, cgtol, maxcgiter,verbose);
	 timeforDavidson2GPU = etime(clock,t3)
      end
   case {'davidson_qe'}
       cgtol = options.cgtol;
       maxcgiter = options.maxcgiter;
       if ~isfield(options,'h_diag')||~isfield(options,'s_diag')
           t1 = tic;
           npol = 1 + mol.noncolin;
           vion_g = fft3(H.vion)/mol.n1/mol.n2/mol.n3;
           h_diag = repmat(H.gkin,npol,1) + vion_g(1);
           [h_diag, s_diag] = addnlc(mol,h_diag,H.vnlsign,H.vnlmat);
           options.h_diag = h_diag;
           options.s_diag = s_diag;
           % fprintf('Precondition time = %20.4e\n', toc(t1));
       end
       [X,ev] = davidson_qe(mol, H,  X, cgtol, maxcgiter,options.iterscf,options.h_diag,options.s_diag);
   case {'davpcg'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      [X, ev] = davpcg(H, X, prec, cgtol, maxcgiter,verbose);
   case {'eigs'}
      eigstol = options.eigstol;
      maxeigsiter = options.maxeigsiter;
      [X, ev] = diagbyeigs(mol, H, ncol, eigstol, maxeigsiter);
   case {'chebyfilt'}
      v0 = genX0(mol,1);
      degree = options.degree;
      if version == 'CPU'
         [T,Q,f]=lanczos_cpu(H,v0,2*ncol);
	 d = sort(real(eig(T)));
	 lb = d(ncol+1);
	 ub = 1.01*d(2*ncol);
	 fprintf('lb = %11.3e, ub = %11.3e, degree = %d\n', ...
	 lb, ub, degree);
	 t4 = clock;
	 Y = chebyfilt_cpu(H,X,degree,lb,ub);
	 timeforChebyfiltCPU = etime(clock,t4)
	 [X,R]=qr(Y,0);
      elseif version == 'GPU'
	 [T,Q,f]=lanczos_gpu(H,v0,2*ncol);
	 d = sort(real(eig(T)));
	 lb = d(ncol+1);
	 ub = 1.01*d(2*ncol);
	 fprintf('lb = %11.3e, ub = %11.3e, degree = %d\n', ...
	 lb, ub, degree);
	 t4 = clock;
 	 Y = chebyfilt_gpu(H,X,degree,lb,ub);
	 timeforChebyfiltGPU = etime(clock,t4)
	 [X,R]=qr(Y,0);
      end
   case {'omm'}
      ommtol = options.cgtol;
      maxommiter = options.maxcgiter;
      %
      % run a few Lanczos iteration to get an upper bound of the spectrum
      %
      v0 = genX0(mol,1);
      maxlan = 10;
      if version == 'CPU'
         [T,V,f]=lanczos_cpu(H,v0,maxlan);
         T = (T+T')/2;
         [S,D]=eig(T);
         z = V*S(:,maxlan);
         eigub = D(maxlan,maxlan);
         eigub = eigub + norm(H*z-eigub*z);
         fprintf('upper bound of the spectrum = %11.3e\n', eigub);
      %
         t5 = clock;
         [X, ev, lvec, rvec] = omm_cpu(H, X, eigub, prec, ommtol, maxommiter,verbose);
         timeforOmmCPU = etime(clock,t5)
      elseif version == 'GPU'
	 [T,V,f]=lanczos_gpu(H,v0,maxlan);
	 T = (T+T')/2;
	 [S,D]=eig(T);
	 z = V*S(:,maxlan);
	 eigub = D(maxlan,maxlan);
	 eigub = eigub + norm(H*z-eigub*z);
	 fprintf('upper bound of the spectrum = %11.3e\n', eigub);
	 %
	 t5 = clock;
	 [X, ev, lvec, rvec] = omm_gpu(H, X, eigub, prec, ommtol, maxommiter,verbose);
	 timeforOmmGPU = etime(clock,t5)
      end
   case {'ppcg'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      t6 = clock;
      [X, ev, lvec, rvec] = ppcg(H, X, prec, cgtol, maxcgiter,verbose);
      timeforPPCG = etime(clock,t6)
   otherwise
      disp('Unknown method for diagonalizing H! Use eigs');
      [X, ev] = diagbyeigs(mol, H, ncol, eigstol, maxeigsiter);
end
%
if ( verbose ==1 )
   HX = H*X;
   G = X'*HX;
   R = HX-X*G;
   ev = sort(real(eig(G)));
   for j = 1:ncol
      resnrm(j) = norm(R(:,j));
      %fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', ...
      %        j, ev(j), resnrm(j));
   end
end
