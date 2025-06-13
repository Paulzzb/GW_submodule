function GWfout(qpenergy)
  
  freq_dep = qpenergy.freq_dep;
  fout = qpenergy.fout;
  
  bandindices = qpenergy.bandindices;
  ev = qpenergy.ev;
  Ex = qpenergy.Ex;
  Esex_x = qpenergy.Esex_x;
  Ecoh = qpenergy.Ecoh;
  Eres = qpenergy.Eres;
  Eint = qpenergy.Eint;
  Sig = qpenergy.Sig;
  Vxc = qpenergy.Vxc;
  Eqp = qpenergy.Eqp;
  
  % open file
  fid = fopen(fout, 'w');
  if fid == -1
      error('Can not open %s to write in.', fout);
  end
  if (freq_dep == 1); freq_dep = 0; end
  switch freq_dep
    case 0
      fprintf(fid, '   n         Eo           X        SX-X          CH         Sig         Vxc        Eqp0\n');
      for i = 1:length(bandindices)
        fprintf(fid, '%6d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n', ...
               bandindices(i), ev(i), Ex(i), Esex_x(i), Ecoh(i), Sig(i), Vxc(i), Eqp(i));
      end
    case 2
      % 打印标题行（可选）
      fprintf(fid, '     n         Eo           X      Re Res      Re Int      Re Sig        Vxc     Re Eqp0  \n');
      fprintf(fid, '                                   Im Res      Im Int      Im Sig                Im Eqp0  \n');

      for i = 1:length(bandindices)
          fprintf(fid, '%6d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f \n', ...
              bandindices(i), ev(i), Ex(i), real(Eres(i)), real(Eint(i)), real(Sig(i)), Vxc(i), real(Eqp(i)));
          
          fprintf(fid, '                              %12.6f%12.6f%12.6f             %12.6f \n', ...
              imag(Eres(i)), imag(Eint(i)), imag(Sig(i)), imag(Eqp(i)));
      end
      fclose(fid);
    otherwise
      error('freq_dep is expected to be 0, 1, or 2.');
  end
        
end % function
