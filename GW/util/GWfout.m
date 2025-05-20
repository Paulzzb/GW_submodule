function GWfout(freq_dep, fout, n, Eo, Ex, Eres, Eint, Sig, Vxc, Eqp)
    % 打开文件
    fid = fopen(fout, 'w');
    if fid == -1
        error('无法打开文件 %s 进行写入。', fout);
    end
    if (freq_dep == 1); freq_dep = 0; end
    switch freq_dep
      case 0
        fprintf(fid, '   n         Eo           X        SX-X          CH         Sig         Vxc        Eqp0\n');
        for i = 1:length(n)
          fprintf(fid, '%6d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n', ...
                 n(i), Eo(i), Ex(i), Eres(i), Eint(i), Sig(i), Vxc(i), Eqp(i));
        end
      case 2
        % 打印标题行（可选）
        fprintf(fid, '     n         Eo           X      Re Res      Re Int      Re Sig        Vxc     Re Eqp0  \n');
        fprintf(fid, '                                   Im Res      Im Int      Im Sig                Im Eqp0  \n');

        % 遍历每个结果
        for i = 1:length(n)
            % 实部数据输出
            fprintf(fid, '%6d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f$\n', ...
                n(i), Eo(i), Ex(i), real(Eres(i)), real(Eint(i)), real(Sig(i)), Vxc(i), real(Eqp(i)));
            
            % 虚部数据输出（右对齐）
            fprintf(fid, '                              %12.6f%12.6f%12.6f             %12.6f$\n', ...
                imag(Eres(i)), imag(Eint(i)), imag(Sig(i)), imag(Eqp(i)));
        end

        % 关闭文件
        fclose(fid);
      otherwise
        error('freq_dep is expected to be 0, 1, or 2.');
    end
        
end % function
