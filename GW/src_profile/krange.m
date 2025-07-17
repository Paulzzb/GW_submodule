function [kpt,gpt]=krange(kpt,tol)
  % Bring kpt into the first BZ zone
  for ii=1:3
    gpt(ii)=0;
    while (kpt(ii) < -tol)
      gpt(ii)=gpt(ii)+1;
      kpt(ii)=kpt(ii)+1;
    end
    while (kpt(ii) > 1-tol)
      gpt(ii)=gpt(ii)-1;
      kpt(ii)=kpt(ii)-1;
    end
  end
end