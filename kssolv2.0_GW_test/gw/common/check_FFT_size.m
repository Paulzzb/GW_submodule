function check = check_FFT_size(Nfft, Nfac)
% Best FFT grid dimension is given by 2^a*3^b*5^c*7^d*11^e*13^f
% where a,b,c,d are arbitrary and e,f are 0 or 1
% Ref: http://www.fftw.org/fftw2_doc/fftw_3.html
% On entry
% Nfft = FFT grid dimension to test
% Nfac = number of factors to test
% On exit
% check_FFT_size = .true. if good FFT grid dimension
maxfac = 6;
fac(1 : maxfac) = [2 3 5 7 11 13];
if ( (Nfft < 1) || (Nfac < 1) || (Nfac > maxfac))
    error('check_FFT_size input')
end

remainder = Nfft;
for ifac = 1 : maxfac
    pow(ifac) = 0;
end
%% factoring
for ifac = 1 : Nfac
    maxpow = fix(log(remainder) / log(fac(ifac))) + 1;
    for ipow = 1 : maxpow
        if (mod(remainder, fac(ifac)) == 0)
            remainder = remainder / fac(ifac);
            pow(ifac) = pow(ifac) + 1;
        end
    end
end
%% check
product = remainder;
for ifac = 1 : Nfac
    for ipow = 1 : pow(ifac)
        product = product * fac(ifac);
    end
end
if (product ~= Nfft)
    error('Internal error in check_FFT_size; factorization failed')
end
%% output: check or not, if need to check return to 0
check = ((remainder == 1) && (pow(5) <= 1) && (pow(5) <=1));
