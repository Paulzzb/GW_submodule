function [Fx_wpbe, d1rfx, d1sfx] = wpbe_analy_erfc_approx_grad(rho,s,omega)
One = 1;
Two = 2;
Three = 3;
Four = 4;
Five = 5;
Six = 6;
Seven = 7;
Eight = 8;
Nine = 9;
Ten = 10;
Fifteen = 15;
Sixteen = 16;

r36=36;
r64=64;
r81=81;
r256=256;
r384=384;

r27=27;
r128=128;
r144=144;
r288=288;
r324=324;
r729=729;

r20=20;
r32=32;
r243=243;
r2187=2187;
r6561=6561;
r40=40;

r12=12;
r25=25;
r30=30;
r54=54;
r75=75;
r105=105;
r135=135;
r1215=1215;
r15309=15309;

f12    = 0.5;
f13    = One/Three;
f14    = 0.25;

f32    = 1.5;
f34    = 0.75;
f94    = 2.25;
f98    = 1.125;
f1516  = Fifteen / Sixteen;

pi     = acos(-One);
pi2    = pi*pi;
srpi   = sqrt(pi);

ea1 = -1.128223946706117;
ea2 = 1.452736265762971;
ea3 = -1.243162299390327;
ea4 = 0.971824836115601;
ea5 = -0.568861079687373;
ea6 = 0.246880514820192;
ea7 = -0.065032363850763;
ea8 = 0.008401793031216;

eb1 = 1.455915450052607;

A      =  1.0161144d0;
B      = -3.7170836d-1;
C      = -7.7215461d-2;
D      =  5.7786348d-1;
E      = -5.1955731d-2;

Ha1    = 9.79681d-3;
Ha2    = 4.10834d-2;
Ha3    = 1.87440d-1;
Ha4    = 1.20824d-3;
Ha5    = 3.47188d-2;

Fc1    = 6.4753871d0;
Fc2    = 4.7965830d-1;

EGa1   = -2.628417880d-2;
EGa2   = -7.117647788d-2;
EGa3   =  8.534541323d-2;

expei1 = 4.03640D0;
expei2 = 1.15198D0;
expei3 = 5.03627D0;
expei4 = 4.19160D0;

EGscut     = 8.0d-2;
wcutoff    = 1.4D1;
expfcutoff = 7.0D2;

xkf    = (Three*pi2.*rho).^f13;

A2 = A*A;
A3 = A2*A;
A12 = sqrt(A);
A32 = A12*A;
A52 = A32*A;

if omega == 0
    w = zeros(size(s));
else
    w = omega ./ xkf;
end
w2    = w .* w;
w3    = w2 .* w;
w4    = w2 .* w2;
w5    = w3 .* w2;
w6    = w5 .* w;
w7    = w6 .* w;
w8    = w7 .* w;

d1rw  = -One/Three./rho.*w;

X      = - Eight/Nine;

s2     = s.*s;
s3     = s2.*s;
s4     = s2.*s2;
s5     = s4.*s;
s6     = s5.*s;

Hnum    = Ha1*s2 + Ha2*s4;
Hden    = One + Ha3*s4 + Ha4*s5 + Ha5*s6;

H       = Hnum./Hden;

d1sHnum = Two*Ha1.*s + Four*Ha2.*s3;
d1sHden = Four*Ha3.*s3 + Five*Ha4.*s4 + Six*Ha5.*s5;

d1sH    = (Hden.*d1sHnum - Hnum.*d1sHden) ./ (Hden.*Hden);

F      = Fc1*H + Fc2;
d1sF   = Fc1*d1sH;

eb1 = eb1*ones(size(rho));
eb1(w>wcutoff) = 2;

Hsbw = s2.*H + eb1.*w2;
Hsbw2 = Hsbw.*Hsbw;
Hsbw3 = Hsbw2.*Hsbw;
Hsbw4 = Hsbw3.*Hsbw;
Hsbw12 = sqrt(Hsbw);
Hsbw32 = Hsbw12.*Hsbw;
Hsbw52 = Hsbw32.*Hsbw;
Hsbw72 = Hsbw52.*Hsbw;

d1sHsbw  = d1sH.*s2 + Two*s.*H;
d1rHsbw  = Two*eb1.*d1rw.*w;

DHsbw = D + s2.*H + eb1.*w2;
DHsbw2 = DHsbw.*DHsbw;
DHsbw3 = DHsbw2.*DHsbw;
DHsbw4 = DHsbw3.*DHsbw;
DHsbw5 = DHsbw4.*DHsbw;
DHsbw12 = sqrt(DHsbw);
DHsbw32 = DHsbw12.*DHsbw;
DHsbw52 = DHsbw32.*DHsbw;
DHsbw72 = DHsbw52.*DHsbw;
DHsbw92 = DHsbw72.*DHsbw;

HsbwA94   = f94/A * Hsbw;
HsbwA942  = HsbwA94.*HsbwA94;
HsbwA943  = HsbwA942.*HsbwA94;
HsbwA945  = HsbwA943.*HsbwA942;
HsbwA9412 = sqrt(HsbwA94);

DHs    = D + s2.*H;
DHs2   = DHs.*DHs;
DHs3   = DHs2.*DHs;
DHs4   = DHs3.*DHs;
DHs72  = DHs3.*sqrt(DHs);
DHs92  = DHs72.*DHs;

d1sDHs = Two*s.*H + s2.*d1sH;

DHsw   = DHs + w2;
DHsw2  = DHsw.*DHsw;
DHsw52 = sqrt(DHsw).*DHsw2;
DHsw72 = DHsw52.*DHsw;

d1rDHsw = Two*d1rw.*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_a    = One/Sixteen*srpi * (Fifteen*E + Six*C*(One+F.*s2).*DHs + Four*B*DHs2 + Eight*A*DHs3)./DHs72...
    - f34*pi*sqrt(A) * exp(f94/A*H.*s2) .* (One - qe_erf(f32/sqrt(A)*s.*sqrt(H)));
d1sG_a = One/r32*srpi * ((r36/A12*sqrt(A) * (Two*H + d1sH.*s) ./ sqrt(H))...
    + (-Eight*A*d1sDHs.*DHs3 - r105*E*d1sDHs - r30*C*d1sDHs.*DHs.*(One+s2.*F)+r12*DHs2.*(-B*d1sDHs + C*s.*(d1sF.*s + Two*F))) ./ DHs92...
    - srpi/A12*r54 * exp(f94/A*H.*s2).*s.*(Two*H + d1sH.*s) .* qe_erfc(f32/sqrt(A)*s.*sqrt(H)));

G_b    = f1516*srpi * s2 ./ DHs72;
d1sG_b = Fifteen*srpi/r32 * s .* (Four*DHs - Seven*d1sDHs.*s) ./ DHs92;

idxsl = s > EGscut;
idxss = ~idxsl;

EG    = zeros(size(s));
d1sEG = zeros(size(s));
EG(idxsl)     = - (f34*pi + G_a(idxsl)) ./ G_b(idxsl);
d1sEG(idxsl)  = (-Four*d1sG_a(idxsl).*G_b(idxsl) + d1sG_b(idxsl).*(Four*G_a(idxsl) + Three*pi)) ./ (Four*G_b(idxsl).^2);
EG(idxss)     = EGa1 + EGa2.*s2(idxss) + EGa3.*s4(idxss);
d1sEG(idxss)  = Two*EGa2.*s(idxss) + Four*EGa3.*s3(idxss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term2 = (DHs2*B + DHs*C + Two*E + DHs.*s2.*F*C + Two*s2.*EG) ./(Two*DHs3);

d1sterm2 = (-Six*d1sDHs.*(EG.*s2 + E)...
    + DHs2 .* (-d1sDHs*B + s*C.*(d1sF.*s + Two*F))...
    + Two*DHs .* (Two*EG.*s - d1sDHs*C + s2 .* (d1sEG - C*d1sDHs.*F)))...
    ./ (Two*DHs4);

term3 = - w  .* (Four*B*DHsw2 + Six*C*DHsw + Fifteen*E...
    + Six*C*DHsw.*s2.*F + Fifteen*s2.*EG) ./...
    (Eight*DHs.*DHsw52);

d1sterm3 = w .* (Two*d1sDHs.*DHsw .* (Four*B*DHsw2...
    + Six*C*DHsw + Fifteen*E...
    + Three*s2.*(Five*EG + Two*C*DHsw.*F))...
    + DHs .* (r75*d1sDHs.*(EG.*s2 + E)...
    + Four*DHsw2.*(d1sDHs*B...
    - Three*C*s.*(d1sF.*s + Two*F))...
    - Six*DHsw.*(-Three*C*d1sDHs...
    + s.*(Ten*EG + Five*d1sEG.*s...
    - Three*C*d1sDHs.*s.*F))))...
    ./ (Sixteen*DHs2.*DHsw72);

d1rterm3 = (-Two*d1rw.*DHsw .* (Four*B*DHsw2...
    + Six*C*DHsw + Fifteen*E...
    + Three*s2.*(Five*EG + Two*C*DHsw.*F))...
    + w .* d1rDHsw .* (r75*(EG.*s2 + E)...
    + Two*DHsw.*(Two*B*DHsw + Nine*C...
    + Nine*C*s2.*F)))...
    ./ (Sixteen*DHs.*DHsw72);

term4 = - w3 .* (DHsw*C + Five*E + C*DHsw.*s2.*F + Five*s2.*EG) ./...
    (Two*DHs2.*DHsw52);

d1sterm4 = (w3 .* (Four*d1sDHs.*DHsw .* (DHsw*C + Five*E...
    + s2 .* (Five*EG + C*DHsw.*F))...
    + DHs .* (r25*d1sDHs.*(EG.*s2 + E)...
    - Two*C*DHsw2.*s.*(d1sF.*s + Two*F)...
    + DHsw .* (Three*C*d1sDHs + s.*(-r20*EG...
    - Ten*d1sEG.*s...
    + Three*C*d1sDHs.*s.*F)))))...
    ./ (Four*DHs3.*DHsw72);

d1rterm4 = (w2 .* (-Six*d1rw.*DHsw .* (DHsw*C + Five*E...
    + s2 .* (Five*EG + C*DHsw.*F))...
    + w .* d1rDHsw .* (r25*(EG.*s2 + E) +...
    Three*C*DHsw.*(One + s2.*F))))...
    ./ (Four*DHs2.*DHsw72);

term5 = - w5 .* (E + s2.*EG) ./...
    (DHs3.*DHsw52);

d1sterm5 = (w5 .* (Six*d1sDHs.*DHsw.*(EG.*s2 + E)...
    + DHs .* (-Two*DHsw.*s .* (Two*EG + d1sEG.*s)...
    + Five*d1sDHs .* (EG.*s2 + E))))...
    ./ (Two*DHs4.*DHsw72);

d1rterm5 = (Five * w4 .* (EG.*s2 + E) .* (-Two*d1rw.*DHsw...
    + d1rDHsw .* w))...
    ./ (Two*DHs3.*DHsw72);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxhs = HsbwA94 < expfcutoff;
idxhl = ~idxhs;

piexperf = zeros(size(s));
expei    = zeros(size(s));
piexperf(idxhs) = pi*exp(HsbwA94(idxhs)).*qe_erfc(HsbwA9412(idxhs));
expei(idxhs)    = exp(HsbwA94(idxhs)).*(-EXPINT(1,HsbwA94(idxhs)));
piexperf(idxhl) = pi*(One/srpi./HsbwA9412(idxhl) - One/Two/sqrt(pi)./sqrt(HsbwA943(idxhl))...
    + Three/Four/sqrt(pi)./sqrt(HsbwA945(idxhl)));
expei(idxhl)    = - (HsbwA942(idxhl) + expei1*HsbwA94(idxhl) + expei2) ./ HsbwA94(idxhl) ./ (HsbwA942(idxhl)...
    + expei3*HsbwA94(idxhl) + expei4);

piexperfd1  = - Three*srpi/sqrt(A)/Two*sqrt(Hsbw)./Hsbw + Nine/(Four*A)*piexperf;
d1spiexperf = d1sHsbw.*piexperfd1;
d1rpiexperf = d1rHsbw.*piexperfd1;

expeid1  = f14*(Four./Hsbw + Nine/A*expei);
d1sexpei = d1sHsbw.*expeid1;
d1rexpei = d1rHsbw.*expeid1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxws   = w == 0;
idxws_iz  = idxws & s == 0;
idxws_nz  = idxws & s > 0;
idxwl   = w > wcutoff;
idxwm   = ~idxws & ~idxwl;

idxt10  = idxws_nz | idxwm;
t10     = zeros(size(s));
d1st10  = zeros(size(s));
d1rt10  = zeros(size(s));
t10(idxt10)    = f12*A*log(Hsbw(idxt10) ./ DHsbw(idxt10));
t10d1          = f12*A./Hsbw(idxt10) - f12*A./DHsbw(idxt10);
d1st10(idxt10) = d1sHsbw(idxt10).*t10d1;
d1rt10(idxt10) = d1rHsbw(idxt10).*t10d1;

Fx_wpbe = zeros(size(s));
d1sfx   = zeros(size(s));
d1rfx   = zeros(size(s));

% w == 0 & s == 0
Fx_wpbe(idxws_iz) = 1;
d1sfx(idxws_iz)   = 0;
d1rfx(idxws_iz)   = 0;

% w == 0 & s > 0
t1       = -f12*A*expei(idxws_nz);
d1st1    = -f12*A*d1sexpei(idxws_nz);
d1rt1    = -f12*A*d1rexpei(idxws_nz);
term1    = t1 + t10(idxws_nz);
d1sterm1 = d1st1 + d1st10(idxws_nz);
d1rterm1 = d1rt1 + d1rt10(idxws_nz);
Fx_wpbe(idxws_nz) = X * (term1 + term2(idxws_nz));
d1sfx(idxws_nz)   = X * (d1sterm1 + d1sterm2(idxws_nz));
d1rfx(idxws_nz)   = X * d1rterm1;

% w > wcutoff
term1 = -f12*A*(expei(idxwl)+log(DHsbw(idxwl)./Hsbw(idxwl)));
term1d1  = - A/Two./DHsbw(idxwl) - f98*expei(idxwl);
d1sterm1 = d1sHsbw(idxwl).*term1d1;
d1rterm1 = d1rHsbw(idxwl).*term1d1;
Fx_wpbe(idxwl)  = X * (term1 + term2(idxwl) + term3(idxwl) + term4(idxwl) + term5(idxwl));
d1sfx(idxwl)    = X * (d1sterm1 + d1sterm2(idxwl) + d1sterm3(idxwl) + d1sterm4(idxwl) + d1sterm5(idxwl));
d1rfx(idxwl)    = X * (d1rterm1 + d1rterm3(idxwl) + d1rterm4(idxwl) + d1rterm5(idxwl));

% 0 < w <= wcutoff
w = w(idxwm);
d1rw = d1rw(idxwm);

np1    = -f32*ea1*A12*w + r27*ea3/(Eight*A12)*w3(idxwm)...
    - r243*ea5/(r32*A32)*w5(idxwm) + r2187*ea7/(r128*A52)*w7(idxwm);
d1rnp1 = (-f32*ea1*A12 + r81*ea3/(Eight*A12)*w2(idxwm)...
    - r1215*ea5/(r32*A32)*w4(idxwm)...
    + r15309*ea7/(r128*A52)*w6(idxwm)).*d1rw;

np2    = -A + f94*ea2*w2(idxwm) - r81*ea4/(Sixteen*A)*w4(idxwm)...
    + r729*ea6/(r64*A2)*w6(idxwm) - r6561*ea8/(r256*A3)*w8(idxwm);
d1rnp2 = (f12*Nine*ea2*w...
    - r81*ea4*w3(idxwm)/(Four*A)...
    + r2187*ea6*w5(idxwm)/(r32*A2)...
    - r6561*ea8*w7(idxwm)/(r32*A3)).*d1rw;

t1    = f12*(np1.*piexperf(idxwm) + np2.*expei(idxwm));
d1st1 = f12*(np1.*d1spiexperf(idxwm) + np2.*d1sexpei(idxwm));
d1rt1 = f12*(d1rnp1.*piexperf(idxwm) + d1rnp2.*expei(idxwm) + np1.*d1rpiexperf(idxwm) + np2.*d1rexpei(idxwm));

d1sHsbw = d1sHsbw(idxwm);
d1rHsbw = d1rHsbw(idxwm);
Hsbw = Hsbw(idxwm);
Hsbw2 = Hsbw2(idxwm);
Hsbw3 = Hsbw3(idxwm);
Hsbw4 = Hsbw4(idxwm);
Hsbw12 = Hsbw12(idxwm);
Hsbw32 = Hsbw32(idxwm);
Hsbw52 = Hsbw52(idxwm);
Hsbw72 = Hsbw72(idxwm);

f2    = f12*ea1*srpi*A ./ DHsbw12(idxwm);
f2d1  = -ea1*srpi*A/Four ./ DHsbw32(idxwm);
d1sf2 = d1sHsbw.*f2d1;
d1rf2 = d1rHsbw.*f2d1;

f3    = f12*ea2*A ./ DHsbw(idxwm);
f3d1  = -ea2*A/Two ./ DHsbw2(idxwm);
d1sf3 = d1sHsbw.*f3d1;
d1rf3 = d1rHsbw.*f3d1;

f4    = ea3*srpi*(-f98 ./ Hsbw12 + f14*A ./ DHsbw32(idxwm));
f4d1  = ea3*srpi*(Nine/Sixteen ./ Hsbw32 - Three*A/Eight ./ DHsbw52(idxwm));
d1sf4 = d1sHsbw.*f4d1;
d1rf4 = d1rHsbw.*f4d1;

f5    = ea4/r128*(-r144./Hsbw + r64*A./DHsbw2(idxwm));
f5d1  = ea4*(f98./Hsbw2 - A./DHsbw3(idxwm));
d1sf5 = d1sHsbw.*f5d1;
d1rf5 = d1rHsbw.*f5d1;

f6    = ea5*Three*srpi/r32/A*(Three*(Nine*Hsbw-Two*A)./Hsbw32 + Four*A2./DHsbw52(idxwm));
f6d1  = ea5*srpi*(r27/r32./Hsbw52 - r81/r64/A./Hsbw32 - Fifteen/Sixteen*A./DHsbw72(idxwm));
d1sf6 = d1sHsbw.*f6d1;
d1rf6 = d1rHsbw.*f6d1;

f7    = ea6/r32*(r32*A./DHsbw3(idxwm) + (-r36 + r81/A*s2(idxwm).*H(idxwm))./Hsbw2);
d1sf7 = ea6/r32/A*(Three*(r27*d1sH(idxwm).*DHsbw4(idxwm).*Hsbw.*s2(idxwm)...
    + Eight*A*d1sHsbw.*(Three*DHsbw4(idxwm) - Four*A*Hsbw3)...
    + r54*DHsbw4(idxwm).*s(idxwm).*(Hsbw - d1sHsbw.*s(idxwm)).*H(idxwm)))...
    ./ (DHsbw4(idxwm).*Hsbw3);
d1rf7 = ea6*d1rHsbw.*(f94./Hsbw3-Three*A./DHsbw4(idxwm) -r81/Sixteen/A*s2(idxwm).*H(idxwm)./Hsbw3);

f8    = ea7/r128/A2*(-Three*srpi*(-r40*A3*Hsbw52+Nine*DHsbw72(idxwm).*(r27*Hsbw2-Six*A*Hsbw+Four*A2))) ./ (DHsbw72(idxwm).*Hsbw52);
f8d1  = ea7*srpi*(r135/r64./Hsbw72 + r729/r256/A2./Hsbw32 -r243/r128/A./Hsbw52 -r105*A/r32./DHsbw92(idxwm));
d1sf8 = d1sHsbw.*f8d1;
d1rf8 = d1rHsbw.*f8d1;

f9    = (r324*ea6*A*eb1(idxwm).*DHsbw4(idxwm).*Hsbw + ea8*(r384*A3*Hsbw3 + DHsbw4(idxwm).*(-r729*Hsbw2 + r324*A*Hsbw - r288*A2)))...
    ./ (r128*A2*DHsbw4(idxwm).*Hsbw3);
f9d1  = -r81*ea6/Sixteen/A*eb1(idxwm)./Hsbw3 + ea8*(r27/Four./Hsbw4+r729/r128/A2./Hsbw2-r81/Sixteen/A./Hsbw3-r12*A./DHsbw5(idxwm));
d1sf9 = d1sHsbw.*f9d1;
d1rf9 = d1rHsbw.*f9d1;

t2t9    = f2.*w + f3.*w2(idxwm) + f4.*w3(idxwm) + f5.*w4(idxwm) + f6.*w5(idxwm)...
    + f7.*w6(idxwm) + f8.*w7(idxwm) + f9.*w8(idxwm);
d1st2t9 = d1sf2.*w + d1sf3.*w2(idxwm) + d1sf4.*w3(idxwm) + d1sf5.*w4(idxwm) + d1sf6.*w5(idxwm)...
    + d1sf7.*w6(idxwm) + d1sf8.*w7(idxwm) + d1sf9.*w8(idxwm);
d1rt2t9 = d1rw.*f2 + (d1rf2 + Two*d1rw.*f3).*w  ...
    + (d1rf3 + Three*d1rw.*f4).*w2(idxwm)       ...
    + (d1rf4 + Four *d1rw.*f5).*w3(idxwm)       ...
    + (d1rf5 + Five *d1rw.*f6).*w4(idxwm)       ...
    + (d1rf6 + Six  *d1rw.*f7).*w5(idxwm)       ...
    + (d1rf7 + Seven*d1rw.*f8).*w6(idxwm)       ...
    + (d1rf8 + Eight*d1rw.*f9).*w7(idxwm) + d1rf9.*w8(idxwm);

term1    = t1 + t2t9 + t10(idxwm);
d1sterm1 = d1st1 + d1st2t9 + d1st10(idxwm);
d1rterm1 = d1rt1 + d1rt2t9 + d1rt10(idxwm);

Fx_wpbe(idxwm)  = X * (term1 + term2(idxwm) + term3(idxwm) + term4(idxwm) + term5(idxwm));
d1sfx(idxwm)    = X * (d1sterm1 + d1sterm2(idxwm) + d1sterm3(idxwm) + d1sterm4(idxwm) + d1sterm5(idxwm));
d1rfx(idxwm)    = X * (d1rterm1 + d1rterm3(idxwm) + d1rterm4(idxwm) + d1rterm5(idxwm));
end