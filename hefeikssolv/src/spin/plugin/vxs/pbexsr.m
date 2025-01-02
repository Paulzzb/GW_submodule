function [sx, V1X, V2X] = pbexsr(RHO,GRHO,OMEGA)

     SMALL=1.D-20;
     SMAL2=1.D-08;
     US=0.161620459673995492;
     AX=-0.738558766382022406;
     UM=0.2195149727645171;
     UK=0.8040;
     UL=UM/UK;
     f1 = -1.10783814957303361;
     alpha = 2.0/3.0;

      RS = RHO^(1.0/3.0);
      VX = (4.0/3.0)*f1*alpha*RS;

      AA    = GRHO;
      RR    = 1.0/(RHO*RS);
      EX    = AX/RR;
      S2    = AA*RR*RR*US*US;

      S = sqrt(S2);
      if S > 8.3
          S = 8.572844D0 - 18.796223D0/S2;
      end
      [FX,D1X,D2X] = wpbe_analy_erfc_approx_grad(RHO,S,OMEGA);

      sx = EX*FX;       
      DSDN = -4.D0/3.D0*S/RHO;
      V1X = VX*FX + (DSDN*D2X+D1X)*EX;
      DSDG = US*RR;
      V2X = EX*1.D0/sqrt(AA)*DSDG*D2X;

end