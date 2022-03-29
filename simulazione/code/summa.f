*==============================================================================================
*     LEPTONIC AND TOP CONTRIBUTION TO VACUUM POLARIZATION
** with correct expansion for |q2| <<< m2!!
      FUNCTION SUMMA(AM,Q2,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NC(9),QF2(9)
      common/ncqf2/nc,qf2
      data nc  /1.d0,1.d0,1.d0,3.d0,
     . 3.d0, !u
     . 3.d0, !d
     . 3.d0, !c
     . 3.d0, !s
     . 3.d0 /!b
      data qf2 /1.d0,1.d0,1.d0,0.44444444444444444444d0,
     . 0.44444444444444444444d0, !u
     . 0.11111111111111111111d0, !d
     . 0.44444444444444444444d0, !c
     . 0.11111111111111111111d0, !s
     . 0.11111111111111111111d0 /!b
* NC AND QF ARE COLOR FACTOR (1 FOR LEPTONS, 
* 3 FOR TOP) AND CHARGE
      AM2=AM*AM

      if (Q2.eq.0.d0) then
         summa = 0.d0
         return
      endif

      z6lim = 1d-8

      IF (Q2.GT.0.D0.AND.Q2.LT.(4.D0*AM2)) THEN
         x = 4.D0*AM2/Q2
         z = sqrt(1.d0/x)
         z2 = z*z
         z4 = z2*z2
         z6 = z2*z4

         if (z6.gt.z6lim) then
            SQ = DSQRT(1d0 - z2)
            s1 = -5.D0/3.D0 - 1.d0/z2
            s2 = (1.d0 + z2- 2.D0*z4)/z/z2*1.d0/SQ*DATAN(z/SQ)
            s1ps2 = s1 + s2
         else
*     with wxmaxima
* s2 = taylor((1+z^2-2*z^4)/sqrt(1-z^2)/z^3*atan(z/sqrt(1-z^2)), z , 0 , 20);
            s1ps2 = -z2*
     .       (4.d0/5.d0+12.d0/35.d0*z2+65.d0/315.d0*z4+32.d0/231.d0*z6)
         endif
         SUMMA=1.d0/3.d0*NC(I)*QF2(I)*s1ps2
      ELSE
           x = -4.D0*AM2/Q2
           if (x.gt.0.d0) then ! space-like
              z  = sqrt(1.d0/x)
              z2 = z*z
              z4 = z2*z2
              z6 = z2*z4
              
              if (z6.gt.z6lim) then
                SQ     = DSQRT(1.D0 + z2)
                ARGLOG = DABS( (z - sq)/(z+sq) )              
                s1 = -5.D0/3.D0 + 1.d0/z2
                s2 = 0.5d0*(1.d0 - z2 - 2.D0*z4)/SQ*DLOG(ARGLOG)/z2/z
                
                s1ps2 = s1 + s2
              else
* with wxmaxima
* s2 = taylor(log( abs ((z - sqrt(1+z^2))/(z + sqrt(1+z^2)) ))/sqrt(1+z^2)*(1-z^2-2*z^4)/2/z^3,z,0,20); 
                 s1ps2 = z2*
     .        (4.d0/5.d0-12.d0/35.d0*z2+65.d0/315.d0*z4-32.d0/231.d0*z6)
              endif
           else ! for x > 0 ie q2 > 4m^2
              x = -x
              z  = sqrt(1.d0/x)
              z2 = z*z
              z4 = z2*z2

              SQ     = DSQRT(z2-1.d0)
              ARGLOG = DABS((z - sq)/(z+sq) )            
              s1 = -5.D0/3.D0 - 1.d0/z2
              s2 =  0.5d0*(1.d0 + z2 - 2.D0*z4)/SQ*DLOG(ARGLOG)/z2/z
              s1ps2 = s1 + s2
           endif
           SUMMA = 1.d0/3.d0 * NC(I)*QF2(I)*s1ps2
       ENDIF
       RETURN
       END
