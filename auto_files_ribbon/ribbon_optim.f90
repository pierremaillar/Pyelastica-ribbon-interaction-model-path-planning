!--------------------------------------------------------
!--------------------------------------------------------
!   Semester_Project: Pierre Maillard
!   File name: ribbon.f90
!   Author: Pierre Maillard pierre.maillard@epfl.ch
!   Date created: 22/08/2024
!   Date last modified: 25/09/2024
!   Auto-07p: .f90 file
!--------------------------------------------------------
!--------------------------------------------------------

      module physical_constants
        implicit none
        DOUBLE PRECISION :: nu = 0.4
        DOUBLE PRECISION :: L = 30      ! Length [mm]
        DOUBLE PRECISION :: t = 0.072*10   ! thikness [mm]
        DOUBLE PRECISION :: a = 5.0     ! Width [mm]
        DOUBLE PRECISION :: Y = 2.77D3  ! Youngs Modulus [N/mm^2]
        DOUBLE PRECISION :: k_mat = 10

      end module physical_constants


!     ---------------
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------------

      use physical_constants

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)



      DOUBLE PRECISION nup2, aSqOvLT, aSqOvLTp2, lambda, lambdap2, rhog
      DOUBLE PRECISION x_pos, y_pos, z_pos, d1x, d1y, d1z, d2x, d2y, d2z, d3x, d3y, d3z
      DOUBLE PRECISION m1, m1p, m2, m2p, m3, m3p, R1, R2, R3, k2, k2p, k2p2, k3, k3p, k3p2, k1, s

      DOUBLE PRECISION w_x_pos, w_y_pos, w_z_pos, w_d1x, w_d1y, w_d1z, w_d2x, w_d2y, w_d2z, w_d3x, w_d3y, w_d3z
      DOUBLE PRECISION w_m1, w_m1p, w_m2, w_m2p, w_m3, w_m3p, w_R1, w_R2, w_R3, w_k2, w_k2p, w_k2p2, w_k3, w_k3p, w_k3p2, w_k1, w_s

      DOUBLE PRECISION f1, f2, f3, me1, me2, me3
      DOUBLE PRECISION phi, phiPrime, phiSecond
      DOUBLE PRECISION mat11, mat12, mat21, mat22, deter
      DOUBLE PRECISION RES(2)

     

      DOUBLE PRECISION X_ini, Y_ini, Z_ini
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
DOUBLE PRECISION, DIMENSION(3) :: coeff_X = (/ &
-5.7179494582421264e-01, &
5.6146045920345822e-01, &
3.6325811827610466e-01 /)
DOUBLE PRECISION, DIMENSION(3) :: coeff_Y = (/ &
-2.5560486181732678e-03, &
4.5172813224856889e-01, &
-6.6868812226108010e-01 /)
DOUBLE PRECISION, DIMENSION(3) :: coeff_Z = (/ &
2.8081132936098291e-01, &
-1.5820334672871192e-01, &
-2.9124526533312312e-01 /)
DOUBLE PRECISION :: X_end = 0.6262083489
DOUBLE PRECISION :: Y_end = 0.59707836478
DOUBLE PRECISION :: Z_end = 0.2673042062
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------


!       "One-dimensional model for elastic ribbons constitutive equation from Basile Audoly and Sébastien Neukirch", in the form to be fed to AUTO.

!       PARAMETERS: 
!    	1: Tip Load  
!	(2: Stiffness)
!	(3: Width)
!       (3: Thickness)
!	(4: Strength of the contact surface)
!	(5: Load orientation)
!	(6: Length of arm)    


          


       nup2 = nu*nu
    
       aSqOvLT   = a*a/(L*t)
       !aSqOvLT = 10 ! reproduire résultats de Audoly


       aSqOvLTp2 = aSqOvLT * aSqOvLT
       lambda   = aSqOvLT * sqrt(12.*(1.-nup2)) 
       lambdap2 = lambda * lambda
 

! x aligne with d3 / y with d1 / z with d2


!----------------------------
!     ODE equation system      
!----------------------------
       x_pos = U(1)
       y_pos = U(2)
       z_pos = U(3)
	
       d3x = U(4)
       d3y = U(5)
       d3z = U(6)
       d2x = U(7)
       d2y = U(8)
       d2z = U(9)
       d1x = U(10)
       d1y = U(11)
       d1z = U(12)

       R1 = U(13) 
       R2 = U(14)
       R3 = U(15)       
       m1 = U(16) 
       m2 = U(17) 
       m3 = U(18) 
       k2 = U(19)
       k3 = U(20)
       s = U(21)


!----------------------------
!     ODE augmentation by Optimization function      
!----------------------------
       w_x_pos = U(22)
       w_y_pos = U(23)
       w_z_pos = U(24)

       w_d3x = U(25)
       w_d3y = U(26)
       w_d3z = U(27)
       w_d2x = U(28)
       w_d2y = U(29)
       w_d2z = U(30)
       w_d1x = U(31)
       w_d1y = U(32)
       w_d1z = U(33)

       w_R1 = U(34) 
       w_R2 = U(35)
       w_R3 = U(36)       
       w_m1 = U(37) 
       w_m2 = U(38) 
       w_m3 = U(39) 
       w_k2 = U(40)
       w_k3 = U(41)
       w_s = U(42)


!--------------------------------  


!     external Forces


       f1 = par(1) 

       f2 = par(2)

       f3 = 0


!     external Moments
       me1 = 0
       me2 = 0
       me3 = 0

!--------------------------------  

!     Constitutive equations


       m1p = k3*m2 - k2*m3 + R2 - me1   !m1'
       m2p = k1*m3 - k3*m1 - R1 - me2  !m2'
       m3p = k2*m1 - k1*m2 - me3       !m3'


!     Call the subroutine for computing phi and its derivatives
       call phiAndDerivatives(lambda * k2, phi, phiPrime, phiSecond)


       k2p2 = k2* k2;
       k3p2 = k3* k3;

!     m2 and m3 expressions

       m2 = k2 + 12.0_8 * aSqOvLTp2 * nu * k2 * (nu * k2p2 + k3p2) * phi + &
        3.0_8 * aSqOvLTp2 * lambda * (nu * k2p2 + k3p2) * (nu * k2p2 + k3p2) * phiPrime

       m3 = 2.0_8 * k3 * (1.0_8 / (1.0_8 + nu) + &
        6.0_8 * aSqOvLTp2 * (nu * k2p2 + k3p2) * phi)


       mat11 = 1. + 24.*aSqOvLTp2*nup2* k2p2*phi + 12.*aSqOvLTp2*nu*(nu* k2p2 + k3p2)*phi + &
        24.*aSqOvLTp2*lambda*nu* k2*(nu* k2p2 + k3p2)*phiPrime + &
        3.*aSqOvLTp2*lambdap2*(nu* k2p2 + k3p2)*(nu* k2p2 + k3p2)*phiSecond
       mat12 = 12.*aSqOvLTp2*k3*(2.*nu* k2*phi + lambda*(nu* k2p2 + k3p2)*phiPrime)
       mat21 = 12.*aSqOvLTp2*k3*(2.*nu* k2*phi + lambda*(nu* k2p2 + k3p2)*phiPrime)
       mat22 = 2.*(1/(1. + nu) + 6.*aSqOvLTp2*(nu* k2p2 + 3.*k3p2)*phi)

       deter = mat11 * mat22 - mat12 * mat21

       k2p = (-1.*m3p*mat12 + m2p*mat22)/deter
       k3p = (m3p*mat11 - 1.*m2p*mat21)/deter



!     Kinematics     
       F(1) = d3x    !x'
       F(2) = d3y    !y'
       F(3) = d3z    !z'

       F(4) = k2*d1x - k1*d2x       !d3x'
       F(5) = k2*d1y - k1*d2y       !d3y'
       F(6) = k2*d1z - k1*d2z       !d3z'

       F(7) = k1*d3x - k3*d1x       !d2x'
       F(8) = k1*d3y - k3*d1y       !d2y'
       F(9) = k1*d3z - k3*d1z       !d2z'

       F(10) = k3*d2x - k2*d3x      !d1x'
       F(11) = k3*d2y - k2*d3y      !d1y'
       F(12) = k3*d2z - k2*d3z      !d1z'


!     local balance
       F(13) = k3*R2 - k2*R3 - f1       !R1'
       F(14) = k1*R3 - k3*R1 - f2       !R2'
       F(15) = k2*R1 - k1*R2 - f3       !R3'
       F(16) = m1p                 !m1'
       F(17) = m2p                 !m2'
       F(18) = m3p                !m3'


!     Constitutive equation
       F(19) = k2p                  !k2'
       F(20) = k3p                  !k3'
       F(21) = 1                    !s'




!     Optimization functions

       X_ini = (coeff_X(1) * s**3 + coeff_X(2) * s**2 + coeff_X(3) * s) * (1 - s) + X_end * s 

       Y_ini = (coeff_Y(1) * s**3 + coeff_Y(2) * s**2 + coeff_Y(3) * s) * (1 - s) + Y_end * s

       Z_ini = (coeff_Z(1) * s**3 + coeff_Z(2) *s**2 + coeff_Z(3) * s) * (1 - s) + Z_end * s

  

       w_k1 = 0

       F(22) = -w_d3x + par(32)*2*(x_pos-X_ini)   !w_x'
       F(23) = -w_d3y + par(32)*2*(y_pos-Y_ini)   !w_y'
       F(24) = -w_d3z  + par(32)*2*(z_pos-Z_ini)  !w_z'

       F(25) = -(w_k2*d1x - w_k1*d2x + k2*w_d1x - k1*w_d2x)      !w_d3x'
       F(26) = -(w_k2*d1y - w_k1*d2y + k2*w_d1y - k1*w_d2y)      !w_d3y'
       F(27) = -(w_k2*d1z - w_k1*d2z + k2*w_d1z - k1*w_d2z)      !w_d3z'

       F(28) = -(w_k1*d3x - w_k3*d1x + k1*w_d3x - k3*w_d1x)      !w_d2x'
       F(29) = -(w_k1*d3y - w_k3*d1y + k1*w_d3y - k3*w_d1y)      !w_d2y'
       F(30) = -(w_k1*d3z - w_k3*d1z + k1*w_d3z - k3*w_d1z)      !w_d2z'

       F(31) = -(w_k3*d2x - w_k2*d3x + k3*w_d2x - k2*w_d3x)     !w_d1x'
       F(32) = -(w_k3*d2y - w_k2*d3y + k3*w_d2y - k2*w_d3y)     !w_d1y'
       F(33) = -(w_k3*d2z - w_k2*d3z + k3*w_d2z - k2*w_d3z)     !w_d1z'


!     local balance
       F(34) = -(w_k3*R2 - w_k2*R3 + k3*w_R2 - k2*w_R3)       !w_R1'
       F(35) = -(w_k1*R3 - w_k3*R1 + k1*w_R3 - k3*w_R1)       !w_R2'
       F(36) = -(w_k2*R1 - w_k1*R2 + k2*w_R1 - k1*w_R2)       !w_R3'
       F(37) = -(w_k3*m2 - w_k2*m3 + k3*w_m2 - k2*w_m3 + w_R2)    !w_m1'
       F(38) = -(w_k1*m3 - w_k3*m1 + k1*w_m3 - k3*w_m1 - w_R1)    !w_m2'
       F(39) = -(w_k2*m1 - w_k1*m2 + w_k2*w_m1 - k1*w_m2)         !w_m3'


!     Constitutive equation
       F(40) = -(-1.*(w_k2*m1 - w_k1*m2 + w_k2*w_m1 - k1*w_m2)*mat12 + &
                (w_k1*m3 - w_k3*m1 + k1*w_m3 - k3*w_m1 - w_R1)*mat22)/deter  !w_k2'
       F(41) = -((w_k2*m1 - w_k1*m2 + w_k2*w_m1 - k1*w_m2)*mat11 - &
             1.*(w_k1*m3 - w_k3*m1 + k1*w_m3 - k3*w_m1 - w_R1)*mat21)/deter  !w_k3'

       F(42) = 0                                                             !w_s'
	 

      END SUBROUTINE FUNC

!----------------------------------------------------------------------

      module solver_storage
        implicit none
        DOUBLE PRECISION :: k2_new = 0.0d0, k3_new = 0.0d0
        DOUBLE PRECISION :: m2 = 0.0d0, m3 = 0.0d0, m2_tresh = 0.0d0, m3_tresh = 0.0d0
        LOGICAL :: first_call = .true.
      end module solver_storage

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------------

      use solver_storage

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
      DOUBLE PRECISION :: f1, f2, f3, m1
      DOUBLE PRECISION :: tolerance
      integer :: max_iter, iter
      LOGICAL :: print_on = .false.


! Fait en sorte que FB soit nuls --> U0 pour la boudary à gauche et U1 pour l'autre + (x) pour la variable. Est lancé de très nombreuse fois pour chaque cas et chaque point
!x_pos = U(1)
!y_pos = U(2)
!z_pos = U(3)

!d3x = U(4)
!d3y = U(5)
!d3z = U(6)
!d2x = U(7)
!d2y = U(8)
!d2z = U(9)
!d1x = U(10)
!d1y = U(11)
!d1z = U(12)

!R1 = U(13)
!R2 = U(14)
!R3 = U(15)
!m1 = U(16)
!m2 = U(17)
!m3 = U(18)

!k2 = U(19)
!k3 = U(20)

!s = U(21)




!--------------------------------  
! Boundary conditions: force and moments at the tip

       f1 = par(1)
       f2 = par(2)
       f3 = 0

       m1 = 0
       m2 = 0
       m3 = 0

!Choose a treshold to recompute k2 an k3 if m2 or m3 changes to much
       if (abs(m2-m2_tresh + m3 - m3_tresh)>1.0e-2) then
          first_call = .true.
          m2_tresh = m2
          m3_tresh = m3
       end if        


 ! Set tolerance and maximum iterations
       tolerance = 1.0e-6
       max_iter = 50000

       if (first_call) then
            call newton_solve_k2_k3(k2_new, k3_new, m2, m3/2, m2, m3, tolerance, max_iter, iter)
            if (iter <= max_iter .AND. print_on)  then
                print *,"__________"
                print *, "Converged in ", iter, " iterations."
                print *, "m2 = ", m2, "m3 = ", m3
                print *, "k2 = ", k2_new, "k3 = ", k3_new
            endif
            first_call = .false.   ! Mark that Newton has been executed once
        endif


!--------------------------------  
! positions at the base

       FB(1)=U0(1)
       FB(2)=U0(2)
       FB(3)=U0(3)
!--------------------------------  
! orientation at the base

       FB(4)=U0(4)-1
       FB(5)=U0(5)
       FB(6)=U0(6)
       FB(7)=U0(7)
       FB(8)=U0(8)
       FB(9)=U0(9)-1
       FB(10)=U0(10)
       FB(11)=U0(11)-1
       FB(12)=U0(12)

!--------------------------------  
! stress fild at the tip
  
       FB(13)=U1(13) - f1
       FB(14)=U1(14) - f2
       FB(15)=U1(15) - f3

       FB(16)=U1(16) - m1
       FB(17)=U1(17) - m2
       FB(18)=U1(18) - m3

!--------------------------------
!curvature at the tip or base
       FB(19)=U1(19) - k2_new
       FB(20)=U1(20) - k3_new
!--------------------------------

       FB(21)=U0(21)



!--------------------------------  
! Optimization function
!--------------------------------  

       FB(22)=U0(22) - par(11)
       FB(23)=U1(22)

       FB(24)=U0(23) - par(12)
       FB(25)=U1(23)

       FB(26)=U0(24) - par(13)
       FB(27)=U1(24)
!--------------------------------  

       FB(28)=U0(25) - par(14)
       FB(29)=U1(25)

       FB(30)=U0(26) - par(15)
       FB(31)=U1(26)

       FB(32)=U0(27) - par(16)
       FB(33)=U1(27)

       FB(34)=U0(28) - par(17)
       FB(35)=U1(28)

       FB(36)=U0(29) - par(18)
       FB(37)=U1(29)

       FB(38)=U0(30) - par(19)
       FB(39)=U1(30)

       FB(40)=U0(31) - par(20)
       FB(41)=U1(31)

       FB(42)=U0(32) - par(21)
       FB(43)=U1(32)

       FB(44)=U0(33) - par(22)
       FB(45)=U1(33) 

!--------------------------------  
  
       FB(46)=U1(34) + par(23)
       FB(47)=U0(34)

       FB(48)=U1(35) + par(24)
       FB(49)=U0(35) 

       FB(50)=U1(36) + par(25)
       FB(51)=U0(36) 

       FB(52)=U1(37) + par(26)
       FB(53)=U0(37)

       FB(54)=U1(38) + par(27)
       FB(55)=U0(38)

       FB(56)=U1(39) + par(28)
       FB(57)=U0(39)

!--------------------------------
       FB(58)=U1(40) + par(29)
       FB(59)=U0(40)

       FB(60)=U1(41) + par(30)
       FB(61)=U0(41)
!--------------------------------

       FB(62)=U0(42) - par(31)
       FB(63)=U1(42)


      END SUBROUTINE BCND
!----------------------------------------------------------------------


      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!--------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

      DOUBLE PRECISION P, E

      DOUBLE PRECISION X_ini, Y_ini, Z_ini
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
DOUBLE PRECISION, DIMENSION(3) :: coeff_X = (/ &
-5.7179494582421264e-01, &
5.6146045920345822e-01, &
3.6325811827610466e-01 /)
DOUBLE PRECISION, DIMENSION(3) :: coeff_Y = (/ &
-2.5560486181732678e-03, &
4.5172813224856889e-01, &
-6.6868812226108010e-01 /)
DOUBLE PRECISION, DIMENSION(3) :: coeff_Z = (/ &
2.8081132936098291e-01, &
-1.5820334672871192e-01, &
-2.9124526533312312e-01 /)
DOUBLE PRECISION :: X_end = 0.6262083489
DOUBLE PRECISION :: Y_end = 0.59707836478
DOUBLE PRECISION :: Z_end = 0.2673042062
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

       X_ini = (coeff_X(1) * U(21)**3 + coeff_X(2) * U(21)**2 + coeff_X(3) * U(21)) * (1 - U(21)) + X_end * U(21) 

       Y_ini = (coeff_Y(1) * U(21)**3 + coeff_Y(2) * U(21)**2 + coeff_Y(3) * U(21)) * (1 - U(21)) + Y_end * U(21)

       Z_ini = (coeff_Z(1) * U(21)**3 + coeff_Z(2) * U(21)**2 + coeff_Z(3) * U(21)) * (1 - U(21)) + Z_end * U(21)


! try adding all the relaxation parameters
       FI(1)=U(22)**2 + U(23)**2 + U(24)**2 + U(25)**2 + U(26)**2 + U(27)**2 + U(28)**2 + U(29)**2 + U(30)**2 + U(31)**2 &
       + U(32)**2 + U(33)**2 + U(34)**2 + U(35)**2 + U(36)**2 + U(37)**2 + U(38)**2 + U(39)**2 + U(40)**2 + U(41)**2 + U(42)**2 &
       + par(11)**2 + par(12)**2 + par(13)**2 + par(14)**2 + par(15)**2+ par(16)**2  + par(17)**2 + par(18)**2 + par(19)**2  &
       + par(20)**2 + par(21)**2 + par(22)**2 + par(23)**2 + par(24)**2 + par(25)**2 + par(26)**2 + par(27)**2 + par(28)**2 &
       + par(29)**2 + par(30)**2 + par(31)**2 + par(32)**2 - par(33) 


! Objective funciton
       !FI(2)=PAR(10)-((U(1)-X_ini)**2 + (U(2)-Y_ini)**2 + (U(3)-Z_ini)**2) &
       !        - 0.01*(PAR(1)**2+PAR(2)**2)

       FI(2)=PAR(10)-(U(1)-U(2))**2 &
               - 0.01*(PAR(1)**2+PAR(2)**2)

       !FI(2)=PAR(10)-((U(1)-0.1)**2 + (U(2)-0.1)**2 + (U(3)-0)**2) &
       !        - 0.01*(PAR(1)**2+PAR(2)**2)

! Integral condition from the co-stat 


       FI(3)=U(34)-PAR(23)-PAR(32)*0.2*PAR(1)

     IF(NINT.EQ.3)RETURN
       FI(4)=U(35)-PAR(24)-PAR(32)*0.2*PAR(2) - PAR(34)

      END SUBROUTINE ICND


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------



subroutine phiAndDerivatives(omega2bar, phi, phiprime, phidouble)
    implicit none
    ! Declare input/output variables
    DOUBLE PRECISION, intent(in)  :: omega2bar
    DOUBLE PRECISION, intent(out) :: phi, phiprime, phidouble

    ! Local variables
    DOUBLE PRECISION :: Omega2bp2, Omega2bp4, absOmega2b, sgnOmega2b
    DOUBLE PRECISION :: c1, c2, c3, sqrtAbsOmega2b, q
    DOUBLE PRECISION :: qp2, qp3, qp5, qp6, qp7
    DOUBLE PRECISION :: cshQ, csQ, snhQ, snQ
    DOUBLE PRECISION :: cshQp2, csQp2, cshQp3, csQp3, snQp2, snhQp2
    DOUBLE PRECISION :: snSum, snSump2, snSump3
    DOUBLE PRECISION :: f0, f1, f2

    ! Calculate squared terms
    Omega2bp2 = omega2bar * omega2bar
    Omega2bp4 = Omega2bp2 * Omega2bp2

    ! Calculate absolute value and sign
    absOmega2b = abs(omega2bar)
    sgnOmega2b = sign(1.0_8, omega2bar)

    ! First case: absOmega2b < 0.3
    if (absOmega2b < 0.3_8) then
        c1 = 0.002777777777777778_8
        c2 = -5.5114638447971785D-6
        c3 = 1.1008092191954626D-8
        phi = c1 + c2 * Omega2bp2 + c3 * Omega2bp4
        phiprime = omega2bar * (2.0_8 * c2 + 4.0_8 * c3 * Omega2bp2)
        phidouble = 2.0_8 * c2 + 12.0_8 * c3 * Omega2bp2

    ! Second case: absOmega2b > 1800
    else if (absOmega2b > 1800.0_8) then
        c1 = -5.656854249492381_8
        sqrtAbsOmega2b = sqrt(absOmega2b)
        phi = c1 / (sqrtAbsOmega2b * Omega2bp2) + 2.0_8 / Omega2bp2
        phiprime = -2.5_8 * sgnOmega2b * c1 / Omega2bp4 * sqrtAbsOmega2b - 4.0_8 / (omega2bar * Omega2bp2)
        phidouble = 2.5_8 * 3.5_8 * c1 / (Omega2bp4 * sqrtAbsOmega2b) + 12.0_8 / Omega2bp4

    ! Third case: 1800 > absOmega2b > 0.3
    else
        q = sqrt(absOmega2b / 2.0_8)
        qp2 = q * q
        qp3 = q * qp2
        qp5 = qp2 * qp3
        qp6 = qp3 * qp3
        qp7 = q * qp6

        cshQ = cosh(q)
        csQ = cos(q)
        snhQ = sinh(q)
        snQ = sin(q)

        cshQp2 = cshQ * cshQ
        csQp2 = csQ * csQ
        cshQp3 = cshQ * cshQp2
        csQp3 = csQ * csQp2
        snQp2 = snQ * snQ
        snhQp2 = snhQ * snhQ

        snSum = snhQ + snQ
        snSump2 = snSum * snSum
        snSump3 = snSum * snSump2

        ! Functions of q = sqrt(w2/2) and their derivatives
        f0 = (0.5_8 * (-2.0_8 * cshQ + 2.0_8 * csQ + q * snSum)) / (qp5 * snSum)
        f1 = (cshQp2 * q - csQp2 * q + 5.0_8 * cshQ * snSum - 5.0_8 * csQ * snSum - 3.0_8 * q * snSump2) / (qp6 * snSump2)
	f2 = (2.0_8 * (-cshQp3 * qp2 + csQp3 * qp2 + csQ * (15.0_8 * snhQp2 + (30.0_8 + qp2) * snhQ * snQ + &
      	 (15.0_8 + qp2) * snQp2) + 5.0_8 * csQp2 * q * snSum - cshQp2 * q * (csQ * q + 5.0_8 * snSum) + &
      	 cshQ * (csQp2 * qp2 + (-15.0_8 + qp2) * snhQ * snSum - 15.0_8 * snQ * snSum) + &
      	 10.0_8 * q * snSump3)) / (qp7 * snSump3)

        ! Original function obtained by composition omega2 -> q -> f
        phi = f0
        phiprime = f1 * sgnOmega2b / (4.0_8 * q)
        phidouble = (f2 * q - f1) / (16.0_8 * qp3)
    end if
end subroutine phiAndDerivatives

!----------------------------------------------------------------------
subroutine newton_solve_k2_k3(k2, k3, k2_ini, k3_ini, m2, m3, tolerance, max_iter, iter)
    
    use physical_constants

    implicit none
    DOUBLE PRECISION, intent(inout) :: k2, k3   
    DOUBLE PRECISION, intent(in) :: m2, m3, k2_ini, k3_ini      
    DOUBLE PRECISION, intent(in) :: tolerance
    integer, intent(in) :: max_iter
    integer, intent(out) :: iter

    DOUBLE PRECISION :: k2p2, k3p2, residual_m2, residual_m3
    DOUBLE PRECISION :: dm2_dk2, dm2_dk3, dm3_dk2, dm3_dk3, detJ
    DOUBLE PRECISION :: delta_k2, delta_k3
    DOUBLE PRECISION :: phi, phiPrime, phiSecond
    DOUBLE PRECISION :: lambda, aSqOvLT, aSqOvLTp2
    
      k2 = k2_ini
      k3 = k3_ini

      aSqOvLT = a*a/(L*t)
      aSqOvLT = aSqOvLT * aSqOvLT
      lambda = aSqOvLT * sqrt(12.*(1.-nu*nu)) 
      

     ! Newton-Raphson method loop
      do iter = 1, max_iter

        ! Compute k2^2 and k3^2
        k2p2 = k2 * k2
        k3p2 = k3 * k3

        ! Call the subroutine to compute phi, phiPrime, and phiSecond
        call phiAndDerivatives(lambda * k2, phi, phiPrime, phiSecond)

        ! Compute residuals
        residual_m2 = k2 + 12.0_8 * aSqOvLTp2 * nu * k2 * (nu * k2p2 + k3p2) * phi + &
                      3.0_8 * aSqOvLTp2 * lambda * (nu * k2p2 + k3p2) * (nu * k2p2 + k3p2) * phiPrime - m2

        residual_m3 = 2.0_8 * k3 * (1.0_8 / (1.0_8 + nu) + &
                      6.0_8 * aSqOvLTp2 * (nu * k2p2 + k3p2) * phi) - m3

        ! Check if residuals are within tolerance
        if (abs(residual_m2) < tolerance .and. abs(residual_m3) < tolerance) then
            return
        endif

        ! Compute the Jacobian matrix elements
        dm2_dk2 = 1.0_8 + 12.0_8 * aSqOvLTp2 * nu * (3.0_8 * nu * k2p2 * phi + &
                   k3p2) + 12.0_8 * aSqOvLTp2 * lambda * nu * k2 * (nu * k2p2 + k3p2) * phiPrime + &
                   12.0_8 * aSqOvLTp2 * lambda * (nu * k2p2 + k3p2) * nu * k2 * phiPrime + &
                   3.0_8 * aSqOvLTp2 * lambda * lambda * (nu * k2p2 + k3p2) * (nu * k2p2 + k3p2) * phiSecond
 
        dm2_dk3 = 12.0_8 * aSqOvLTp2 * nu * k2 * 2.0_8 * k3 * phi + &
                  6.0_8 * aSqOvLTp2 * lambda * 2.0_8 * k3 * (nu * k2p2 + k3p2) * phiPrime

        dm3_dk2 = 2.0_8 * k3 * 6.0_8 * aSqOvLTp2 * (nu * 2.0_8 * k2 * phi + (nu * k2p2 + k3p2) * phiPrime) 
                   
        dm3_dk3 = 2.0_8 * (1.0_8 / (1.0_8 + nu) + 6.0_8 * aSqOvLTp2 * (nu * k2p2 + k3p2) * phi) + &
                  12.0_8 * aSqOvLTp2 * 2.0_8 * k3 * phi

        ! Compute the determinant of the Jacobian
        detJ = dm2_dk2 * dm3_dk3 - dm2_dk3 * dm3_dk2

        ! Solve for delta_k2 and delta_k3 using the inverse of the Jacobian
        delta_k2 = -(residual_m2 * dm3_dk3 - residual_m3 * dm2_dk3) / detJ
        delta_k3 = -(residual_m3 * dm2_dk2 - residual_m2 * dm3_dk2) / detJ

        ! Update k2 and k3
        k2 = k2 + 1.5*delta_k2
        k3 = k3 + 1.5*delta_k3

      end do

      print*, "Did not converge, M2 and M3 residuals are", residual_m2, residual_m3
      print*, "M2 = ", m2, "M3 =", m3

      end subroutine newton_solve_k2_k3

!     -----------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      REAL(8), PARAMETER :: PI = 4 * atan(1.0)
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      INTEGER :: I


! Est lancé nombreuses fois avant le début de la résolution       
       DO I=1,34
       PAR(I)=0.0
       ENDDO

       PAR(1)=0.001	
       PAR(10)=1.0


       U(1)= 0.0
       U(2)= 0.0
       U(3)= 0.0

       U(4)= 1.0
       U(5)=0.0
       U(6)=0.0
       U(7)=0.0
       U(8)=0.0
       U(9)= 1.0
       U(10)=0.0
       U(11)= 1.0
       U(12)=0.0

       U(13)=0.0
       U(14)=0.0
       U(15)=0.0
       U(16)=0.0
       U(17)=0.0
       U(18)=0.0

       U(19)=0.0
       U(20)=0.0
       U(21)= 0.0


       U(22)= 0.0
       U(23)= 0.0
       U(24)= 0.0

       U(25)= 0.0
       U(26)=0.0
       U(27)=0.0
       U(28)=0.0
       U(29)=0.0
       U(30)= 0.0
       U(31)=0.0
       U(32)= 0.0
       U(33)=0.0

       U(34)=0.0
       U(35)=0.0
       U(36)=0.0
       U(37)=0.0
       U(38)=0.0
       U(39)=0.0

       U(40)=0.0
       U(41)=0.0
       U(42)= 0.0



      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


      SUBROUTINE PVLS(NDIM,U,PAR)
      END SUBROUTINE PVLS
