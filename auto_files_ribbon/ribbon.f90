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
        DOUBLE PRECISION :: nu = 0.34
        DOUBLE PRECISION :: L = 50      ! Length [mm]
        DOUBLE PRECISION :: t = 0.5    !t = 0.1     ! thickness [mm]
        DOUBLE PRECISION :: a = 5.0     ! Width [mm]
        DOUBLE PRECISION :: Y = 2.77e3  ! Youngs Modulus of the ribbon [N/mm^2]
	
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
      DOUBLE PRECISION fe1, fe2, fe3, me1, me2, me3
      DOUBLE PRECISION phi, phiPrime, phiSecond
      DOUBLE PRECISION mat11, mat12, mat21, mat22, deter
      DOUBLE PRECISION RES(2)     

!---!
DOUBLE PRECISION :: w1 = 65.9781449733566
DOUBLE PRECISION :: w2 = 61.53227236821222
DOUBLE PRECISION :: w3 = 9.4691987778594
DOUBLE PRECISION :: w4 = -59.64790742879249
DOUBLE PRECISION :: w5 = 92.47116690320811
DOUBLE PRECISION :: w6 = 5.508010461615811
DOUBLE PRECISION :: w7 = 37.92694606317317
DOUBLE PRECISION :: w8 = 56.90369584173426
!---!

 
!       "One-dimensional model for elastic ribbons constitutive equation from Basile Audoly and Sébastien Neukirch", in the form to be fed to AUTO.

!       PARAMETERS: 
!    	1: Convergence parameter 
!       2: epsilon 1
!       3: epsilon 2
!       2: epsilon 3
!       3: epsilon 4
!       2: epsilon 5
!       3: epsilon 6
!       2: epsilon 7
!       3: epsilon 8
 

       nup2 = nu*nu
    
       aSqOvLT   = a*a/(L*t)


       aSqOvLTp2 = aSqOvLT * aSqOvLT
       lambda   = aSqOvLT * sqrt(12.*(1.-nup2)) 
       lambdap2 = lambda * lambda
 

! axis convension at the base with no deformations:
! x aligne with d3 / y with d1 / z with d2
! d3 aligne to the length of the ribbon/ d2 with the width/ d1 with the thickness


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

!--------------------------------  


!!!---external Forces---!!!

!       w_s = (/ w1 , w2 , w3 , w4 , w5, w6 , w7 , w8 , w9 , w10 /)
!       call bezier(w_s, s, fe1) ! bezier's curve


       !fe1 = par(1)*(w1*s + w2*s**2 + w3*s**3 + w4*s**4) !4 order polynom
       !fe1 = par(1)*(10*s -30*s**2 + 15*s**3)            !3 order polynom
       !fe1 = par(1)*((-w1 + tanh(w2 * (s - w3)))*w4)     !Hyperbolic tangent funciton (step function)

       !! sum of gaussian shapes
       !fe1 = par(1)*(w1* exp(-((s - w2/100)**2) / (2.0 * (w3/100)**2)) + w4 * exp(-((s - w5/100)**2) / (2.0 * (w6/100)**2)))
       !!

       fe1 = 0

       fe2 = 0

       fe3 = 0

!!!---external Moments---!!!
       me1 = 0
       me2 = 0
       me3 = 0
       !me3 = par(1)**2*((s-w8/100)*w7) !first order polynom

!--------------------------------  

!!!     Constitutive equations

!     Call the subroutine for computing phi and its derivatives
       call phiAndDerivatives(lambda * k2, phi, phiPrime, phiSecond)

       m1p = k3*m2 - k2*m3 + R2 - me1   !m1'
       m2p = k1*m3 - k3*m1 - R1 - me2  !m2'
       m3p = k2*m1 - k1*m2 - me3       !m3'


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
       F(13) = k3*R2 - k2*R3 - fe1   !R1'
       F(14) = k1*R3 - k3*R1 - fe2   !R2'
       F(15) = k2*R1 - k1*R2 - fe3   !R3'
       F(16) = m1p                  !m1'
       F(17) = m2p                  !m2'
       F(18) = m3p                  !m3'


!     Constitutive equation
       F(19) = k2p                  !k2'
       F(20) = k3p                  !k3'
       F(21) = 1                    !s'

	 

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
      DOUBLE PRECISION :: fe1, fe2, fe3, m1
      DOUBLE PRECISION :: tolerance
      integer :: max_iter, iter
      LOGICAL :: print_on = .false.


!---!
DOUBLE PRECISION :: w1 = 1
DOUBLE PRECISION :: w2 = 30
DOUBLE PRECISION :: w3 = 0.4
DOUBLE PRECISION :: w4 = 100
DOUBLE PRECISION :: w7 = 0.16
DOUBLE PRECISION :: w8 = 2
!---!


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


!!!---external Forces at the tip---!!!

       fe1 = +15.0*U1(12)*par(1)+0.15*U1(11)*par(1)
       fe2 = +15.0*U1(9)*par(1)+0.15*U1(8)*par(1)
       fe3 = +15.0*U1(6)*par(1)+0.15*U1(5)*par(1)



!!!---external Moments at the tip---!!!
       m1 = 0
       m2 = 0
       m3 = 0

!Choose a threshold of m2 and m3 variation from which to recalculate k2 and k3 at tip.
       if (abs(m2-m2_tresh + m3 - m3_tresh)>1.0e-6) then
          first_call = .true.
          m2_tresh = m2
          m3_tresh = m3
       end if        


         


 ! Set tolerance and maximum iterations
       tolerance = 1.0e-6
       max_iter = 50000

 ! solve the consitutive equation at the tip to solve issue when imposing moment (m2 != 0, m3 != 0) boundary conditions
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
  
       FB(13)=U1(13) - fe1
       FB(14)=U1(14) - fe2
       FB(15)=U1(15) - fe3


       FB(16)=U1(16) - m1
       FB(17)=U1(17) - m2
       FB(18)=U1(18) - m3

!--------------------------------
!curvature at the tip or base
       FB(19)=U1(19) - k2_new
       FB(20)=U1(20) - k3_new
!--------------------------------

       FB(21)=U0(21)



      END SUBROUTINE BCND

!----------------------------------------------------------------------

subroutine bezier(control_y, s, y_s)
    implicit none
    DOUBLE PRECISION, intent(in), dimension(5) :: control_y  ! Y-coordinates of control points
    DOUBLE PRECISION, intent(in) :: s                        ! Parameter value (0 <= s <= 1)
    DOUBLE PRECISION, intent(out) :: y_s                     ! Y-coordinate on the Bézier curve at s
    integer :: i, n
    DOUBLE PRECISION :: bernstein

    n = 4  ! Degree of the Bézier curve (5 points -> degree 4)

    

    ! Sum up contributions from all Bernstein polynomials
    y_s = 0.0_8
    do i = 0, n
        bernstein = comb(n, i) * s**i * (1.0_8 - s)**(n - i)  ! Bernstein polynomial B_i^n(s)
        y_s = y_s + bernstein * control_y(i + 1)
    end do

contains

    ! Function to compute binomial coefficient "n choose i"
    pure function comb(n, k) result(c)
        integer, intent(in) :: n, k
        integer :: j
        DOUBLE PRECISION :: c

        if (k == 0 .or. k == n) then
            c = 1.0_8
        else
            c = 1.0_8
            do j = 1, k
                c = c * real(n - j + 1, 8) / real(j, 8)
            end do
        end if
    end function comb

end subroutine bezier

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


subroutine exforce(error_x, error_y, error_z, k, fx, fy, fz)
    implicit none
    ! Declare input/output variables
    DOUBLE PRECISION, intent(in)  :: error_x, error_y, error_z, k
    DOUBLE PRECISION, intent(out) :: fx,fy,fz


    if (error_x > 0.0) then
    fx = k*error_x - 0.0*k
    else if (error_x < -0.0) then
    fx = k*error_x + 0.0*k
    else
    fx = 0
    end if

    if (error_y > 0.0) then
    fy = k*error_y - 0.0
    else if (error_y < -0.0) then
    fy = k*error_y + 0.0*k
    else
    fy = 0
    end if

    if (error_z > 0.0) then
    fz = k*error_z - 0.0*k
    else if (error_z < -0.0) then
    fz = k*error_z + 0.0*k
    else
    fz = 0
    end if

end subroutine exforce


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
    

       PAR(1)=0.0


       U(1)= T
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



      END SUBROUTINE STPNT
!----------------------------------------------------------------------
      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      !----------------------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      !     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
      !----------------------------------------------------------------------
      ! NOTE :
      ! Parameters set in this subroutine should be considered as ``solution
      ! measures'' and be used for output purposes only.
      !
      ! They should never be used as `true'' continuation parameters.
      !
      ! They may, however, be added as ``over-specified parameters'' in the
      ! parameter list associated with the AUTO-Constant NICP, in order to
      ! print their values on the screen and in the ``p.xxx file.
      !
      ! They may also appear in the list associated with AUTO-constant NUZR.
      !
      !----------------------------------------------------------------------
      ! For algebraic problems the argument U is, as usual, the state vector.
      ! For differential equations the argument U represents the approximate
      ! solution on the entire interval [0,1]. In this case its values can
      ! be accessed indirectly by calls to GETP, as illustrated below, or
      ! by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
      ! U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS.
      !----------------------------------------------------------------------

      ! Set PAR(4) equal to the value of U(2) at the left boundary.
       ! Tip position
       PAR(2)=GETP('BV1',1,U)
       PAR(3)=GETP('BV1',2,U)
       PAR(4)=GETP('BV1',3,U)

       !director vector at tip
       !d3
       PAR(5)=GETP('BV1',4,U)
       PAR(6)=GETP('BV1',5,U)
       PAR(7)=GETP('BV1',6,U)
       !d2
       PAR(8)=GETP('BV1',7,U)
       PAR(9)=GETP('BV1',8,U)
       PAR(10)=GETP('BV1',9,U)
       !d1
       PAR(11)=GETP('BV1',10,U)
       PAR(12)=GETP('BV1',11,U)
       PAR(13)=GETP('BV1',12,U)

       !stresses at base
       !stress
       PAR(14)=GETP('BV0',13,U)
       PAR(15)=GETP('BV0',14,U)
       PAR(16)=GETP('BV0',15,U)
       !Moments
       PAR(17)=GETP('BV0',16,U)
       PAR(18)=GETP('BV0',17,U)
       PAR(19)=GETP('BV0',18,U)

       !Curvatures at the tip
       PAR(20)=GETP('BV1',19,U)
       PAR(21)=GETP('BV1',20,U)




      !----------------------------------------------------------------------
      ! The first argument of GETP may be one of the following:
      !        'NRM' (L2-norm),     'MAX' (maximum),
      !        'INT' (integral),    'BV0 (left boundary value),
      !        'MIN' (minimum),     'BV1' (right boundary value).
      !        'MNT' (t value for minimum)
      !        'MXT' (t value for maximum)
      !        'NDIM', 'NDX' (effective (active) number of dimensions)
      !        'NTST' (NTST from constant file)
      !        'NCOL' (NCOL from constant file)
      !        'NBC'  (active NBC)
      !        'NINT' (active NINT)
      !        'DTM'  (delta t for all t values, I=1...NTST)
      !        'WINT' (integration weights used for interpolation, I=0...NCOL)
      !
      ! Also available are
      !   'STP' (Pseudo-arclength step size used).
      !   'FLD' (`Fold function', which vanishes at folds).
      !   'BIF' (`Bifurcation function', which vanishes at singular points).
      !   'HBF' (`Hopf function'; which vanishes at Hopf points).
      !   'SPB' ( Function which vanishes at secondary periodic bifurcations).
      !   'EIG' ( Eigenvalues/multipliers, I=1...2*NDIM, alternates real/imag parts).
      !   'STA' ( Number of stable eigenvalues/multipliers).
      !----------------------------------------------------------------------

      END SUBROUTINE PVLS

!----------------------------------------------------------------------
