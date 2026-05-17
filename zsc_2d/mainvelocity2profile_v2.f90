!-------------------------------------------------------------------------------
! Subroutine  : Common
!! Modules that contain common variables
!-------------------------------------------------------------------------------

!> @brief Definition of single and double data types
 MODULE NTypes
 IMPLICIT NONE
 SAVE

 INTEGER, PARAMETER :: SGL = KIND(1.0)
 INTEGER, PARAMETER :: DBL = KIND(1.D0)

 END MODULE NTypes

!> @brief Parameters related to the geometry and the time intervals
 MODULE Domain
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! Maximum time of steps and data dump times
 INTEGER :: MaxStep, RelaxStep, iCount, iStep, tCall, tDump, tStat 

! Domain size
 INTEGER :: xmin, xmax, ymin, ymax, xmin_f, xmax_f, ymin_f, ymax_f
 INTEGER :: now, nxt, NX, NY, NX_f, NY_f

! ndim = spatial dimension
 INTEGER, PARAMETER :: ndim = 2

! Neighbour arrays
 INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ni, ni_f

! Initial volume of the discrete phase
 REAL(KIND = DBL) :: invInitVol

! Define constants used in the code (inv6 = 1/6, inv12 = 1/12, invPi = 1/pi)
 REAL(KIND = DBL), PARAMETER :: inv6  = 0.16666666666666666667D0
 REAL(KIND = DBL), PARAMETER :: inv12 = 0.08333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: invPi = 0.31830988618379067154D0

 END MODULE Domain

!> @brief Parameters related to hydrodynamics quantities
 MODULE FluidParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

 INTEGER :: nBubbles
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:) :: bubbles
 REAL(KIND = DBL), DIMENSION(1:11) :: Convergence
 REAL(KIND = DBL) :: kappa, kappaG, kappa_6, alpha, alpha4
 REAL(KIND = DBL) :: D, sigma, IntWidth, delt
 REAL(KIND = DBL) :: rhoL, rhoH, tauRho, invTauRho, invTauRho2
 REAL(KIND = DBL) :: tauPhi, invTauPhi
 REAL(KIND = DBL) :: phistar, phistar2, phistar4
 REAL(KIND = DBL) :: Gamma, eps, pConv
 REAL(KIND = DBL) :: eta, eta2, invEta2

! Hydrodynamics arrays
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:)   :: p, p_f, phi_f
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:) :: u, u_f

! Differential term arrays
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:) :: lapPhi_f
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:) :: gradPhiX_f, gradPhiY_f

 END MODULE FluidParams

!> @brief Parameters related to the LBM discretization
 MODULE LBMParams
 USE NTypes, ONLY : DBL
 USE Domain, ONLY : ndim
 IMPLICIT NONE
 SAVE

! fdim = order parameter distribution dimension - 1 (D2Q5)
! gdim = momentum distribution dimension - 1 (D2Q9)
 INTEGER, PARAMETER :: fdim = 4
 INTEGER, PARAMETER :: gdim = 8

! Distribution functions
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:)   :: f, fcol
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:,:) :: g

! D2Q9 Lattice speed of sound (Cs = 1/DSQRT(3), Cs_sq = 1/3, invCs_sq = 3)
 REAL(KIND = DBL), PARAMETER :: Cs       = 0.57735026918962576451D0
 REAL(KIND = DBL), PARAMETER :: Cs_sq    = 0.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: invCs_sq = 3.00000000000000000000D0

 !REAL(KIND = DBL), PARAMETER :: Cs       = 17.0D0
 !REAL(KIND = DBL), PARAMETER :: Cs_sq    = 289.0D0
 !REAL(KIND = DBL), PARAMETER :: invCs_sq = 3.40D-3


! Distributions weights (D2Q9: Eg0 = 4/9, Eg1 = 1/9, Eg2 = 1/36)
 REAL(KIND = DBL), PARAMETER :: Eg0 = 0.44444444444444444444D0
 REAL(KIND = DBL), PARAMETER :: Eg1 = 0.11111111111111111111D0
 REAL(KIND = DBL), PARAMETER :: Eg2 = 0.02777777777777777778D0

! Modified distribution weight (to avoid operations: EgiC = Egi*invCs_sq)
 REAL(KIND = DBL), PARAMETER :: Eg0C = 1.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: Eg1C = 0.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: Eg2C = 0.08333333333333333333D0
 REAL(KIND = DBL) :: Eg0T, Eg1T, Eg2T

 END MODULE LBMParams

!-------------------------------------------------------------------------------
! Program  : ZSC-2D-DGR
!-------------------------------------------------------------------------------
!! Driver for the dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM
!> @details
!! Driver for the dual grid implementation of the Zheng-Shu-Chew multiphase LBM
!! using D2Q5/D2Q9 discretization and periodic boundary conditions. For details:
!!
!! Journal of Computational Physics 218: 353-371, 2006.
!!
!! Serial Implementation. Convergence test and post-processing functions
!! designed for single drop/bubble cases.
!!
!! In this dual grid implementation the order parameter grid is twice larger
!! than the coarser momentum grid. This provides improved numerical
!! differentials at the interface and leads to better mass conservation and
!! greater stability. All differential terms are calculated in the order
!! parameter grid and then copied over to the momentum grid (the grids overlap).
!! The velocity and pressure are calculated in the momentum grid and
!! interpolated in the order parameter grid. Two collision-stream steps are
!! are taken in the order parameter grid for every step taken in the momentum
!! grid to ensure synchronization.
!!
!! The average velocity, mass conservation factor, effective radius of the drop,
!! pressure difference between the inside and the outside of the drop and the
!! error with respect to the analytical value given by Laplace's equation are
!! written to file "stats.out"

!-------------------------------------------------------------------------------

	PROGRAM main

! Common Variables
	USE Domain,      ONLY : iStep, MaxStep, now, nxt, RelaxStep, tDump, tStat
	USE FluidParams, ONLY : eps, pConv
	IMPLICIT NONE


! Read input parameters, allocate memory for common arrays and initialize
	CALL Parameters
	CALL MemAlloc(1)
	CALL Init

! Main iteration loop
	DO iStep = 1, MaxStep
	print*,iStep,MaxStep
! Full Step for momentum
	CALL CollisionG
	CALL Update

! Fractional Step 1 for order parameter
	CALL CollisionF
	CALL Stream
	CALL Differentials

! Fractional Step 2 for order parameter
	CALL CollisionF
	CALL Stream
	CALL Differentials

! Update momentum mesh index
	   now = 1 - now
	   nxt = 1 - nxt

! Save data if required
	   IF ( MOD(istep,tStat) == 0 ) CALL Stats
	   IF ( MOD(istep,tDump) == 0 ) CALL VtkPlane

! Test pressure convergence
	   IF ( eps < pConv ) EXIT

	 END DO
	 RelaxStep = iStep - 1

! Save final data
	 CALL VtkPlane
	 CALL FinalDump

! Free memory
	 CALL MemAlloc(2)


END PROGRAM
!*************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : Parameters
!-------------------------------------------------------------------------------
!! Read input parameters and define simulation constants
!> @details
!! Read input parameters from files "properties.in" and "discrete.in" and define
!! constants for the simulation in the dual grid D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM.

!-------------------------------------------------------------------------------

	 SUBROUTINE Parameters

! Common Variables
	 USE Domain
	 USE FluidParams
	 USE LBMParams
	 IMPLICIT NONE

! Local Variables
	 INTEGER :: i
	 INTEGER :: IO_ERR


! Read parameter data
	 OPEN(UNIT = 10, FILE = "properties.in", STATUS = "OLD", ACTION ="READ",IOSTAT = IO_ERR)
	 IF ( IO_ERR == 0 ) THEN
	     READ(10,*)
	     READ(10,*) Maxstep
	     READ(10,*)
	     READ(10,*) tStat, tDump, delt
	     READ(10,*)
	     READ(10,*) xmin, xmax, ymin, ymax
	     READ(10,*)
	     READ(10,*) rhoL, rhoH
	     READ(10,*)
	     READ(10,*) tauRho, tauPhi
	     READ(10,*)
	     READ(10,*) IntWidth, sigma, Gamma
	     READ(10,*)
	     READ(10,*) pConv
	     CLOSE(UNIT = 10)
	 ELSE
	   STOP "Error: Unable to open input file 'properties.in'."
	 END IF

!  Read bubble positions
	 OPEN(UNIT = 10, FILE = "discrete.in", STATUS = "OLD", ACTION = "READ", IOSTAT = IO_ERR)
	 IF ( IO_ERR == 0 ) THEN
	     READ(10,*)
	     READ(10,*) nBubbles
	     READ(10,*)
	     ALLOCATE( bubbles(1:nBubbles,1:3) )
	     DO i = 1, nBubbles
       READ(10,*) bubbles(i,1), bubbles(i,2), bubbles(i,3)
	     END DO
	     CLOSE(UNIT = 10)
	 ELSE
	   STOP "Error: Unable to open input file 'discrete.in'."
	 END IF

! Scale parameters for order parameter mesh
	 sigma    = 2.D0*sigma
	 IntWidth = 2.D0*IntWidth
	 DO i = 1, nBubbles
	  bubbles(i,1) = 2*bubbles(i,1) - 1
	   bubbles(i,2) = 2*bubbles(i,2) - 1
	   bubbles(i,3) = 2*bubbles(i,3)
	
	 END DO

! Fluid properties
	 phistar    = 0.5D0*(rhoH - rhoL)
	 phistar2   = phistar*phistar
	 phistar4   = phistar2*phistar2
	 eta        = 1.D0/(tauPhi + 0.5D0)
	 eta2       = 1.D0 - eta
	 invEta2    = 0.5D0/eta
	 invTauPhi  = 1.D0/tauPhi
	 invTauRho  = 1.D0/tauRho
	 invTauRho2 = 1.D0 - 0.5D0*invTauRho

! Chemical potential parameters
	 alpha   = 0.75D0*sigma/(IntWidth*phistar4)
	 alpha4  = 4.D0*alpha
	 kappa   = 0.50D0*alpha*(IntWidth*phistar)**2.D0
	 kappaG  = 0.25D0*kappa
	 kappa_6 = kappa*inv6

! Set domain limits for order parameter mesh
	 xmin_f = xmin
	 xmax_f = 2*xmax - 1
	 ymin_f = ymin
	 ymax_f = 2*ymax - 1

! Modified distribution weights for use in CollisionG
	 Eg0T = Eg0C*invTauRho2
	 Eg1T = Eg1C*invTauRho2
	 Eg2T = Eg2C*invTauRho2

! Redefine MaxStep (iterate between 0 and MaxStep-1)
	 MaxStep = MaxStep - 1

	 RETURN
	 END SUBROUTINE Parameters
!*********************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : memAlloc
!-------------------------------------------------------------------------------
!> @file
!! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
!> @details
!! If input FLAG is 1 then memory is allocated for the following common arrays:
!!
!! - \b p   (pressure)
!! - \b u   (velocity)
!! - \b phi (order parameter)
!! - \b ni  (near neighbors)
!! - \b f   (order parameter distribution function)
!! - \b g   (momentum distribution function)
!!
!! These are the arrays required by the dual grid D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM, and are deallocated when FLAG is different from one.

!-------------------------------------------------------------------------------

 SUBROUTINE MemAlloc(FLAG)

! Common Variables
 USE Domain,      ONLY : ni, ni_f, xmax, xmax_f, xmin, xmin_f, ymax, ymax_f, ymin, ymin_f
 USE FluidParams, ONLY : gradPhiX_f, gradPhiY_f, lapPhi_f, p, p_f, phi_f, u, u_f
 USE LBMParams,   ONLY : f, fcol, fdim, g, gdim
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Hydrodynamics and neighbor arrays - Momentum/Pressure mesh
 ALLOCATE(   p(xmin:xmax,ymin:ymax) )
 ALLOCATE(   u(xmin:xmax,ymin:ymax,2) )
 ALLOCATE(  ni(xmin:xmax,ymin:ymax,1:4) )

! Hydrodynamics and neighbor arrays - Order parameter mesh
 ALLOCATE(   p_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(   u_f(xmin_f:xmax_f,ymin_f:ymax_f,2) )
 ALLOCATE( phi_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(  ni_f(xmin_f:xmax_f,ymin_f:ymax_f,1:4) )

! Differential term arrays - Order parameter mesh
 ALLOCATE( gradPhiX_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( gradPhiY_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(   lapPhi_f(xmin_f:xmax_f,ymin_f:ymax_f) )

! Memory allocation for distribution functions
 ALLOCATE(    f(xmin_f:xmax_f,ymin_f:ymax_f,0:fdim) )
 ALLOCATE( fcol(xmin_f:xmax_f,ymin_f:ymax_f,0:fdim) )
 ALLOCATE(    g(xmin:xmax,ymin:ymax,0:gdim,0:1) )

 ELSE

! Free momentum/pressure mesh arrays
 DEALLOCATE( p )
 DEALLOCATE( u )
 DEALLOCATE( ni )
 DEALLOCATE( g )

! Free order parameter mesh arrays
 DEALLOCATE( f )
 DEALLOCATE( fcol )
 DEALLOCATE( ni_f )
 DEALLOCATE( phi_f )
 DEALLOCATE( p_f )
 DEALLOCATE( u_f )
 DEALLOCATE( gradPhiX_f )
 DEALLOCATE( gradPhiY_f )
 DEALLOCATE( lapPhi_f )

 END IF

 RETURN
 END SUBROUTINE MemAlloc

!******************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : Init
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase
!! LBM. Smeared interface initialized using equilibrium order parameter function
!! for each drop defined in the input (in the range [R-IntWidth,R+IntWidth]).
!! The distribution functions f and g are initialized to their equilibrium
!! values for zero velocity.

!-------------------------------------------------------------------------------

 SUBROUTINE Init

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams, ONLY : Eg0, Eg1, Eg2, f, g
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js, m, xs, ys
 REAL(KIND = DBL) :: R, rhon, phin, muPhi, lapPhi
 REAL(KIND = DBL) :: c(0:8)
 REAL(KIND = DBL) :: Af0, Af1
 REAL(KIND = DBL) :: Ag0, Ag1, Eg1A, Eg2A
 integer, parameter:: opposite(0:8) = (/0,3,4,1,2,7,8,5,6/)

! Initialize counters
 now   = 0
 nxt   = 1
 iStep = 0
 tCall = 1
 eps   = 1.0D0

! Set vtk limits
 NX = xmax - xmin + 1
 NY = ymax - ymin + 1

 NX_f = xmax_f - xmin_f + 1
 NY_f = ymax_f - ymin_f + 1

! Set initial average density
 rhon = 0.5D0*(rhoH + rhoL)

!--------- Initialize the order parameter and near neigbor list ----------------
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Initialize phi assuming heavy fluid is continuous phase

      phi_f(i,j) = -phistar
      DO m = 1, nBubbles
       R =  DSQRT ( ((DBLE(i)-bubbles(m,1))**2) + ( (DBLE(j)-bubbles(m,2))**2))
        IF ( R <= (DBLE(bubbles(m,3)) + IntWidth)) THEN
        phi_f(i,j) = +phistar*(TANH( 2.D0*( DBLE(bubbles(m,3)) - R )/IntWidth ))
        END IF
      END DO
      
! Assign near neigbors
      ni_f(i,j,1) = i + 1
      ni_f(i,j,2) = j + 1
      ni_f(i,j,3) = i - 1
      ni_f(i,j,4) = j - 1

    END DO
 END DO
 
! Correct near neighbours at domain edges (periodic boundary conditions)
 ni_f(xmin_f,:,3) = xmax_f
 ni_f(xmax_f,:,1) = xmin_f
 ni_f(:,ymin_f,4) = ymax_f
 ni_f(:,ymax_f,2) = ymin_f
! DO j = ymin_f, ymax_f
!DO i = 1,4
!c(i)= ni_f(xmin_f,j,i)

!ENDDO
!DO i = 1,4
!ni_f(xmin_f,j,i)=c(opposite(i))

!ENDDO

!ENDDO

!DO j = ymin_f, ymax_f
!DO i = 1,4
!c(i)= ni_f(xmax_f,j,i)

!ENDDO
!DO i = 1,4
!ni_f(xmax_f,j,i)=c(opposite(i))

!ENDDO

!ENDDO

 	  DO j = ymin_f, ymax_f
	  DO i = xmin_f, xmax_f
	    u_f(i,j,1:2) = 0.0D0
        u_f(i,j,2) = -0.0D-5 * ((abs(1 - dble(i-1) / (dble(xmax_f-1) * 0.5D0)))**2 - 1.0D0)

		END DO
      END DO
!---------- Initialize the order parameter distribution function f -------------
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Local value of the order parameter
     phin  = phi_f(i,j)

! Identify near neigbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Laplacian of the order parameter
     lapPhi_f(i,j) = ( phi_f(ie,jn) + phi_f(ie,js) + phi_f(iw,jn)    &
                   + phi_f(iw,js) + 4.D0*( phi_f(ie,j) + phi_f(iw,j) &
                   + phi_f(i,jn) + phi_f(i,js) ) - 20.D0*phin )*inv6

! Local value of the chemical potential
     muPhi = alpha4*phin*( phin*phin - phistar2 ) - kappa*lapPhi_f(i,j)

! Equilibrium coefficients
     Af1 = 0.5D0*Gamma*muPhi
     Af0 = -2.D0*Gamma*muPhi

! Equilibrium distribution for zero initial velocity
     f(i,j,0)   = Af0 + phin
!	 f(i,j,1:4) = Af1
     f(i,j,1) = Af1+invEta2*phin*u_f(i,j,1)
	 f(i,j,2) = Af1+invEta2*phin*u_f(i,j,2)
	 f(i,j,3) = Af1-invEta2*phin*u_f(i,j,1)
	 f(i,j,4) = Af1-invEta2*phin*u_f(i,j,2)

   END DO
 END DO
 DO i = xmin_f+1,xmax_f-1

f(i,ymin_f,2) = f(i,ymin_f,4)+(phin/eta)*u_f(i,ymin_f,2)
f(i,ymax_f,4) = f(i,ymax_f,2)+(phin/eta)*u_f(i,ymax_f,2)
END DO

!---------- Initialize g, neighbors in momentum mesh, pressure and velocity ----
 	  DO j = ymin, ymax
	  DO i = xmin, xmax
	    u(i,j,1) = 0.0D0
        u(i,j,2) = -0.0D-6 * ((abs(1 - dble(i-1) / (dble(xmax-1) * 0.5D0)))**2 - 1.0D0)

		END DO
      END DO

 DO j = ymin, ymax
   DO i = xmin, xmax

! Identify near neighbors
      ni(i,j,1) = i + 1
      ni(i,j,2) = j + 1
      ni(i,j,3) = i - 1
      ni(i,j,4) = j - 1

! Define source nodes in the order parameter mesh
     xs = 2*i - 1
     ys = 2*j - 1

! Local values of the phase and the chemical potential
     phin   = phi_f(xs,ys)
     lapPhi = 4.D0*lapPhi_f(xs,ys)
     muPhi  = alpha4*phin*( phin*phin - phistar2 ) - kappaG*lapPhi

! Coefficients
      Ag1  = rhon + 3.D0*phin*muPhi
      Ag0  = 2.25D0*rhon - 1.25D0*Ag1
      Eg1A = Eg1*Ag1
      Eg2A = Eg2*Ag1
	  
! Equilibrium distribution for zero initial velocity
      g(i,j,0,now) = Eg0*Ag0+Eg0*rhon*(-3.d0*(u(i,j,1)*u(i,j,1)+u(i,j,2)*u(i,j,2)))
      g(i,j,1,now) = Eg1A+Eg1*rhon*(1.5d0*1.0d0*u(i,j,1)+4.5d0*(u(i,j,1)*u(i,j,1))-1.5d0*(u(i,j,1)*u(i,j,1)+u(i,j,2)*u(i,j,2)))
      g(i,j,2,now) = Eg1A+Eg1*rhon*(1.5d0*1.0d0*u(i,j,2)+4.5d0*(u(i,j,2)*u(i,j,2))-1.5d0*(u(i,j,1)*u(i,j,1)+u(i,j,2)*u(i,j,2)))
      g(i,j,3,now) = Eg1A+Eg1*rhon*(-1.5d0*1.0d0*u(i,j,1)+4.5d0*(u(i,j,1)*u(i,j,1))-1.5d0*(u(i,j,1)*u(i,j,1)+u(i,j,2)*u(i,j,2)))
      g(i,j,4,now) = Eg1A+Eg1*rhon*(-1.5d0*1.0d0*u(i,j,2)+4.5d0*(u(i,j,2)*u(i,j,2))-1.5d0*(u(i,j,1)*u(i,j,1)+u(i,j,2)*u(i,j,2)))
      g(i,j,5,now) = Eg2A+Eg2*rhon*(1.5d0*(u(i,j,1)+u(i,j,2))+4.5d0*((u(i,j,1)+u(i,j,2))**2)-1.5d0*(u(i,j,1)*u(i,j,1) &
                        +u(i,j,2)*u(i,j,2)))
      g(i,j,6,now) = Eg2A+Eg2*rhon*(1.5d0*(-u(i,j,1)+u(i,j,2))+4.5d0*((-u(i,j,1)+u(i,j,2))**2)-1.5d0*(u(i,j,1)*u(i,j,1) &
                        +u(i,j,2)*u(i,j,2)))
      g(i,j,7,now) = Eg2A+Eg2*rhon*(-1.5d0*(u(i,j,1)+u(i,j,2))+4.5d0*((u(i,j,1)+u(i,j,2))**2)-1.5d0*(u(i,j,1)*u(i,j,1) &
                        +u(i,j,2)*u(i,j,2)))
      g(i,j,8,now) = Eg2A+Eg2*rhon*(+1.5d0*(u(i,j,1)-u(i,j,2))+4.5d0*((u(i,j,1)-u(i,j,2))**2)-1.5d0*(u(i,j,1)*u(i,j,1) &
                        +u(i,j,2)*u(i,j,2)))

! Initialize pressure and velocity arrays
      p(i,j)   = 1.0D0
!      u(i,j,:) = 0.0D0


   END DO
 END DO
DO i = xmin+1,xmax-1
g(i,ymin,2,now) = g(i,ymin,4,now)+(2.d0/3.d0)*u(i,ymin,2)
g(i,ymin,6,now)= g(i,ymin,8,now)+0.5*(g(i,ymin,1,now)-g(i,ymin,3,now))-0.5d0*rhon*(u(i,ymin,1))+0.167d0*rhon*u(i,ymin,2)
g(i,ymin,5,now)= g(i,ymin,7,now)+0.5*(g(i,ymin,3,now)-g(i,ymin,1,now))+0.5d0*rhon*(u(i,ymin,1))+0.167d0*rhon*u(i,ymin,2)

g(i,ymax,4,now) = g(i,ymax,2,now)-(2.d0/3.d0)*u(i,ymax,2)
g(i,ymax,8,now)= g(i,ymax,6,now)+0.5*(g(i,ymax,3,now)-g(i,ymax,1,now))+0.5d0*rhon*u(i,ymax,1)-0.167d0*rhon*u(i,ymax,2)
g(i,ymax,7,now)= g(i,ymax,5,now)+0.5*(g(i,ymax,1,now)-g(i,ymax,3,now))-0.5d0*rhon*(u(i,ymin,1))-0.167d0*rhon*u(i,ymin,2)
END DO
! Correct near neighbours at domain edges (periodic boundary conditions)
 ni(xmin,:,3) = xmax
 ni(xmax,:,1) = xmin
 ni(:,ymin,4) = ymax
 ni(:,ymax,2) = ymin
!DO j = ymin, ymax
!DO i = 1,8
!c(i)= g(xmin,j,i,now)

!ENDDO
!DO i = 1,8
!g(xmin,j,i,now)=c(opposite(i))

!ENDDO

!ENDDO

!DO j = ymin, ymax
!DO i = 1,8
!c(i)= g(xmax,j,i,now)

!ENDDO
!DO i = 1,8
!g(xmax,j,i,now)=c(opposite(i))

!ENDDO

!ENDDO
! Calculate differential terms needed in Collision
 CALL Differentials

! Save initialized data
 CALL Stats
 CALL VtkPlane

 RETURN
 END SUBROUTINE Init
!*****************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : Stats
!-------------------------------------------------------------------------------
!> @file
!! Calculate intermediate simulation results and save to file 'stats.out'
!> @details
!! Calculate average velocity, mass conservation factor, effective radius of the
!! drop, pressure difference between the inside and the outside of the drop and
!! the error with respect to the analytical value given by Laplace's equation,
!! and write them to file "stats.out" for the dual grid D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM.
!!
!! The effective radius is calculated assuming the drop is a perfect circle with
!! area given by A = Pi*R*R.
!!
!! The pressure inside and the pressure outside of the drop are calculated as
!! the average pressures inside and outside the drop, excluding the interface
!! area

!-------------------------------------------------------------------------------

 SUBROUTINE Stats

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain,      ONLY : invInitVol, invPi, iStep, tCall, xmax_f, xmin_f, ymax_f, ymin_f
 USE FluidParams, ONLY : bubbles, Convergence, eps, IntWidth, p_f, phi_f, sigma, u_f,delt
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, nodesIn, nodesOut
 INTEGER :: IO_ERR
 REAL(KIND = DBL) :: Pin, Pout, Pdif, Perr, Ro, R1, Ref, Ux, Uy, Vef, Vol

! Initialize
 Ux   = 0.D0
 Uy   = 0.D0
 Vol  = 0.D0
 Pin  = 0.D0
 Pout = 0.D0
 nodesIn  = 0
 nodesOut = 0
 Ro = bubbles(1,3)

! Loop through all nodes inside the bubble
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

     IF ( phi_f(i,j) >= 0.D0 ) THEN
       Ux  = Ux  + u_f(i,j,1)
       Uy  = Uy  + u_f(i,j,2)
       Vol = Vol + 1.D0
     END IF

     R1 =  DSQRT((DBLE(i)-bubbles(1,1))**2 + (DBLE(j)-bubbles(1,2))**2)
     IF (R1 < (Ro - IntWidth)) THEN
       Pin = Pin + p_f(i,j)
       nodesIn = nodesIn + 1
     ELSE IF (R1 > (Ro + IntWidth)) THEN
       Pout = Pout + p_f(i,j)
       nodesOut = nodesOut + 1
     END IF

   END DO
 END DO

! Save the initial volume during the first time step
 IF( iStep == 0 ) invInitVol = 1.D0/Vol

! Average velocity and pahse conservation (re-scale Ref to momentum grid)
 Ux  = Ux/Vol
 Uy  = Uy/Vol
 Ref = 0.5D0*DSQRT( Vol*invPi )
 Vef = Vol*invInitVol

! Compliance with Laplace's Law
 Pin  = Pin/DBLE(nodesIn)
 Pout = Pout/DBLE(nodesOut)
 Pdif = Pin - Pout
 Perr = (sigma/Ro - Pdif)*Ro/sigma

 OPEN(UNIT = 10, FILE = "stats.out", STATUS = "UNKNOWN", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(10,'(I9,6ES19.9)')istep,Ux,Uy,Vef,Ref,Pdif,Perr
   CLOSE(UNIT = 10)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output file 'stats.out'."
 END IF

! Analyze convergence - Done over a certain period of time ( 10 x tStats ) to 
! ensure that long range fluctuations do not give false positives
 Convergence(tCall) = Pdif
 IF ( MOD(tCall,11) == 0 ) THEN
   eps   = 0.D0
   tCall = 1
   DO i = 2, 11
     eps = eps + DABS( Convergence(i) - Convergence(i-1) )
   END DO
   eps = eps*0.1D0/Pdif
 ELSE
   tCall = tCall + 1
 END IF

 RETURN
 END SUBROUTINE Stats
!******************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : VtkPlane
!-------------------------------------------------------------------------------
!> @file
!! Save simulation data in VTK format
!> @details
!! Save order parameter, pressure, and velocity data in VTK format for the
!! dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM. Filenames are
!! 'DATA_XXXXXXX.VTK' where XXXXXXX is the time step. These files are understood
!! by Paraview.

!-------------------------------------------------------------------------------

 SUBROUTINE VtkPlane

! Common Variables
 USE Domain,      ONLY : NX, NY, xmax, xmin, ymax, ymin
 USE FluidParams, ONLY : p, phi_f, u
 IMPLICIT NONE

! Required functions
 CHARACTER(8) :: filer

! Local Variables
 INTEGER :: i, j
 INTEGER :: IO_ERR
 CHARACTER(8)  :: filer1
 CHARACTER(16) :: fmtd1, fmtd3
 CHARACTER(17) :: filer2

!-------- Set output files formats and names -----------------------------------
 fmtd1 = '(1(ES21.11E3))'
 fmtd3 = '(3(ES21.11E3))'
 filer1 = filer()
 filer2 = "DATA"//filer1(1:8)//".vtk"

 OPEN(UNIT = 12, FILE = filer2, STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(12,'(A)')"# vtk DataFile Version 2.0"
   WRITE(12,'(A)')"MP-LABS v1.0"
   WRITE(12,'(A)')"ASCII"
   WRITE(12,*)
   WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
   WRITE(12,*)"DIMENSIONS",NX,NY,1
   WRITE(12,*)"ORIGIN",xmin,ymin,1
   WRITE(12,*)"SPACING",1,1,1
   WRITE(12,*)
   WRITE(12,*)"POINT_DATA",NX*NY
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Phi double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO j = ymin, ymax
       DO i = xmin, xmax
          WRITE(UNIT = 12,FMT = fmtd1)phi_f(2*i-1,2*j-1)
       END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Pressure double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO j = ymin, ymax
       DO i = xmin, xmax
          WRITE(UNIT = 12,FMT = fmtd1)p(i,j)
       END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"VECTORS Velocity double"
   DO j = ymin, ymax
       DO i = xmin, xmax
         WRITE(UNIT = 12,FMT = fmtd3)u(i,j,1),u(i,j,2),0.D0
      END DO
   END DO
   CLOSE(UNIT = 12)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output vtk file."
 END IF

 RETURN
 END SUBROUTINE VtkPlane


!-------------------------------------------------------------------------------
! Function : filer
!-------------------------------------------------------------------------------
!> @brief Set vtk file name as a function of processor number and time step.
 FUNCTION filer()

 USE DOMAIN, ONLY: iStep
 IMPLICIT NONE

! Function type
 CHARACTER(8) :: filer


!-------- Decide how to name the files depending on the timestep (iStep) -------
   IF (iStep < 10) THEN
     WRITE(filer,9001)iStep
   ELSE IF (iStep < 100) THEN
     WRITE(filer,9002)iStep
   ELSE IF (iStep < 1000) THEN
     WRITE(filer,9003)iStep
   ELSE IF (iStep < 10000) THEN
     WRITE(filer,9004)iStep
   ELSE IF (iStep < 100000) THEN
     WRITE(filer,9005)iStep
   ELSE IF (iStep < 1000000) THEN
     WRITE(filer,9006)iStep
   ELSE IF (iStep < 10000000) THEN
     WRITE(filer,9007)iStep
   END IF

!-------- Define the file name format depending on the timestep (iStep) --------
 9001 FORMAT('_000000',i1)
 9002 FORMAT('_00000',i2)
 9003 FORMAT('_0000',i3)
 9004 FORMAT('_000',i4)
 9005 FORMAT('_00',i5)
 9006 FORMAT('_0',i6)
 9007 FORMAT('_',i7)

 RETURN
 END FUNCTION filer
!*********************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : CollisionG
!-------------------------------------------------------------------------------
!> @file
!! Collision and streaming steps for distribution function g
!> @details
!! Collision and streaming steps for the momentum distribution function g and
!! calculation of updated macroscopic quantities (velocity, pressure and order
!! parameter) in the dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

!-------------------------------------------------------------------------------

 SUBROUTINE CollisionG

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : ni, now, nxt, xmax, xmin, ymax, ymin
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js, xs, ys
 integer, parameter:: opposite(0:8) = (/0,3,4,1,2,7,8,5,6/)
 REAL(KIND = DBL) :: sFx, sFy, iFs
 REAL(KIND = DBL) :: Vg, Vsq, Vg2, ux, uy,c(0:8)
 REAL(KIND = DBL) :: rhon, invrho, phin, phin2, muPhi
 REAL(KIND = DBL) :: gradPhiX, gradPhiY, gradPhiSq, lapPhi
 REAL(KIND = DBL) :: Ag0, Ag1, Eg1A, Eg2A, Eg1R, Eg2R, geq


!---------- Hydrodynamics, collision and stream for g in the same loop ---------
 DO j = ymin, ymax
   DO i = xmin, xmax

! Identify near neighbors (for the streaming step)
     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

! Define source nodes in the order parameter mesh
     xs = 2*i - 1
     ys = 2*j - 1

! Copy and re-scale the values of the order parameter and its differentials
     phin      = phi_f(xs,ys)
     phin2     = phin*phin
     lapPhi    = 4.D0*lapPhi_f(xs,ys)
     gradPhiX  = 2.D0*gradPhiX_f(xs,ys)
     gradPhiY  = 2.D0*gradPhiY_f(xs,ys)
     gradPhiSq = gradPhiX*gradPhiX + gradPhiY*gradPhiY

! Local value of the chemical potential
     muPhi = alpha4*phin*( phin2 - phistar2 ) - kappaG*lapPhi

!---------- Hydrodynamics ------------------------------------------------------
!  Local density value
     rhon   = g(i,j,0,now) + g(i,j,1,now) + g(i,j,2,now) + g(i,j,3,now) &
            + g(i,j,4,now) + g(i,j,5,now) + g(i,j,6,now) + g(i,j,7,now) &
            + g(i,j,8,now)
     invRho = 1.D0/rhon

! Interfacial force
     sFx = muPhi*gradPhiX
	 if (phin.gt.0.0d0)then
     sFy = muPhi*gradPhiY+2.0d0*phistar*3.8D-5
	 else
	 sFy = muPhi*gradPhiY
	 endif

! Velocity field at each node
     ux = (g(i,j,1,now) - g(i,j,3,now) + g(i,j,5,now) - g(i,j,6,now) &
        - g(i,j,7,now) + g(i,j,8,now) + 0.5D0*sFx)*invRho

     uy = (g(i,j,2,now) - g(i,j,4,now) + g(i,j,5,now) + g(i,j,6,now) &
        - g(i,j,7,now) - g(i,j,8,now) + 0.5D0*sFy)*invRho

     u(i,j,1) = ux
     u(i,j,2) = uy

! Pressure at each node
     p(i,j) = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 )&
            - kappaG*( phin*lapPhi + 0.5D0*gradPhiSq ) + Cs_sq*rhon

!---------- Collision and Stream for momentum distribution function ------------
     Ag1  = rhon + 3.D0*phin*muPhi
     Ag0  = 2.25D0*rhon - 1.25D0*Ag1
     Eg1A = Eg1*Ag1
     Eg2A = Eg2*Ag1
     Eg1R = Eg1C*rhon
     Eg2R = Eg2C*rhon
     Vsq  = 0.5D0*( ux*ux + uy*uy )

! Direction 0
     geq = Eg0*( Ag0 - rhon*3.D0*Vsq )
     iFs = Eg0T*( (-ux)*sFx - uy*sFy )
     g(i,j,0,nxt) = g(i,j,0,now) + invTauRho*delt*( geq - g(i,j,0,now) ) + iFs

! Direction 1
     geq = Eg1A + Eg1R*( ux + 1.5D0*ux*ux - Vsq )
     iFs = Eg1T*( (1.D0 + 2.D0*ux)*sFx - uy*sFy )
     g(ie,j,1,nxt) = g(i,j,1,now) + invTauRho*delt*( geq - g(i,j,1,now) ) + iFs

! Direction 2
     geq = Eg1A + Eg1R*( uy + 1.5D0*uy*uy - Vsq )
     iFs = Eg1T*( (-ux)*sFx + (1.D0 + 2.D0*uy)*sFy )
     g(i,jn,2,nxt) = g(i,j,2,now) + invTauRho*delt*( geq - g(i,j,2,now) ) + iFs

! Direction 3
     geq = Eg1A + Eg1R*( (-ux) + 1.5D0*ux*ux - Vsq )
     iFs = Eg1T*( (2.D0*ux - 1.D0)*sFx - uy*sFy )
     g(iw,j,3,nxt) = g(i,j,3,now) + invTauRho*delt*( geq - g(i,j,3,now) ) + iFs

! Direction 4
     geq = Eg1A + Eg1R*( (-uy) + 1.5D0*uy*uy - Vsq )
     iFs = Eg1T*( (-ux)*sFx + (2.D0*uy - 1)*sFy )
     g(i,js,4,nxt) = g(i,j,4,now) + invTauRho*delt*( geq - g(i,j,4,now) ) + iFs

! Direction 5
     Vg = ux + uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (1.D0 - ux) + Vg2 )*sFx + ( (1.D0 - uy) + Vg2 )*sFy )
     g(ie,jn,5,nxt) = g(i,j,5,now) + invTauRho*delt*( geq - g(i,j,5,now) ) + iFs

! Direction 6
     Vg = -ux + uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (-1.D0 - ux) - Vg2 )*sFx + ( (1.D0 - uy) + Vg2 )*sFy )
     g(iw,jn,6,nxt) = g(i,j,6,now) + invTauRho*delt*( geq - g(i,j,6,now) ) + iFs

! Direction 7
     Vg = -ux - uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (-1.D0 - ux) - Vg2 )*sFx + ( (-1.D0 - uy) - Vg2 )*sFy )
     g(iw,js,7,nxt) = g(i,j,7,now) + invTauRho*delt*( geq - g(i,j,7,now) ) + iFs

! Direction 8
     Vg = ux - uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (1.D0 - ux) + Vg2 )*sFx + ( (-1.D0 - uy) - Vg2 )*sFy )
     g(ie,js,8,nxt) = g(i,j,8,now) + invTauRho*delt*( geq - g(i,j,8,now) ) + iFs

   END DO
 END DO


 RETURN
 END SUBROUTINE CollisionG

!**********************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : Update
!-------------------------------------------------------------------------------
!> @file
!! Update velocity and pressure in the order parameter grid
!> @details
!! Update velocity and pressure in the order parameter grid for the dual grid
!! D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM. This function uses a bilinear
!! interpolation of the values in the momentum grid to update the order
!! parameter grid values.

!-------------------------------------------------------------------------------

 SUBROUTINE Update

! Common Variables
 USE Domain,      ONLY : ni, ni_f, xmax, xmin, ymax, ymin
 USE FluidParams, ONLY : p, p_f, u, u_f
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, xs1, xs2, ys1, ys2, xt1, xt2, xt3, yt1, yt2, yt3


 DO j = ymin, ymax-1
   DO i = xmin, xmax-1

! Coordinates of source nodes for the interpolation (xs,ys)
     xs1 = i
     ys1 = j
     xs2 = ni(i,j,1)
     ys2 = ni(i,j,2)

! Coordinates of target nodes for the interpolation (xt,yt)
     xt1 = 2*i - 1
     yt1 = 2*j - 1
     xt2 = ni_f(xt1,yt1,1)
     yt2 = ni_f(xt1,yt1,2)
     xt3 = ni_f(xt2,yt2,1)
     yt3 = ni_f(xt2,yt2,2)
 
! Target 1
     u_f(xt1,yt1,1) = u(xs1,ys1,1)
     u_f(xt1,yt1,2) = u(xs1,ys1,2)
     p_f(xt1,yt1)   = p(xs1,ys1)

! Target 2
     u_f(xt2,yt1,1) = 0.5D0*( u(xs1,ys1,1) + u(xs2,ys1,1) )
     u_f(xt2,yt1,2) = 0.5D0*( u(xs1,ys1,2) + u(xs2,ys1,2) )
     p_f(xt2,yt1)   = 0.5D0*( p(xs1,ys1)   + p(xs2,ys1)   )

! Target 3
     u_f(xt3,yt1,1) = u(xs2,ys1,1)
     u_f(xt3,yt1,2) = u(xs2,ys1,2)
     p_f(xt3,yt1)   = p(xs2,ys1)

! Target 4
    u_f(xt1,yt2,1) = 0.5D0*( u(xs1,ys1,1) + u(xs1,ys2,1) )
    u_f(xt1,yt2,2) = 0.5D0*( u(xs1,ys1,2) + u(xs1,ys2,2) )
    p_f(xt1,yt2)   = 0.5D0*( p(xs1,ys1)   + p(xs1,ys2)   )

! Target 5
    u_f(xt2,yt2,1) = 0.25D0*( u(xs1,ys1,1) + u(xs2,ys1,1) + u(xs1,ys2,1) &
                   + u(xs2,ys2,1) )
    u_f(xt2,yt2,2) = 0.25D0*( u(xs1,ys1,2) + u(xs2,ys1,2) + u(xs1,ys2,2) &
                   + u(xs2,ys2,2) )
    p_f(xt2,yt2)   = 0.25D0*( p(xs1,ys1)   + p(xs2,ys1)   + p(xs1,ys2)   &
                   + p(xs2,ys2)   )

! Target 6
    u_f(xt3,yt2,1) = 0.5D0*( u(xs2,ys1,1) + u(xs2,ys2,1) )
    u_f(xt3,yt2,2) = 0.5D0*( u(xs2,ys1,2) + u(xs2,ys2,2) )
    p_f(xt3,yt2)   = 0.5D0*( p(xs2,ys1)   + p(xs2,ys2)   )

! Target 7
    u_f(xt1,yt3,1) = u(xs1,ys2,1)
    u_f(xt1,yt3,2) = u(xs1,ys2,2)
    p_f(xt1,yt3)   = p(xs1,ys2)

! Target 8
    u_f(xt2,yt3,1) = 0.5D0*( u(xs1,ys2,1) + u(xs2,ys2,1) )
    u_f(xt2,yt3,2) = 0.5D0*( u(xs1,ys2,2) + u(xs2,ys2,2) )
    p_f(xt2,yt3)   = 0.5D0*( p(xs1,ys2)   + p(xs2,ys2)   )

! Target 9
    u_f(xt3,yt3,1) = u(xs2,ys2,1)
    u_f(xt3,yt3,2) = u(xs2,ys2,2)
    p_f(xt3,yt3)   = p(xs2,ys2)

   END DO
 END DO 

 RETURN
 END SUBROUTINE Update
!*********************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : CollisionF
!-------------------------------------------------------------------------------
!> @file
!! Collision step for distribution function f
!> @details
!! Collision step for the order parameter distribution function f in the Dual
!! Grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

!-------------------------------------------------------------------------------

 SUBROUTINE CollisionF

! Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain,      ONLY : xmax_f, xmin_f, ymax_f, ymin_f
 USE FluidParams
 USE LBMParams,   ONLY : f, fcol
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j
 !integer, parameter:: opposite(0:4) = (/0,3,4,1,2/)
 REAL(KIND = DBL) :: ux, uy,c(0:4)
 REAL(KIND = DBL) :: phin, muPhi
 REAL(KIND = DBL) :: Af0, Af1, Cfp


 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Local values of phi and the chemical potential
     ux    = u_f(i,j,1)
     uy    = u_f(i,j,2)
     phin  = phi_f(i,j)
     muPhi = alpha4*phin*( phin*phin - phistar2 ) - kappa*lapPhi_f(i,j)

! Collision
     Af1 = 0.5D0*Gamma*muPhi
     Af0 = -2.D0*Gamma*muPhi
     Cfp = invEta2*phin
     f(i,j,0)    = f(i,j,0) + invTauPhi*( Af0 + phin   - f(i,j,0) )
     fcol(i,j,1) = f(i,j,1) + invTauPhi*( Af1 + Cfp*ux - f(i,j,1) )
     fcol(i,j,2) = f(i,j,2) + invTauPhi*( Af1 + Cfp*uy - f(i,j,2) )
     fcol(i,j,3) = f(i,j,3) + invTauPhi*( Af1 - Cfp*ux - f(i,j,3) )
     fcol(i,j,4) = f(i,j,4) + invTauPhi*( Af1 - Cfp*uy - f(i,j,4) )

   END DO
 END DO



 RETURN
 END SUBROUTINE CollisionF

!******************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : Stream
!-------------------------------------------------------------------------------
!> @file
!! Relaxation and streaming step for distribution function f.
!> @details
!! Relaxation and streaming step for distribution function f in the dual grid
!! D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

!-------------------------------------------------------------------------------

 SUBROUTINE Stream

! Common Variables
 USE Domain,      ONLY : ni_f, xmax_f, xmin_f, ymax_f, ymin_f
 USE LBMParams,   ONLY : f, fcol
 USE FluidParams, ONLY : eta, eta2
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js

 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

     f(ie, j,1) = eta*fcol(i,j,1) + eta2*fcol(ie,j,1)
     f( i,jn,2) = eta*fcol(i,j,2) + eta2*fcol(i,jn,2)
     f(iw, j,3) = eta*fcol(i,j,3) + eta2*fcol(iw,j,3)
     f( i,js,4) = eta*fcol(i,j,4) + eta2*fcol(i,js,4)
   END DO
 END DO

 RETURN
 END SUBROUTINE Stream
!**********************************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : Differentials
!-------------------------------------------------------------------------------
!> @file
!! Order parameter calculation and differential terms for the interfacial force.
!> @details
!! Calculate the order parameter phi_f using updated values of the distribution
!! function f and then obtain the differential terms necessary for the
!! interfacial force calculation (the Laplacian and the gradient of the order
!! parameter phi_f) in the dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.
!! The calculation is done directly on the order parameter grid.

!-------------------------------------------------------------------------------

 SUBROUTINE Differentials

! Common Variables
 USE Domain,      ONLY : inv12, inv6, ni_f, xmax_f, xmin_f, ymax_f, ymin_f
 USE FluidParams, ONLY : gradPhiX_f, gradPhiY_f, lapPhi_f, phi_f
 USE LBMParams,   ONLY : f
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js

! Order parameter
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f
     phi_f(i,j) = f(i,j,0) + f(i,j,1) + f(i,j,2) + f(i,j,3) + f(i,j,4)
   END DO
 END DO

! Differentials needed for the interfacial terms 
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Laplacian
     lapPhi_f(i,j) = ( phi_f(ie,jn) + phi_f(ie,js) + phi_f(iw,jn)          &
                   + phi_f(iw,js) + 4.D0*( phi_f(ie,j) + phi_f(iw,j)       &
                   + phi_f(i,jn) + phi_f(i,js) ) - 20.D0*phi_f(i,j) )*inv6
! Gradients
     gradPhiX_f(i,j)  = ( 4.D0*( phi_f(ie,j) - phi_f(iw,j) ) + phi_f(ie,jn) &
                      - phi_f(iw,js) + phi_f(ie,js) - phi_f(iw,jn) )*inv12

     gradPhiY_f(i,j)  = ( 4.D0*( phi_f(i,jn) - phi_f(i,js) ) + phi_f(ie,jn) &
                      - phi_f(iw,js) - phi_f(ie,js) + phi_f(iw,jn) )*inv12

   END DO
 END DO

 RETURN
 END SUBROUTINE Differentials

!*********************************************************************************
!-------------------------------------------------------------------------------
! Subroutine : FinalDump
!-------------------------------------------------------------------------------
!> @file
!! Save relevant data at the end of the simulation run to file 'final.out'
!> @details
!! Generates the final data file for the simulation in the dual grid D2Q5/D2Q9
!! Zheng-Shu-Chew multiphase LBM, which contains:
!!
!! - Input parameters
!! - Estimated memory usage
!! - Pressure difference between the inside and the outside of the drop
!! - Error in the verification of Laplace's Law for the pressure
!! - Mass conservation factor
!! - Effective drop radius
!! - Maximum velocity in the domain

!-------------------------------------------------------------------------------

 SUBROUTINE FinalDump

!  Common Variables
 USE NTypes, ONLY : DBL
 USE Domain
 USE FluidParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, nodesIn, nodesOut
 INTEGER :: IO_ERR
 REAL(KIND = DBL) :: Pin, Pout, Ro, R1, Pdif, Perr, Uloc, Umax, Ref, Vef, Vol
 REAL(KIND = DBL) :: distroMem, auxMem, totalMem, memUnitF, memUnitG

! Initialize
 Pin  =  0.D0
 Pout =  0.D0
 Vol  =  0.D0
 Umax = -1.D0
 nodesIn  = 0
 nodesOut = 0
 Ro = bubbles(1,3)

! Pressure inside and outside the bubble, maximum velocity and effective radius
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

     R1 =  DSQRT( ( DBLE(i)-bubbles(1,1) )**2 + ( DBLE(j)-bubbles(1,2) )**2 )
     IF ( R1 < (Ro - IntWidth) ) THEN
       Pin = Pin + p_f(i,j)
       nodesIn = nodesIn + 1
     ELSE IF ( R1 > (Ro + IntWidth) ) THEN
       Pout = Pout + p_f(i,j)
       nodesOut = nodesOut + 1
     END IF

     Uloc = DSQRT( u_f(i,j,1)*u_f(i,j,1) + u_f(i,j,2)*u_f(i,j,2) )
     IF ( Uloc > Umax ) Umax = Uloc
     IF ( phi_f(i,j) >= 0.D0 ) Vol = Vol + 1.D0

   END DO
 END DO

! Calculate compliance with Laplace Law
 Pin  = Pin/DBLE(nodesIn)
 Pout = Pout/DBLE(nodesOut)
 Pdif = Pin - Pout
 Perr = (sigma/Ro - Pdif)*Ro/sigma

! Calculate phase conservation
 Ref = DSQRT( Vol*invPi )
 Vef = Vol*invInitVol

! Estimate memory usage (Mb)
 memUnitG = NX*NY/( 1024.D0*1024.D0 )
 memUnitF = NX_f*NY_f/( 1024.D0*1024.D0 )
 distroMem = 8.D0*( 18.D0*memUnitG + 10.D0*memUnitF )
 auxMem    = 8.D0*( 3.D0*( memUnitG + memUnitF ) + 4.D0*memUnitF ) &
           + 4.D0*( 4.D0*( memUnitG + memUnitF ) )
 totalMem  = distroMem + auxMem

! The effective radius, sigma and IntWidth are re-scaled to the momentum grid
 OPEN(UNIT = 10, FILE = "final.out", STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(10,'(A)')'*** Multiphase Zheng-Shu-Chew LBM 2D Simulation ***'
   WRITE(10,'(A)')'*** Dual Grid Implementation (Serial)           ***'
   WRITE(10,*)
   WRITE(10,'(A)')'INPUT PARAMETERS'
   WRITE(10,'(A,I9)')'Total Iterations      = ',MaxStep+1
   WRITE(10,'(A,I9)')'Relaxation Iterations = ',RelaxStep+1
   WRITE(10,'(A,I9)')'Length in X Direction = ',xmax
   WRITE(10,'(A,I9)')'Length in Y Direction = ',ymax
   WRITE(10,'(A,ES15.5)')'Interface Width   = ',0.5D0*IntWidth
   WRITE(10,'(A,ES15.5)')'Interface Tension = ',0.5D0*sigma
   WRITE(10,'(A,ES15.5)')'Interface Mobility= ',Gamma
   WRITE(10,'(A,ES15.5)')'RhoL    = ',rhoL
   WRITE(10,'(A,ES15.5)')'RhoH    = ',rhoH
   WRITE(10,'(A,ES15.5)')'TauRho  = ',tauRho
   WRITE(10,'(A,ES15.5)')'TauPhi  = ',tauPhi
   WRITE(10,*)
   WRITE(10,'(A)')'MEMORY USAGE (Mb)'
   WRITE(10,'(A,ES15.5)')'Distributions     = ',distroMem
   WRITE(10,'(A,ES15.5)')'Auxiliary Arrays  = ',auxMem
   WRITE(10,'(A,ES15.5)')'Total Memory Used = ',totalMem
   WRITE(10,*)
   WRITE(10,'(A)')'OUTPUT RESULTS (RELAXATION)'
   WRITE(10,'(A,ES19.9)')'Effective Radius   = ',0.5D0*Ref
   WRITE(10,'(A,ES19.9)')'Phase Conservation = ',Vef
   WRITE(10,'(A,ES19.9)')'(Pin - Pout)       = ',Pdif
   WRITE(10,'(A,ES19.9)')'Laplace Error      = ',Perr
   WRITE(10,'(A,ES19.9)')'Parasitic Velocity = ',Umax
   WRITE(10,*)
   WRITE(10,'(A)')'***       Simulation Finished Succesfully       ***'
   CLOSE(UNIT = 10)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: unable to open output file 'final.out'."
 END IF

 RETURN
 END SUBROUTINE FinalDump
!*******************************************************************************************
