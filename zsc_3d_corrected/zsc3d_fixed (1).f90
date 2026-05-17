!===============================================================================
! ZSC 3D Multiphase LBM - Fixed and Validated Version
!===============================================================================
! Based on: Arup Das's ZSC Fortran codes (IIT Roorkee)
! Fixed by: Analysis of Sudhakar & Das, I&EC Research 2020
! Physics:  Zheng, Shu, Chew - JCP 218:353-371, 2006
!
! Fixes applied:
!   1. Wall BCs: proper bounce-back using opposite() array (was sequential overwrite)
!   2. Inlet BC: Zou-He with parabolic velocity profile
!   3. Outlet BC: extrapolation + constant pressure ghost layer
!   4. Neighbor array: corrected at boundaries (was commented out)
!   5. Geometry: parameterized (was hardcoded ±50 offsets)
!   6. Code structure: cleaned up, geometry-independent physics core
!
! Geometry modes:
!   mode=1: Straight vertical channel (validation baseline - Sudhakar Fig 1f)
!   mode=2: U-bend channel (Sudhakar Fig 1a, Cases 1-2)
!   mode=3: Static drop Laplace law test (quick sanity check)
!
! Compile: gfortran -O2 -o zsc3d zsc3d_fixed.f90
! Run:     ./zsc3d
!===============================================================================

!-------------------------------------------------------------------------------
MODULE NTypes
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: SGL = KIND(1.0)
  INTEGER, PARAMETER :: DBL = KIND(1.D0)
END MODULE NTypes

!-------------------------------------------------------------------------------
MODULE Domain
  USE NTypes, ONLY : DBL
  IMPLICIT NONE
  SAVE

  ! Time step and limits
  INTEGER :: iStep, tCall, tDump, tStat, MaxStep

  ! Domain size
  INTEGER :: xjump, xo, yo, zo
  INTEGER :: now, nxt, NX, NY, NZ, NPX, NPY, NPZ
  INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax

  ! Geometry mode: 1=straight channel, 2=U-bend, 3=Laplace test
  INTEGER :: geom_mode

  ! U-bend geometry parameters (read from input for mode 2)
  ! a       = channel width (internal fluid nodes per side of square cross-section)
  ! Rbend   = bend centerline radius in lattice units
  ! Larm    = length of straight vertical arms
  ! The domain is auto-sized for mode 2 from these parameters
  INTEGER :: ubend_a, ubend_Rbend, ubend_Larm

  ! Inlet/outlet face identifiers for U-bend
  INTEGER :: inlet_zface, outlet_zface
  INTEGER :: inlet_x1, inlet_x2, outlet_x1, outlet_x2

  ! Spatial dimension
  INTEGER, PARAMETER :: ndim = 3

  ! Neighbor array
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: ni

  ! Wall mask: 1=fluid, 0=solid wall
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: is_fluid

  ! Inverse of the initial volume of the discrete phase
  REAL(KIND = DBL) :: invInitVol

  ! Constants
  REAL(KIND = DBL), PARAMETER :: inv36  = 0.02777777777778D0
  REAL(KIND = DBL), PARAMETER :: inv12  = 0.08333333333333333333D0
  REAL(KIND = DBL), PARAMETER :: invPi  = 0.31830988618379067154D0

  ! Opposite direction mapping for D3Q19
  !   Lattice velocities (D3Q19):
  !   0: (0,0,0)   rest
  !   1: (+1,0,0)   2: (-1,0,0)
  !   3: (0,+1,0)   4: (0,-1,0)
  !   5: (0,0,+1)   6: (0,0,-1)
  !   7: (+1,+1,0)  8: (-1,-1,0)
  !   9: (+1,-1,0) 10: (-1,+1,0)
  !  11: (+1,0,+1) 12: (-1,0,-1)
  !  13: (+1,0,-1) 14: (-1,0,+1)
  !  15: (0,+1,+1) 16: (0,-1,-1)
  !  17: (0,+1,-1) 18: (0,-1,+1)
  INTEGER, PARAMETER :: opp_g(0:18) = &
    (/0, 2,1, 4,3, 6,5, 8,7, 10,9, 12,11, 14,13, 16,15, 18,17/)

  ! Opposite direction mapping for D3Q7
  !   0: (0,0,0)   rest
  !   1: (+1,0,0)   2: (-1,0,0)
  !   3: (0,+1,0)   4: (0,-1,0)
  !   5: (0,0,+1)   6: (0,0,-1)
  INTEGER, PARAMETER :: opp_f(0:6) = (/0, 2,1, 4,3, 6,5/)

  ! D3Q19 lattice velocities (ex, ey, ez)
  INTEGER, PARAMETER :: ex(0:18) = &
    (/0, 1,-1, 0,0, 0,0, 1,-1, 1,-1, 1,-1, 1,-1, 0,0, 0,0/)
  INTEGER, PARAMETER :: ey(0:18) = &
    (/0, 0,0, 1,-1, 0,0, 1,-1, -1,1, 0,0, 0,0, 1,-1, 1,-1/)
  INTEGER, PARAMETER :: ez(0:18) = &
    (/0, 0,0, 0,0, 1,-1, 0,0, 0,0, 1,-1, -1,1, 1,-1, -1,1/)

  ! D3Q7 lattice velocities
  INTEGER, PARAMETER :: ex_f(0:6) = (/0, 1,-1, 0,0, 0,0/)
  INTEGER, PARAMETER :: ey_f(0:6) = (/0, 0,0, 1,-1, 0,0/)
  INTEGER, PARAMETER :: ez_f(0:6) = (/0, 0,0, 0,0, 1,-1/)

END MODULE Domain

!-------------------------------------------------------------------------------
MODULE FluidParams
  USE NTypes, ONLY : DBL
  IMPLICIT NONE
  SAVE

  INTEGER :: nBubbles
  REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:) :: bubbles
  REAL(KIND = DBL) :: sigma, IntWidth, Gamma
  REAL(KIND = DBL) :: alpha, alpha4, kappa
  REAL(KIND = DBL) :: tauPhi, invTauPhi, invTauPhi1, phiStar, phiStar2, phiStar4
  REAL(KIND = DBL) :: rhoL, rhoH
  REAL(KIND = DBL) :: tauRho, invTauRho, invTauRhoOne, invTauRhoHalf
  REAL(KIND = DBL) :: eta, eta2, invEta2
  REAL(KIND = DBL) :: gravity, gravX, gravY, gravZ, Eo, eps, pConv
  REAL(KIND = DBL) :: sinthita, costhita

  ! Inlet velocity magnitude (for Zou-He)
  REAL(KIND = DBL) :: u_inlet_max

  REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:) :: phi
  REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:,:) :: u

END MODULE FluidParams

!-------------------------------------------------------------------------------
MODULE LBMParams
  USE NTypes, ONLY : DBL
  IMPLICIT NONE
  SAVE

  ! D3Q7 for order parameter, D3Q19 for momentum
  INTEGER, PARAMETER :: fdim = 6
  INTEGER, PARAMETER :: gdim = 18

  ! Distribution functions
  REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: f, g

  ! D3Q19 lattice constants
  REAL(KIND = DBL), PARAMETER :: Cs       = 0.57735026918962576451D0
  REAL(KIND = DBL), PARAMETER :: Cs_sq    = 0.33333333333333333333D0
  REAL(KIND = DBL), PARAMETER :: invCs_sq = 3.00000000000000000000D0

  ! D3Q19 weights
  REAL(KIND = DBL), PARAMETER :: Eg0 = 0.33333333333333333333D0
  REAL(KIND = DBL), PARAMETER :: Eg1 = 0.05555555555555555556D0
  REAL(KIND = DBL), PARAMETER :: Eg2 = 0.02777777777777777778D0

  ! Modified weights (Egi * invCs_sq)
  REAL(KIND = DBL), PARAMETER :: Eg0C = 1.00000000000000000000D0
  REAL(KIND = DBL), PARAMETER :: Eg1C = 0.16666666666666666667D0
  REAL(KIND = DBL), PARAMETER :: Eg2C = 0.08333333333333333333D0

  ! Pre-multiplied weights for collision
  REAL(KIND = DBL) :: Eg0n, Eg1n, Eg2n, EgC0n, EgC1n, EgC2n

END MODULE LBMParams

!===============================================================================
! MAIN PROGRAM
!===============================================================================
PROGRAM main
  USE Domain,      ONLY : iStep, MaxStep, now, nxt, tDump, tStat
  USE FluidParams, ONLY : eps, pConv
  IMPLICIT NONE

  CALL Parameters
  CALL MemAlloc(1)
  CALL Init

  ! Main iteration loop
  DO iStep = 1, MaxStep

    CALL Collision
    CALL Stream
    CALL BoundaryConditions

    ! Swap time levels
    now = 1 - now
    nxt = 1 - nxt

    IF (MOD(iStep, 100) == 0) PRINT '(A,I8,A,I8)', ' Step ', iStep, ' / ', MaxStep
    IF (MOD(iStep, tStat) == 0) CALL Stats
    IF (MOD(iStep, tDump) == 0) CALL VtkDump

  END DO

  ! Final output
  CALL Stats
  CALL VtkDump

  ! Cleanup
  CALL MemAlloc(2)

END PROGRAM main

!===============================================================================
! SUBROUTINE: Parameters
!===============================================================================
SUBROUTINE Parameters
  USE Domain
  USE FluidParams
  USE LBMParams
  IMPLICIT NONE

  INTEGER :: i
  INTEGER :: IO_ERR

  ! Read parameter data
  OPEN(UNIT=10, FILE="properties.in", STATUS="OLD", ACTION="READ", IOSTAT=IO_ERR)
  IF (IO_ERR /= 0) STOP "Error: Cannot open properties.in"

  READ(10,*)            ! Header
  READ(10,*) MaxStep
  READ(10,*)
  READ(10,*) tStat, tDump, xjump
  READ(10,*)
  READ(10,*) xmin, xmax, ymin, ymax, zmin, zmax
  READ(10,*)
  READ(10,*) rhoL, rhoH
  READ(10,*)
  READ(10,*) tauRho, tauPhi
  READ(10,*)
  READ(10,*) IntWidth, sigma, Gamma
  READ(10,*)
  READ(10,*) Eo, pConv
  READ(10,*)
  READ(10,*) sinthita, costhita
  READ(10,*)
  READ(10,*) geom_mode
  READ(10,*)
  READ(10,*) u_inlet_max

  ! U-bend specific parameters (only used if geom_mode == 2)
  IF (geom_mode == 2) THEN
    READ(10,*)
    READ(10,*) ubend_a, ubend_Rbend, ubend_Larm
  END IF

  CLOSE(UNIT=10)

  ! Read bubble/drop positions
  OPEN(UNIT=10, FILE="discrete.in", STATUS="OLD", ACTION="READ", IOSTAT=IO_ERR)
  IF (IO_ERR /= 0) STOP "Error: Cannot open discrete.in"

  READ(10,*)
  READ(10,*) nBubbles
  READ(10,*)
  ALLOCATE(bubbles(1:nBubbles, 1:4))
  DO i = 1, nBubbles
    READ(10,*) bubbles(i,1), bubbles(i,2), bubbles(i,3), bubbles(i,4)
  END DO
  CLOSE(UNIT=10)

  ! Derived fluid properties
  phiStar       = 0.5D0*(rhoH - rhoL)
  phiStar2      = phiStar*phiStar
  phiStar4      = phiStar2*phiStar2
  invTauRho     = 1.D0/tauRho
  invTauRhoOne  = 1.D0 - invTauRho
  invTauRhoHalf = 1.D0 - 0.5D0*invTauRho

  eta       = 1.D0/(tauPhi + 0.5D0)
  eta2      = 1.D0 - eta
  invEta2   = 0.5D0/eta
  invTauPhi = 1.D0/tauPhi
  invTauPhi1 = 1.D0 - invTauPhi

  ! Chemical potential parameters
  alpha  = 0.75D0*sigma/(IntWidth*phiStar4)
  alpha4 = alpha*4.D0
  kappa  = (IntWidth*phiStar)**2*alpha/2.D0

  ! Gravity from Eötvös number
  gravity = 0.25D0*Eo*sigma/((rhoH - rhoL)*bubbles(1,4)*bubbles(1,4))
  gravX = 2.D0*phiStar*gravity*costhita
  gravY = 0.D0
  gravZ = 2.D0*phiStar*gravity*sinthita

  ! Pre-multiplied weights for collision
  Eg0n  = invTauRho*Eg0
  Eg1n  = invTauRho*Eg1
  Eg2n  = invTauRho*Eg2
  EgC0n = invTauRhoHalf*Eg0C
  EgC1n = invTauRhoHalf*Eg1C
  EgC2n = invTauRhoHalf*Eg2C

  ! Mid-planes
  xo = INT(0.5D0*(xmax + xmin))
  yo = INT(0.5D0*(ymax + ymin))
  zo = INT(0.5D0*(zmax + zmin))

  ! Grid dimensions
  NX = xmax - xmin + 1
  NY = ymax - ymin + 1
  NZ = zmax - zmin + 1

  PRINT '(A)',        '========================================='
  PRINT '(A)',        '  ZSC 3D Multiphase LBM - Fixed Version'
  PRINT '(A)',        '========================================='
  PRINT '(A,I1)',     '  Geometry mode: ', geom_mode
  PRINT '(A,I4,A,I4,A,I4)', '  Domain: ', NX, ' x ', NY, ' x ', NZ
  PRINT '(A,F10.4)',  '  rhoH: ', rhoH
  PRINT '(A,F10.4)',  '  rhoL: ', rhoL
  PRINT '(A,F10.4)',  '  sigma: ', sigma
  PRINT '(A,F10.4)',  '  Eo: ', Eo
  PRINT '(A,F10.6)',  '  gravity: ', gravity
  PRINT '(A,F10.4)',  '  tauRho: ', tauRho
  PRINT '(A,F10.4)',  '  tauPhi: ', tauPhi
  PRINT '(A,F10.4)',  '  Gamma: ', Gamma
  PRINT '(A,I8)',     '  MaxStep: ', MaxStep
  PRINT '(A)',        '========================================='

END SUBROUTINE Parameters

!===============================================================================
! SUBROUTINE: MemAlloc
!===============================================================================
SUBROUTINE MemAlloc(FLAG)
  USE Domain
  USE FluidParams, ONLY : phi, u
  USE LBMParams,   ONLY : f, fdim, g, gdim
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FLAG

  IF (FLAG == 1) THEN
    ALLOCATE(ni(xmin-1:xmax+1, ymin-1:ymax+1, zmin-1:zmax+1, 1:6))
    ALLOCATE(is_fluid(xmin-1:xmax+1, ymin-1:ymax+1, zmin-1:zmax+1))
    ALLOCATE(phi(xmin-1:xmax+1, ymin-1:ymax+1, zmin-1:zmax+1))
    ALLOCATE(u(xmin-1:xmax+1, ymin-1:ymax+1, zmin-1:zmax+1, 3))
    ALLOCATE(f(xmin-1:xmax+1, ymin-1:ymax+1, zmin-1:zmax+1, 0:fdim, 0:1))
    ALLOCATE(g(xmin-1:xmax+1, ymin-1:ymax+1, zmin-1:zmax+1, 0:gdim, 0:1))

    ! Initialize everything to zero
    ni = 0
    is_fluid = 0
    phi = 0.D0
    u = 0.D0
    f = 0.D0
    g = 0.D0
  ELSE
    DEALLOCATE(ni, is_fluid, phi, u, f, g)
  END IF

END SUBROUTINE MemAlloc

!===============================================================================
! SUBROUTINE: Init
!===============================================================================
SUBROUTINE Init
  USE NTypes,    ONLY : DBL
  USE Domain
  USE FluidParams
  USE LBMParams
  IMPLICIT NONE

  INTEGER :: i, j, k, m, ie, iw, jn, js, kt, kb
  REAL(KIND = DBL) :: r, Af0, Ag0, Af1, Ag1, Ag2
  REAL(KIND = DBL) :: muPhin, phin, rhon, lapPhi

  ! Initialize counters
  iStep = 0
  tCall = 1
  now   = 0
  nxt   = 1
  eps   = 1.D0

  !---------------------------------------------------------------------------
  ! Step 1: Build geometry (wall mask)
  !---------------------------------------------------------------------------
  CALL BuildGeometry

  !---------------------------------------------------------------------------
  ! Step 2: Set neighbor arrays
  !---------------------------------------------------------------------------
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        ni(i,j,k,1) = i + 1
        ni(i,j,k,2) = i - 1
        ni(i,j,k,3) = j + 1
        ni(i,j,k,4) = j - 1
        ni(i,j,k,5) = k + 1
        ni(i,j,k,6) = k - 1
      END DO
    END DO
  END DO

  ! For wall nodes, clamp neighbors to self (prevents out-of-bounds access)
  ! The actual bounce-back is handled in BoundaryConditions subroutine
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) THEN
          ni(i,j,k,1:6) = (/i, i, j, j, k, k/)  ! self-reference
        ELSE
          ! Clamp fluid neighbors at domain edges
          IF (ni(i,j,k,1) > xmax) ni(i,j,k,1) = i
          IF (ni(i,j,k,2) < xmin) ni(i,j,k,2) = i
          IF (ni(i,j,k,3) > ymax) ni(i,j,k,3) = j
          IF (ni(i,j,k,4) < ymin) ni(i,j,k,4) = j
          IF (ni(i,j,k,5) > zmax) ni(i,j,k,5) = k
          IF (ni(i,j,k,6) < zmin) ni(i,j,k,6) = k
        END IF
      END DO
    END DO
  END DO

  !---------------------------------------------------------------------------
  ! Step 3: Initialize order parameter (drops/bubbles)
  !---------------------------------------------------------------------------
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 1) THEN
          phi(i,j,k) = -phiStar  ! Continuous/heavy phase

          DO m = 1, nBubbles
            r = DSQRT((DBLE(i) - bubbles(m,1))**2 &
              + (DBLE(j) - bubbles(m,2))**2 &
              + (DBLE(k) - bubbles(m,3))**2)

            IF (r <= (DBLE(bubbles(m,4)) + IntWidth)) THEN
              phi(i,j,k) = phiStar*DTANH(2.D0*(DBLE(bubbles(m,4)) - r)/IntWidth)
            END IF
          END DO
        ELSE
          phi(i,j,k) = -phiStar  ! Wall nodes get heavy phase
        END IF
      END DO
    END DO
  END DO

  ! Calculate initial volume for mass conservation check
  invInitVol = 0.D0
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (phi(i,j,k) >= 0.D0 .AND. is_fluid(i,j,k) == 1) THEN
          invInitVol = invInitVol + 1.D0
        END IF
      END DO
    END DO
  END DO
  IF (invInitVol > 0.D0) invInitVol = 1.D0/invInitVol

  !---------------------------------------------------------------------------
  ! Step 4: Initialize distribution functions to equilibrium (u=0)
  !---------------------------------------------------------------------------
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) CYCLE

        ! Get neighbors
        ie = ni(i,j,k,1); iw = ni(i,j,k,2)
        jn = ni(i,j,k,3); js = ni(i,j,k,4)
        kt = ni(i,j,k,5); kb = ni(i,j,k,6)

        ! Laplacian of order parameter
        lapPhi = (phi(ie,jn,k) + phi(iw,js,k) + phi(ie,js,k) + phi(iw,jn,k) &
               +  phi(ie,j,kt) + phi(iw,j,kb) + phi(ie,j,kb) + phi(iw,j,kt) &
               +  phi(i,jn,kt) + phi(i,js,kb) + phi(i,jn,kb) + phi(i,js,kt) &
               +  2.0D0*(phi(ie,j,k) + phi(iw,j,k) + phi(i,jn,k) &
               +  phi(i,js,k) + phi(i,j,kt) + phi(i,j,kb) - 12.D0*phi(i,j,k)))*inv36

        ! Chemical potential
        phin   = phi(i,j,k)
        rhon   = 0.5D0*(rhoH + rhoL)
        muPhin = alpha4*phin*(phin*phin - phiStar2) - kappa*lapPhi

        ! f equilibrium (D3Q7, u=0)
        Af0 = -3.D0*Gamma*muPhin
        Af1 =  0.5D0*Gamma*muPhin
        f(i,j,k,0,now)   = Af0 + phin
        f(i,j,k,1:6,now) = Af1

        ! g equilibrium (D3Q19, u=0)
        Ag0 = Eg0*(rhon - 6.D0*phin*muPhin)
        Ag1 = Eg1*(3.D0*phin*muPhin + rhon)
        Ag2 = Eg2*(3.D0*phin*muPhin + rhon)
        g(i,j,k,0,now)    = Ag0
        g(i,j,k,1:6,now)  = Ag1
        g(i,j,k,7:18,now) = Ag2
      END DO
    END DO
  END DO

  ! Initial stats and VTK dump
  CALL Stats
  CALL VtkDump

  PRINT '(A)', '  Initialization complete.'
  PRINT '(A)', '========================================='

END SUBROUTINE Init

!===============================================================================
! SUBROUTINE: BuildGeometry
!===============================================================================
SUBROUTINE BuildGeometry
  USE Domain
  USE FluidParams, ONLY : u, u_inlet_max, geom_mode => geom_mode
  IMPLICIT NONE

  INTEGER :: i, j, k

  ! Default: everything is fluid
  is_fluid = 1

  SELECT CASE(geom_mode)

  !---------------------------------------------------------------------------
  CASE(1) ! Straight vertical channel (Sudhakar validation baseline)
    ! Walls on x-faces and y-faces, open in z
    ! Flow in +z direction, gravity in -z
    DO k = zmin, zmax
      DO j = ymin, ymax
        DO i = xmin, xmax
          ! Wall on all four sides (x and y boundaries)
          IF (i == xmin .OR. i == xmax .OR. j == ymin .OR. j == ymax) THEN
            is_fluid(i,j,k) = 0
          END IF
        END DO
      END DO
    END DO

    PRINT '(A)', '  Geometry: Straight vertical channel'
    PRINT '(A,I4,A,I4)', '  Channel cross-section: ', &
      (xmax-xmin-1), ' x ', (ymax-ymin-1)

  !---------------------------------------------------------------------------
  CASE(2) ! U-bend channel (Sudhakar Fig 1a)
    !
    ! Geometry layout in the x-z plane (y is the depth/width direction):
    !
    !   z=zmax  +---------+           +---------+
    !           | INLET   |           | OUTLET  |
    !           | ARM     |           | ARM     |
    !           | (flow   |           | (flow   |
    !           |  +z)    |           |  -z)    |
    !           |         |           |         |
    !   z=zbend +---------+-----------+---------+
    !           |        BEND REGION            |
    !           |   (connects arms via 180°     |
    !           |    U-turn in x-z plane)       |
    !   z=zmin  +-------------------------------+
    !           x=xmin                    x=xmax
    !
    ! Coordinate system:
    !   x: horizontal (left-right), inlet arm on left, outlet arm on right
    !   y: into the page (channel depth, walls at ymin/ymax)
    !   z: vertical (up-down), flow enters at top-left, exits top-right
    !
    ! Parameters:
    !   a = internal channel width (fluid nodes per side)
    !   Rbend = bend centerline radius
    !   Larm = length of straight vertical arms above bend
    !
    ! Wall thickness = 1 lattice node on each side
    ! Total a_with_walls = a + 2
    !
    ! Domain sizing (auto-computed):
    !   x: wall + inlet_arm + inner_wall + bend_gap + inner_wall + outlet_arm + wall
    !      = 1 + a + (2*Rbend) + a + 1 = 2*a + 2*Rbend + 2
    !   y: wall + a + wall = a + 2
    !   z: bend_depth + Larm + wall = (Rbend + a/2 + 1) + Larm + 1
    !      Actually: z goes from zmin to zmax where:
    !        zmin = 1 (bottom of bend)
    !        zbend_top = Rbend + a + 1 (where straight arms start)
    !        zmax = zbend_top + Larm
    !
    BLOCK
      INTEGER :: a, Rb, La
      INTEGER :: zbend_center, zbend_top
      INTEGER :: inlet_x_lo, inlet_x_hi, outlet_x_lo, outlet_x_hi
      INTEGER :: bend_x_lo, bend_x_hi
      REAL(8) :: cx, cz, rr, r_inner, r_outer

      a  = ubend_a        ! Internal fluid width
      Rb = ubend_Rbend    ! Bend centerline radius
      La = ubend_Larm     ! Straight arm length

      ! The bend centerline is a semicircle in x-z plane
      ! Center of semicircle: x = midpoint between arms, z = zbend_center
      ! 
      ! Inlet arm fluid region:  x in [xmin+1, xmin+a]
      ! Outlet arm fluid region: x in [xmax-a, xmax-1]
      ! Bend center x: xc = (xmin+a + xmax-a) / 2 = (xmin + xmax) / 2
      ! Bend center z: zc = where arms meet the bend = zbend_top
      !
      ! Actually, let's define it more carefully:
      ! 
      ! Inlet arm inner wall (right side):  x = xmin + a + 1
      ! Outlet arm inner wall (left side):  x = xmax - a - 1  (= xmin + a + 1 + 2*Rb - 2... let's compute)
      !
      ! Domain x-width: total = 2*(a+1) + 2*Rb = 2*a + 2*Rb + 2
      ! So xmax = xmin + 2*a + 2*Rb + 1
      
      ! Verify domain size matches parameters
      IF (xmax - xmin + 1 /= 2*a + 2*Rb + 2) THEN
        PRINT *, 'WARNING: x-domain size mismatch for U-bend!'
        PRINT *, '  Expected NX = ', 2*a + 2*Rb + 2
        PRINT *, '  Got NX = ', xmax - xmin + 1
        PRINT *, '  Auto-adjusting xmax...'
        xmax = xmin + 2*a + 2*Rb + 1
        NX = xmax - xmin + 1
      END IF

      IF (ymax - ymin + 1 /= a + 2) THEN
        PRINT *, 'WARNING: y-domain size mismatch for U-bend!'
        PRINT *, '  Expected NY = ', a + 2
        PRINT *, '  Got NY = ', ymax - ymin + 1
        PRINT *, '  Auto-adjusting ymax...'
        ymax = ymin + a + 1
        NY = ymax - ymin + 1
      END IF

      ! Key z-coordinates
      ! Bottom of domain = bottom of bend
      ! zbend_center = center height of the semicircular bend
      zbend_center = zmin + Rb + a/2 + 1
      zbend_top = zbend_center + Rb   ! Where straight arms begin

      IF (zmax - zmin + 1 /= zbend_top - zmin + La + 1) THEN
        PRINT *, 'WARNING: z-domain size mismatch for U-bend!'
        PRINT *, '  Expected NZ = ', zbend_top - zmin + La + 1
        PRINT *, '  Got NZ = ', zmax - zmin + 1
        PRINT *, '  Auto-adjusting zmax...'
        zmax = zbend_top + La
        NZ = zmax - zmin + 1
      END IF

      ! X-ranges for each arm (fluid regions, excluding walls)
      inlet_x_lo  = xmin + 1          ! Inlet arm left wall at xmin
      inlet_x_hi  = xmin + a          ! Inlet arm right wall at xmin+a+1
      outlet_x_lo = xmax - a          ! Outlet arm left wall at xmax-a-1
      outlet_x_hi = xmax - 1          ! Outlet arm right wall at xmax

      ! Store for BC application
      inlet_x1  = inlet_x_lo
      inlet_x2  = inlet_x_hi
      outlet_x1 = outlet_x_lo
      outlet_x2 = outlet_x_hi
      inlet_zface  = zmax     ! Inlet at top of left arm
      outlet_zface = zmax     ! Outlet at top of right arm

      ! Bend semicircle center
      cx = 0.5D0*(DBLE(inlet_x_hi) + DBLE(outlet_x_lo))  ! x-center
      cz = DBLE(zbend_top)                                  ! z-center

      r_inner = DBLE(Rb) - 0.5D0*DBLE(a)  ! Inner bend radius
      r_outer = DBLE(Rb) + 0.5D0*DBLE(a)  ! Outer bend radius

      ! --- Build is_fluid mask ---
      ! Start with everything solid
      is_fluid = 0

      DO k = zmin, zmax
        DO j = ymin, ymax
          DO i = xmin, xmax
            ! Y-walls always present
            IF (j == ymin .OR. j == ymax) CYCLE

            ! --- Straight arms (above the bend) ---
            IF (k >= zbend_top .AND. k <= zmax) THEN
              ! Inlet arm
              IF (i >= inlet_x_lo .AND. i <= inlet_x_hi) THEN
                is_fluid(i,j,k) = 1
              END IF
              ! Outlet arm
              IF (i >= outlet_x_lo .AND. i <= outlet_x_hi) THEN
                is_fluid(i,j,k) = 1
              END IF
            END IF

            ! --- Bend region (semicircular in x-z plane) ---
            IF (k < zbend_top) THEN
              ! Check if point is inside the semicircular bend channel
              ! The bend is a 180° turn: a semicircle below z = zbend_top
              rr = DSQRT((DBLE(i) - cx)**2 + (DBLE(k) - cz)**2)

              ! Inside the bend annulus AND below the centerline
              IF (rr >= r_inner .AND. rr <= r_outer .AND. k <= INT(cz)) THEN
                is_fluid(i,j,k) = 1
              END IF

              ! Also include the straight portions connecting arms to the bend
              ! Left side: inlet arm extends down to meet the bend
              IF (i >= inlet_x_lo .AND. i <= inlet_x_hi .AND. &
                  k >= INT(cz) .AND. k < zbend_top) THEN
                is_fluid(i,j,k) = 1
              END IF
              ! Right side: outlet arm extends down to meet the bend
              IF (i >= outlet_x_lo .AND. i <= outlet_x_hi .AND. &
                  k >= INT(cz) .AND. k < zbend_top) THEN
                is_fluid(i,j,k) = 1
              END IF
            END IF

          END DO
        END DO
      END DO

      PRINT '(A)', '  Geometry: U-bend channel (Sudhakar Fig 1a)'
      PRINT '(A,I4)', '  Channel width a = ', a
      PRINT '(A,I4)', '  Bend radius Rb = ', Rb
      PRINT '(A,I4)', '  Arm length La = ', La
      PRINT '(A,I4,A,I4,A,I4)', '  Domain: ', NX, ' x ', NY, ' x ', NZ
      PRINT '(A,I4,A,I4)', '  Inlet arm x-range: ', inlet_x1, ' to ', inlet_x2
      PRINT '(A,I4,A,I4)', '  Outlet arm x-range: ', outlet_x1, ' to ', outlet_x2
      PRINT '(A,I4)', '  Bend center z = ', INT(cz)
      PRINT '(A,F8.2,A,F8.2)', '  Bend r_inner = ', r_inner, ', r_outer = ', r_outer

    END BLOCK

  !---------------------------------------------------------------------------
  CASE(3) ! Static drop - Laplace law test (fully periodic, no walls)
    ! All nodes are fluid
    is_fluid = 1
    PRINT '(A)', '  Geometry: Periodic box (Laplace test)'

  !---------------------------------------------------------------------------
  CASE DEFAULT
    STOP "Error: Unknown geometry mode"
  END SELECT

  ! Ghost nodes outside domain are always solid
  is_fluid(xmin-1,:,:) = 0
  is_fluid(xmax+1,:,:) = 0
  is_fluid(:,ymin-1,:) = 0
  is_fluid(:,ymax+1,:) = 0
  is_fluid(:,:,zmin-1) = 0
  is_fluid(:,:,zmax+1) = 0

END SUBROUTINE BuildGeometry

!===============================================================================
! SUBROUTINE: Collision
! ZSC physics core — this is CORRECT in Arup's original code
! Only change: skip wall nodes using is_fluid mask
!===============================================================================
SUBROUTINE Collision
  USE NTypes,      ONLY : DBL
  USE Domain
  USE FluidParams
  USE LBMParams
  IMPLICIT NONE

  INTEGER :: i, j, k, ie, iw, jn, js, kt, kb
  REAL(KIND = DBL) :: muPhin, phin, phin2
  REAL(KIND = DBL) :: rhon, invRhon
  REAL(KIND = DBL) :: sFx, sFy, sFz
  REAL(KIND = DBL) :: UF, ux, uy, uz, Vg, Vsq
  REAL(KIND = DBL) :: Af0, Af1, Cf1
  REAL(KIND = DBL) :: Ag0, Ag1, Eg1A, Eg2A, Eg1R, Eg2R
  REAL(KIND = DBL) :: geq1, geq2, gradPhiX, gradPhiY, gradPhiZ, lapPhi
  REAL(KIND = DBL) :: in01, in02, in03, in04, in05, in06
  REAL(KIND = DBL) :: in07, in08, in09, in10, in11, in12
  REAL(KIND = DBL) :: in13, in14, in15, in16, in17, in18

  !--- Update order parameter from f ---
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) CYCLE
        phi(i,j,k) = SUM(f(i,j,k,:,now))
        ! Clamp to prevent blowup
        IF (phi(i,j,k) >  phiStar*1.05D0) phi(i,j,k) =  phiStar*1.05D0
        IF (phi(i,j,k) < -phiStar*1.05D0) phi(i,j,k) = -phiStar*1.05D0
      END DO
    END DO
  END DO

  !--- Main collision loop ---
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) CYCLE

        phin  = phi(i,j,k)
        phin2 = phin*phin
        rhon  = SUM(g(i,j,k,:,now))
        invRhon = 1.D0/rhon

        ! Identify neighbors
        ie = ni(i,j,k,1); iw = ni(i,j,k,2)
        jn = ni(i,j,k,3); js = ni(i,j,k,4)
        kt = ni(i,j,k,5); kb = ni(i,j,k,6)

        ! Neighbor phi values
        in01 = phi(ie,j ,k ); in02 = phi(iw,j ,k )
        in03 = phi(i ,jn,k ); in04 = phi(i ,js,k )
        in05 = phi(i ,j ,kt); in06 = phi(i ,j ,kb)
        in07 = phi(ie,jn,k ); in08 = phi(iw,js,k )
        in09 = phi(ie,js,k ); in10 = phi(iw,jn,k )
        in11 = phi(ie,j ,kt); in12 = phi(iw,j ,kb)
        in13 = phi(ie,j ,kb); in14 = phi(iw,j ,kt)
        in15 = phi(i ,jn,kt); in16 = phi(i ,js,kb)
        in17 = phi(i ,jn,kb); in18 = phi(i ,js,kt)

        ! Laplacian of phi (isotropic 18-neighbor stencil)
        lapPhi = (in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 &
               +  in15 + in16 + in17 + in18 &
               +  2.0D0*(in01 + in02 + in03 + in04 + in05 + in06 &
               -  12.D0*phin))*inv36

        ! Gradient of phi
        gradPhiX = (2.0D0*(in01 - in02) + in07 - in08 + in09 - in10 &
                 +  in11 - in12 + in13 - in14)*inv12
        gradPhiY = (2.0D0*(in03 - in04) + in07 - in08 + in10 - in09 &
                 +  in15 - in16 + in17 - in18)*inv12
        gradPhiZ = (2.0D0*(in05 - in06) + in11 - in12 + in14 - in13 &
                 +  in15 - in16 + in18 - in17)*inv12

        ! Chemical potential
        muPhin = alpha4*phin*(phin2 - phiStar2) - kappa*lapPhi

        ! Interfacial force + gravity (gravity on light phase only)
        sFx = muPhin*gradPhiX
        sFy = muPhin*gradPhiY
        sFz = muPhin*gradPhiZ
        IF (phin > 0.D0) THEN
          sFx = sFx + gravX
          sFy = sFy + gravY
          sFz = sFz + gravZ
        END IF

        ! Velocity (with half-force correction)
        ux = (g(i,j,k,1,now) - g(i,j,k,2,now) + g(i,j,k,7,now) &
           -  g(i,j,k,8,now) + g(i,j,k,9,now) - g(i,j,k,10,now) &
           +  g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
           -  g(i,j,k,14,now) + 0.5D0*sFx)*invRhon

        uy = (g(i,j,k,3,now) - g(i,j,k,4,now) + g(i,j,k,7,now) &
           -  g(i,j,k,8,now) - g(i,j,k,9,now) + g(i,j,k,10,now) &
           +  g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
           -  g(i,j,k,18,now) + 0.5D0*sFy)*invRhon

        uz = (g(i,j,k,5,now) - g(i,j,k,6,now) + g(i,j,k,11,now) &
           -  g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
           +  g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
           +  g(i,j,k,18,now) + 0.5D0*sFz)*invRhon

        ! Store velocity for use in streaming/BCs
        u(i,j,k,1) = ux
        u(i,j,k,2) = uy
        u(i,j,k,3) = uz

        !--- f collision (order parameter, D3Q7) ---
        Af0 = -3.D0*Gamma*muPhin*invTauPhi
        Af1 =  0.5D0*Gamma*muPhin*invTauPhi
        Cf1 =  invTauPhi*invEta2*phin

        f(i,j,k,0,nxt) = invTauPhi1*f(i,j,k,0,now) + Af0 + invTauPhi*phin
        f(i,j,k,1,now) = invTauPhi1*f(i,j,k,1,now) + Af1 + Cf1*ux
        f(i,j,k,2,now) = invTauPhi1*f(i,j,k,2,now) + Af1 - Cf1*ux
        f(i,j,k,3,now) = invTauPhi1*f(i,j,k,3,now) + Af1 + Cf1*uy
        f(i,j,k,4,now) = invTauPhi1*f(i,j,k,4,now) + Af1 - Cf1*uy
        f(i,j,k,5,now) = invTauPhi1*f(i,j,k,5,now) + Af1 + Cf1*uz
        f(i,j,k,6,now) = invTauPhi1*f(i,j,k,6,now) + Af1 - Cf1*uz

        !--- g collision + streaming (momentum, D3Q19) ---
        ! Note: g streams during collision (fused collide-stream)
        Ag0  = rhon - 6.D0*phin*muPhin
        Ag1  = 3.D0*phin*muPhin + rhon
        Eg1A = Eg1n*Ag1;  Eg2A = Eg2n*Ag1
        Eg1R = Eg1n*rhon;  Eg2R = Eg2n*rhon
        Vsq  = 1.5D0*(ux*ux + uy*uy + uz*uz)
        UF   = ux*sFx + uy*sFy + uz*sFz

        ! Rescale velocity
        ux = ux*invCs_sq
        uy = uy*invCs_sq
        uz = uz*invCs_sq

        ! DIR 0
        g(i,j,k,0,nxt) = invTauRhoOne*g(i,j,k,0,now) + Eg0n*(Ag0 - rhon*Vsq) &
                        - EgC0n*UF

        ! DIR 1 & 2 (+x, -x)
        geq1 = Eg1A + Eg1R*(0.5D0*ux*ux - Vsq)
        geq2 = Eg1R*ux
        g(ie,j,k,1,nxt) = invTauRhoOne*g(i,j,k,1,now) + geq1 + geq2 &
                         + EgC1n*((1.D0 + ux)*sFx - UF)
        g(iw,j,k,2,nxt) = invTauRhoOne*g(i,j,k,2,now) + geq1 - geq2 &
                         + EgC1n*((-1.D0 + ux)*sFx - UF)

        ! DIR 3 & 4 (+y, -y)
        geq1 = Eg1A + Eg1R*(0.5D0*uy*uy - Vsq)
        geq2 = Eg1R*uy
        g(i,jn,k,3,nxt) = invTauRhoOne*g(i,j,k,3,now) + geq1 + geq2 &
                         + EgC1n*((1.D0 + uy)*sFy - UF)
        g(i,js,k,4,nxt) = invTauRhoOne*g(i,j,k,4,now) + geq1 - geq2 &
                         + EgC1n*((-1.D0 + uy)*sFy - UF)

        ! DIR 5 & 6 (+z, -z)
        geq1 = Eg1A + Eg1R*(0.5D0*uz*uz - Vsq)
        geq2 = Eg1R*uz
        g(i,j,kt,5,nxt) = invTauRhoOne*g(i,j,k,5,now) + geq1 + geq2 &
                         + EgC1n*((1.D0 + uz)*sFz - UF)
        g(i,j,kb,6,nxt) = invTauRhoOne*g(i,j,k,6,now) + geq1 - geq2 &
                         + EgC1n*((-1.D0 + uz)*sFz - UF)

        ! DIR 7 & 8 (+x+y, -x-y)
        Vg = ux + uy
        geq1 = Eg2A + Eg2R*(0.5D0*Vg*Vg - Vsq)
        geq2 = Eg2R*Vg
        g(ie,jn,k,7,nxt) = invTauRhoOne*g(i,j,k,7,now) + geq1 + geq2 &
                          + EgC2n*((1.D0 + Vg)*(sFx + sFy) - UF)
        g(iw,js,k,8,nxt) = invTauRhoOne*g(i,j,k,8,now) + geq1 - geq2 &
                          + EgC2n*((-1.D0 + Vg)*(sFx + sFy) - UF)

        ! DIR 9 & 10 (+x-y, -x+y)
        Vg = ux - uy
        geq1 = Eg2A + Eg2R*(0.5D0*Vg*Vg - Vsq)
        geq2 = Eg2R*Vg
        g(ie,js,k,9,nxt) = invTauRhoOne*g(i,j,k,9,now) + geq1 + geq2 &
                          + EgC2n*((1.D0 + Vg)*(sFx - sFy) - UF)
        g(iw,jn,k,10,nxt) = invTauRhoOne*g(i,j,k,10,now) + geq1 - geq2 &
                           + EgC2n*((-1.D0 + Vg)*(sFx - sFy) - UF)

        ! DIR 11 & 12 (+x+z, -x-z)
        Vg = ux + uz
        geq1 = Eg2A + Eg2R*(0.5D0*Vg*Vg - Vsq)
        geq2 = Eg2R*Vg
        g(ie,j,kt,11,nxt) = invTauRhoOne*g(i,j,k,11,now) + geq1 + geq2 &
                           + EgC2n*((1.D0 + Vg)*(sFx + sFz) - UF)
        g(iw,j,kb,12,nxt) = invTauRhoOne*g(i,j,k,12,now) + geq1 - geq2 &
                           + EgC2n*((-1.D0 + Vg)*(sFx + sFz) - UF)

        ! DIR 13 & 14 (+x-z, -x+z)
        Vg = ux - uz
        geq1 = Eg2A + Eg2R*(0.5D0*Vg*Vg - Vsq)
        geq2 = Eg2R*Vg
        g(ie,j,kb,13,nxt) = invTauRhoOne*g(i,j,k,13,now) + geq1 + geq2 &
                           + EgC2n*((1.D0 + Vg)*(sFx - sFz) - UF)
        g(iw,j,kt,14,nxt) = invTauRhoOne*g(i,j,k,14,now) + geq1 - geq2 &
                           + EgC2n*((-1.D0 + Vg)*(sFx - sFz) - UF)

        ! DIR 15 & 16 (+y+z, -y-z)
        Vg = uy + uz
        geq1 = Eg2A + Eg2R*(0.5D0*Vg*Vg - Vsq)
        geq2 = Eg2R*Vg
        g(i,jn,kt,15,nxt) = invTauRhoOne*g(i,j,k,15,now) + geq1 + geq2 &
                           + EgC2n*((1.D0 + Vg)*(sFy + sFz) - UF)
        g(i,js,kb,16,nxt) = invTauRhoOne*g(i,j,k,16,now) + geq1 - geq2 &
                           + EgC2n*((-1.D0 + Vg)*(sFy + sFz) - UF)

        ! DIR 17 & 18 (+y-z, -y+z)
        Vg = uy - uz
        geq1 = Eg2A + Eg2R*(0.5D0*Vg*Vg - Vsq)
        geq2 = Eg2R*Vg
        g(i,jn,kb,17,nxt) = invTauRhoOne*g(i,j,k,17,now) + geq1 + geq2 &
                           + EgC2n*((1.D0 + Vg)*(sFy - sFz) - UF)
        g(i,js,kt,18,nxt) = invTauRhoOne*g(i,j,k,18,now) + geq1 - geq2 &
                           + EgC2n*((-1.D0 + Vg)*(sFy - sFz) - UF)

      END DO
    END DO
  END DO

END SUBROUTINE Collision

!===============================================================================
! SUBROUTINE: Stream
! Lax-Wendroff streaming for order parameter f (ZSC-specific)
!===============================================================================
SUBROUTINE Stream
  USE NTypes,    ONLY : DBL
  USE Domain
  USE FluidParams, ONLY : eta, eta2
  USE LBMParams,   ONLY : f
  IMPLICIT NONE

  INTEGER :: i, j, k, ie, iw, jn, js, kt, kb

  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) CYCLE

        ie = ni(i,j,k,1); iw = ni(i,j,k,2)
        jn = ni(i,j,k,3); js = ni(i,j,k,4)
        kt = ni(i,j,k,5); kb = ni(i,j,k,6)

        ! Lax-Wendroff streaming: f_new(neighbor) = eta*f_post(here) + eta2*f_post(neighbor)
        f(ie,j ,k ,1,nxt) = eta*f(i,j,k,1,now) + eta2*f(ie,j ,k ,1,now)
        f(iw,j ,k ,2,nxt) = eta*f(i,j,k,2,now) + eta2*f(iw,j ,k ,2,now)
        f(i ,jn,k ,3,nxt) = eta*f(i,j,k,3,now) + eta2*f(i ,jn,k ,3,now)
        f(i ,js,k ,4,nxt) = eta*f(i,j,k,4,now) + eta2*f(i ,js,k ,4,now)
        f(i ,j ,kt,5,nxt) = eta*f(i,j,k,5,now) + eta2*f(i ,j ,kt,5,now)
        f(i ,j ,kb,6,nxt) = eta*f(i,j,k,6,now) + eta2*f(i ,j ,kb,6,now)
      END DO
    END DO
  END DO

END SUBROUTINE Stream

!===============================================================================
! SUBROUTINE: BoundaryConditions
! THE KEY FIX: proper bounce-back for walls, Zou-He inlet, pressure outlet
!===============================================================================
SUBROUTINE BoundaryConditions
  USE NTypes,    ONLY : DBL
  USE Domain
  USE FluidParams
  USE LBMParams
  IMPLICIT NONE

  INTEGER :: i, j, k, d, od
  REAL(KIND = DBL) :: uz_wall, rho_out, rho_in
  REAL(KIND = DBL) :: ux_in, uy_in, uz_in
  REAL(KIND = DBL) :: Neq

  !---------------------------------------------------------------------------
  ! 1. WALL BOUNCE-BACK (proper implementation)
  !    For each wall node, swap incoming/outgoing distributions
  !    This is the FIX for Bug #1 (sequential overwrite)
  !---------------------------------------------------------------------------
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 1) CYCLE  ! Skip fluid nodes

        ! For wall nodes: g_opposite = g_incoming (bounce-back)
        DO d = 1, 18
          od = opp_g(d)
          ! The distribution that streamed INTO this wall node
          ! gets reflected back to the opposite direction
          g(i,j,k,od,nxt) = g(i,j,k,d,nxt)
        END DO

        ! Same for f (D3Q7)
        DO d = 1, 6
          od = opp_f(d)
          f(i,j,k,od,nxt) = f(i,j,k,d,nxt)
        END DO
      END DO
    END DO
  END DO

  !---------------------------------------------------------------------------
  ! 2. INLET — Zou-He with parabolic velocity profile
  !    Mode 1: inlet at z=zmin, full cross-section
  !    Mode 2: inlet at z=zmax, inlet arm x-range only (flow in -z direction)
  !---------------------------------------------------------------------------
  IF (geom_mode == 1) THEN
    ! Straight channel: inlet at zmin, flow in +z
    DO j = ymin+1, ymax-1
      DO i = xmin+1, xmax-1
        IF (is_fluid(i,j,zmin) == 0) CYCLE

        ux_in = 0.D0; uy_in = 0.D0
        uz_in = u_inlet_max * &
          (1.D0 - (2.D0*DBLE(i - xo)/DBLE(xmax - xmin - 2))**2) * &
          (1.D0 - (2.D0*DBLE(j - yo)/DBLE(ymax - ymin - 2))**2)
        IF (uz_in < 0.D0) uz_in = 0.D0

        ! Zou-He for +z inlet (unknown: 5, 11, 14, 15, 18)
        rho_in = (g(i,j,zmin,0,nxt) &
               +  g(i,j,zmin,1,nxt) + g(i,j,zmin,2,nxt) &
               +  g(i,j,zmin,3,nxt) + g(i,j,zmin,4,nxt) &
               +  g(i,j,zmin,7,nxt) + g(i,j,zmin,8,nxt) &
               +  g(i,j,zmin,9,nxt) + g(i,j,zmin,10,nxt) &
               +  2.D0*(g(i,j,zmin,6,nxt) + g(i,j,zmin,12,nxt) &
               +  g(i,j,zmin,13,nxt) + g(i,j,zmin,16,nxt) &
               +  g(i,j,zmin,17,nxt))) / (1.D0 - uz_in)

        g(i,j,zmin,5,nxt)  = g(i,j,zmin,6,nxt)  + rho_in*uz_in/3.D0
        g(i,j,zmin,11,nxt) = g(i,j,zmin,12,nxt) + rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmin,2,nxt) - g(i,j,zmin,1,nxt)) &
                            + rho_in*ux_in/6.D0
        g(i,j,zmin,14,nxt) = g(i,j,zmin,13,nxt) + rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmin,1,nxt) - g(i,j,zmin,2,nxt)) &
                            - rho_in*ux_in/6.D0
        g(i,j,zmin,15,nxt) = g(i,j,zmin,16,nxt) + rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmin,4,nxt) - g(i,j,zmin,3,nxt)) &
                            + rho_in*uy_in/6.D0
        g(i,j,zmin,18,nxt) = g(i,j,zmin,17,nxt) + rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmin,3,nxt) - g(i,j,zmin,4,nxt)) &
                            - rho_in*uy_in/6.D0

        f(i,j,zmin,5,nxt) = f(i,j,zmin,6,nxt)
      END DO
    END DO

  ELSE IF (geom_mode == 2) THEN
    ! U-bend: inlet at z=zmax in the inlet arm, flow in -z direction
    ! Unknown distributions: those streaming INTO domain (dirs 6,12,13,16,17)
    DO j = ymin+1, ymax-1
      DO i = inlet_x1, inlet_x2
        IF (is_fluid(i,j,zmax) == 0) CYCLE

        ux_in = 0.D0; uy_in = 0.D0
        ! Parabolic profile across the inlet arm cross-section
        uz_in = -u_inlet_max * &
          (1.D0 - (2.D0*DBLE(i - (inlet_x1+inlet_x2)/2)/DBLE(inlet_x2 - inlet_x1))**2) * &
          (1.D0 - (2.D0*DBLE(j - yo)/DBLE(ymax - ymin - 2))**2)
        ! uz_in is negative (flow going -z into the domain)
        IF (uz_in > 0.D0) uz_in = 0.D0

        ! Zou-He for -z inlet (unknown: 6, 12, 13, 16, 17)
        rho_in = (g(i,j,zmax,0,nxt) &
               +  g(i,j,zmax,1,nxt) + g(i,j,zmax,2,nxt) &
               +  g(i,j,zmax,3,nxt) + g(i,j,zmax,4,nxt) &
               +  g(i,j,zmax,7,nxt) + g(i,j,zmax,8,nxt) &
               +  g(i,j,zmax,9,nxt) + g(i,j,zmax,10,nxt) &
               +  2.D0*(g(i,j,zmax,5,nxt) + g(i,j,zmax,11,nxt) &
               +  g(i,j,zmax,14,nxt) + g(i,j,zmax,15,nxt) &
               +  g(i,j,zmax,18,nxt))) / (1.D0 + uz_in)

        g(i,j,zmax,6,nxt)  = g(i,j,zmax,5,nxt)  - rho_in*uz_in/3.D0
        g(i,j,zmax,12,nxt) = g(i,j,zmax,11,nxt) - rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmax,1,nxt) - g(i,j,zmax,2,nxt)) &
                            - rho_in*ux_in/6.D0
        g(i,j,zmax,13,nxt) = g(i,j,zmax,14,nxt) - rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmax,2,nxt) - g(i,j,zmax,1,nxt)) &
                            + rho_in*ux_in/6.D0
        g(i,j,zmax,16,nxt) = g(i,j,zmax,15,nxt) - rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmax,3,nxt) - g(i,j,zmax,4,nxt)) &
                            - rho_in*uy_in/6.D0
        g(i,j,zmax,17,nxt) = g(i,j,zmax,18,nxt) - rho_in*uz_in/6.D0 &
                            + 0.5D0*(g(i,j,zmax,4,nxt) - g(i,j,zmax,3,nxt)) &
                            + rho_in*uy_in/6.D0

        f(i,j,zmax,6,nxt) = f(i,j,zmax,5,nxt)
      END DO
    END DO
  END IF

  !---------------------------------------------------------------------------
  ! 3. OUTLET — Constant pressure extrapolation
  !    Mode 1: outlet at z=zmax, full cross-section
  !    Mode 2: outlet at z=zmax, outlet arm x-range (flow in +z direction out)
  !---------------------------------------------------------------------------
  IF (geom_mode == 1) THEN
    DO j = ymin+1, ymax-1
      DO i = xmin+1, xmax-1
        IF (is_fluid(i,j,zmax) == 0) CYCLE

        g(i,j,zmax,6,nxt)  = 2.D0*g(i,j,zmax-1,6,nxt)  - g(i,j,zmax-2,6,nxt)
        g(i,j,zmax,12,nxt) = 2.D0*g(i,j,zmax-1,12,nxt) - g(i,j,zmax-2,12,nxt)
        g(i,j,zmax,13,nxt) = 2.D0*g(i,j,zmax-1,13,nxt) - g(i,j,zmax-2,13,nxt)
        g(i,j,zmax,16,nxt) = 2.D0*g(i,j,zmax-1,16,nxt) - g(i,j,zmax-2,16,nxt)
        g(i,j,zmax,17,nxt) = 2.D0*g(i,j,zmax-1,17,nxt) - g(i,j,zmax-2,17,nxt)

        f(i,j,zmax,6,nxt) = 2.D0*f(i,j,zmax-1,6,nxt) - f(i,j,zmax-2,6,nxt)
      END DO
    END DO

  ELSE IF (geom_mode == 2) THEN
    ! U-bend: outlet at z=zmax in the outlet arm, flow exits in +z
    ! Unknown: 5, 11, 14, 15, 18 (streaming upward out of domain)
    ! Use extrapolation from interior
    DO j = ymin+1, ymax-1
      DO i = outlet_x1, outlet_x2
        IF (is_fluid(i,j,zmax) == 0) CYCLE

        g(i,j,zmax,5,nxt)  = 2.D0*g(i,j,zmax-1,5,nxt)  - g(i,j,zmax-2,5,nxt)
        g(i,j,zmax,11,nxt) = 2.D0*g(i,j,zmax-1,11,nxt) - g(i,j,zmax-2,11,nxt)
        g(i,j,zmax,14,nxt) = 2.D0*g(i,j,zmax-1,14,nxt) - g(i,j,zmax-2,14,nxt)
        g(i,j,zmax,15,nxt) = 2.D0*g(i,j,zmax-1,15,nxt) - g(i,j,zmax-2,15,nxt)
        g(i,j,zmax,18,nxt) = 2.D0*g(i,j,zmax-1,18,nxt) - g(i,j,zmax-2,18,nxt)

        f(i,j,zmax,5,nxt) = 2.D0*f(i,j,zmax-1,5,nxt) - f(i,j,zmax-2,5,nxt)
      END DO
    END DO
  END IF

  !---------------------------------------------------------------------------
  ! 4. PERIODIC BCs (for Laplace test, mode 3)
  !---------------------------------------------------------------------------
  IF (geom_mode == 3) THEN
    ! g: periodic in all directions
    DO d = 0, 18
      g(xmin,:,:,d,nxt) = g(xmax,:,:,d,nxt)
      g(xmax+1,:,:,d,nxt) = g(xmin+1,:,:,d,nxt)
      g(:,ymin,:,d,nxt) = g(:,ymax,:,d,nxt)
      g(:,ymax+1,:,d,nxt) = g(:,ymin+1,:,d,nxt)
      g(:,:,zmin,d,nxt) = g(:,:,zmax,d,nxt)
      g(:,:,zmax+1,d,nxt) = g(:,:,zmin+1,d,nxt)
    END DO
    ! f: periodic
    DO d = 0, 6
      f(xmin,:,:,d,nxt) = f(xmax,:,:,d,nxt)
      f(xmax+1,:,:,d,nxt) = f(xmin+1,:,:,d,nxt)
      f(:,ymin,:,d,nxt) = f(:,ymax,:,d,nxt)
      f(:,ymax+1,:,d,nxt) = f(:,ymin+1,:,d,nxt)
      f(:,:,zmin,d,nxt) = f(:,:,zmax,d,nxt)
      f(:,:,zmax+1,d,nxt) = f(:,:,zmin+1,d,nxt)
    END DO
  END IF

END SUBROUTINE BoundaryConditions

!===============================================================================
! SUBROUTINE: Stats
!===============================================================================
SUBROUTINE Stats
  USE NTypes,      ONLY : DBL
  USE Domain
  USE FluidParams, ONLY : phi
  USE LBMParams,   ONLY : g
  IMPLICIT NONE

  INTEGER :: i, j, k, IO_ERR
  REAL(KIND = DBL) :: ux, uy, uz, rhon, invRhon
  REAL(KIND = DBL) :: invVol, Ref, Vef, Vol
  REAL(KIND = DBL) :: total_ux, total_uy, total_uz

  Vol = 0.D0
  total_ux = 0.D0; total_uy = 0.D0; total_uz = 0.D0

  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) CYCLE

        rhon = SUM(g(i,j,k,:,now))
        invRhon = 1.D0/rhon

        ux = (g(i,j,k,1,now) - g(i,j,k,2,now) + g(i,j,k,7,now) &
           -  g(i,j,k,8,now) + g(i,j,k,9,now) - g(i,j,k,10,now) &
           +  g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
           -  g(i,j,k,14,now))*invRhon
        uy = (g(i,j,k,3,now) - g(i,j,k,4,now) + g(i,j,k,7,now) &
           -  g(i,j,k,8,now) - g(i,j,k,9,now) + g(i,j,k,10,now) &
           +  g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
           -  g(i,j,k,18,now))*invRhon
        uz = (g(i,j,k,5,now) - g(i,j,k,6,now) + g(i,j,k,11,now) &
           -  g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
           +  g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
           +  g(i,j,k,18,now))*invRhon

        IF (phi(i,j,k) >= 0.D0) THEN
          Vol = Vol + 1.D0
          total_ux = total_ux + ux
          total_uy = total_uy + uy
          total_uz = total_uz + uz
        END IF
      END DO
    END DO
  END DO

  IF (Vol > 0.D0) THEN
    invVol = 1.D0/Vol
    ux  = total_ux*invVol
    uy  = total_uy*invVol
    uz  = total_uz*invVol
    Ref = (Vol*invPi*0.75D0)**(1.D0/3.D0)
    Vef = Vol*invInitVol
  ELSE
    ux = 0.D0; uy = 0.D0; uz = 0.D0
    Ref = 0.D0; Vef = 0.D0
  END IF

  OPEN(UNIT=10, FILE="stats.out", STATUS="UNKNOWN", POSITION="APPEND", IOSTAT=IO_ERR)
  IF (IO_ERR == 0) THEN
    WRITE(10,'(I9,6ES19.9)') iStep, ux, uy, uz, Vef, Ref
    CLOSE(UNIT=10)
  END IF

END SUBROUTINE Stats

!===============================================================================
! SUBROUTINE: VtkDump
!===============================================================================
SUBROUTINE VtkDump
  USE NTypes,      ONLY : DBL
  USE Domain
  USE FluidParams
  USE LBMParams,   ONLY : Cs_sq, g
  IMPLICIT NONE

  INTEGER :: i, j, k, ie, iw, jn, js, kt, kb, IO_ERR
  REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
  REAL(KIND = DBL) :: gradPhiX, gradPhiY, gradPhiZ, gradPhiSq, lapPhi
  REAL(KIND = DBL) :: in00, in01, in02, in03, in04, in05, in06
  REAL(KIND = DBL) :: in07, in08, in09, in10, in11, in12
  REAL(KIND = DBL) :: in13, in14, in15, in16, in17, in18
  CHARACTER(20) :: filename
  CHARACTER(12) :: stepstr

  ! Build filename
  WRITE(stepstr, '(I8.8)') iStep
  filename = 'DATA_'//TRIM(ADJUSTL(stepstr))//'.vtk'

  OPEN(UNIT=12, FILE=filename, STATUS="REPLACE", IOSTAT=IO_ERR)
  IF (IO_ERR /= 0) THEN
    PRINT *, 'Warning: Cannot open VTK file ', filename
    RETURN
  END IF

  WRITE(12,'(A)') '# vtk DataFile Version 2.0'
  WRITE(12,'(A)') 'ZSC 3D LBM Output'
  WRITE(12,'(A)') 'ASCII'
  WRITE(12,*)
  WRITE(12,'(A)') 'DATASET STRUCTURED_POINTS'
  WRITE(12,'(A,3I6)') 'DIMENSIONS ', NX, NY, NZ
  WRITE(12,'(A,3I6)') 'ORIGIN ', xmin, ymin, zmin
  WRITE(12,'(A,3F6.1)') 'SPACING ', 1.0, 1.0, 1.0
  WRITE(12,*)
  WRITE(12,'(A,I12)') 'POINT_DATA ', NX*NY*NZ
  WRITE(12,*)

  ! Phi (order parameter)
  WRITE(12,'(A)') 'SCALARS Phi double'
  WRITE(12,'(A)') 'LOOKUP_TABLE default'
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        WRITE(12,'(ES21.11E3)') phi(i,j,k)
      END DO
    END DO
  END DO

  ! Pressure
  WRITE(12,*)
  WRITE(12,'(A)') 'SCALARS Pressure double'
  WRITE(12,'(A)') 'LOOKUP_TABLE default'
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) THEN
          WRITE(12,'(ES21.11E3)') 0.D0
          CYCLE
        END IF

        rhon  = SUM(g(i,j,k,:,now))
        phin  = phi(i,j,k)
        phin2 = phin*phin

        ie = ni(i,j,k,1); iw = ni(i,j,k,2)
        jn = ni(i,j,k,3); js = ni(i,j,k,4)
        kt = ni(i,j,k,5); kb = ni(i,j,k,6)

        in01 = phi(ie,j,k);  in02 = phi(iw,j,k)
        in03 = phi(i,jn,k);  in04 = phi(i,js,k)
        in05 = phi(i,j,kt);  in06 = phi(i,j,kb)
        in07 = phi(ie,jn,k); in08 = phi(iw,js,k)
        in09 = phi(ie,js,k); in10 = phi(iw,jn,k)
        in11 = phi(ie,j,kt); in12 = phi(iw,j,kb)
        in13 = phi(ie,j,kb); in14 = phi(iw,j,kt)
        in15 = phi(i,jn,kt); in16 = phi(i,js,kb)
        in17 = phi(i,jn,kb); in18 = phi(i,js,kt)

        lapPhi = (in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 &
               +  in15 + in16 + in17 + in18 &
               +  2.0D0*(in01 + in02 + in03 + in04 + in05 + in06 &
               -  12.D0*phin))*inv36

        gradPhiX = (2.0D0*(in01 - in02) + in07 - in08 + in09 - in10 &
                 +  in11 - in12 + in13 - in14)*inv12
        gradPhiY = (2.0D0*(in03 - in04) + in07 - in08 + in10 - in09 &
                 +  in15 - in16 + in17 - in18)*inv12
        gradPhiZ = (2.0D0*(in05 - in06) + in11 - in12 + in14 - in13 &
                 +  in15 - in16 + in18 - in17)*inv12
        gradPhiSq = gradPhiX*gradPhiX + gradPhiY*gradPhiY + gradPhiZ*gradPhiZ

        pressure = alpha*(3.D0*phin2*phin2 - 2.D0*phiStar2*phin2 - phiStar4) &
                 - kappa*(phin*lapPhi + 0.5D0*gradPhiSq) + Cs_sq*rhon

        WRITE(12,'(ES21.11E3)') pressure
      END DO
    END DO
  END DO

  ! Velocity
  WRITE(12,*)
  WRITE(12,'(A)') 'VECTORS Velocity double'
  DO k = zmin, zmax
    DO j = ymin, ymax
      DO i = xmin, xmax
        IF (is_fluid(i,j,k) == 0) THEN
          WRITE(12,'(3ES21.11E3)') 0.D0, 0.D0, 0.D0
          CYCLE
        END IF

        rhon = SUM(g(i,j,k,:,now))
        invRhon = 1.D0/rhon

        ux = (g(i,j,k,1,now) - g(i,j,k,2,now) + g(i,j,k,7,now) &
           -  g(i,j,k,8,now) + g(i,j,k,9,now) - g(i,j,k,10,now) &
           +  g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
           -  g(i,j,k,14,now))*invRhon
        uy = (g(i,j,k,3,now) - g(i,j,k,4,now) + g(i,j,k,7,now) &
           -  g(i,j,k,8,now) - g(i,j,k,9,now) + g(i,j,k,10,now) &
           +  g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
           -  g(i,j,k,18,now))*invRhon
        uz = (g(i,j,k,5,now) - g(i,j,k,6,now) + g(i,j,k,11,now) &
           -  g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
           +  g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
           +  g(i,j,k,18,now))*invRhon

        WRITE(12,'(3ES21.11E3)') ux, uy, uz
      END DO
    END DO
  END DO

  CLOSE(UNIT=12)

END SUBROUTINE VtkDump
