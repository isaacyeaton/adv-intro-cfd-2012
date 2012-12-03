
!##############################################################################
! This code solves for the viscous flow in a lid-driven cavity
!##############################################################################
module set_precision ! Sets the default precision for floating point numbers

  implicit none
  
  integer, parameter :: Sngl = Selected_Real_Kind(6 , 37)  
  integer, parameter :: Dbl  = Selected_Real_Kind(15, 307)
  integer, parameter :: Quad = Selected_Real_Kind(33, 4931)
  integer, parameter :: Prec = Dbl ! Set default precision to double precision

end module set_precision

  
!##############################################################################
module set_constants

  use set_precision, only : Prec

  implicit none

  real(Prec), parameter :: zero   = 0.0_Prec
  real(Prec), parameter :: tenth  = 0.1_Prec
  real(Prec), parameter :: sixth  = 1.0_Prec/6.0_Prec
  real(Prec), parameter :: fifth  = 0.2_Prec
  real(Prec), parameter :: fourth = 0.25_Prec
  real(Prec), parameter :: third  = 1.0_Prec/3.0_Prec
  real(Prec), parameter :: half   = 0.5_Prec
  real(Prec), parameter :: one    = 1.0_Prec
  real(Prec), parameter :: two    = 2.0_Prec
  real(Prec), parameter :: three  = 3.0_Prec
  real(Prec), parameter :: four   = 4.0_Prec
  real(Prec), parameter :: six    = 6.0_Prec
!$$$$$$   real(Prec), parameter :: big    = huge(zero) ! Largest resolvable #
!$$$$$$   real(Prec), parameter :: small  = tiny(zero) ! Smallest resolvable #

end module set_constants

  
!##############################################################################
module set_inputs ! Sets the input variables for the code

  use set_precision, only : Prec
  use set_constants, only : zero, one
  
  implicit none
  
  ! Set input values
  Integer, Parameter ::   imax = 65             ! Number of points in the x-direction (use odd numbers only)
  Integer, Parameter ::   jmax = 65             ! Number of points in the y-direction (use odd numbers only)
  Integer, Parameter ::   neq = 3               ! Number of equation to be solved ( = 3: mass, x-mtm, y-mtm)
  Integer ::              nmax = 500000         ! Maximum number of iterations
  Integer ::              iterout = 5000        ! Number of time steps between solution output
  Integer ::              imms = 0              ! Manufactured solution flag: = 1 for manuf. sol., = 0 otherwise
  Integer ::              isgs = 1              ! Symmetric Gauss-Seidel  flag: = 1 for SGS, = 0 for point Jacobi
  Integer ::              irstr = 0             ! Restart flag: = 1 for restart (file 'restart.in', = 0 for initial run
  Integer ::              ipgorder = 0          ! Order of pressure gradient: 0 = 2nd, 1 = 3rd (not needed)
  Integer ::              lim = 1               ! variable to be used as the limiter sensor (= 1 for pressure)

  Real(kind=Prec) ::      cfl  = 0.9_Prec       ! CFL number used to determine time step
  Real(kind=Prec) ::      Cx = 0.01_Prec        ! Parameter for 4th order artificial viscosity in x
  Real(kind=Prec) ::      Cy = 0.01_Prec        ! Parameter for 4th order artificial viscosity in y
  Real(kind=Prec) ::      toler = 1.e-10_Prec   ! Tolerance for iterative residual convergence
  Real(kind=Prec) ::      rkappa = 0.1_Prec     ! Time derivative preconditioning constant
  Real(kind=Prec) ::      Re = 100.0_Prec       ! Reynolds number = rho*Uinf*L/rmu
  Real(kind=Prec) ::      pinf = 0.801333844662_Prec ! Initial pressure (N/m^2) -> from MMS value at cavity center
  Real(kind=Prec) ::      uinf = one            ! Lid velocity (m/s)
  Real(kind=Prec) ::      rho = one             ! Density (kg/m^3)
  Real(kind=Prec) ::      xmin = zero           ! Cavity dimensions...: minimum x location (m)
  Real(kind=Prec) ::      xmax = 0.05_Prec      !                       maximum x location (m)
  Real(kind=Prec) ::      ymin = zero           !                       maximum y location (m)
  Real(kind=Prec) ::      ymax = 0.05_Prec      !                       maximum y location (m)
  Real(kind=Prec) ::      Cx2 = 0.0_Prec        ! Coefficient for 2nd order damping (not required)
  Real(kind=Prec) ::      Cy2 = 0.0_Prec        ! Coefficient for 2nd order damping (not required)
  Real(kind=Prec) ::      fsmall = 1.e-20_Prec  ! small parameter

  ! Derived input quantities (set in subroutine 'set_derived_inputs' below)
  Real(kind=Prec) ::      rhoinv =  -99.9_Prec  ! Inverse density, 1/rho (m^3/kg)
  Real(kind=Prec) ::      rlength = -99.9_Prec  ! Characteristic length (m) [cavity width]
  Real(kind=Prec) ::      rmu =     -99.9_Prec  ! Viscosity (N*s/m^2)
  Real(kind=Prec) ::      vel2ref = -99.9_Prec  ! Reference velocity squared (m^2/s^2)
  Real(kind=Prec) ::      dx =      -99.9_Prec  ! Delta x (m)
  Real(kind=Prec) ::      dy =      -99.9_Prec  ! Delta y (m)
  real(Prec) ::           rpi =     -99.9_Prec  ! Pi = 3.14159... (defined below)
  
  contains

  !#######
  subroutine set_derived_inputs

    implicit none
  
    rhoinv = one/rho                            ! Inverse density, 1/rho (m^3/kg)
    rlength = xmax - xmin                       ! Characteristic length (m) [cavity width]
    rmu = rho*uinf*rlength/Re                   ! Viscosity (N*s/m^2)
    vel2ref = uinf*uinf                         ! Reference velocity squared (m^2/s^2)
    dx = (xmax - xmin)/float(imax - 1)          ! Delta x (m)
    dy = (ymax - ymin)/float(jmax - 1)          ! Delta y (m)
    rpi = acos(-one)                            ! Pi = 3.14159... 
  
    write(*,*) 'rho,V,L,mu,Re: ',rho,uinf,rlength,rmu,Re

  end subroutine set_derived_inputs

end module set_inputs


!##############################################################################
module variables

  use set_precision, only : Prec
  use set_inputs, only : imax, jmax, neq
  
  implicit none

  Integer ::  ninit                                                  ! Initial iteration number (used for restart file)

  Real(kind=Prec),Dimension(imax,jmax,neq) :: u = -99.9_Prec         ! Solution vector [p, u, v]^T at each node
  Real(kind=Prec),Dimension(imax,jmax,neq) :: uold = -99.9_Prec      ! Previous (old) solution vector
  Real(kind=Prec),Dimension(imax,jmax,neq) :: s = -99.9_Prec         ! Source term
  Real(kind=Prec),Dimension(imax,jmax) ::     dt = -99.9_Prec        ! Local time step at each node
  Real(kind=Prec),Dimension(imax,jmax) ::     artviscx = -99.9_Prec  ! Artificial viscosity in x-direction
  Real(kind=Prec),Dimension(imax,jmax) ::     artviscy = -99.9_Prec  ! Artificial viscosity in y-direction
  Real(kind=Prec),Dimension(neq) ::           res = -99.9_Prec       ! Iterative residual for each equation
  Real(kind=Prec),Dimension(neq) ::           resinit = -99.9_Prec   ! Initial iterative residual for each equation (from iteration 1)
  Real(kind=Prec),Dimension(neq) ::           rL1norm = -99.9_Prec   ! L1 norm of discretization error for each equation
  Real(kind=Prec),Dimension(neq) ::           rL2norm = -99.9_Prec   ! L2 norm of discretization error for each equation
  Real(kind=Prec),Dimension(neq) ::           rLinfnorm = -99.9_Prec ! Linfinity norm of discretization error for each equation
  Real(kind=Prec) ::                          rtime = -99.9_Prec     ! Variable to estimate simulation time
  Real(kind=Prec) ::                          dtmin = +1.e99_Prec    ! Minimum time step for a given iteration (initialized large)

end module variables

  
!##############################################################################
module manufactured_solutions

  use set_precision, only : Prec
  use set_constants, only : zero, one, two, half, third, fourth, fifth, sixth, tenth
  use set_inputs, only : neq

  implicit none

  Real(kind=Prec),Dimension(neq) :: phi0=(/ fourth, 0.3_Prec, fifth /)          ! MMS constant 
  Real(kind=Prec),Dimension(neq) :: phix=(/ half, 0.15_Prec, sixth /)           ! MMS amplitude constant
  Real(kind=Prec),Dimension(neq) :: phiy=(/ 0.4_Prec, fifth, fourth /)          ! MMS amplitude constant
  Real(kind=Prec),Dimension(neq) :: phixy=(/ third, fourth, tenth /)            ! MMS amplitude constant
  Real(kind=Prec),Dimension(neq) :: apx=(/ half, third, 7._Prec/17._Prec /)     ! MMS frequency constant
  Real(kind=Prec),Dimension(neq) :: apy=(/ fifth, fourth, sixth /)              ! MMS frequency constant
  Real(kind=Prec),Dimension(neq) :: apxy=(/ 2._Prec/7._Prec, 0.4_Prec, third /) ! MMS frequency constant
  Real(kind=Prec),Dimension(neq) :: fsinx=(/ zero, one, zero /)                 ! MMS constant to determine sine vs. cosine
  Real(kind=Prec),Dimension(neq) :: fsiny=(/ one, zero, zero /)                 ! MMS constant to determine sine vs. cosine
  Real(kind=Prec),Dimension(neq) :: fsinxy=(/ one, one, zero /)                 ! MMS constant to determine sine vs. cosine
                                                                                ! Note: fsin = 1 means the sine function
                                                                                ! Note: fsin = 0 means the cosine function
                                                                                ! Note: arrays here refer to the 3 variables:
                                                                                !                                   [p, u, v]
  contains

!#######
! Function umms
  function umms(x,y,k)

  use set_precision, only : Prec
  use set_constants, only : one
  use set_inputs, only : rpi, rlength

  implicit none  
  
  Integer :: k                                ! Solution variable being evaluated (p = 1, u = 2, v = 3)

  Real(kind=Prec) ::  umms                    ! Define function as double precision
  Real(kind=Prec) ::  x                       ! x location input to MMS function
  Real(kind=Prec) ::  y                       ! y location input to MMS function
  Real(kind=Prec) ::  termx = -99.9_Prec      ! Temp variable
  Real(kind=Prec) ::  termy = -99.9_Prec      ! Temp variable
  Real(kind=Prec) ::  termxy = -99.9_Prec     ! Temp variable
  Real(kind=Prec) ::  argx = -99.9_Prec       ! Temp variable
  Real(kind=Prec) ::  argy = -99.9_Prec       ! Temp variable
  Real(kind=Prec) ::  argxy = -99.9_Prec      ! Temp variable

     
  ! This function returns the MMS exact solution
  
  argx = apx(k)*rpi*x/rlength
  argy = apy(k)*rpi*y/rlength
  argxy = apxy(k)*rpi*x*y/rlength/rlength
  termx = phix(k)*(fsinx(k)*sin(argx)+(one-fsinx(k))*cos(argx))
  termy = phiy(k)*(fsiny(k)*sin(argy)+(one-fsiny(k))*cos(argy))
  termxy = phixy(k)*(fsinxy(k)*sin(argxy)  &
    &        +(one-fsinxy(k))*cos(argxy))
  
  umms = phi0(k) + termx + termy + termxy
 
  return
 
  end function umms

!#######
! Function srcmms_mass
  
  function srcmms_mass(x,y)

  use set_precision, only : Prec
  use set_inputs, only : rho, rpi, rlength
    
  implicit none  
  
  Real(kind=Prec) ::  srcmms_mass  ! Define function as double precision
  Real(kind=Prec) ::  x            ! x location input to MMS function
  Real(kind=Prec) ::  y            ! y location input to MMS function
  Real(kind=Prec) ::  dudx         ! Temp variable: u velocity gradient in x direction
  Real(kind=Prec) ::  dvdy         ! Temp variable: v velocity gradient in y direction

! This function returns the MMS mass source term

  dudx = phix(2)*apx(2)*rpi/rlength*cos(apx(2)*rpi*x/rlength)  &
 &       + phixy(2)*apxy(2)*rpi*y/rlength/rlength  &
 &       * cos(apxy(2)*rpi*x*y/rlength/rlength)
  
  dvdy = -phiy(3)*apy(3)*rpi/rlength*sin(apy(3)*rpi*y/rlength)  &
 &       - phixy(3)*apxy(3)*rpi*x/rlength/rlength  &
 &       * sin(apxy(3)*rpi*x*y/rlength/rlength)

  srcmms_mass = rho*dudx + rho*dvdy

  return

  end function srcmms_mass

!#######
! Function srcmms_xmtm
  function srcmms_xmtm(x,y)

  use set_precision, only : Prec
  use set_inputs, only : rho, rpi, rmu, rlength
    
  implicit none  
  
  Real(kind=Prec) ::  srcmms_xmtm  ! Define function as double precision
  Real(kind=Prec) ::  x            ! x location input to MMS function
  Real(kind=Prec) ::  y            ! y location input to MMS function
  Real(kind=Prec) ::  dudx         ! Temp variable: u velocity gradient in x direction
  Real(kind=Prec) ::  dudy         ! Temp variable: u velocity gradient in y direction
  Real(kind=Prec) ::  termx        ! Temp variable
  Real(kind=Prec) ::  termy        ! Temp variable
  Real(kind=Prec) ::  termxy       ! Temp variable
  Real(kind=Prec) ::  uvel         ! Temp variable: u velocity 
  Real(kind=Prec) ::  vvel         ! Temp variable: v velocity 
  Real(kind=Prec) ::  dpdx         ! Temp variable: pressure gradient in x direction
  Real(kind=Prec) ::  d2udx2       ! Temp variable: 2nd derivative of u velocity in x direction
  Real(kind=Prec) ::  d2udy2       ! Temp variable: 2nd derivative of u velocity in y direction

! This function returns the MMS x-momentum source term

  termx = phix(2)*dsin(apx(2)*rpi*x/rlength)
  termy = phiy(2)*dcos(apy(2)*rpi*y/rlength)
  termxy = phixy(2)*dsin(apxy(2)*rpi*x*y/rlength/rlength)
  uvel = phi0(2) + termx + termy + termxy
  
  termx = phix(3)*dcos(apx(3)*rpi*x/rlength)
  termy = phiy(3)*dcos(apy(3)*rpi*y/rlength)
  termxy = phixy(3)*dcos(apxy(3)*rpi*x*y/rlength/rlength)
  vvel = phi0(3) + termx + termy + termxy
  
  dudx = phix(2)*apx(2)*rpi/rlength*dcos(apx(2)*rpi*x/rlength)  &
 &       + phixy(2)*apxy(2)*rpi*y/rlength/rlength  &
 &       * dcos(apxy(2)*rpi*x*y/rlength/rlength)
  
  dudy = -phiy(2)*apy(2)*rpi/rlength*dsin(apy(2)*rpi*y/rlength)  &
 &       + phixy(2)*apxy(2)*rpi*x/rlength/rlength  &
 &       * dcos(apxy(2)*rpi*x*y/rlength/rlength)
  
  dpdx = -phix(1)*apx(1)*rpi/rlength*dsin(apx(1)*rpi*x/rlength)  &
 &       + phixy(1)*apxy(1)*rpi*y/rlength/rlength  &
 &       * dcos(apxy(1)*rpi*x*y/rlength/rlength)

  d2udx2 = -phix(2)*(apx(2)*rpi/rlength)**2  &
 &         * dsin(apx(2)*rpi*x/rlength)  &
 &         - phixy(2)*(apxy(2)*rpi*y/rlength/rlength)**2  &
 &         * dsin(apxy(2)*rpi*x*y/rlength/rlength)  
 
  d2udy2 = -phiy(2)*(apy(2)*rpi/rlength)**2  &
 &         * dcos(apy(2)*rpi*y/rlength)  &
 &         - phixy(2)*(apxy(2)*rpi*x/rlength/rlength)**2  &
 &         * dsin(apxy(2)*rpi*x*y/rlength/rlength)
  
  srcmms_xmtm = rho*uvel*dudx + rho*vvel*dudy + dpdx  &
 &              - rmu*( d2udx2 + d2udy2 )

  return

  end function srcmms_xmtm

!#######
!  Function srcmms_ymtm

  function srcmms_ymtm(x,y)

  use set_precision, only : Prec
  use set_inputs, only : rho, rpi, rmu, rlength
    
  implicit none  
  
  Real(kind=Prec) ::  srcmms_ymtm  ! Define function as double precision
  Real(kind=Prec) ::  x            ! x location input to MMS function
  Real(kind=Prec) ::  y            ! y location input to MMS function
  Real(kind=Prec) ::  dvdx         ! Temp variable: v velocity gradient in x direction
  Real(kind=Prec) ::  dvdy         ! Temp variable: v velocity gradient in y direction
  Real(kind=Prec) ::  termx        ! Temp variable
  Real(kind=Prec) ::  termy        ! Temp variable
  Real(kind=Prec) ::  termxy       ! Temp variable
  Real(kind=Prec) ::  uvel         ! Temp variable: u velocity 
  Real(kind=Prec) ::  vvel         ! Temp variable: v velocity 
  Real(kind=Prec) ::  dpdy         ! Temp variable: pressure gradient in y direction
  Real(kind=Prec) ::  d2vdx2       ! Temp variable: 2nd derivative of v velocity in x direction
  Real(kind=Prec) ::  d2vdy2       ! Temp variable: 2nd derivative of v velocity in y direction

! This function returns the MMS y-momentum source term

  termx = phix(2)*dsin(apx(2)*rpi*x/rlength)
  termy = phiy(2)*dcos(apy(2)*rpi*y/rlength)
  termxy = phixy(2)*dsin(apxy(2)*rpi*x*y/rlength/rlength)
  uvel = phi0(2) + termx + termy + termxy
  
  termx = phix(3)*dcos(apx(3)*rpi*x/rlength)
  termy = phiy(3)*dcos(apy(3)*rpi*y/rlength)
  termxy = phixy(3)*dcos(apxy(3)*rpi*x*y/rlength/rlength)
  vvel = phi0(3) + termx + termy + termxy
  
  dvdx = -phix(3)*apx(3)*rpi/rlength*dsin(apx(3)*rpi*x/rlength)  &
  &       - phixy(3)*apxy(3)*rpi*y/rlength/rlength  &
  &       * dsin(apxy(3)*rpi*x*y/rlength/rlength)
  
  dvdy = -phiy(3)*apy(3)*rpi/rlength*dsin(apy(3)*rpi*y/rlength)  &
  &       - phixy(3)*apxy(3)*rpi*x/rlength/rlength  &
  &       * dsin(apxy(3)*rpi*x*y/rlength/rlength)
  
  dpdy = phiy(1)*apy(1)*rpi/rlength*dcos(apy(1)*rpi*y/rlength)  &
  &       + phixy(1)*apxy(1)*rpi*x/rlength/rlength  &
  &       * dcos(apxy(1)*rpi*x*y/rlength/rlength)
  
  d2vdx2 = -phix(3)*(apx(3)*rpi/rlength)**2  &
  &         * dcos(apx(3)*rpi*x/rlength)  &
  &         - phixy(3)*(apxy(3)*rpi*y/rlength/rlength)**2  &
  &         * dcos(apxy(3)*rpi*x*y/rlength/rlength)
  
  d2vdy2 = -phiy(3)*(apy(3)*rpi/rlength)**2  &
  &         * dcos(apy(3)*rpi*y/rlength)  &
  &         - phixy(3)*(apxy(3)*rpi*x/rlength/rlength)**2  &
  &         * dcos(apxy(3)*rpi*x*y/rlength/rlength)
  
  srcmms_ymtm = rho*uvel*dvdx + rho*vvel*dvdy + dpdy  &
  &              - rmu*( d2vdx2 + d2vdy2 )
  
  return
  
  end function srcmms_ymtm

!#######

end module manufactured_solutions


!##############################################################################
module functions

  implicit none

  contains

  !#######
  function first_derivative(delta, var_left, var_right)
  ! NOTE: This is an example of an easily testable function to compute a 1st derivative
  !       It is not currently used in the code
  
    use set_precision, only : Prec
    use set_constants, only : half
  
    real(Prec), intent(in) :: delta, var_left, var_right
    real(Prec)             :: first_derivative
  
  continue
  
    first_derivative = half*(var_right - var_left)/delta
  
  end function first_derivative
  
end module functions


!##############################################################################
module subroutines

  implicit none

  contains

  !#######
  subroutine output_file_headers 

  use set_inputs, only : imms

  ! Note: The vector of primitive variables is: 
  !               u = [p, u, v]^T
  ! Set up output files (history and solution)      
  open(30,file='history.dat',status='unknown')
  write(30,*) 'TITLE = "Cavity Iterative Residual History"'
  write(30,*) 'variables="Iteration""Time(s)""Res1""Res2""Res3"'
  
  open(40,file='cavity.dat',status='unknown')
  write(40,*) 'TITLE = "Cavity Field Data"'
  if(imms.eq.1) then
    write(40,*) 'variables="x(m)""y(m)""p(N/m^2)""u(m/s)""v(m/s)" &
              & "p-exact""u-exact""v-exact""DE-p""DE-u""DE-v"'
  elseif(imms.eq.0) then
    write(40,*) 'variables="x(m)""y(m)""p(N/m^2)""u(m/s)""v(m/s)"'
  else
    write(*,*) 'ERROR! imms must equal 0 or 1!!!'
    stop
  endif
  
  ! Header for Screen Output
  write(*,*) '   Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum'

  end subroutine output_file_headers

!#######
!  Subroutine initial

  subroutine initial
  
  use set_precision, only : Prec
  use set_constants, only : zero, one
  use set_inputs, only : irstr, imax, jmax, neq, uinf, pinf
  use variables, only : ninit, rtime, u, s, resinit

  implicit none

  Integer ::              i                    ! i index (x direction)
  Integer ::              j                    ! j index (y direction)
  Integer ::              k                    ! k index (# of equations)

  Real(kind=Prec) ::      x = -99.9_Prec       ! Temporary variable for x location
  Real(kind=Prec) ::      y = -99.9_Prec       ! Temporary variable for y location

  ! This subroutine sets inital conditions in the cavity 

  ! Note: The vector of primitive variables is: 
  !               u = [p, u, v]^T

  if(irstr.eq.0) then  ! Starting run from scratch
    
    ninit = 1          ! set initial iteration to one
    rtime = zero       ! set initial time to zero
    do k = 1, neq
      resinit(k) = one
    enddo
    do i = 1, imax
    do j = 1, jmax
      u(i,j,1) = pinf
      u(i,j,2) = zero
      u(i,j,3) = zero
      s(i,j,1) = zero
      s(i,j,2) = zero
      s(i,j,3) = zero
    enddo
    u(i,jmax,2) = uinf ! Initialize lid (top) to freestream velocity
    enddo
    
  elseif(irstr.eq.1) then  ! Restarting from previous run (file 'restart.in')
  
    open(60,file='restart.in',status='old')  ! Note: 'restart.in' must exist!
    read(60,*) ninit,rtime           ! Need to known current iteration # and time value
    read(60,*) (resinit(k),k=1,neq)  ! Needs initial iterative residuals for scaling
    do j = 1, jmax
    do i = 1, imax
      read(60,*)x,y,(u(i,j,k),k=1,neq)
    enddo
    enddo
    ninit = ninit + 1
    write(*,*) 'Restarting at iteration ',ninit
    close(60)
    
  else
    
    write(*,*) 'ERROR: irstr must equal 0 or 1!'
    stop
    
  endif

  return

  end subroutine initial
 
!#######
!  Subroutine set_boundary_conditions

  subroutine set_boundary_conditions
  
  use set_inputs, only : imms
      
  implicit none

  ! This subroutine determines the appropriate BC routines to call
    if(imms.eq.0) then
      call bndry
    elseif(imms.eq.1) then
      call bndrymms  
    else
      write(*,*) 'ERROR: imms must equal 0 or 1!'
      stop
    endif

  end subroutine set_boundary_conditions
     
!#######
!  Subroutine bndry

  subroutine bndry
  
  use set_precision, only : Prec
  use set_constants, only : zero, one, two, half
  use set_inputs, only : imax, jmax, uinf
  use variables, only : u
      
  implicit none

  Integer ::              i                    ! i index (x direction)
  Integer ::              j                    ! j index (y direction)

! This applies the cavity boundary conditions
 

!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************


  return

  end subroutine bndry

!#######
! Subroutine bndrymms
  subroutine bndrymms

  use set_precision, only : Prec
  use set_constants, only : two
  use set_inputs, only : imax, jmax, neq, xmax, xmin, ymax, ymin, rlength
  use variables, only : u
  use manufactured_solutions, only : umms

  implicit none
  
  Integer ::              i                    ! i index (x direction)
  Integer ::              j                    ! j index (y direction)
  Integer ::              k                    ! k index (# of equations)

  Real(kind=Prec) ::      x = -99.9_Prec       ! Temporary variable for x location
  Real(kind=Prec) ::      y = -99.9_Prec       ! Temporary variable for y location

! This applies the cavity boundary conditions for the manufactured solution

! Side Walls
  do j = 2, jmax-1
    y = (ymax - ymin)*float(j - 1)/float(jmax - 1)
    i = 1
    x = xmin
    do k = 1,neq
      u(i,j,k) = umms(x,y,k)
    enddo
    u(1,j,1) = two*u(2,j,1) - u(3,j,1)   ! 2nd Order BC
!$$$$$$     u(1,j,1) = u(2,j,1)                  ! 1st Order BC

    i=imax
    x = xmax
    do k = 1,neq
      u(i,j,k) = umms(x,y,k)
    enddo
    u(imax,j,1) = two*u(imax-1,j,1) - u(imax-2,j,1)   ! 2nd Order BC
!$$$$$$     u(imax,j,1) = u(imax-1,j,1)                       ! 1st Order BC
  enddo

! Top/Bottom Walls
  do i = 1, imax 
    x = (xmax - xmin)*float(i - 1)/float(imax - 1)
    j = 1
    y = ymin
    do k = 1,neq
      u(i,j,k) = umms(x,y,k)
    enddo
    u(i,1,1) = two*u(i,2,1) - u(i,3,1)   ! 2nd Order BC
!$$$$$$     u(i,1,1) = u(i,2,1)                  ! 1st Order BC

    j = jmax
    y = ymax
    do k = 1,neq
      u(i,j,k) = umms(x,y,k)
    enddo
    u(i,jmax,1) = two*u(i,jmax-1,1) - u(i,jmax-2,1)   ! 2nd Order BC
!$$$$$$     u(i,jmax,1) = u(i,jmax-1,1)                       ! 1st Order BC
  enddo

  return

  end subroutine bndrymms

!#######
  subroutine compute_time_step

  use set_precision, only : Prec
  use set_constants, only : one, two, four, half, fourth
  use set_inputs, only : vel2ref, rmu, rho, dx, dy, cfl, rkappa, imax, jmax
  use variables, only : u, dt, dtmin

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      dtvisc = -99.9_Prec      ! Viscous time step stability criteria (constant over domain)
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Local velocity squared
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared paramete for time derivative preconditioning
  Real(kind=Prec) ::      lambda_x = -99.9_Prec    ! Max absolute value eigenvalue in (x,t)
  Real(kind=Prec) ::      lambda_y = -99.9_Prec    ! Max absolute value eigenvalue in (y,t)
  Real(kind=Prec) ::      lambda_max = -99.9_Prec  ! Max absolute value eigenvalue (used in convective time step computation)
  Real(kind=Prec) ::      dtconv = -99.9_Prec      ! Local convective time step restriction


!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************


  end subroutine compute_time_step

  !#######
  subroutine compute_source_terms

  use set_precision, only : Prec
!  use set_constants, only : zero, one, two, four, half, fourth
  use set_inputs, only : imax, jmax, imms, rlength, xmax, xmin, ymax, ymin
  use variables, only : s
  use manufactured_solutions

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location

  ! Evaluate Source Terms Once at Beginning (only interior points; will be zero for standard cavity)
  
  do j = 2, jmax-1
  do i = 2, imax-1
    x = (xmax - xmin)*float(i - 1)/float(imax - 1)
    y = (ymax - ymin)*float(j - 1)/float(jmax - 1)
    s(i,j,1) = float(imms)*srcmms_mass(x,y)
    s(i,j,2) = float(imms)*srcmms_xmtm(x,y)
    s(i,j,3) = float(imms)*srcmms_ymtm(x,y)
  enddo
  enddo

  end subroutine compute_source_terms

  !#######
  subroutine Compute_Artificial_Viscosity

  use set_precision, only : Prec
  use set_constants, only : zero, one, two, four, six, half, fourth
  use set_inputs, only : imax, jmax, lim, rho, dx, dy, Cx, Cy, Cx2, Cy2, &
           &             fsmall, vel2ref, rkappa
  use variables, only : u, artviscx, artviscy

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Local velocity squared
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared paramete for time derivative preconditioning
  Real(kind=Prec) ::      lambda_x = -99.9_Prec    ! Max absolute value e-value in (x,t)
  Real(kind=Prec) ::      lambda_y = -99.9_Prec    ! Max absolute value e-value in (y,t)
  Real(kind=Prec) ::      d4pdx4 = -99.9_Prec      ! 4th derivative of pressure w.r.t. x
  Real(kind=Prec) ::      d4pdy4 = -99.9_Prec      ! 4th derivative of pressure w.r.t. y



!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************



  end subroutine Compute_Artificial_Viscosity

  !#######


  subroutine point_Jacobi

  use set_precision, only : Prec
  use set_constants, only :  two, three, six, half
  use set_inputs, only : imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, &
                       & xmax, xmin, ymax, ymin, rmu, vel2ref
  use variables, only : u, uold, artviscx, artviscy, dt, s

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      dpdx = -99.9_Prec        ! First derivative of pressure w.r.t. x
  Real(kind=Prec) ::      dudx = -99.9_Prec        ! First derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      dvdx = -99.9_Prec        ! First derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      dpdy = -99.9_Prec        ! First derivative of pressure w.r.t. y
  Real(kind=Prec) ::      dudy = -99.9_Prec        ! First derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      dvdy = -99.9_Prec        ! First derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      d2udx2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      d2vdx2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      d2udy2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      d2vdy2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared parameter for time derivative preconditioning
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Velocity squared 


  ! Point Jacobi method

!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************



  end subroutine point_Jacobi

  !#######

  
  subroutine SGS_forward_sweep

  use set_precision, only : Prec
  use set_constants, only :  two, three, six, half
  use set_inputs, only : imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, &
                       & xmax, xmin, ymax, ymin, rmu, vel2ref
  use variables, only : u, artviscx, artviscy, dt, s

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      dpdx = -99.9_Prec        ! First derivative of pressure w.r.t. x
  Real(kind=Prec) ::      dudx = -99.9_Prec        ! First derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      dvdx = -99.9_Prec        ! First derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      dpdy = -99.9_Prec        ! First derivative of pressure w.r.t. y
  Real(kind=Prec) ::      dudy = -99.9_Prec        ! First derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      dvdy = -99.9_Prec        ! First derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      d2udx2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      d2vdx2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      d2udy2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      d2vdy2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared parameter for time derivative preconditioning
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Velocity squared 


  ! Symmetric Gauss-Siedel: Forward Sweep

!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************



  end subroutine SGS_forward_sweep

  !#######
  subroutine SGS_backward_sweep

  use set_precision, only : Prec
  use set_constants, only :  two, three, six, half
  use set_inputs, only : imax, jmax, ipgorder, rho, rhoinv, dx, dy, rkappa, &
                       &  xmax, xmin, ymax, ymin, rmu, vel2ref
  use variables, only : u, artviscx, artviscy, dt, s

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)

  Real(kind=Prec) ::      dpdx = -99.9_Prec        ! First derivative of pressure w.r.t. x
  Real(kind=Prec) ::      dudx = -99.9_Prec        ! First derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      dvdx = -99.9_Prec        ! First derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      dpdy = -99.9_Prec        ! First derivative of pressure w.r.t. y
  Real(kind=Prec) ::      dudy = -99.9_Prec        ! First derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      dvdy = -99.9_Prec        ! First derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      d2udx2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. x
  Real(kind=Prec) ::      d2vdx2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. x
  Real(kind=Prec) ::      d2udy2 = -99.9_Prec      ! Second derivative of x velocity w.r.t. y
  Real(kind=Prec) ::      d2vdy2 = -99.9_Prec      ! Second derivative of y velocity w.r.t. y
  Real(kind=Prec) ::      beta2 = -99.9_Prec       ! Beta squared parameter for time derivative preconditioning
  Real(kind=Prec) ::      uvel2 = -99.9_Prec       ! Velocity squared 

  ! Symmetric Gauss-Siedel: Backward Sweep

!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************



  end subroutine SGS_backward_sweep

  !#######
  subroutine pressure_rescaling

  use set_precision, only : Prec
  use set_inputs, only : imax, jmax, imms, xmax, xmin, ymax, ymin, rlength, pinf
  use variables, only : u 
  use manufactured_solutions, only : umms

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              iref                     ! i index location of pressure rescaling point
  Integer ::              jref                     ! j index location of pressure rescaling point

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location
  Real(kind=Prec) ::      deltap = -99.9_Prec      ! Temporary variable: y location

  iref = (imax-1)/2+1     ! Set reference pressure to center of cavity
  jref = (jmax-1)/2+1
  if(imms.eq.1) then
    x = (xmax - xmin)*float(iref - 1)/float(imax - 1)
    y = (ymax - ymin)*float(jref - 1)/float(jmax - 1)
    deltap = u(iref,jref,1) - umms(x,y,1) ! Constant in MMS
  else
    deltap = u(iref,jref,1) - pinf ! Reference pressure 
  endif
  do i = 1, imax
  do j = 1, jmax
    u(i,j,1) = u(i,j,1) - deltap
  enddo
  enddo

  end subroutine pressure_rescaling

  !#######
  subroutine check_iterative_convergence(n,conv)

  use set_precision, only : Prec
  use set_constants, only : zero
  use set_inputs, only : imax, jmax, neq, fsmall
  use variables, only : u, uold, dt, res, resinit, ninit, rtime, dtmin

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                        ! k index (over # of equations)
  Integer ::              n                        ! Current iteration number

  Real(kind=Prec) ::      conv                     ! Minimum of iterative residual norms from three equations
  
  ! Compute iterative residuals to monitor iterative convergence

!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************


  conv = amax1(res(1),res(2),res(3)) ! Place current maximum scaled iterative residual into variable 'conv'
  
  ! Write iterative residuals every 10 iterations
  if( (mod(n,10).eq.0).or.n.eq.ninit ) then
    write (30,*) n,rtime,(res(k),k=1,3)
    write(*,300) n,rtime,dtmin,(res(k),k=1,3)
300   format(i8,2(e11.3),3(e14.5))
  endif   
       
  ! Write header for iterative residuals every 200 iterations
  if( (mod(n,200).eq.0).or.n.eq.ninit ) then
    write(*,*) '   Iter. Time (s)   dt (s)      Continuity    x-Momentum    y-Momentum'
  endif        

  end subroutine check_iterative_convergence

  !#######
  subroutine write_output(n)

  use set_precision, only : Prec
  use set_inputs, only : imax, jmax, neq, xmax, xmin, ymax, ymin, &
      &                  rlength, imms
  use variables, only : u, dt, resinit, ninit, rtime
  use manufactured_solutions, only : umms

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                        ! k index (over # of equations)
  Integer ::              n                        ! Current iteration number

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location
   
  ! Field output
  write(40,*) 'zone T="n = ',n,' " '
  write(40,*) 'I=',imax,' J=',jmax
  write(40,*) 'DATAPACKING=POINT'
  if(imms.eq.1) then
    do j = 1, jmax
    do i = 1, imax
      x = (xmax - xmin)*float(i - 1)/float(imax - 1)
      y = (ymax - ymin)*float(j - 1)/float(jmax - 1)
      write(40,*)x,y,(u(i,j,k),k=1,neq),  &
         &    (umms(x,y,k),k=1,neq),  &
         &    ((u(i,j,k) - umms(x,y,k)),k=1,neq)
    enddo
    enddo
  elseif(imms.eq.0) then
    do j = 1, jmax
    do i = 1, imax
      x = (xmax - xmin)*float(i - 1)/float(imax - 1)
      y = (ymax - ymin)*float(j - 1)/float(jmax - 1)
      write(40,*)x,y,(u(i,j,k),k=1,neq)
    enddo
    enddo
  else
    write(*,*) 'ERROR! imms must equal 0 or 1!!!'
    stop
  endif
  
  ! Restart file: overwrites every 'iterout' iteration
  open(50,file='restart.out',status='unknown')  
  write(50,*) n, rtime
  write(50,*) (resinit(k),k=1,neq)
  do j = 1, jmax
  do i = 1, imax
    x = (xmax - xmin)*float(i - 1)/float(imax - 1)
    y = (ymax - ymin)*float(j - 1)/float(jmax - 1)
    write(50,*)x,y,(u(i,j,k),k=1,neq)
  enddo
  enddo
  close(50)

  end subroutine write_output

  !#######
  subroutine Discretization_Error_Norms

  use set_precision, only : Prec
  use set_constants, only : zero
  use set_inputs, only : imax, jmax, neq, imms, xmax, xmin, ymax, ymin, rlength
  use variables, only : u, rL1norm, rL2norm, rLinfnorm 
  use manufactured_solutions, only : umms

  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                        ! k index (over # of equations)

  Real(kind=Prec) ::      x = -99.9_Prec           ! Temporary variable: x location
  Real(kind=Prec) ::      y = -99.9_Prec           ! Temporary variable: y location
  Real(kind=Prec) ::      DE = -99.9_Prec          ! Discretization error (absolute value)

  if(imms.eq.1) then ! Only compute discretization error norms for manufactured solution

!**************************************************************
!************ADD CODING HERE FOR INTRO CFD STUDENTS************
!**************************************************************


  endif

  end subroutine Discretization_Error_Norms

end module subroutines


!##############################################################################
!####################### P R O G R A M  C A V I T Y ###########################
!##############################################################################
Program Cavity

  ! Note: The vector of primitive variables is: 
  !               u = [p, u, v]^T
 
  use set_precision, only : Prec
  use set_constants, only : zero
  use set_inputs, only : imax, jmax, xmax, xmin, ymax, ymin, rlength, nmax, &
                  &      neq, toler, iterout, isgs, set_derived_inputs
  use variables, only : ninit, u, uold, res, resinit, artviscx, artviscy, &
                  &      rtime, dtmin
  use subroutines  
      
  implicit none
  
  Integer ::              i                        ! i index (x direction)
  Integer ::              j                        ! j index (y direction)
  Integer ::              k                    ! k index (# of equations)
  Integer ::              n                    ! Iteration number index
 
  Real(kind=Prec) ::      conv = -99.9_Prec    ! Minimum of iterative residual norms from three equations
  
  ! Set derived input quantities
  call set_derived_inputs
 
  ! Set up headers for output files
  call output_file_headers
 
  ! Set Initial Profile for u vector
  call initial
 
  ! Set Boundary Conditions for u
  call set_boundary_conditions
 
  ! Write out inital conditions to solution file
  call write_output(ninit)
 
  ! Initialize Artificial Viscosity arrays to zero (note: artviscx(i,j) and artviscy(i,j)
  do j = 1, jmax
  do i = 1, imax
    artviscx(i,j) = zero
    artviscy(i,j) = zero
  enddo
  enddo
      
  ! Evaluate Source Terms Once at Beginning (only interior points; will be zero for standard cavity)
  call compute_source_terms
  
  ! ========== Main Loop ==========
  do n = ninit, nmax
 
    ! Calculate time step   
    call compute_time_step

    ! Save u values at time level n (u and uold are 2D arrays)
    do j = 1, jmax
    do i = 1, imax
    do k = 1, neq
      uold(i,j,k) = u(i,j,k)
    enddo
    enddo
    enddo

    if (isgs.eq.1) then ! ==Symmetric Gauss Seidel==
      
      ! Artificial Viscosity
      call Compute_Artificial_Viscosity
      
      ! Symmetric Gauss-Siedel: Forward Sweep
      call SGS_forward_sweep
   
      ! Set Boundary Conditions for u
      call set_boundary_conditions
   
      ! Artificial Viscosity
      call Compute_Artificial_Viscosity
   
      ! Symmetric Gauss-Siedel: Backward Sweep
      call SGS_backward_sweep

      ! Set Boundary Conditions for u
      call set_boundary_conditions
 
    elseif (isgs.eq.0) then ! ==Point Jacobi==

      ! Artificial Viscosity
      call Compute_Artificial_Viscosity
      
      ! Symmetric Gauss-Siedel: Forward Sweep
      call point_Jacobi
   
      ! Set Boundary Conditions for u
      call set_boundary_conditions

    else
      write(*,*) 'ERROR: isgs must be 0 or 1!!!'
      stop
    endif
 
    ! Pressure Rescaling (based on center point)
    call pressure_rescaling
 
    ! Update the time
    rtime = rtime + dtmin
                 
    ! Check iterative convergence using L2 norms of iterative residuals
    call check_iterative_convergence(n, conv)
    
    if(conv.le.toler) then
      write(30,*) n,rtime,(res(k),k=1,3)
      goto 200  !exit iteration loop if converged
    else
      continue
    endif
                
    ! Output solution and restart file every 'iterout' steps
    if( (mod(n,iterout).eq.0) ) then
      call write_output(n)
    endif
 
  enddo  ! ========== End Main Loop ==========
  
    
  write(*,*) 'Solution failed to converge in ',nmax,' iterations!!!'
    
  goto 201
    
  200  continue  ! go here once solution is converged

  write(*,*) 'Solution converged in ',n,' iterations!!!'

  201  continue

  ! Calculate and Write Out Discretization Error Norms (will do this for MMS only)
  call Discretization_Error_Norms 

  ! Output solution and restart file 
  call write_output(n)

  ! Close open files
  close(50)
  close(30)
  close(40)

end Program Cavity
