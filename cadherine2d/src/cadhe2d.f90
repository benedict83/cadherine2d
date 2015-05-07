module prec
!===========================================
! PREC
!===========================================
  use, intrinsic :: iso_fortran_env, only: int32, int64, sp => real32, dp => real64, &
                                           stdin => input_unit, &
                                           stdout => output_unit, &
                                           stderr => error_unit
  implicit none
  private
  public :: int32, int64, sp, dp
  public :: stdin, stdout, stderr
end module prec

module parametres
!===========================================
! PARAMETRES
!===========================================
  use prec
  implicit none
  character(6) :: nom_fichier
  character(12) :: nom_fichier_param
  character(10) :: nom_fichier_out
  integer ( int32 ) :: u_out
  real ( dp ) :: xmin,xmax
  real ( dp ) :: ymin,ymax
  real ( dp ) :: tbeg,tend
  real ( dp ) :: rho,l1,l2,eps
  real ( dp ) :: a0,a1
  real ( dp ) :: sigma
  real ( dp ) :: h1,h2,x0,y0,l
  real ( dp ) :: dx,dy,dt
  real ( dp ) :: oodx2
  integer ( int32 ) :: ndim
  integer ( int32 ) :: ns
  integer ( int32 ) :: iinit
end module parametres

module dimensions
!===========================================
! DIMENSIONS
!===========================================
  use prec
  implicit none
  integer ( int32 ), parameter :: nip = 19            ! Number of interior points
  integer ( int32 ) :: nup                            ! Number of unknown points
  integer ( int32 ) :: neqn                           ! Number of equations
contains
subroutine calc_neqn ( )
  use parametres
  implicit none

  nup = (nip+2)**ndim
  neqn = 2*nup
end subroutine calc_neqn
end module dimensions

module mymod
!===========================================
! MYMOD
!===========================================
  implicit none
  private
  public :: phi
contains
function phi (x,y)
!===========================================
! PHI computes the kernel phi
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in) :: x,y
  real ( dp ) :: phi
! Data dictionary: declare local variable types & definitions
  real ( dp ) :: distance
  real ( dp ) :: l3

! partie dependante de la dimension
  if (ndim == 2) then
    l3 = (2.0_dp**(0.5_dp))*l1
    distance = (x**2+y**2)**(0.5_dp)
  else
! dimension 1
    l3 = 2.0_dp*l1
    distance = abs(x)
  end if
  if (distance > l3) then
    phi = 0
  else if (distance > l1) then
    phi = -1
  else
    phi = -1
  end if
!!  write(stdout,*) l3,distance,phi
end function phi
end module mymod

module mymod2
!===========================================
! MYMOD2
!===========================================
  implicit none
  private
  public :: convol
contains
function convol (x,y,vec)
!===========================================
! CONVOL computes the convolution term
!===========================================
  use prec
  use parametres
  use dimensions
  use mymod, only: phi
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in) :: x,y
  real ( dp ), dimension(0:neqn-1), intent(in) :: vec
  real ( dp ) :: convol
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: ix,iy
  real ( dp ) :: position_x,position_y
  real ( dp ) :: termephi

! partie dependante de la dimension
  if(ndim == 2)then
! bord gauche i.e. {0}*[0,1]
  position_x = x
  position_y = y
  termephi = phi (position_x,position_y)
  convol = 0.25_dp * vec(nup) * termephi * dx * dy
  do iy=1,nip
    position_y = position_y-dy
    termephi = phi (position_x,position_y)
    convol = convol + 0.5_dp * vec(nup+iy*(nip+2)) * termephi * dx * dy
  end do
  position_y = position_y-dy
  termephi = phi (position_x,position_y)
  convol = convol + 0.25_dp * vec(nup+(nip+1)*(nip+2)) * termephi * dx * dy
! interieur i.e. ]0,1[*[0,1]
  do ix=1,nip
    position_x = position_x-dx
    position_y = y
    termephi = phi (position_x,position_y)
    convol = convol + 0.5_dp * vec(nup+ix) * termephi * dx * dy
    do iy=1,nip
      position_y = position_y-dy
      termephi = phi (position_x,position_y)
      convol = convol + vec(nup+iy*(nip+2)+ix) * termephi * dx * dy
    end do
    position_y = position_y-dy
    termephi = phi (position_x,position_y)
    convol = convol + 0.5_dp * vec(nup+(nip+1)*(nip+2)+ix) * termephi * dx * dy
  end do
! bord droit i.e. {1}*[0,1]
  position_x = position_x-dx
  position_y = y
  termephi = phi (position_x,position_y)
  convol = convol + 0.25_dp * vec(nup+nip+1) * termephi * dx * dy
  do iy=1,nip
    position_y = position_y-dy
    termephi = phi (position_x,position_y)
    convol = convol + 0.5_dp * vec(nup+iy*(nip+2)+nip+1) * termephi * dx * dy
  end do
  position_y = position_y-dy
  termephi = phi (position_x,position_y)
  convol = convol + 0.25_dp * vec(nup+(nip+1)*(nip+2)+nip+1) * termephi * dx * dy
  else
! dimension 1
  position_x = x
  position_y = 0.0_dp ! peu importe
  termephi = phi (position_x,position_y)
  convol = 0.5_dp * vec(nup) * termephi * dx
  do ix=1,nip
    position_x = position_x-dx
    termephi = phi (position_x,position_y)
    convol = convol + vec(nup+ix) * termephi * dx
  end do
  position_x = position_x-dx
  termephi = phi (position_x,position_y)
  convol = convol + 0.5_dp * vec(nup+(nip+1)) * termephi * dx
  end if

end function convol
end module mymod2

module mymod3
!===========================================
! MYMOD3
!===========================================
  implicit none
  private
  public :: fcadhe
contains
subroutine fcadhe (n,t,vec,dvec)
!===========================================
! FCADHE compute the value of
! dvec=(du,dv)
!===========================================
  use prec
  use parametres
  use dimensions
  use mymod2, only: convol
  implicit none
! Data dictionary: declare calling parameter types & definitions
  integer ( int32 ), intent(in) :: n     ! Moralement n == neqn
  real ( dp ), intent(in) :: t
  real ( dp ), dimension(0:n-1), intent(in) :: vec
  real ( dp ), dimension(0:n-1), intent(out) :: dvec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i,ix,iy
  real ( dp ) :: uleft,vleft,uright,vright
  real ( dp ) :: ulow,vlow,uup,vup
  real ( dp ) :: ui,vi
  real ( dp ) :: x,y,termeconvol
  real ( dp ) :: fixation,liberation
  real ( dp ) :: reacplus,reacmoins,reac

! big loop
  do i=0,nup-1
! left neighbour
    if(mod(i,nip+2) == 0)then
! Neumann homogeneous boundary condition
      uleft=vec(i+1)
      vleft=vec(nup+i+1)
    else
     uleft=vec(i-1)
     vleft=vec(nup+i-1)
    end if
! right neighbour
    if(mod(i,nip+2) == nip+1)then
! Neumann homogeneous boundary condition
      uright=vec(i-1)
      vright=vec(nup+i-1)
    else
      uright=vec(i+1)
      vright=vec(nup+i+1)
    end if
! partie dependante de la dimension
    if(ndim == 2)then
  ! lower neighbour
      if(i <= nip+1)then
  ! Neumann homogeneous boundary condition
        ulow=vec(i+nip+2)
        vlow=vec(nup+i+nip+2)
      else
        ulow=vec(i-nip-2)
        vlow=vec(nup+i-nip-2)
      end if
  ! upper neighbour
      if(i >= (nip+1)*(nip+2))then
  ! Neumann homogeneous boundary condition
        uup=vec(i-nip-2)
        vup=vec(nup+i-nip-2)
      else
        uup=vec(i+nip+2)
        vup=vec(nup+i+nip+2)
      end if
    else
! dimension 1
    end if
! the derivative
    ui=vec(i)
    vi=vec(nup+i)
! partie dependante de la dimension
    if(ndim == 2)then
      dvec(i)=sigma*oodx2*(uleft+uright+ulow+uup-4.0_dp*ui)
      dvec(nup+i)=0.0_dp
    else
! dimension 1
      dvec(i)=sigma*oodx2*(uleft+uright-2.0_dp*ui)
      dvec(nup+i)=0.0_dp
    end if
! appel a convol
    ix = mod(i,nip+2)
    iy = i/(nip+2)
    x = ix*dx
    y = iy*dy
    termeconvol=convol(x,y,vec)
    fixation = min(1.0_dp,(max(0.0_dp,1.0_dp-termeconvol)))
    liberation = eps * (1 - fixation)
! reaction 
    reacplus = fixation * ui * (rho-vi)
    reacmoins = -1.0_dp * liberation * vi
    reac = reacplus + reacmoins

    dvec(i) = dvec(i) - reac
    dvec(nup+i) = dvec(nup+i) + reac
  end do
end subroutine fcadhe
end module mymod3

program main
!===========================================
! MAIN solves a system of ODEs resulting from the 2-dimensional space
! discretization of the system of reaction-diffusion equations for the cadherin model
!
! u_t = sigma (u_{xx}+u_{yy}) -r(u,v)
! or
! u_t = sigma u_{xx} -r(u,v)
! v_t = r(u,v)
!
! and periodic boundary conditions
!
! We discretize the space variables with
! x_i=i/(nip+1) for i=1,...,nip
! or
! x_i=i/(nip+1) and y_i=i/(nip+1) for i=1,...,nip
! We obtain a system of neqn equations
! The spectral radius of the Jacobian can
! be estimated with the Gershgorin theorem
! Thus we provide an external function RADIUS,
! giving the spectral radius of the Jacobian
! matrix
!===========================================
  use prec
  use parametres
  use dimensions
  use mymod3, only: fcadhe
  implicit none
! Data dictionary: declare local variable types & definitions
!!  integer ( int32 ) :: u
  integer ( int32 ) :: k
  integer ( int32 ) :: idid
  integer ( int32 ), dimension(12) :: iwork
  real ( dp ) :: t0,t1
  real ( dp ), dimension(:), allocatable :: vec
  real ( dp ) :: rtol,atol
  real ( dp ) :: h
! If integrating with ROCK2 define work of length 4neqn
  real ( dp ), dimension(:), allocatable :: work

! parameters
  call param ( )
! systematic computations
  call sys ( )
  allocate (vec(0:neqn-1))
  allocate (work(1:7*neqn)) ! Work is of length 7neqn because the radius is computed externally (otherwise should be 8neqn)

  open(newunit=u_out,file=nom_fichier_out,access="sequential",status ="replace",form="formatted")

! condition initiale
  call init (vec)
! premier intervalle de temps
  t0 = tbeg
  t1 = tbeg + dt
! 0-eme sortie
  write(stdout,*) ''
! faire un header pour le paragraphe dans le fichier
  write(u_out,*) '# Time is',t0
! appeler out pour ecrire dans le fichier de sortie
  call out (t0,vec)

! boucle d'integration
  do k=1,ns
! required tolerance
    rtol=0.1_dp**2
    atol=rtol
! initial step size
    h=1.0e-4_dp
! Initialize iwork:
    iwork(1)=1  ! RADIUS returns an upper bound for the spectral radius
    iwork(2)=1  ! The Jacobian is constant (RADIUS is called once)
    iwork(3)=0  ! Return and solution at t1
    iwork(4)=0  ! Atol and rtol are scalars
    write(stdout,*) ''
    write(stdout,"(1x,a,f7.2,a,f7.2)") 'Integration du systeme de ',t0,' a ',t1
    call rock4 (neqn,t0,t1,h,vec,fcadhe,atol,rtol,work,iwork,idid)
! k-eme sortie
    write(stdout,*) '...'
    write(stdout,"(1x,a,f7.2)") 'Fait, et maintenant le temps est ',t1
! print statistics
    write(stdout,*) ''
!!    write(stdout,*) '--Solution is tabulated in file ',nom_fichier_out
!!    write(stdout,*) 'The value of dt is',dt
!!    write(stdout,*) 'The value of dx is',dx
!!    write(stdout,*) 'The value of oodx2 is',oodx2
    write(stdout,*) 'The value of IDID is',idid
    write(stdout,*) 'Max estimation of the spectral radius=',iwork(11)
    write(stdout,*) 'Min estimation of the spectral radius=',iwork(12)
    write(stdout,*) 'Max number of stages used=',iwork(10)
    write(stdout,*) 'Number of f evaluations for the spectral radius=',iwork(9)
    write(stdout,"(1x,a,i4,a,i4,a,i4,a,i3)") 'Number of f evaluations=', &
                                             iwork(5),' steps=',iwork(6),' accpt=',iwork(7),' rejct=',iwork(8)
! passer deux lignes dans le fichier IMPORTANT POUR GNUPLOT
    write(u_out,*) ''
    write(u_out,*) ''
! faire un header pour le paragraphe dans le fichier
    write(u_out,*) '# Time is',t1
! appeler out pour ecrire dans le fichier de sortie
    call out (t1,vec)
! prochain intervalle de temps
    t0 = t1
    t1 = t1 + dt
  end do
end program main

function radius ( )
!===========================================
! RADIUS gives an estimation of the spectral
! radius of the Jacobian matrix of the problem
! This is a bound for the whole interval and
! thus RADIUS is called once
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ) :: radius

  radius = 8.0_dp*nup*sigma + 2.0_dp
end function radius

subroutine init (vec)
!===========================================
! INIT computes the initial condition
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), dimension (0:neqn-1), intent(out) :: vec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i,ix,iy
  integer ( int32 ) :: u
  real ( dp ) :: x,r
  real ( dp ) :: xlaouonest,ylaouonest
  real ( dp ) :: gauss

  if (iinit == 0) then
! read initial condition from file
!!    open(newunit=u,file="fort.10",access="sequential",form ="formatted",status="old")
    open(newunit=u,file="fort.10",status="old")
    do i=0,nup-1
      read(u,*) x,vec(i)
    end do
    close(u)
  else if (iinit == 1) then
! petit plateau pour v autour de x0,y0, v nulle ailleurs

  else if (iinit == 2) then
! petite bosse pour v autour de x0,y0, v nulle ailleurs
! big loop
    do i=0,nup-1
      ix = mod(i,nip+2)
      iy = i/(nip+2)
      xlaouonest = ix*dx
      ylaouonest = iy*dy
      gauss = h2 * exp(-((xlaouonest-x0)**2+(ylaouonest-y0)**2)/(2.0_dp*l))
      vec(i) = 1.0_dp - gauss
      vec(nup+i) = gauss
!!      write(stdout,*) ix,iy,xlaouonest,ylaouonest,gauss
    end do
  else if (iinit == 3) then
! condition initiale aleatoire
    call init_random_seed ( )
! big loop
    do i=0,nup-1
      call random_number (r)
      vec(i) = 1-h2*r
      vec(nup+i) = h2*r
    end do
  else if (iinit == 4) then
! condition initiale aleatoire en une dimension
    call init_random_seed ( )
    do ix=0,nip+1
      call random_number (r)
      do iy=0,nip+1
        vec(iy*(nip+2)+ix) = 1-h2*r
        vec(nup+iy*(nip+2)+ix) = h2*r
      end do
    end do
  end if
end subroutine init

subroutine out (t,vec)
!===========================================
! 7.OUT writes t,x,y,u(x,y),v(x,y) in the file with label u_out
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in) :: t
  real ( dp ), dimension(0:neqn-1), intent(in) :: vec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: ix,iy
  real ( dp ) :: x,y,u,v
  real ( dp ) :: umin,vmin = 1000
  real ( dp ) :: umax,vmax = 0

  do ix=0,nip+1
    x = ix*dx
! partie dependante de la dimension
    if (ndim == 2) then
      do iy=0,nip+1
        y = iy*dy
        u = real(vec(iy*(nip+2)+ix))
        v = real(vec(nup+iy*(nip+2)+ix))
        write(u_out,"(5(1f9.4,1x))") t,x,y,u,v
        if (u < umin) then
          umin = u
        end if
        if (u > umax) then
          umax = u
        end if
        if (v < vmin) then
          vmin = v
        end if
        if (v > vmax) then
          vmax = v
        end if
        if (u < - 1.0e-3_dp) then
          write(stdout,*) 'Attention'
          write(stdout,*) x,y
          write(stdout,*) 'est un endroit ou u est negatif et vaut'
          write(stdout,*) u
        end if
        if (v < - 1.0e-3_dp) then
          write(stdout,*) 'Attention'
          write(stdout,*) x,y
          write(stdout,*) 'est un endroit ou v est negatif et vaut'
          write(stdout,*) v
        end if
      end do
    else
! dimension 1
      y = 0
      u = real(vec(ix))
      v = real(vec(nup+ix))
      write(u_out,"(5(1f9.4,1x))") t,x,y,u,v
    end if
  end do
  write(stdout,*) 'u est compris entre'
  write(stdout,*) umin,umax
  write(stdout,*) 'v est compris entre'
  write(stdout,*) vmin,vmax
end subroutine out

subroutine param ( )
!===========================================
! 8.PARAM reads parameters from a data file
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: ierror
  integer ( int32 ) :: u

  write(stdout,"(1x,a)",advance='no') 'nom du run (6 caracteres) ? '
  read (stdin,"(a6)") nom_fichier
  nom_fichier_param =  nom_fichier//'.param'
  nom_fichier_out =  nom_fichier//'.out'
  open(newunit=u,file=nom_fichier_param,status="old",form="formatted",action="read",iostat=ierror)
  if (ierror == 0) then
! Open was ok, read values
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror) ndim,xmin,xmax,ymin,ymax
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror) tbeg,tend,ns
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror) rho,l1,l2,eps
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror) sigma
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror) iinit,h1,h2,x0,y0,l
    close(u)
  end if
end subroutine param

subroutine sys ( )
!===========================================
! 9.SYS makes the systematic computations
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare local variable types & definitions
!!  real ( dp ) :: tiny = 1.0e-10_dp

! compute time step
!!  dt = (tend-tbeg+tiny)/ns
  dt = (tend-tbeg)/ns
! compute space step
  dx = (xmax-xmin)/(nip+1)
  dy = (ymax-ymin)/(nip+1)
! compute useful quantities
  oodx2 = 1.0_dp/(dx*dx)

end subroutine sys

subroutine init_random_seed ( )
!===========================================
! 10.INIT_RANDOM_SEED initializes a pseudo-random number sequence
!===========================================
  use prec
  implicit none
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i, n, clock
  integer ( int32 ), dimension(:), allocatable :: seed

  call random_seed (size = n)
  allocate(seed(n))
  call system_clock (count = clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed (put = seed)
  deallocate(seed)
end subroutine init_random_seed

