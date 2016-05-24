module prec
!===========================================
! PREC
!===========================================
  use, intrinsic :: iso_fortran_env, only: int32                , int64,        &
                                           sp => real32         , dp => real64, &
                                           stdin => input_unit  ,               &
                                           stdout => output_unit,               &
                                           stderr => error_unit
  implicit none
  private
  public :: int32, int64, sp, dp
  public :: stdin, stdout, stderr
end module prec

module dimensions
!===========================================
! DIMENSIONS AND PARAMETERS
!===========================================
  use prec
  implicit none
  integer   ( int32 ), parameter :: nip = 99           ! Number of interior points
  integer   ( int32 ) :: nup1d                          ! Number of unknown points in 1d
  integer   ( int32 ) :: nuptotal                       ! Number of total unknown points
  integer   ( int32 ) :: neqn                           ! Number of equations
  character ( 3 )     :: bc
  character ( 8 )     :: nom_fichier
  character ( 14 )    :: nom_fichier_param
  character ( 12 )    :: nom_fichier_out
  integer   ( int32 ) :: u_out
  real      ( dp )    :: xmin,xmax
  real      ( dp )    :: ymin,ymax
  real      ( dp )    :: tbeg,tend
  real      ( dp )    :: rho,d1,d2,integrale,eps1
  real      ( dp )    :: sigma,eps2
  real      ( dp )    :: h1,h2,x0,y0,l
  real      ( dp )    :: dx,dy,dt
  real      ( dp )    :: oodx2,soodx2,e2soodx2
  integer   ( int32 ) :: ndim
  integer   ( int32 ) :: ns
  integer   ( int32 ) :: iinit
contains
subroutine param ( )
!===========================================
! PARAM reads parameters from a data file
!===========================================
  implicit none
! Data dictionary: declare local variable types & definitions
  integer   ( int32 ) :: ierror
  integer   ( int32 ) :: u_param

  write (stdout,"(1x,a)",advance='no') 'nom du run (8 caracteres) ? '
  read (stdin,"(a8)") nom_fichier
  nom_fichier_param =  nom_fichier//'.param'
  nom_fichier_out   =  nom_fichier//'.out'
  open (newunit=u_param,file=nom_fichier_param,status="old",form="formatted",action="read",iostat=ierror)
  if (ierror == 0) then
! Open was ok, read values
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror) ndim,xmin,xmax,ymin,ymax
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror) ns,tbeg,tend
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror) rho,d1,d2,integrale,eps1
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror) sigma,eps2
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror) iinit,h1,h2,x0,y0,l
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror)
    read (u_param,*,iostat=ierror) bc
    close (u_param)
  end if
end subroutine param

subroutine calc_neqn ( )
  implicit none

  select case (bc)
  case ("per")
    nup1d = nip+1 ! if periodic BC, then nip+1 unknown points
  case ("neu")
    nup1d = nip+2 ! if Neumann  BC, then nip+2 unknown points
  case default
  end select
  nuptotal = nup1d**ndim
  neqn = 2*nuptotal
end subroutine calc_neqn
end module dimensions

module mod_functions
!===========================================
! MODULE FOR FUNCTIONS PHI F AND G
!===========================================
  implicit none
  private
  public :: phi, fg
contains
function phi (x, y)
!===========================================
! PHI computes the kernel phi
!===========================================
  use prec
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in) :: x, y
  real ( dp ) :: phi
! Data dictionary: declare local variable types & definitions
  real ( dp ) :: norme
  real ( dp ) :: k1, k2

  select case (ndim)
  case (1)
    norme = abs(x)
  case (2)
    norme = (x**2+y**2)**(0.5_dp)
  case default
  end select
!!  select case ()
!!  case (gaussienne)
!!    k1  = exp(-(module/d1)**2)
!!    k2  = 2.0_dp * (exp(-(module/d2)**2))
!!    phi = (k1 - k2) / integrale
!!  case (creneau)
  if (norme > d1) then
    phi = 0.0_dp
  else if (norme > d2) then
    phi = 0.25_dp/d2   ! Il faut que l'integrale vaille  1/4 (soit  0.25) sur [d2, d1]
  else
    phi = - 0.25_dp/d2 ! Il faut que l'integrale vaille -1/4 (soit -0.25) sur [0, d2]
  end if
!!  end select
end function phi
function fg (h)
!===========================================
! FG computes the functions F and G
!===========================================
  use prec
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in)   :: h
  real ( dp ), dimension(2) :: fg
!!  fg (1) = 1.0_dp - tanh( 4.0_dp * ( h - rho * 0.5_dp ) )
!!  fg (2) = 1.0_dp
!! ajout pour un modele sans convolution, ici on calcule fg (v) donc il faut imaginer que h=v
  fg (1) = 0.5_dp + tanh( h )
  fg (2) = 1.0_dp - tanh( 2.0_dp * h )
end function fg
end module mod_functions

module mod_convol
!===========================================
! MODULE FOR CONVOL
!===========================================
  implicit none
  private
  public :: convol
contains
function convol (x, y, vec)
!===========================================
! CONVOL computes the convolution term
!===========================================
  use prec
  use dimensions
  use mod_functions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real      ( dp ),                      intent(in) :: x, y
  real      ( dp ), dimension(0:neqn-1), intent(in) :: vec
  real      ( dp )                                  :: convol
! Data dictionary: declare local variable types & definitions
  integer   ( int32 ) :: ix, iy
  real      ( dp )    :: position_x, position_y
  real      ( dp )    :: termephi

  select case (ndim)
  case (1)
    position_x = x
    position_y = 0.0_dp ! peu importe
    termephi   = phi (position_x, position_y)
    convol     = 0.5_dp * vec(nuptotal) * termephi * dx
    xloop_1d: do ix = 1, nip
      position_x = position_x-dx ! position_x = x - ix*dx
      termephi   = phi (position_x, position_y)
      convol     = convol + vec(nuptotal+ix) * termephi * dx
    end do xloop_1d
    position_x = position_x-dx ! position_x = x - (nip+1)*dx
    termephi   = phi (position_x, position_y)
    select case (bc)
    case ("per")
      convol = convol + 0.5_dp * vec(nuptotal) * termephi * dx
    case ("neu")
      convol = convol + 0.5_dp * vec(nuptotal+(nip+1)) * termephi * dx
    case default
    end select
  case (2)
! bord gauche i.e. {0}*[0, 1]
    position_x = x
    position_y = y
    termephi = phi (position_x, position_y)
    convol   = 0.25_dp * vec(nuptotal) * termephi * dx * dy
    yloop_gauche: do iy = 1, nip
      position_y = position_y-dy ! position_y = y - iy*dy
      termephi   = phi (position_x, position_y)
      convol     = convol + 0.5_dp * vec(nuptotal+iy*nup1d) * termephi * dx * dy
    end do yloop_gauche
    position_y = position_y-dy ! position_y = y - (nip+1)*dy
    termephi   = phi (position_x, position_y)
    select case (bc)
    case ("per")
      convol = convol + 0.25_dp * vec(nuptotal) * termephi * dx * dy
      convol = 2 * convol ! in case of periodic BC, we need to take two times the left slice
    case ("neu")
      convol = convol + 0.25_dp * vec(nuptotal+(nip+1)*nup1d) * termephi * dx * dy
    case default
    end select
! interieur i.e. ]0, 1[*[0, 1]
    xloop_2d: do ix = 1, nip
      position_x = position_x-dx ! position_x = x - ix*dx
      position_y = y
      termephi   = phi (position_x, position_y)
      convol     = convol + 0.5_dp * vec(nuptotal+ix) * termephi * dx * dy
      yloop_interieur: do iy = 1, nip
        position_y = position_y-dy ! position_y = y - iy*dy
        termephi   = phi (position_x, position_y)
        convol     = convol + vec(nuptotal+iy*nup1d+ix) * termephi * dx * dy
      end do yloop_interieur
      position_y = position_y-dy ! position_y = y - (nip+1)*dy
      termephi   = phi (position_x, position_y)
      select case (bc)
      case ("per")
        convol = convol + 0.5_dp * vec(nuptotal+ix) * termephi * dx * dy
      case ("neu")
        convol = convol + 0.5_dp * vec(nuptotal+(nip+1)*nup1d+ix) * termephi * dx * dy
      case default
      end select
    end do xloop_2d
! bord droit i.e. {1}*[0, 1]
    select case (bc)
    case ("per")
    ! do nothing : the right slice was already integrated with the left one  
    case ("neu")
    position_x = position_x-dx ! position_x = x - (nip+1)*dx
    position_y = y
    termephi   = phi (position_x, position_y)
    convol     = convol + 0.25_dp * vec(nuptotal+nip+1) * termephi * dx * dy
    yloop_droit: do iy = 1, nip
      position_y = position_y-dy ! position_y = y - iy*dy
      termephi   = phi (position_x, position_y)
      convol     = convol + 0.5_dp * vec(nuptotal+iy*nup1d+nip+1) * termephi * dx * dy
    end do yloop_droit
    position_y = position_y-dy ! position_y = y - (nip+1)*dy
    termephi   = phi (position_x, position_y)
    convol = convol + 0.25_dp * vec(nuptotal+(nip+1)*nup1d+nip+1) * termephi * dx * dy
    case default
    end select
  case default
  end select
end function convol
end module mod_convol

module mod_fcadhe
!===========================================
! MODULE FOR FCADHE
!===========================================
  implicit none
  private
  public :: fcadhe
contains
subroutine fcadhe (n, t, vec, dvec)
!===========================================
! FCADHE compute the value of
! dvec=(du, dv)
!===========================================
  use prec
  use dimensions
  use mod_functions
  use mod_convol, only: convol
  implicit none
! Data dictionary: declare calling parameter types & definitions
  integer ( int32 ), intent(in) :: n     ! Moralement n == neqn
  real ( dp )      , intent(in) :: t
  real ( dp )      , dimension(0:n-1), intent(in)  :: vec
  real ( dp )      , dimension(0:n-1), intent(out) :: dvec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i, ix, iy
  real ( dp ) :: uleft, uright
  real ( dp ) :: ulow, uup
  real ( dp ) :: vleft, vright
  real ( dp ) :: vlow, vup
  real ( dp ) :: ui, vi
  real ( dp ) :: x, y
  real ( dp ) :: termeconvol
  real ( dp ) :: h
  real ( dp ), dimension (2) :: fixlib
  real ( dp ) :: fix, lib
  real ( dp ) :: reacplus, reacmoins
  real ( dp ) :: termereac

  bigloop: do i = 0, nuptotal-1
! left neighbour
    if(mod(i, nup1d) == 0)then
      select case (bc)
      case ("per") ! periodic boundary condition
        uleft = vec(nip+i)
        vleft = vec(nuptotal+nip+i)
      case ("neu") ! Neumann homogeneous boundary condition
        uleft = vec(i+1)
        vleft = vec(nuptotal+i+1)
      case default
      end select
    else
      uleft = vec(i-1)
      vleft = vec(nuptotal+i-1)
    end if
! right neighbour
    if(mod(i, nup1d) == nup1d-1)then
      select case (bc)
      case ("per")    ! periodic boundary condition
        uright = vec(-nip+i)
        vright = vec(nuptotal-nip+i)
!!        dvec (i)          =  dvec (-nip+i-1)         ! on force la valeur de du
!!        dvec (nuptotal+i) =  dvec (nuptotal-nip+i-1) ! on force la valeur de dv
!!        cycle bigloop ! on passe a l'iteration suivante de la boucle
      case ("neu")    ! Neumann homogeneous boundary condition
        uright = vec(i-1)
        vright = vec(nuptotal+i-1)
      case default
      end select
    else
      uright = vec(i+1)
      vright = vec(nuptotal+i+1)
    end if
    select case (ndim)
    case (1)
! no lower and upper neighbours in one dimension -> nothing to do
    case (2)
! lower neighbour
      if(i <= nup1d-1)then
        select case (bc)
        case ("per") ! periodic boundary condition
          ulow = vec(nip*nup1d+i)
          vlow = vec(nuptotal+nip*nup1d+i)
        case ("neu") ! Neumann homogeneous boundary condition
          ulow = vec(i+nup1d)
          vlow = vec(nuptotal+i+nup1d)
        case default
        end select
      else
        ulow = vec(i-nup1d)
        vlow = vec(nuptotal+i-nup1d)
      end if
! upper neighbour
      if(i >= (nup1d-1)*nup1d)then
        select case (bc)
        case ("per") ! periodic boundary condition
          ulow = vec(i-nip*nup1d)
          vlow = vec(nuptotal+i-nip*nup1d)
        case ("neu") ! Neumann homogeneous boundary condition
          ulow = vec(i-nup1d)
          vlow = vec(nuptotal+i-nup1d)
        case default
        end select
      else
        uup = vec(i+nup1d)
        vup = vec(nuptotal+i+nup1d)
      end if
    case default
    end select
! calcul de la derivee en trois etapes
! 1. partie diffusive
    ui = vec (i)
    vi = vec (nuptotal+i)
    select case (ndim)
    case (1)
      dvec (i)          = soodx2   * (uleft + uright - 2.0_dp * ui)
      dvec (nuptotal+i) = 0.0_dp
!!      dvec(nuptotal+i) = e2soodx2 * (vleft + vright - 2.0_dp * vi)
    case (2)
      dvec (i)          = soodx2   * (uleft + uright + ulow + uup - 4.0_dp * ui)
      dvec (nuptotal+i) = e2soodx2 * (vleft + vright + vlow + vup - 4.0_dp * vi)
    case default
    end select
! 2. partie reactive
    ! calcul du terme de convolution
    ix = mod(i, nup1d) ! ces 4 lignes peuvent etre changees
    iy = i/(nup1d)     ! ces 4 lignes peuvent etre changees
    x  = ix*dx         ! ces 4 lignes peuvent etre changees
    y  = iy*dy         ! ces 4 lignes peuvent etre changees
!!    termeconvol = convol (x, y, vec)
!!    h           = termeconvol + rho * 0.5_dp
    fixlib      = fg (vi)
    fix         = fixlib (1)
    lib         = fixlib (2)
    ! calcul du terme de reaction 
    reacplus  =     ui * (rho-vi) * fix
    reacmoins = - eps1 * vi * lib
    termereac = reacplus + reacmoins
! 3. ajout de la partie diffusive et de la partie reactive
    dvec (i)          = dvec (i)          - termereac
    dvec (nuptotal+i) = dvec (nuptotal+i) + termereac
  end do bigloop
end subroutine fcadhe
end module mod_fcadhe

module mod_sys
!===========================================
! MODULE FOR SYS
!===========================================
  implicit none
  private
  public :: sys

contains
subroutine sys ( )
!===========================================
! SYS makes the systematic computations
!===========================================
  use prec
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
  soodx2 = sigma * oodx2
  e2soodx2 = eps2 * soodx2
!!  write(stdout, *) ''
!!  write(stdout, *) soodx2, e2soodx2

end subroutine sys

end module mod_sys

module mod_init_random_seed
!===========================================
! MODULE FOR INIT RANDOM SEED
!===========================================
  implicit none
  private
  public :: init_random_seed

contains
subroutine init_random_seed ( )
!===========================================
! INIT_RANDOM_SEED initializes a pseudo-random number sequence
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

end module mod_init_random_seed

module mod_init
!===========================================
! MODULE FOR INIT
!===========================================
  implicit none
  private
  public :: init

contains
subroutine init (vec)
!===========================================
! INIT computes the initial condition
!===========================================
  use prec
  use dimensions
  use mod_init_random_seed, only: init_random_seed
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), dimension (0:neqn-1), intent(out) :: vec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i, ix, iy
  integer ( int32 ) :: u
  real ( dp ) :: x, r
  real ( dp ) :: xlaouonest, ylaouonest
  real ( dp ) :: gauss

  select case (iinit)
  case (0) ! read initial condition from file
    open(newunit=u, file="fort.10", access="sequential", form ="formatted", status="old")
    readloop: do i=0, nuptotal-1
      read(u, *) x, vec(i)
    end do readloop
    close(u)

  case (1) ! petit plateau pour v autour de (x0, y0), v nulle ailleurs

  case (2) ! petite bosse  pour v autour de (x0, y0), v nulle ailleurs
    xlaouonest = 0.0_dp
    ylaouonest = 0.0_dp
    bigloop2: do i=0, nuptotal-1
      if (xlaouonest > xmax) then
        xlaouonest = xlaouonest - (xmax - xmin)
        ylaouonest = ylaouonest + dy
      end if
      gauss = h2 * exp(-((xlaouonest-x0)**2+(ylaouonest-y0)**2)/(2.0_dp*l))
      vec(i) = 1.0_dp - gauss
      vec(nuptotal+i) = gauss
!!      write(stdout, *) ix, iy, xlaouonest, ylaouonest, gauss
      xlaouonest = xlaouonest + dx
    end do bigloop2

  case (3) ! condition initiale aleatoire
    call init_random_seed ( )
    bigloop3: do i=0, nuptotal-1
      call random_number (r)
      vec(i) = 1-h2*r
      vec(nuptotal+i) = h2*r
    end do bigloop3

  case (4) ! condition initiale aleatoire en une dimension
    call init_random_seed ( )
    select case (ndim)
    case (1)
      xloop: do ix=0, nup1d-1
        call random_number (r)
          vec(ix) = 1-h2*r
          vec(nuptotal+ix) = h2*r
      end do xloop
    case (2)
      xloop_2d: do ix=0, nup1d-1
        call random_number (r)
        yloop: do iy=0, nup1d-1
          vec(iy*nup1d+ix) = 1-h2*r
          vec(nuptotal+iy*nup1d+ix) = h2*r
        end do yloop
      end do xloop_2d
    case default
    end select
  case default
  end select
end subroutine init

end module mod_init

module mod_out
!===========================================
! MODULE FOR OUT
!===========================================
  implicit none
  private
  public :: out

contains
subroutine out (t, vec)
!===========================================
! OUT writes t, x, y, u(x, y), v(x, y) in the file with label u_out
!===========================================
  use prec
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in) :: t
  real ( dp ), dimension(0:neqn-1), intent(in) :: vec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: ix               , iy
  integer ( int32 ) :: i
  real ( dp )       :: x = 0.0_dp       , y = 0.0_dp
  real ( dp )       :: u                , v
  real ( dp )       :: umin = 1000.0_dp , vmin = 1000.0_dp
  real ( dp )       :: umax = 0.0_dp    , vmax = 0.0_dp

  x = 0.0_dp     ! necessaire mais je ne comprends pas encore bien pourquoi
  xloop: do ix = 0, nup1d-1
    i = ix
    y = 0.0_dp   ! necessaire et la je comprends pourquoi
    select case (ndim)
    case (1)
      u = real(vec(i))
      v = real(vec(nuptotal+i))
      write(u_out, "(5(1f10.4, 1x))") t, x, y, u, v
    case (2)
      yloop: do iy = 0, nup1d-1
!!        i = iy*(nup1d)+ix
        u = real(vec(i))
        v = real(vec(nuptotal+i))
        write(u_out, "(5(1f10.4, 1x))") t, x, y, u, v
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
          write(stdout, *) 'Attention'
          write(stdout, *) x, y
          write(stdout, *) 'est un endroit ou u est negatif et vaut'
          write(stdout, *) u
        end if
        if (v < - 1.0e-3_dp) then
          write(stdout, *) 'Attention'
          write(stdout, *) x, y
          write(stdout, *) 'est un endroit ou v est negatif et vaut'
          write(stdout, *) v
        end if
        y = y + dy ! y = (iy+1) * dy
        i = i + nup1d
      end do yloop
    case default
    end select
    x = x + dx ! x = (ix+1) * dx
  end do xloop
  write(stdout, "(1x, a, f10.4, a, f10.4)") 'u est compris entre', umin , ' et ', umax
  write(stdout, "(1x, a, f10.4, a, f10.4)") 'v est compris entre', vmin , ' et ', vmax
end subroutine out

end module mod_out

program cadhe2d
!===========================================
! MAIN solves a system of ODEs resulting from the 2-dimensional space
! discretization of the system of reaction-diffusion equations for the cadherin model
!
! u_t = sigma (u_{xx}+u_{yy}) -r(u, v)
! or
! u_t = sigma u_{xx} -r(u, v)
! v_t = r(u, v)
!
! and periodic boundary conditions
!
! We discretize the space variables with
! x_i=i/(nip+1) for i=1, ..., nip
! or
! x_i=i/(nip+1) and y_i=i/(nip+1) for i=1, ..., nip
! We obtain a system of neqn equations
! The spectral radius of the Jacobian can
! be estimated with the Gershgorin theorem
! Thus we provide an external function RADIUS,
! giving the spectral radius of the Jacobian
! matrix
!===========================================
  use prec
  use dimensions
  use mod_fcadhe, only: fcadhe
  use mod_sys,    only: sys
  use mod_init,   only: init
  use mod_out,    only: out
  implicit none
! Data dictionary: declare local variable types & definitions
  integer ( int32 )                      :: k
  integer ( int32 )                      :: idid
  integer ( int32 ), dimension(12)       :: iwork
  real ( dp )                            :: t0, t1
  real ( dp ), dimension(:), allocatable :: vec
  real ( dp )                            :: rtol, atol
  real ( dp )                            :: h
  real ( dp ), dimension(:), allocatable :: work

! parameters
  call param ( )
! dimensions
  call calc_neqn ( )
! systematic computations
  call sys ( )
  allocate (vec(0:neqn-1))
! If integrating with ROCK2 allocate work of length 4*neqn
  allocate (work(7*neqn)) ! Allocate work of length 7*neqn because the radius is computed externally
! Otherwise allocate work of length 8*neqn

  open(newunit=u_out, file=nom_fichier_out, access="sequential", status ="replace", form="formatted")

! condition initiale
  call init (vec)
! premier intervalle de temps
  t0 = tbeg
  t1 = tbeg + dt
! 0-eme sortie
  write(stdout, *) ''
! faire un header pour le paragraphe dans le fichier
  write(u_out, *) '# Time is', t0
! appeler out pour ecrire dans le fichier de sortie
  call out (t0, vec)

! boucle d'integration
  timeloop: do k=1, ns
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
    write(stdout, *) ''
    write(stdout, *) ''
    write(stdout, "(1x, a, f10.4, a, f10.4)") 'Integration du systeme de ', t0, ' a ', t1
    write(stdout, *) '...'
    call rock4 (neqn, t0, t1, h, vec, fcadhe, atol, rtol, work, iwork, idid)
! k-eme sortie
    write(stdout, *) '...'
    write(stdout, "(1x, a, f10.4)") 'Fait, et maintenant le temps est ', t1
! print statistics
!!    write(stdout, *) '--Solution is tabulated in file ', nom_fichier_out
!!    write(stdout, *) 'The value of dt is', dt
!!    write(stdout, *) 'The value of dx is', dx
!!    write(stdout, *) 'The value of oodx2 is', oodx2
    write(stdout, *) 'The value of IDID is', idid
    write(stdout, *) 'Max estimation of the spectral radius=', iwork(11)
    write(stdout, *) 'Min estimation of the spectral radius=', iwork(12)
    write(stdout, *) 'Max number of stages used=', iwork(10)
    write(stdout, *) 'Number of f evaluations for the spectral radius=', iwork(9)
    write(stdout, "(1x, a, i4, a, i4, a, i4, a, i3)") 'Number of f evaluations=', &
                                             iwork(5), ' steps=', iwork(6), ' accpt=', iwork(7), ' rejct=', iwork(8)
! passer deux lignes dans le fichier IMPORTANT POUR GNUPLOT
    write(u_out, *) ''
    write(u_out, *) ''
! faire un header pour le paragraphe dans le fichier
    write(u_out, *) '# Time is', t1
! appeler out pour ecrire dans le fichier de sortie
    call out (t1, vec)
! prochain intervalle de temps
    t0 = t1
    t1 = t1 + dt
  end do timeloop
  deallocate(vec)
  deallocate(work)
end program cadhe2d

function radius ( )
!===========================================
! RADIUS gives an estimation of the spectral
! radius of the Jacobian matrix of the problem
! This is a bound for the whole interval and
! thus RADIUS is called once
!===========================================
  use prec
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ) :: radius

  radius = 8.0_dp*nuptotal*sigma + 2.0_dp
end function radius
