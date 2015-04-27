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
  integer ( int32 ) :: iinit,ns
end module parametres

module dimensions
!===========================================
! DIMENSIONS
!===========================================
  use prec
  implicit none
  integer ( int32 ), parameter :: ndim = 1
  integer ( int32 ), parameter :: nsd = 49
  integer ( int32 ), parameter :: nssq = nsd**ndim
  integer ( int32 ), parameter :: nsnsm1 = nsd*(nsd-1)
  integer ( int32 ), parameter :: neqn = 2*(nssq)
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
    l3 = (2.0_dp**(0.5_dp)-1.0_dp)*l1
    distance = (x**2+y**2)**(0.5_dp)
  else
! dimension 1
    l3 = 2.0_dp*l1
    distance = abs(x)
  end if
  if (distance > l3) then
    phi = 0
  else if (distance > l1) then
    phi = 1
  else
    phi = -1
  end if
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
  real ( dp ), dimension(neqn), intent(in) :: vec
  real ( dp ) :: convol
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: ix
  real ( dp ) :: ans
  real ( dp ) :: termephi

  ans = real(nsd)      ! necessaire pour faire de nsd un dp

  termephi = phi (x,y)
  convol = 0.5_dp * vec(nssq+1) * termephi
  do ix=1,nsd
    termephi= phi (x-ix/ans,y)
    convol = convol + vec(nssq+ix) * termephi
  end do
  termephi= phi (x-1,y)
  convol = convol + 0.5_dp * vec(nssq+nsd) * termephi

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
  integer ( int32 ), intent(in) :: n
  real ( dp ), intent(in) :: t
  real ( dp ), dimension(n), intent(in) :: vec
  real ( dp ), dimension(n), intent(out) :: dvec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i,ix,iy
  real ( dp ) :: uleft,vleft,uright,vright
  real ( dp ) :: ulow,vlow,uup,vup
  real ( dp ) :: uij,vij
  real ( dp ) :: x,y,termeconvol
  real ( dp ) :: fixation,liberation
  real ( dp ) :: reacplus,reacmoins,reac
  real ( dp ) :: ans

  ans = real(nsd)      ! necessaire pour faire de nsd un dp
! big loop
  do i=1,nssq
! left neighbour
    if(mod(i,nsd) == 1)then
! periodic boundary condition
!!      uleft=vec(i+nsd-1)
!!      vleft=vec(nssq+i+nsd-1)
! Neumann homogeneous boundary condition
      uleft=vec(i+1)
      vleft=vec(nssq+i+1)
    else
     uleft=vec(i-1)
     vleft=vec(nssq+i-1)
    end if
! right neighbour
    if(mod(i,nsd) == 0)then
! periodic boundary condition
!!      uright=vec(i-nsd+1)
!!      vright=vec(nssq+i-nsd+1)
! Neumann homogeneous boundary condition
      uright=vec(i-1)
      vright=vec(nssq+i-1)
    else
      uright=vec(i+1)
      vright=vec(nssq+i+1)
    end if
! partie dependante de la dimension
    if(ndim == 2)then
  ! lower neighbour
      if(i <= nsd)then
  ! periodic boundary condition
!!        ulow=vec(i+nsnsm1)
!!        vlow=vec(nssq+i+nsnsm1)
  ! Neumann homogeneous boundary condition
        ulow=vec(i+nsd)
        vlow=vec(nssq+i+nsd)
      else
        ulow=vec(i-nsd)
        vlow=vec(nssq+i-nsd)
      end if
  ! upper neighbour
      if(i > nsnsm1)then
  ! periodic boundary condition
!!        uup=vec(i-nsnsm1)
!!        vup=vec(nssq+i-nsnsm1)
  ! Neumann homogeneous boundary condition
        uup=vec(i-nsd)
        vup=vec(nssq+i-nsd)
      else
        uup=vec(i+nsd)
        vup=vec(nssq+i+nsd)
      end if
    else
! dimension 1
    end if
! the derivative
    uij=vec(i)
    vij=vec(i+nssq)
! partie dependante de la dimension
    if(ndim == 2)then
      dvec(i)=sigma*nssq*(uleft+uright+ulow+uup-4.0_dp*uij)
      dvec(nssq+i)=0.0_dp
    else
! dimension 1
      dvec(i)=sigma*nssq*(uleft+uright-2.0_dp*uij)
      dvec(nssq+i)=0.0_dp
    end if
! appel a convol
    iy=(i-1)/nsd+1
    ix=i-(iy-1)*nsd
    x=ix*dx
    y=iy*dy
    termeconvol=convol(x,y,vec)
! reaction 
    fixation=max(0.0_dp,1.0_dp-termeconvol)
    liberation=eps*min(1.0_dp,termeconvol)
    reacplus=fixation*uij*(rho-vij)
    reacmoins=-liberation*vij
    reac=reacplus+reacmoins

    dvec(i)=dvec(i) - reac
    dvec(i+nssq)=dvec(i+nssq) + reac
!!if (i == 1) then
!!  write(stdout,*) 'du est',dvec(i)
!!  write(stdout,*) 'dv est',dvec(i+nssq)
!!  write(stdout,*) 'Le temps est',t0
!!else
!!end if
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
! x_i=i/(nsd+1) for i=1,...,nsd
! or
! x_i=i/(nsd+1) and y_i=i/(nsd+1) for i=1,...,nsd
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
  real ( dp ), dimension(neqn) :: vec
  real ( dp ) :: rtol,atol
  real ( dp ) :: h
! If integrating with ROCK2 define work of length 4neqn
  real ( dp ), dimension(7*neqn) :: work ! Work is of length 7neqn because the radius is computed externally (otherwise should be 8neqn)

! parameters
  call param ( )
! systematic computations
  call sys ( )
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
    write(*,"(1x,a,f5.2,a,f5.2)") 'Integration du systeme de ',t0,' a ',t1
    call rock4 (neqn,t0,t1,h,vec,fcadhe,atol,rtol,work,iwork,idid)
! k-eme sortie
    write(stdout,*) '...'
    write(*,"(1x,a,f5.2)") 'Fait, et maintenant le temps est ',t1
! print statistics
    write(stdout,*) ''
!!    write(stdout,*) '--Solution is tabulated in file ',nom_fichier_out
    write(stdout,*) 'The value of IDID is',idid
    write(stdout,*) 'Max estimation of the spectral radius=',iwork(11)
    write(stdout,*) 'Min estimation of the spectral radius=',iwork(12)
    write(stdout,*) 'Max number of stages used=',iwork(10)
    write(stdout,*) 'Number of f evaluations for the spectral radius=',iwork(9)
    write(*,"(1x,a,i4,a,i4,a,i4,a,i3)") 'Number of f evaluations=',iwork(5),' steps=',iwork(6),' accpt=',iwork(7),' rejct=',iwork(8)
! passer une ligne dans le fichier IMPORTANT POUR GNUPLOT
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

  radius = 8.0_dp*nssq*sigma + 2.0_dp
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
  real ( dp ), dimension (neqn), intent(out) :: vec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: i,ix,iy
  integer ( int32 ) :: u
  real ( dp ) :: x,r
  real ( dp ) :: ans
  real ( dp ) :: xlaouonest,ylaouonest

  ans = real(nsd)      ! necessaire pour faire de nsd un dp

  if (iinit == 0) then
! read initial condition from file
!!    open(newunit=u,file="fort.10",access="sequential",form ="formatted",status="old")
    open(newunit=u,file="fort.10",status="old")
    do i=1,nssq
      read(u,*) x,vec(i)
    end do
    close(u)
  else if (iinit == 1) then
! petit plateau pour g autour de x0, g nulle ailleurs

  else if (iinit == 2) then
! petite bosse pour g autour de x0, g nulle ailleurs

  else if (iinit == 3) then
! condition initiale aleatoire
    call init_random_seed ( )
! big loop
    do i=1,nssq
      call random_number (r)
      vec(i) = 1-h2*r
      vec(nssq+i) = h2*r
    end do
  else if (iinit == 4) then
! gaussienne
! big loop
    do i=1,nssq
      iy=(i-1)/nsd+1
      ix=i-(iy-1)*nsd
      xlaouonest = (ix-1)/(ans-1)
      ylaouonest = (iy-1)/(ans-1)
      vec(i) = exp(-((xlaouonest-x0)**2+(ylaouonest-y0)**2)/(2.0_dp*l))
      vec(nssq+i) = 0
    end do
  else if (iinit == 5) then
! condition initiale aleatoire en une dimension
    call init_random_seed ( )
    do ix=1,nsd
      call random_number (r)
      do iy=1,nsd
        vec((iy-1)*nsd+ix) = 1-h2*r
        vec(nssq+(iy-1)*nsd+ix) = h2*r
      end do
    end do
  end if
end subroutine init

subroutine out (t,vec)
!===========================================
! 7.OUT writes t,x,y,u(x,y),v(x,y) in the file with label u
!===========================================
  use prec
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real ( dp ), intent(in) :: t
  real ( dp ), dimension(neqn), intent(in) :: vec
! Data dictionary: declare local variable types & definitions
  integer ( int32 ) :: ix,iy
  real ( dp ) :: x,y,u,v
  real ( dp ) :: umin,vmin = 1000
  real ( dp ) :: umax,vmax = 0

  do ix=1,nsd
    x = ix*dx
! partie dependante de la dimension
    if (ndim == 2) then
      do iy=1,nsd
        y = iy*dy
        u = real(vec((iy-1)*nsd+ix))
        v = real(vec(nssq+(iy-1)*nsd+ix))
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
        if (u < 0) then
          write(stdout,*) 'Attention'
          write(stdout,*) x,y
          write(stdout,*) 'est un endroit ou u est negatif et vaut'
          write(stdout,*) u
        end if
        if (v < 0) then
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
      v = real(vec(nssq+ix))
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

  write(*,"(1x,a)",advance='no') 'nom du run (6 caracteres) ? '
  read (*,"(a6)") nom_fichier
  nom_fichier_param =  nom_fichier//'.param'
  nom_fichier_out =  nom_fichier//'.out'
  open(newunit=u,file=nom_fichier_param,status="old",form="formatted",action="read",iostat=ierror)
  if (ierror == 0) then
! Open was ok, read values
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror)
    read(u,*,iostat=ierror) xmin,xmax,ymin,ymax
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
  dx = (xmax-xmin)/(nsd+1)
  dy = (ymax-ymin)/(nsd+1)

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
