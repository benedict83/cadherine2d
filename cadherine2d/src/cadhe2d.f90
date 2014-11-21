module parametres
!===========================================
! 1.PARAMETRES
!===========================================
implicit none
save
  character(6) :: nom_fichier
  character(12) :: nom_fichier_param
  character(10) :: nom_fichier_out
  integer, parameter :: double = 8 ! Compiler dependent value
  real ( kind = double ) :: xmin,xmax,tbeg,tend
  real ( kind = double ) :: rho,eps
  real ( kind = double ) :: a0,a1
  real ( kind = double ) :: sigma
  real ( kind = double ) :: h1,h2,x0,y0,l
  real ( kind = double ) :: dt,dx
  real ( kind = double ) :: l1,l2
  integer :: iinit,ngx,npri,ns
end module parametres

module dimensions
!===========================================
! 2.DIMENSIONS
!===========================================
implicit none
save
  integer, parameter :: nsd = 32
  integer, parameter :: nssq = nsd*nsd
  integer, parameter :: nsnsm1 = nsd*(nsd-1)
  integer, parameter :: neqn = nsd*nsd*2
end module dimensions

program main
!===========================================
! 3.MAIN
!  solves a system of ODEs resulting from the 2-dimensional space
!  discretization of the system of reaction-diffusion equations for the cadherin model

!  u_t = sigma (u_{xx}+u_{yy}) -r(u,v)
!  v_t = r(u,v)

!  and periodic boundary conditions

!  u(x+1,y,t)=u(x,y,t),  v(x,y+1,t)=v(x,y,t).

!  We discretize the space variables with
!  x_i=i/(N+1), y_i=i/(N+1) for i=0,1,...,N,
!  with N=128. We obtain a system of 32768
!  equations. The spectral radius of the Jacobian
!  can be estimated with the Gershgorin theorem
!  (13200 is an estimation for it). Thus we
!  provide an external function RADIUS, giving
!  the spectral radius of the Jacobian matrix.

!===========================================
  use parametres
  use dimensions
  implicit none
  integer :: k
  integer :: idid
  integer, dimension(12) :: iwork
  real ( kind = double ) :: t0,t1
  real ( kind = double ), dimension(neqn) :: vec
  real ( kind = double ) :: rtol,atol
  real ( kind = double ) :: h
  real ( kind = double ), dimension(7*neqn) :: work
!--------------------------------------------------------
!    Work is of length 7*neqn because the spectral radius
!    is computed externally (otherwise should be 8*neqn).
!    If integrating with ROCK2 define work(4*neqn).
!--------------------------------------------------------
  external fcadhe
! parameters
  call param
! systematic computations
  call sys(neqn)

! condition initiale
  call init(nssq,vec)
! premier intervalle de temps
  t0 = tbeg
  t1 = tbeg + dt
! 0-ème sortie
  call out(t0,neqn,vec)
! boucle d'integration
  do k=1,ns
! ----- required tolerance -----
    rtol=0.1d0**2
    atol=rtol
! ----- initial step size -----
    h=1.0d-4
!--------------------------------------------------------
!   Initialize iwork:
!    iwork(1)=1  RADIUS returns an upper bound for the spectral radius.
!    iwork(2)=1  The Jacobian is constant (RADIUS is called once).
!    iwork(3)=0  Return and solution at t1.
!    iwork(4)=0  Atol and rtol are scalars.
!--------------------------------------------------------
    iwork(1)=1
    iwork(2)=1
    iwork(3)=0
    iwork(4)=0
    write(*,*) ''
    write(*,"(1x,a,f5.2,a,f5.2)") 'Integration du systeme de ',t0,' a ',t1
! ----- to integrate with rock2.f90
!!    call rock2(neqn,t0,t1,h,vec,fcadhe,atol,rtol,work,iwork,idid)
    call rock4(neqn,t0,t1,h,vec,fcadhe,atol,rtol,work,iwork,idid)
! k-ème sortie
    write(*,*) '...'
    write(*,"(1x,a,f5.2)") 'Fait, et maintenant le temps est ',t1
! ----- print statistics -----
    write(*,*) ''
!!    write(*,*) '--Solution is tabulated in file ',nom_fichier_out
    write(*,*) 'The value of IDID is',idid
    write(*,*) 'Max estimation of the spectral radius=',iwork(11)
    write(*,*) 'Min estimation of the spectral radius=',iwork(12)
    write(*,*) 'Max number of stages used=',iwork(10)
    write(*,*) 'Number of f eval. for the spectr. radius=',iwork(9)
    write(*,"(1x,a,i4,a,i4,a,i4,a,i3)") 'Number of f evaluations=',iwork(5),' steps=',iwork(6),' accpt=',iwork(7),' rejct=',iwork(8)
! appeler out pour écrire dans le fichier de sortie
    call out(t1,neqn,vec)
! prochain intervalle de temps
    t0 = t1
    t1 = t1 + dt
! fin de boucle
  end do
end program main

function radius ( )
!===========================================
! 4.RADIUS gives an estimation of the spectral
!   radius of the Jacobian matrix of the problem.
!   This is a bound for the whole interval and
!   thus RADIUS is called once.
!===========================================
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare variable types & definitions
  real ( kind = double ) :: radius

  radius = 8.0d0*nssq*sigma + 2.d0
end function radius

subroutine fcadhe (n,t,vec,dvec)
!===========================================
! 5.FCADHE compute the value of
!   dvec(x,y)=(du(x,y),dv(x,y)).
!===========================================
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  integer, intent(in) :: n
  real ( kind = double ), intent(in) :: t
  real ( kind = double ), dimension(n), intent(in) :: vec
  real ( kind = double ), dimension(n), intent(out) :: dvec
! Data dictionary: declare local variable types & definitions
  integer :: i,ix,iy,j,k
  integer :: denom
  real ( kind = double ) :: uleft,vleft,uright,vright,ulow,vlow,uup,vup
  real ( kind = double ) :: uij,vij
  real ( kind = double ) :: xlaouonest,ylaouonest
  real ( kind = double ) :: xlaouonconvole,ylaouonconvole
  real ( kind = double ) :: distanceenx,distanceeny
  real ( kind = double ) :: somme
  real ( kind = double ) :: convol,fixation,liberation
  real ( kind = double ) :: reacplus,reacmoins,reac
  real ( kind = double ) :: ans

  ans=nsd      ! ------ necessaire ?
! ----- big loop -----
  do i=1,nssq
! ----- left neighbour -----
    if(mod(i,nsd) == 1)then
! periodic boundary condition
      uleft=vec(i+nsd-1)
      vleft=vec(nssq+i+nsd-1)
! Neumann homogeneous boundary condition
!!      uleft=vec(i+1)
!!      vleft=vec(nssq+i+1)
    else
     uleft=vec(i-1)
     vleft=vec(nssq+i-1)
    end if
! ----- right neighbour -----
    if(mod(i,nsd) == 0)then
! periodic boundary condition
      uright=vec(i-nsd+1)
      vright=vec(nssq+i-nsd+1)
! Neumann homogeneous boundary condition
!!      uright=vec(i-1)
!!      vright=vec(nssq+i-1)
    else
      uright=vec(i+1)
      vright=vec(nssq+i+1)
    end if
! ----- lower neighbour -----
    if(i <= nsd)then
! periodic boundary condition
      ulow=vec(i+nsnsm1)
      vlow=vec(nssq+i+nsnsm1)
! Neumann homogeneous boundary condition
!!      ulow=vec(i+nsd)
!!      vlow=vec(nssq+i+nsd)
    else
      ulow=vec(i-nsd)
      vlow=vec(nssq+i-nsd)
    end if
! ----- upper neighbour -----
    if(i > nsnsm1)then
! periodic boundary condition
      uup=vec(i-nsnsm1)
      vup=vec(nssq+i-nsnsm1)
! Neumann homogeneous boundary condition
!!      uup=vec(i-nsd)
!!      vup=vec(nssq+i-nsd)
    else
      uup=vec(i+nsd)
      vup=vec(nssq+i+nsd)
    end if
! ----- the derivative -----
    uij=vec(i)
    vij=vec(i+nssq)
! ----- partie diffusive -----
    dvec(i)=sigma*nssq*(uleft+uright+ulow+uup-4.d0*uij)
    dvec(i+nssq)=0.d0
! ----- convolution -----
    somme=0.d0
    denom=0
    iy=(i-1)/nsd+1
    ix=i-(iy-1)*nsd
    xlaouonest=(ix-1)/(ans-1)
    ylaouonest=(iy-1)/(ans-1)
    do j=1,nsd
      xlaouonconvole=(j-1)/(ans-1)
!! en cas de periodic boundary conditions il faut s'assurer de bien prendre la bonne distance
      distanceenx=min(abs(xlaouonest-xlaouonconvole),abs(xlaouonest-xlaouonconvole-1),abs(xlaouonest-xlaouonconvole+1))
      do k=1,nsd
        ylaouonconvole=(k-1)/(ans-1)
!! en cas de periodic boundary conditions il faut s'assurer de bien prendre la bonne distance
        distanceeny=min(abs(ylaouonest-ylaouonconvole),abs(ylaouonest-ylaouonconvole-1),abs(ylaouonest-ylaouonconvole+1))
        if (((distanceenx**2 +distanceeny**2)**(0.5d0)) <= l1) then
          somme=somme-vec((k-1)*nsd+j+nssq)
          denom=denom+1
        else if (((distanceenx**2 +distanceeny**2)**(0.5d0)) <= ((2**(0.5d0)-1)*l1)) then
          somme=somme+vec((k-1)*nsd+j+nssq)
          denom=denom+1
        end if
      end do
    end do
    convol=(somme+rho*denom)/(rho*denom)
! ----- reaction -----
    fixation=max(0.d0,1.d0-convol)
    liberation=eps*min(1.d0,convol)
    reacplus=fixation*uij*(rho-vij)
    reacmoins=-liberation*vij
    reac=reacplus+reacmoins
! derivee uij_t
    dvec(i)=dvec(i) - reac
! derivee vij_t
    dvec(i+nssq)=dvec(i+nssq) + reac
!!if (i==1) then
!!  write(*,*) 'du est',dvec(i)
!!  write(*,*) 'dv est',dvec(i+nssq)
!!  write(*,*) 'Le temps est',t0
!!else
!!end if
  end do
end subroutine fcadhe

subroutine init (m,vec)
!===========================================
! 6.INIT computes the initial condition
!===========================================
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  integer, intent(in) :: m
  real ( kind = double ), dimension (2*m), intent(out) :: vec
! Data dictionary: declare local variable types & definitions
  integer :: i,ix,iy
  real ( kind = double ) :: x,r
  real ( kind = double ) :: ans
  real ( kind = double ) :: xlaouonest,ylaouonest

  ans=nsd      ! ------ necessaire ?

  if (iinit == 0) then
! read initial condition from file
    open(10,file='fort.10',access='sequential',form ='formatted')
    do i=1,m
      read(10,*) x,vec(i)
    end do
    close(10)
  else if (iinit == 1) then
! petit plateau pour g autour de x0, g nulle ailleurs

 else if (iinit == 2) then
! petite bosse pour g autour de x0, g nulle ailleurs

  else if (iinit == 3) then
! condition initiale aleatoire
    call init_random_seed ( )
    do i=1,m
      call random_number (r)
      vec(i) = 1-h2*r
      vec(m+i) = h2*r
    end do
  else if (iinit == 4) then
! gaussienne
    do i=1,m         ! ----- big loop -----
      iy=(i-1)/nsd+1
      ix=i-(iy-1)*nsd
      xlaouonest = (ix-1)/(ans-1)
      ylaouonest = (iy-1)/(ans-1)
      vec(i) = exp(-((xlaouonest-x0)**2+(ylaouonest-y0)**2)/(2.d0*l))
      vec(m+i) = 0
    end do
  else if (iinit == 5) then
! condition initiale aleatoire en une dimension
    call init_random_seed ( )
    do ix=1,nsd
!!      xlaouonest = (ix-1)/(ans-1)
      call random_number (r)
      do iy=1,nsd
        vec((iy-1)*nsd+ix) = 1-h2*r
        vec(m+(iy-1)*nsd+ix) = h2*r
      end do
    end do
  end if
end subroutine init

subroutine out (t,m,vec)
!===========================================
! 7.OUT writes t,x,y,u(x,y),v(x,y) in the file with label 16
!===========================================
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  integer, intent(in) :: m
  real ( kind = double ), intent(in) :: t
  real ( kind = double ), dimension(m), intent(in) :: vec
! Data dictionary: declare local variable types & definitions
  integer :: i,j
  real :: x,y,u,v
  real :: umin,vmin = 1000
  real :: umax,vmax = 0

!!  if (npri > 0) then
  if (t == 0) then
! passer une ligne dans la sortie standard
    write(*,*) ''
  endif
  if (t > 0) then
! passer une ligne dans le fichier
    write(16,*) ''
  endif
  write(16,*) '# Time is',t
  do i=1,nsd
    x = (i-1)/(real(nsd-1))
    do j=1,nsd
      y = (j-1)/(real(nsd-1))
      u = real(vec((j-1)*nsd+i))
      v = real(vec((j-1)*nsd+i+nssq))
      write(16,"(5(f9.4,x))") t,x,y,u,v
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
        write(*,*) 'Attention'
        write(*,*) x,y
        write(*,*) 'est un endroit ou u est negatif et vaut'
        write(*,*) u
      end if
      if (v < 0) then
        write(*,*) 'Attention'
        write(*,*) x,y
        write(*,*) 'est un endroit ou v est negatif et vaut'
        write(*,*) v
      end if
    end do
  end do
  write(*,*) 'u est compris entre'
  write(*,*) umin,umax
  write(*,*) 'v est compris entre'
  write(*,*) vmin,vmax
!!  end if
end subroutine out

subroutine param
!===========================================
! 8.PARAM reads parameters from a data file
!===========================================
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare local variable types & definitions
  integer :: ierror

  write(*,"(1x,a)",advance='no') 'nom du run (6 caracteres) ? '
  read (*,"(a6)") nom_fichier
  nom_fichier_param =  nom_fichier//'.param'
  nom_fichier_out =  nom_fichier//'.out'
  open(33,file=nom_fichier_param,status='old',form='formatted',action='read',iostat=ierror)
  if (ierror == 0) then
! Open was ok. Read values.
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror) xmin,xmax,ngx
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror) tbeg,tend,ns
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror) rho,l1,l2,eps
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror) sigma
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror)
    read(33,*,iostat=ierror) iinit,h1,h2,x0,y0,l
    close(33)
  end if
end subroutine param

subroutine sys (n)
!===========================================
! 9.SYS makes the systematic computations
!===========================================
  use parametres
  use dimensions
  implicit none
! Data dictionary: declare calling parameter types & definitions
  integer, intent(in) :: n
! Data dictionary: declare local variable types & definitions
  real ( kind = double ) :: tiny = 1.d-10

! compute time step
!!  dt = (tend-tbeg+tiny)/ns
  dt = (tend-tbeg)/ns
! compute space step
  dx = (xmax-xmin)/n
! compute n print
!!  npri = (n-1)/(ngx-1)
  npri = 1
! open nom_fichier.out if npri>0
!!  if (npri > 0) then
    open(16,file=nom_fichier_out,access='sequential',status ='replace',form='formatted')
!!  end if
end subroutine sys

subroutine init_random_seed ( )
!===========================================
! 10.INIT_RANDOM_SEED initializes a pseudo-random number sequence
!===========================================
  implicit none
! Data dictionary: declare local variable types & definitions
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  deallocate(seed)
end subroutine init_random_seed
