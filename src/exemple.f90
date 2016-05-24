module integral
  implicit none
  private
  public :: f
contains
function f( x, params ) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
! Data dictionary: declare calling parameter types & definitions
  real( kind = c_double )        :: f
  real( kind = c_double ), value :: x
  type( c_ptr ), value           :: params
 
  f = sin( x ) / x
end function f 
end module integral

program exemple
  use fgsl
  use integral
  use, intrinsic :: iso_c_binding
  implicit none
! Data dictionary: declare calling parameter types & definitions
! Data dictionary: declare local variable types & definitions
  real( kind = fgsl_double )   :: result, error
  integer( kind = fgsl_size_t) :: neval
  integer( kind = fgsl_int)    :: i
  type( fgsl_function )        :: func

  func = fgsl_function_init( f, c_null_ptr )
     
  i = fgsl_integration_qng ( func,              &
                         0.0_fgsl_double,       &
                         1.0_fgsl_double,       &
                         1e-9_fgsl_double,      & 
                         1e-9_fgsl_double,      &   
                         result, error, neval ) 
  if (i.ne.0) then
    write( *, * ) "There was a problem with integration: code ", i 
  else
    write( *, "(1x,a,f10.7,a,f10.7,a,i3,a)" ) "Result ", result, " +/-", error, " from", neval, " evaluations"
    write( *, * ) "Result =", result
  end if
! Ausgabe:
!   Ergebnis = 0.946083070367183
end program exemple
