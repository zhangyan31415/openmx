program stack_reproducer
 implicit none
 character(32) :: arg
 integer :: n
 call get_command_argument(1,arg)
 read(arg,*)n
 call workspace(n)
contains
 subroutine workspace(k)
  integer,intent(in)::k
  real(8)::a(k)
  call random_number(a)
  if(sum(a*a)<=0.d0)stop 1
  print *, 'automatic-array workspace PASS',k
 end subroutine
end program
