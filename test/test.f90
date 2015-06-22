program test
  integer, parameter :: num = 10
  real, dimension(num,num) :: A
  integer :: i, j, ij

do i=1,num
  if (i .eq. 4) cycle
  do j=1,num
    write(*,*) i,'.',j
  enddo
enddo

end program test
