program test
  integer, parameter :: num = 4
  real, dimension(num,num) :: A
  integer :: i, j, ij

  ij=0
  do i=1,num
    do j=1,num
      ij=ij+1
      write(*,*) 'i=', i, ' j=',j, ' ij=',ij
    end do
  end do
end program test
