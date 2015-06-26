program test
  integer, parameter :: npart = 100
  integer :: i, j, ij

ij=0
do i=1,npart-1
  if (i.ge.4 .or. i.eq.30) then
    do j=i+1,npart
      ij=ij+1
    enddo
  endif
  if (i.ge.4 .or. i.eq.30) cycle
  do j=i+1,npart
    if (j.eq.2 .or. j.eq.13) ij=ij+1
    if (j.eq.2 .or. j.eq.13) cycle
    ij=ij+1
  enddo
enddo
write(*,*) ((npart-1)*npart)/2, ij

end program test
