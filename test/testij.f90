program testij
implicit none

integer, parameter :: npart=4
integer :: i,j,ij,kl,calij

ij=0
do i=1,npart-1
   do j=i+1,npart
      ij=ij+1
      kl=calij(i,j)
      write(*,*) i,j,ij, kl
   enddo
enddo

end program

function calij(i,j)
integer, parameter :: npart=4
integer :: i,j,calij,n,mysum

mysum=0
do n=1,i-1
   mysum=mysum-n
enddo
calij=(j-i)+(i-1)*npart+mysum

end function
