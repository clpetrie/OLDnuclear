module kshellold
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

   integer(kind=i4), private, parameter :: nshmax = 23
   integer(kind=i4), private, parameter :: nkmax = 246
   integer(kind=i4), private :: nshnum(nshmax)= &
      (/1,6,12,8,6,24,24,12,30,24,24,8,24,48,6,48,12,24,24,48,24,24,30/)
   integer(kind=i4), private :: ik(3,nkmax) = reshape((/ &
      0,0,0, 0,0,1, 0,1,0, 1,0,0, 0,1,1, 1,0,1, 1,1,0, 0,-1,1, -1,0,1, &
     -1,1,0, 1,1,1, 1,-1,1, -1,-1,1, -1,1,1, 0,0,2, 0,2,0, 2,0,0, 0,1,2, &
      1,0,2, 1,2,0, 0,2,1, 2,0,1, 2,1,0, 0,-1,2, -1,0,2, -1,2,0, 0,-2,1, &
     -2,0,1, -2,1,0, 1,1,2, 1,2,1, 2,1,1, 1,-1,2, 1,-2,1, 2,-1,1, -1,1,2, &
     -1,2,1, -2,1,1, -1,-1,2, -1,-2,1, -2,-1,1, 0,2,2, 2,0,2, 2,2,0, 0,-2,2, &
     -2,0,2, -2,2,0, 1,2,2, 2,1,2, 2,2,1, 1,-2,2, 2,-1,2, 2,-2,1, -1,2,2, &
     -2,1,2, -2,2,1, -1,-2,2, -2,-1,2, -2,-2,1, 0,0,3, 0,3,0, 3,0,0, 0,1,3, &
      1,0,3, 1,3,0, 0,3,1, 3,0,1, 3,1,0, 0,-1,3, -1,0,3, -1,3,0, 0,-3,1, &
     -3,0,1, -3,1,0, 1,1,3, 1,3,1, 3,1,1, 1,-1,3, 1,-3,1, 3,-1,1, -1,1,3, &
     -1,3,1, -3,1,1, -1,-1,3, -1,-3,1, -3,-1,1, 2,2,2, 2,-2,2, -2,-2,2, &
     -2,2,2, 0,2,3, 2,0,3, 2,3,0, 0,3,2, 3,0,2, 3,2,0, 0,-2,3, -2,0,3, &
     -2,3,0, 0,-3,2, -3,0,2, -3,2,0, 1,2,3, 2,3,1, 3,1,2, 2,1,3, 3,2,1, &
      1,3,2, 1,-2,3, 2,-3,1, 3,-1,2, 2,-1,3, 3,-2,1, 1,-3,2, -1,2,3, -2,3,1, &
     -3,1,2, -2,1,3, -3,2,1, -1,3,2, -1,-2,3, -2,-3,1, -3,-1,2, -2,-1,3, &
     -3,-2,1, -1,-3,2, 0,0,4, 0,4,0, 4,0,0, 2,2,3, 2,3,2, 3,2,2, 2,-2,3, &
      2,-3,2, 3,-2,2, -2,2,3, -2,3,2, -3,2,2, -2,-2,3, -2,-3,2, -3,-2,2, &
      0,1,4, 1,0,4, 1,4,0, 0,4,1, 4,0,1, 4,1,0, 0,-1,4, -1,0,4, -1,4,0, &
      0,-4,1, -4,0,1, -4,1,0, 0,3,3, 3,0,3, 3,3,0, 0,-3,3, -3,0,3, -3,3,0, &
      1,3,3, 3,1,3, 3,3,1, 1,-3,3, 3,-1,3, 3,-3,1, -1,3,3, -3,1,3, -3,3,1, &
     -1,-3,3, -3,-1,3, -3,-3,1, 0,2,4, 2,0,4, 2,4,0, 0,4,2, 4,0,2, 4,2,0, &
      0,-2,4, -2,0,4, -2,4,0, 0,-4,2, -4,0,2, -4,2,0, 1,2,4, 2,4,1, 4,1,2, &
      2,1,4, 4,2,1, 1,4,2, 1,-2,4, 2,-4,1, 4,-1,2, 2,-1,4, 4,-2,1, 1,-4,2, &
     -1,2,4, -2,4,1, -4,1,2, -2,1,4, -4,2,1, -1,4,2, -1,-2,4, -2,-4,1, &
     -4,-1,2, -2,-1,4, -4,-2,1, -1,-4,2, 2,3,3, 3,2,3, 3,3,2, 2,-3,3, 3,-2,3, &
      3,-3,2, -2,3,3, -3,2,3, -3,3,2, -2,-3,3, -3,-2,3, -3,-3,2, 2,2,4, &
      2,4,2, 4,2,2, 2,-2,4, 2,-4,2, 4,-2,2, -2,2,4, -2,4,2, -4,2,2, &
      -2,-2,4, -2,-4,2, -4,-2,2, 0,3,4, 3,0,4, 3,4,0, 0,4,3, 4,0,3, 4,3,0, &
      0,-3,4, -3,0,4, -3,4,0, 0,-4,3, -4,0,3, -4,3,0, 0,0,5, 0,5,0, 5,0,0 &
      /),(/3,nkmax/))
!
! nshnum = number of states in corresponding k shell
! ik = k vectors for the k states. The rightmost value is always positive.
!
   contains
      subroutine setupk(el,nsh,vkin,ak,ak2,vk,nk)
!
! given el the side of the cell, and the number of shell amplitudes, nsh
! along with their amplitudes vkin, return a the k vectors in ak, their
! squares in ak2, their amplitudes in vk, and the total number in nk.
!
      integer(kind=i4) :: nsh,nk
      real(kind=r8) :: vkin(nsh),el
      real(kind=r8), pointer :: ak(:,:),ak2(:),vk(:)
      real(kind=r8) :: pi,piel
      integer(kind=i4) :: i,j,k
      if (nsh.gt.nshmax) then
         write (6,'(1x,''nsh too large'',2i5)') nsh,nshmax
         stop
         endif
!     if((associated(ak).or.associated(ak2)).or.associated(vk)) then
!        write (6,'(1x,''ak, ak2, or vk already associated'')') 
!        stop
!        endif
      pi=4.0_r8*atan(1.0_r8)
      piel=2.0_r8*pi/el
      nk=1
      do i=2,nsh
         nk=nk+nshnum(i)/2
      enddo
      allocate(ak(3,nk))
      allocate(vk(nk))
      allocate(ak2(nk))
      k=0
      do i=1,nsh
         do j=1,max(1,nshnum(i)/2)
            k=k+1
            ak(1,k)=piel*ik(1,k)
            ak(2,k)=piel*ik(2,k)
            ak(3,k)=piel*ik(3,k)
            ak2(k)=ak(1,k)**2+ak(2,k)**2+ak(3,k)**2
            vk(k)=vkin(i)
         enddo
      enddo
      end subroutine setupk

      subroutine shell(nbin,nshell,filled)
!
! given nbin, return the number of shells, nshell so that there are
! at least nbin states. If the shells are filled, filled = .true.
!
      integer(kind=i4) :: nbin,nshell,nsum,i
      logical :: filled
      nsum=0
      do i=1,nshmax
         nsum=nsum+nshnum(i)
         nshell=i
         filled=nsum.eq.nbin
         if (nsum.ge.nbin) return
      enddo
      stop
      end subroutine shell
end module kshellold
