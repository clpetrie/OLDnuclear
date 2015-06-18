module euclidean
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   integer(kind=i4), private, save :: npart,nspop,nqtab
   real(kind=r8), private, save :: dq
contains
   subroutine initeuc(neucopin,nspopin,npartin)
   use mympi
   integer(kind=i4) :: neucopin,nspopin,npartin
   real(kind=r8) :: qmax
   qmax=2.5_r8
   nqtab=6
   dq=qmax/nqtab
   npart=npartin
   nspop=1
   neucopin=npart
   nspopin=nspop
   if (myrank().eq.0) then
      write(6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write(6,'(''nqtab ='',t40,i10)') nqtab
      write(6,'(''qmax = '',t40,f10.5)') qmax
      write(6,'(''dq = '',t40,f10.5)') dq
      write(6,'(''number of spin operators ='',t40,i10)') nspop
      write(6,'(''total number of operators ='',t40,i10)') npart
      write(6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
   endif
   end subroutine initeuc
 
   subroutine setupeuc(w)
   use stack
   type (walker) :: w
   w%x0=w%x
   w%sp0=w%sp
print*,'sp=',w%sp
print*,'!!!!'
   call ospin(w%sp0,w%osp)
   end subroutine setupeuc

   subroutine ospin(spin,osp)
   complex(kind=r8) :: spin(:,:),osp(:,:,:,:)
   integer(kind=i4) :: i
! setup each spin stack as the original
! and then apply tau_z to the spinner of particle i
   do i=1,npart
      osp(i,:,1,:)=spin(:,:)
      call tauz(spin(:,i),osp(i,:,1,i))
   enddo
print*,'op1=',osp(1,:,1,:)
   end subroutine ospin

   subroutine overlap(w)
   use stack
   use wavefunction
   type (walker) :: w
   complex(kind=r8) :: psi,cop,csum1,csum2
   complex(kind=r8) :: spsave(4,npart),spsave1(4),tzsp(4)
   integer(kind=i4) :: i,j
   spsave=w%sp
   call hpsi(w,.false.)
print*,'psi=',w%psi
   psi=w%psi
   csum1=czero
   csum2=czero
   do i=1,npart
      w%sp=w%osp(i,:,1,:)
      call hpsi(w,.false.)
print*,'i, psi1=',i,w%psi
      csum1=csum1+w%psi/psi
      w%sp=spsave
      call tauz(w%sp(:,i),tzsp)
      w%sp(:,i)=tzsp
      call hpsi(w,.false.)
print*,'i, psi2=',i,w%psi
      csum1=csum1+w%psi/psi ! the sum of this is zero
      w%sp=spsave
      do j=1,npart
         spsave1=w%osp(i,:,1,j)
         call tauz(spsave1,tzsp)
         w%osp(i,:,1,j)=tzsp
         w%sp=w%osp(i,:,1,:)
         call hpsi(w,.false.)
print*,'j, psi=',j,w%psi
         csum2=csum2+w%psi/psi
         w%osp(i,:,1,j)=spsave1
      enddo
      w%sp=spsave
   enddo
   cop=1.0_r8+(csum1+csum2)
print*,'csum1=',csum1
print*,'csum2=',csum2
stop
   end subroutine overlap
   
   subroutine tauz(spin,tspin)
   complex(kind=r8) :: spin(:),tspin(:)
   tspin(1)=spin(1)
   tspin(2)=spin(2)
   tspin(3)=-spin(3)
   tspin(4)=-spin(4)
   end subroutine tauz

end module euclidean
