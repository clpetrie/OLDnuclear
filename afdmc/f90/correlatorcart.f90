module correlator
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer, private, parameter :: i8=selected_int_kind(15)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   real(kind=r8), private, parameter :: tiny=1e-5_r8
   integer(kind=i4), private, save :: npart,npair,ntrip
   complex(kind=r8), private, save, allocatable :: sxz0(:,:,:)
   real(kind=r8), private, save, allocatable :: fsvec(:,:,:),fsval(:,:) &
      ,ft(:),fs(:,:,:),fst(:,:,:),f3(:,:,:,:)
   complex(kind=r8), private, allocatable :: f3b(:,:,:,:)
   real(kind=r8), private, save, allocatable :: fstl(:,:,:),fstr(:,:,:)
   real(kind=r8), private, save, allocatable :: fstval(:,:)
!! PP and NN corelations
   real(kind=r8), private, save, allocatable :: ftpp(:),ftnn(:),ftau1(:)
   real(kind=r8), private, save :: fctau
   logical, private, save, allocatable :: dofst(:),doft(:),dofs(:),doftpp(:)
   logical, private, save, allocatable :: doftnn(:)
   integer(kind=i4), private, parameter :: levi(2,3) = &
      reshape((/2,3, 3,1, 1,2/),(/2,3/))
   private :: opmult
   complex(kind=r8), private, save, allocatable :: tau(:,:,:),sigma(:,:,:)
   complex(kind=r8), private, save, allocatable :: sigtau(:,:,:,:,:)
   complex(kind=r8), private, save, allocatable :: sigma1(:,:),tau1(:,:)
   complex(kind=r8), private, save, allocatable :: sigtau1(:,:,:)
   complex(kind=r8), private, save, allocatable :: np0(:),np1(:),pp(:),nn(:)
   complex(kind=r8), private, save, allocatable :: f1b(:,:),f2b(:,:,:)
   real(kind=r8), private, save :: el,eli
   integer(kind=i4), private, save :: dov3
   real(kind=r8), private, save, allocatable :: probinvijk(:)
   logical, private, save, allocatable :: dotrip(:)
!  logical, private, save :: dof3 = .true.
   logical, private, save :: dof3
   logical, private, save :: doindpair1 = .true. !CODY
   logical, private, save :: doindpair2 = .false. !CODY
contains
   subroutine initcormod(npartin,elin)
   integer(kind=i4) :: npartin
   real(kind=r8) :: elin
   npart=npartin
   npair=(npart*(npart-1))/2
   ntrip=(npart*(npart-1)*(npart-2))/6
   if (allocated(sxz0)) then
      deallocate(sxz0,fsvec,fsval,ft,fs,fst,fstl,fstr,fstval,ftpp,ftnn,ftau1)
      deallocate(dofst,doft,dofs,doftpp,doftnn)
      deallocate(tau,sigma,sigtau,sigma1,tau1,sigtau1,np0,np1,pp,nn)
      deallocate(f1b,f2b,f3,f3b)
      deallocate(probinvijk)
   endif
   allocate(sxz0(4,npart,npart))
   allocate(fstl(3,3,npair),fstr(3,3,npair),fstval(3,npair))
   allocate(fsvec(3,3,npair),fsval(3,npair))
   allocate(fs(3,3,npair),fst(3,3,npair))
   allocate(ft(npair),dofst(npair),doft(npair),dofs(npair))
   allocate(doftpp(npair),doftnn(npair))
   allocate(ftpp(npair),ftnn(npair),ftau1(npart))
   allocate(np0(npair),np1(npair),pp(npair),nn(npair))
   allocate(tau(3,3,npair),sigma(3,3,npair),sigtau(3,3,3,3,npair))
   allocate(sigma1(3,npart),tau1(3,npart),sigtau1(3,3,npart))
   allocate(probinvijk(ntrip),dotrip(ntrip))
   allocate(f1b(4,npart),f2b(4,4,npair),f3(3,3,3,ntrip),f3b(4,4,4,ntrip))
   if (elin.gt.0.0_r8) then
      el=elin
      eli=1.0_r8/el
   else 
      el=0.0_r8
      eli=0.0_r8
   endif
   end subroutine initcormod

   subroutine setdov3(dov3in)
   integer(kind=i4) :: dov3in
   dov3=dov3in
   end subroutine setdov3

   subroutine setf3(dof3in)
   logical :: dof3in
   dof3=dof3in
   end subroutine setf3

   subroutine calfop(ftauin,fsigin,fsigtauin,ftauppin,ftaunnin,ff3x,ff3d,cf3x,cf3xd,cf3dd,sp,cut,diag)
   use matrixmod
   real(kind=r8) :: ftauin(:,:),fsigin(:,:,:,:),fsigtauin(:,:,:,:),cut
   real(kind=r8) :: ftauppin(:,:),ftaunnin(:,:)
   complex(kind=r8) :: sp(:,:),spx(4,15,npart)
   integer(kind=i4) :: i,j,k,ij,ijk,jz,is,js,ic,jc,kc,ks,it
   real(kind=r8) :: ff3x(3,npart,3,npart),ff3d(3,npart,3,npart),cf3x,cf3xd,cf3dd
   real(kind=r8) :: f3fac
   logical :: diag
   spx=opmult(sp)
   ij=0
   ftau1=0.0_r8
   fctau=0.0_r8
   f1b=czero
   f2b=czero
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         dofst(ij)=maxval(abs(fsigtauin(:,i,:,j))).gt.cut
         dofs(ij)=maxval(abs(fsigin(:,i,:,j))).gt.cut
         doft(ij)=abs(ftauin(i,j)).gt.cut
         doftpp(ij)=abs(ftauppin(i,j)).gt.cut
         doftnn(ij)=abs(ftaunnin(i,j)).gt.cut
         if (dofst(ij)) then
            fst(:,:,ij)=fsigtauin(:,i,:,j)
            do is=1,3
               do js=1,3
                  do jz=1,4
                     f2b(:,jz,ij)=f2b(:,jz,ij)+fst(is,js,ij)* &
                        (spx(:,3*is+4,i)*spx(jz,3*js+4,j) &
                        +spx(:,3*is+5,i)*spx(jz,3*js+5,j) &
                        +spx(:,3*is+6,i)*spx(jz,3*js+6,j))
                  enddo
               enddo
            enddo
            if (diag) then
                call rsqsvd(fsigtauin(:,i,:,j),fstr(:,:,ij),fstval(:,ij) &
                   ,fstl(:,:,ij),3) 
            endif
         else
            fst(:,:,ij)=0.0_r8
            fstval(:,ij)=0.0_r8
         endif
         if (dofs(ij)) then
            fs(:,:,ij)=fsigin(:,i,:,j)
            do is=1,3
               do js=1,3
                  do jz=1,4
                     f2b(:,jz,ij)=f2b(:,jz,ij)+fs(is,js,ij) &
                        *spx(:,is,i)*spx(jz,js,j)
                  enddo
               enddo
            enddo
            if (diag) then
               fsvec(:,:,ij)=fsigin(:,i,:,j)
               call eigenrs(fsvec(:,:,ij),fsval(:,ij),3)
            endif
         else
            fs(:,:,ij)=0.0_r8
            fsval(:,ij)=0.0_r8
         endif
         if (doft(ij)) then
            ft(ij)=ftauin(i,j)
            do jz=1,4
               f2b(:,jz,ij)=f2b(:,jz,ij)+ft(ij)* &
                  (spx(:,4,i)*spx(jz,4,j)+spx(:,5,i)*spx(jz,5,j) &
                  +spx(:,6,i)*spx(jz,6,j))
            enddo
         else
            ft(ij)=0.0_r8
         endif
         if (doftpp(ij)) then
            ftpp(ij)=ftauppin(i,j)
            ftau1(i)=ftau1(i)+ftpp(ij)
            ftau1(j)=ftau1(j)+ftpp(ij)
            fctau=fctau+0.25_r8*ftpp(ij)
            f1b(:,i)=f1b(:,i)+0.25_r8*ftpp(ij)*spx(:,6,i)
            f1b(:,j)=f1b(:,j)+0.25_r8*ftpp(ij)*spx(:,6,j)
            do jz=1,4
               f2b(:,jz,ij)=f2b(:,jz,ij)+0.25_r8*ftpp(ij)*spx(:,6,i)*spx(jz,6,j)
            enddo
         else
            ftpp(ij)=0.0_r8
         endif
         if (doftnn(ij)) then
            ftnn(ij)=ftaunnin(i,j)
            ftau1(i)=ftau1(i)-ftnn(ij)
            ftau1(j)=ftau1(j)-ftnn(ij)
            fctau=fctau+0.25_r8*ftnn(ij)
            f1b(:,i)=f1b(:,i)-0.25_r8*ftnn(ij)*spx(:,6,i)
            f1b(:,j)=f1b(:,j)-0.25_r8*ftnn(ij)*spx(:,6,j)
            do jz=1,4
               f2b(:,jz,ij)=f2b(:,jz,ij)+0.25_r8*ftnn(ij)*spx(:,6,i)*spx(jz,6,j)
            enddo
         else
            ftnn(ij)=0.0_r8
         endif
      enddo
   enddo
   f3b=0.0_r8
   f3=0.0_r8
   if (.not.dof3) return !skip 3-body correlation
   ijk=0
   do i=1,npart-2
      do j=i+1,npart-1
         do k=j+1,npart
            ijk=ijk+1
            do it=1,3
               do ic=1,3
                  do jc=1,3
                     do kc=1,3
                        f3fac=ff3x(ic,i,levi(1,jc),j)*ff3x(levi(2,jc),j,kc,k) &
                             -ff3x(ic,i,levi(2,jc),j)*ff3x(levi(1,jc),j,kc,k) &
                             +ff3x(jc,j,levi(1,kc),k)*ff3x(levi(2,kc),k,ic,i) &
                             -ff3x(jc,j,levi(2,kc),k)*ff3x(levi(1,kc),k,ic,i) &
                             +ff3x(kc,k,levi(1,ic),i)*ff3x(levi(2,ic),i,jc,j) &
                             -ff3x(kc,k,levi(2,ic),i)*ff3x(levi(1,ic),i,jc,j)
                        f3(ic,jc,kc,ijk)=cf3x*f3fac
                        f3fac=ff3x(ic,i,levi(1,jc),j)*ff3d(levi(2,jc),j,kc,k) &
                             -ff3x(ic,i,levi(2,jc),j)*ff3d(levi(1,jc),j,kc,k) &
                             +ff3x(jc,j,levi(1,kc),k)*ff3d(levi(2,kc),k,ic,i) &
                             -ff3x(jc,j,levi(2,kc),k)*ff3d(levi(1,kc),k,ic,i) &
                             +ff3x(kc,k,levi(1,ic),i)*ff3d(levi(2,ic),i,jc,j) &
                             -ff3x(kc,k,levi(2,ic),i)*ff3d(levi(1,ic),i,jc,j) &
                             +ff3d(ic,i,levi(1,jc),j)*ff3x(levi(2,jc),j,kc,k) &
                             -ff3d(ic,i,levi(2,jc),j)*ff3x(levi(1,jc),j,kc,k) &
                             +ff3d(jc,j,levi(1,kc),k)*ff3x(levi(2,kc),k,ic,i) &
                             -ff3d(jc,j,levi(2,kc),k)*ff3x(levi(1,kc),k,ic,i) &
                             +ff3d(kc,k,levi(1,ic),i)*ff3x(levi(2,ic),i,jc,j) &
                             -ff3d(kc,k,levi(2,ic),i)*ff3x(levi(1,ic),i,jc,j)
                        f3(ic,jc,kc,ijk)=f3(ic,jc,kc,ijk)+cf3xd*f3fac
                        f3fac=ff3d(ic,i,levi(1,jc),j)*ff3d(levi(2,jc),j,kc,k) &
                             -ff3d(ic,i,levi(2,jc),j)*ff3d(levi(1,jc),j,kc,k) &
                             +ff3d(jc,j,levi(1,kc),k)*ff3d(levi(2,kc),k,ic,i) &
                             -ff3d(jc,j,levi(2,kc),k)*ff3d(levi(1,kc),k,ic,i) &
                             +ff3d(kc,k,levi(1,ic),i)*ff3d(levi(2,ic),i,jc,j) &
                             -ff3d(kc,k,levi(2,ic),i)*ff3d(levi(1,ic),i,jc,j)
                        f3(ic,jc,kc,ijk)=f3(ic,jc,kc,ijk)+cf3dd*f3fac
                        do js=1,4
                           do ks=1,4
                              f3b(:,js,ks,ijk)=f3b(:,js,ks,ijk)+f3(ic,jc,kc,ijk) &
                                 *(spx(:,3*ic+it+3,i)*spx(js,3*jc+levi(1,it)+3,j)*spx(ks,3*kc+levi(2,it)+3,k) &
                                  -spx(:,3*ic+it+3,i)*spx(js,3*jc+levi(2,it)+3,j)*spx(ks,3*kc+levi(1,it)+3,k))
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   end subroutine calfop

   subroutine cordet(detrat,sxzin,sp) !CODY added sp here, make sure to remove all if you remove this
! calfop must be called first
   complex(kind=r8) :: detrat,sxzin(:,:,:)
   complex(kind=r8) :: sxmall(npart,15,npart)
   complex(kind=r8) :: ctmp1,d1,d2,d3,d4
   integer(kind=i4) :: i,j,ij,is,js,it
   complex(kind=r8) :: d1b(4,npart),d2b(4,4,npair),d3b(4,4,4,ntrip)
   integer(kind=i4) :: k,l,kl,kop,ks,ls,kt !Added variables start here CODY
   complex(kind=r8) :: sx15(4,15,npart,npart),sx15l(4,15,npart)
   complex(kind=r8) :: d15(15),fkl,sxzk(4,npart,npart,15),sxzl(4,npart,npart)
   complex(kind=r8), intent(in) :: sp(:,:)
   sxz0=sxzin
   d1b=czero
   d2b=czero
   d3b=czero
   detrat=cone+fctau
   call g1bval(d1b,sxz0,cone)
   detrat=detrat+sum(d1b*f1b)
   !call g2bval(d2b,sxz0,cone)
! Here are the explicit loops from g2bval and independent pair stuff, CODY
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         do js=1,4
            d2b(:,js,ij)=d2b(:,js,ij) &
               +cone*(sxz0(:,i,i)*sxz0(js,j,j)-sxz0(:,i,j)*sxz0(js,j,i))
            detrat=detrat+sum(d2b(:,js,ij)*f2b(:,js,ij))
            if(doindpair1) then !Start independant pair terms here
               do k=1,npart
                  sx15(:,:,:,k)=conjg(opmult(conjg(sxz0(:,k,:))))
               enddo
               do k=1,npart-1
                  if(k.le.i .or. k.eq.j) cycle !independent pairs
!???                  if(k.lt.i) cycle !all second order terms
                  do kop=1,15
                     call sxzupdate(sxzk(:,:,:,kop),d15(kop),sxz0,k,sx15(:,kop,:,k),sp(:,k))
                  enddo
                  do l=k+1,npart
                     if(l.eq.i .or. l.eq.j) cycle !independent pairs
!???                     if((k.eq.i .and. l.eq.j) .or. ((i+j)>(k+l))) cycle ! all second order terms
!Mathematica ij = FullSimplify[(j - i) + (i - 1)*npart + Sum[-n, {n, 1, i - 1}]]
                     kl=l-k*(1+k-2*npart)/2-npart
                     if(doft(kl) .or. doftpp(kl) .or. doftnn(kl)) then
                        do kt=1,3
                           sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,3+kt))))
                           call sxzupdate(sxzl,d2,sxzk(:,:,:,3+kt),l,sx15l(:,3+kt,:),sp(:,l))
                           fkl=d15(3+kt)*d2*ft(kl)
                           if (kt==3 .and. doftpp(kl)) fkl=fkl+0.25_r8*ftpp(kl)
                           if (kt==3 .and. doftnn(kl)) fkl=fkl+0.25_r8*ftnn(kl)
!???                           fkl=ft(kl)
                           call g2bval(d2b,sxzl,fkl)
                        enddo
                     endif
                     if(dofs(kl)) then
                        do ks=1,3
                           sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,ks))))
                           do ls=1,3
                              call sxzupdate(sxzl,d2,sxzk(:,:,:,ks),l,sx15l(:,ls,:),sp(:,l))
                              fkl=d15(ks)*d2*fs(ks,ls,kl)
!???                              fkl=fs(ks,ls,kl)
                              call g2bval(d2b,sxzl,fkl)
                           enddo
                        enddo
                     endif
                     if(dofst(kl)) then
                        do kt=1,3
                           do ks=1,3
                              sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,3*ks+kt+3))))
                              do ls=1,3
                                 call sxzupdate(sxzl,d2,sxzk(:,:,:,3*ks+kt+3),l,sx15l(:,3*ls+kt+3,:),sp(:,l))
                                 fkl=d15(3*ks+kt+3)*d2*fst(ks,ls,kl)
!???                                 fkl=fst(ks,ls,kl)
                                 call g2bval(d2b,sxzl,fkl)
                              enddo
                           enddo
                        enddo
                     endif
                     do ls=1,4
                        do is=1,4
!???                           detrat=detrat+sum(d2b(:,ls,kl)*f2b(:,ls,kl))
!???                           detrat=detrat+sum(d2b(:,ls,kl)*f2b(:,ls,kl)*f2b(:,js,ij))
                           detrat=detrat+sum(d2b(:,ls,kl)*f2b(:,ls,kl)*f2b(is,js,ij))!*detrat
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo
!
!!!   detrat=cone+fctau+sum(d1b*f1b)+sum(d2b*f2b)
   if (dof3) then
      call g3bval(d3b,sxz0,cone,.true.)
      detrat=detrat+sum(d3b*f3b)
   endif
   end subroutine cordet

   subroutine corindpair(sp,sxzin,detratin,i,j,d1b,d2b,d3b) !CODY
   complex(kind=r8), intent(in) :: sp(:,:)
   complex(kind=r8), intent(inout) :: d1b(:,:),d2b(:,:,:),d3b(:,:,:,:)
   complex(kind=r8), intent(in) :: sxzin(:,:,:)
   integer(kind=i4), intent(in) :: i,j
   complex(kind=r8) :: sxzk(4,npart,npart,15)
   complex(kind=r8) :: fkl
   complex(kind=r8) :: detrat,detratin
   complex(kind=r8) :: sxzl(4,npart,npart),d2ind,d15ind(15) !d2ind is the same as d2 but for corindpair
   complex(kind=r8) :: sx15ind(4,15,npart,npart),sx15l(4,15,npart)
   integer(kind=i4) :: k,l,kl,kop,ks,kt,ls
   do k=1,npart
     sx15ind(:,:,:,k)=conjg(opmult(conjg(sxzin(:,k,:))))
   enddo
   do k=1,npart-1
      if (k.le.i .or. k.eq.j) cycle
      do kop=1,15
         call sxzupdate(sxzk(:,:,:,kop),d15ind(kop),sxzin,k,sx15ind(:,kop,:,k),sp(:,k))
      enddo
      do l=k+1,npart
         if (l.eq.i .or. l.eq.j) cycle
!Mathematica ij = FullSimplify[(j - i) + (i - 1)*npart + Sum[-n, {n, 1, i - 1}]]
         kl=l-k*(1+k-2*npart)/2-npart
         if (doft(kl)) then
            do kt=1,2
               fkl=ft(kl)
               sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,3+kt))))
               call sxzupdate(sxzl,d2ind,sxzk(:,:,:,3+kt),l,sx15l(:,3+kt,:),sp(:,l))
               detrat=detratin*d15ind(3+kt)*d2ind
               fkl=detrat*fkl
               call g1bval(d1b,sxzl,fkl)
               call g2bval(d2b,sxzl,fkl)
               call g3bval(d3b,sxzl,fkl,.false.)
            enddo
         endif
         if (doft(kl).or.doftpp(kl).or.doftnn(kl)) then
            kt=3
            fkl=ft(kl)
            if (doftpp(kl)) fkl=fkl+0.25_r8*ftpp(kl)
            if (doftnn(kl)) fkl=fkl+0.25_r8*ftnn(kl)
            sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,3+kt))))
            call sxzupdate(sxzl,d2ind,sxzk(:,:,:,3+kt),l,sx15l(:,3+kt,:),sp(:,l))
            detrat=detratin*d15ind(3+kt)*d2ind
            fkl=detrat*fkl
            call g1bval(d1b,sxzl,fkl)
            call g2bval(d2b,sxzl,fkl)
            call g3bval(d3b,sxzl,fkl,.false.)
         endif
         if (dofs(kl)) then
            do ks=1,3
               sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,ks))))
               do ls=1,3
                  call sxzupdate(sxzl,d2ind,sxzk(:,:,:,ks),l,sx15l(:,ls,:),sp(:,l))
                  detrat=detratin*d15ind(ks)*d2ind
                  fkl=detrat*fs(ks,ls,kl)
                  call g1bval(d1b,sxzl,fkl)
                  call g2bval(d2b,sxzl,fkl)
                  call g3bval(d3b,sxzl,fkl,.false.)
               enddo
            enddo
         endif
         if (dofst(kl)) then
            do kt=1,3
               do ks=1,3
                  sx15l(:,:,:)=conjg(opmult(conjg(sxzk(:,l,:,3*ks+kt+3))))
                  do ls=1,3
                     call sxzupdate(sxzl,d2ind,sxzk(:,:,:,3*ks+kt+3),l,sx15l(:,3*ls+kt+3,:),sp(:,l))
                     detrat=detratin*d15ind(3*ks+kt+3)*d2ind
                     fkl=detrat*fst(ks,ls,kl)
                     call g1bval(d1b,sxzl,fkl)
                     call g2bval(d2b,sxzl,fkl)
                     call g3bval(d3b,sxzl,fkl,.false.)
                  enddo
               enddo
            enddo
         endif
      enddo
   enddo
   end subroutine corindpair

   subroutine corpsi(sp,d1b,d2b,d3b)
   complex(kind=r8), intent(in) :: sp(:,:)
   complex(kind=r8), intent(inout):: d1b(:,:),d2b(:,:,:),d3b(:,:,:,:)
   complex(kind=r8) :: fij,f1,fijk
   complex(kind=r8) :: detrat,sxzi(4,npart,npart,15)
   complex(kind=r8) :: sxzj(4,npart,npart),d1,d2,d15(15)
   complex(kind=r8) :: sx15(4,15,npart,npart),sx15j(4,15,npart)
   complex(kind=r8) :: sx15j1(4,15,npart),sx15j2(4,15,npart)
   complex(kind=r8) :: sx15k1(4,15,npart),sx15k2(4,15,npart)
   complex(kind=r8) :: sxzi1(4,npart,npart),di1
   complex(kind=r8) :: sxzj1(4,npart,npart),sxzj2(4,npart,npart)
   complex(kind=r8) :: sxzk(4,npart,npart),dj1,dj2,dk
   integer(kind=i4) :: i,j,ij,iop,is,it,js,k,ks,ijk
   d1b=czero
   d2b=czero
   d3b=czero
   call g1bval(d1b,sxz0,cone+fctau)
   call g2bval(d2b,sxz0,cone+fctau)
   call g3bval(d3b,sxz0,cone+fctau,.false.)
   do i=1,npart
      sx15(:,:,:,i)=conjg(opmult(conjg(sxz0(:,i,:))))
   enddo
   do i=1,npart
      if (abs(ftau1(i)).le.0.0_r8) cycle
      call sxzupdate(sxzj,d1,sxz0,i,sx15(:,6,:,i),sp(:,i))
      f1=d1*0.25_r8*ftau1(i)
      call g1bval(d1b,sxzj,f1)
      call g2bval(d2b,sxzj,f1)
      call g3bval(d3b,sxzj,f1,.false.)
   enddo
   ij=0
   do i=1,npart-1
      do iop=1,15
         call sxzupdate(sxzi(:,:,:,iop),d15(iop),sxz0,i,sx15(:,iop,:,i),sp(:,i))
      enddo
      do j=i+1,npart
         ij=ij+1
         if (doft(ij)) then
            do it=1,2
               fij=ft(ij)
               sx15j(:,:,:)=conjg(opmult(conjg(sxzi(:,j,:,3+it))))
               call sxzupdate(sxzj,d2,sxzi(:,:,:,3+it),j,sx15j(:,3+it,:),sp(:,j))
               detrat=d15(3+it)*d2
               fij=detrat*fij
               call g1bval(d1b,sxzj,fij)
               call g2bval(d2b,sxzj,fij)
               call g3bval(d3b,sxzj,fij,.false.)
               if (doindpair2) then
                  call corindpair(sp,sxzj,detrat,i,j,d1b,d2b,d3b) !CODY
               endif
            enddo
         endif
         if (doft(ij).or.doftpp(ij).or.doftnn(ij)) then
            it=3
            fij=ft(ij)
            if (doftpp(ij)) fij=fij+0.25_r8*ftpp(ij)
            if (doftnn(ij)) fij=fij+0.25_r8*ftnn(ij)
            sx15j(:,:,:)=conjg(opmult(conjg(sxzi(:,j,:,3+it))))
            call sxzupdate(sxzj,d2,sxzi(:,:,:,3+it),j,sx15j(:,3+it,:),sp(:,j))
            detrat=d15(3+it)*d2
            fij=detrat*fij
            call g1bval(d1b,sxzj,fij)
            call g2bval(d2b,sxzj,fij)
            call g3bval(d3b,sxzj,fij,.false.)
            if (doindpair2) then
               call corindpair(sp,sxzj,detrat,i,j,d1b,d2b,d3b) !CODY
            endif
         endif
         if (dofs(ij)) then
            do is=1,3
               sx15j(:,:,:)=conjg(opmult(conjg(sxzi(:,j,:,is))))
               do js=1,3
                  call sxzupdate(sxzj,d2,sxzi(:,:,:,is),j,sx15j(:,js,:),sp(:,j))
                  detrat=d15(is)*d2
                  fij=detrat*fs(is,js,ij)
                  call g1bval(d1b,sxzj,fij)
                  call g2bval(d2b,sxzj,fij)
                  call g3bval(d3b,sxzj,fij,.false.)
                  if (doindpair2) then
                     call corindpair(sp,sxzj,detrat,i,j,d1b,d2b,d3b) !CODY
                  endif
               enddo
            enddo
         endif         
         if (dofst(ij)) then 
            do it=1,3
               do is=1,3
                  sx15j(:,:,:)=conjg(opmult(conjg(sxzi(:,j,:,3*is+it+3))))
                  do js=1,3
                     call sxzupdate(sxzj,d2,sxzi(:,:,:,3*is+it+3),j &
                        ,sx15j(:,3*js+it+3,:),sp(:,j))
                     detrat=d15(3*is+it+3)*d2
                     fij=detrat*fst(is,js,ij)
                     call g1bval(d1b,sxzj,fij)
                     call g2bval(d2b,sxzj,fij)
                     call g3bval(d3b,sxzj,fij,.false.)
                     if (doindpair2) then
                        call corindpair(sp,sxzj,detrat,i,j,d1b,d2b,d3b) !CODY
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo
   enddo
   if (.not.dof3) return !skip 3-body correlation
   do i=1,npart-2
      do is=1,3
         do it=1,3
            call sxzupdate(sxzi1,di1,sxz0,i,sx15(:,3*is+it+3,:,i),sp(:,i))
            do j=i+1,npart-1
               sx15j(:,:,:)=conjg(opmult(conjg(sxzi1(:,j,:))))
               do js=1,3
                  call sxzupdate(sxzj1,dj1,sxzi1,j,sx15j(:,3*js+3+levi(1,it),:),sp(:,j))
                  call sxzupdate(sxzj2,dj2,sxzi1,j,sx15j(:,3*js+3+levi(2,it),:),sp(:,j))
                  do k=j+1,npart
!maple ijk := simplify(sum((n-l)*(n-l-1)/2,l=1..i-1)+sum(n-l,l=i+1..j-1)+k-j);
                     ijk=(i*(i-1)*(i-3*npart+4))/6 &
                        +((npart-2)*(npart-1)*(i-1))/2-2 &
                        +((2*npart-4-j+1)*(j-2))/2+k
                     sx15k1(:,:,:)=conjg(opmult(conjg(sxzj1(:,k,:))))
                     sx15k2(:,:,:)=conjg(opmult(conjg(sxzj2(:,k,:))))
                     do ks=1,3
                        call sxzupdate(sxzk,dk,sxzj1,k,sx15k1(:,3*ks+3+levi(2,it),:),sp(:,k))
                        fijk=f3(is,js,ks,ijk)*di1*dj1*dk
                        call g1bval(d1b,sxzk,fijk)
                        call g2bval(d2b,sxzk,fijk)
                        call g3bval(d3b,sxzk,fijk,.false.)
                        call sxzupdate(sxzk,dk,sxzj2,k,sx15k2(:,3*ks+3+levi(1,it),:),sp(:,k))
                        fijk=-f3(is,js,ks,ijk)*di1*dj2*dk
                        call g1bval(d1b,sxzk,fijk)
                        call g2bval(d2b,sxzk,fijk)
                        call g3bval(d3b,sxzk,fijk,.false.)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   end subroutine corpsi

   subroutine v6tot(x,sp,v2,v3,v4,v5,v6,cvs,tnic,tni2pia,tni2pitm,tni2pic,&
      tni2picxd,tni2picdd,tnivd1,tnivd2,tnive,tni2piaxd,tni2piadd,tni2piapr,&
      tni2piaxdpr,tni2piaddpr)
   use v3bpot
! calfop and cordet must be called first.
   real(kind=r8) :: x(:,:)
   complex(kind=r8) :: sp(:,:),spx(4,15,npart)
   complex(kind=r8) :: cvs(:)
   complex(kind=r8) :: d1b(4,npart),d2b(4,4,npair),d3b(4,4,4,ntrip)
   integer(kind=i4) :: i,j,ij,iop,is,it,js
   real(kind=r8) :: tnic
   complex(kind=r8) :: tni2pia,tni2pitm,tni2pic,tni2picxd,tni2picdd,tnivd1
   complex(kind=r8) :: tnivd2,tnive,tni2piaxd,tni2piadd
   complex(kind=r8) :: tni2piapr,tni2piaxdpr,tni2piaddpr
   real(kind=r8) :: v2(:,:),v3(:,:),v4(:,:),v5(:,:,:,:),v6(:,:,:,:)
   call corpsi(sp,d1b,d2b,d3b)
   if (dov3.ne.0) then
      do i=1,ntrip
        if (dotrip(i)) d3b(:,:,:,i)=d3b(:,:,:,i)*probinvijk(i)
      enddo
   endif
   spx=opmult(sp)
   call op1val(d1b,spx)
   call op2val(d2b,sp,spx)
   cvs=czero
   call calpot(cvs,v2,v3,v4,v5,v6)
   call v3bval(x,spx,d2b,d3b,tnic,tni2pia,tni2pitm,tni2pic,tni2picxd,tni2picdd,&
      tnivd1,tnivd2,tnive,tni2piaxd,tni2piadd,tni2piapr,tni2piaxdpr,tni2piaddpr)
   end subroutine v6tot

   subroutine g1bval(d1b,sxz,fij)
   complex(kind=r8), intent(inout) :: d1b(:,:)
   complex(kind=r8), intent(in) :: sxz(:,:,:),fij
   integer(kind=i4) :: i
   do i=1,npart
      d1b(:,i)=d1b(:,i)+fij*sxz(:,i,i)
   enddo
   end subroutine g1bval

   subroutine g2bval(d2b,sxz,fij)
   complex(kind=r8), intent(inout) :: d2b(:,:,:)
   complex(kind=r8), intent(in) :: sxz(:,:,:),fij
   integer(kind=i4) :: i,j,ij,js
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         do js=1,4
            d2b(:,js,ij)=d2b(:,js,ij) &
               +fij*(sxz(:,i,i)*sxz(js,j,j)-sxz(:,i,j)*sxz(js,j,i))
         enddo
      enddo
   enddo
   end subroutine g2bval

   subroutine g3bval(d3b,sxz,fij,doall)
   complex(kind=r8), intent(inout):: d3b(:,:,:,:)
   complex(kind=r8), intent(in):: sxz(:,:,:),fij
   integer(kind=i4) :: i,j,k,js,ks,ijk
   logical :: doall
   if (dov3.eq.0) return
   ijk=0
   do i=1,npart-2
      do j=i+1,npart-1
         do k=j+1,npart
            ijk=ijk+1
            if (dotrip(ijk).or.doall) then
               do ks=1,4
                  do js=1,4
                     d3b(:,js,ks,ijk)=d3b(:,js,ks,ijk)+fij*( &
                         sxz(:,i,i)*(sxz(js,j,j)*sxz(ks,k,k) &
                        -sxz(js,j,k)*sxz(ks,k,j))+sxz(:,i,j) &
                        *(sxz(js,j,k)*sxz(ks,k,i)-sxz(js,j,i)*sxz(ks,k,k)) &
                        +sxz(:,i,k)*(sxz(js,j,i)*sxz(ks,k,j) &
                         -sxz(js,j,j)*sxz(ks,k,i)))
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo
   end subroutine g3bval

   subroutine sxzupdate(sxznew,detrat,sxzold,i,opi,sp)
   complex(kind=r8), intent(out) :: sxznew(4,npart,npart)
   complex(kind=r8), intent(out) :: detrat
   complex(kind=r8), intent(in) :: sxzold(4,npart,npart)
   complex(kind=r8), intent(in) :: sp(4),opi(4,npart)
   integer(kind=i4), intent(in) :: i
   integer(kind=i4) :: m
   complex(kind=r8) :: di(npart),detinv
   di=opi(1,:)*sp(1)+opi(2,:)*sp(2)+opi(3,:)*sp(3)+opi(4,:)*sp(4)
   detrat=di(i)
   detinv=cone/detrat
   di=di*detinv
   do m=1,npart
      sxznew(:,:,m)=sxzold(:,:,m)-di(m)*sxzold(:,:,i)
   enddo
   do m=1,npart
      sxznew(:,i,m)=opi(:,m)-di(m)*opi(:,i)
   enddo
   sxznew(:,:,i)=sxzold(:,:,i)*detinv
   sxznew(:,i,i)=opi(:,i)*detinv
   end subroutine sxzupdate

   function opmult(sp)
   complex(kind=r8) :: sp(:,:),opmult(4,15,npart)
!
! The order is 1-3 sx,sy,sz, 4-6 tx,ty,tx, 7-9 sx*(tx,ty,tz)
! 10-12 sy*(tx,ty,tz), 13-15 sz*(tx,ty,tz)
!
! multiply by sigma
!
   opmult(1,1,:)=sp(2,:)
   opmult(2,1,:)=sp(1,:)
   opmult(3,1,:)=sp(4,:)
   opmult(4,1,:)=sp(3,:)
   opmult(1,2,:)=-ci*sp(2,:)
   opmult(2,2,:)=ci*sp(1,:)
   opmult(3,2,:)=-ci*sp(4,:)
   opmult(4,2,:)=ci*sp(3,:)
   opmult(1,3,:)=sp(1,:)
   opmult(2,3,:)=-sp(2,:)
   opmult(3,3,:)=sp(3,:)
   opmult(4,3,:)=-sp(4,:)
!
! multiply by tau
!
   opmult(1,4,:)=sp(3,:)
   opmult(2,4,:)=sp(4,:)
   opmult(3,4,:)=sp(1,:)
   opmult(4,4,:)=sp(2,:)
   opmult(1,5,:)=-ci*sp(3,:)
   opmult(2,5,:)=-ci*sp(4,:)
   opmult(3,5,:)=ci*sp(1,:)
   opmult(4,5,:)=ci*sp(2,:)
   opmult(1,6,:)=sp(1,:)
   opmult(2,6,:)=sp(2,:)
   opmult(3,6,:)=-sp(3,:)
   opmult(4,6,:)=-sp(4,:)
!
! multiply by sigma tau
!
   opmult(1,7:13:3,:)=opmult(3,1:3:1,:)
   opmult(2,7:13:3,:)=opmult(4,1:3:1,:)
   opmult(3,7:13:3,:)=opmult(1,1:3:1,:)
   opmult(4,7:13:3,:)=opmult(2,1:3:1,:)
   opmult(1,8:14:3,:)=-ci*opmult(3,1:3:1,:)
   opmult(2,8:14:3,:)=-ci*opmult(4,1:3:1,:)
   opmult(3,8:14:3,:)=ci*opmult(1,1:3:1,:)
   opmult(4,8:14:3,:)=ci*opmult(2,1:3:1,:)
   opmult(1,9:15:3,:)=opmult(1,1:3:1,:)
   opmult(2,9:15:3,:)=opmult(2,1:3:1,:)
   opmult(3,9:15:3,:)=-opmult(3,1:3:1,:)
   opmult(4,9:15:3,:)=-opmult(4,1:3:1,:)
   end function opmult

   subroutine op1val(d1b,spx)
   complex(kind=r8), intent(in) :: d1b(:,:),spx(:,:,:)
   integer(kind=i4) :: i,ic,it
   do i=1,npart
      do ic=1,3
         sigma1(ic,i)=sum(d1b(:,i)*spx(:,ic,i))
         tau1(ic,i)=sum(d1b(:,i)*spx(:,ic+3,i))
         do it=1,3
            sigtau1(ic,it,i)=sum(d1b(:,i)*spx(:,3*(ic-1)+it+6,i))
         enddo
      enddo
   enddo
   end subroutine op1val

   subroutine op2val(d2b,sp,spx)
   complex(kind=r8), intent(in) :: d2b(:,:,:),sp(:,:),spx(:,:,:)
   integer(kind=i4) :: i,j,ic,jc,it,ij,k,js
   complex(kind=r8) :: tz(4,4),sz(4,4),stz(4,4),ppz(4,4),nnz(4,4),np0z(4,4),np1z(4,4)
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         do ic=1,3
            do k=1,4
               tz(:,k)=spx(:,ic+3,i)*spx(k,ic+3,j)
            enddo
            do jc=1,3
               tau(ic,jc,ij)=sum(d2b(:,:,ij)*tz(:,:))
               do k=1,4
                  sz(:,k)=spx(:,ic,i)*spx(k,jc,j)
               enddo
               sigma(ic,jc,ij)=sum(d2b(:,:,ij)*sz(:,:))
               do it=1,3
                  do k=1,4
                     stz(:,k)=spx(:,3*(ic-1)+it+6,i)*spx(k,3*(jc-1)+it+6,j)
                  enddo
                  sigtau(ic,jc,it,it,ij)=sum(d2b(:,:,ij)*stz(:,:))
               enddo
            enddo
         enddo
         ppz=0.0_r8
         nnz=0.0_r8
         np0z=0.0_r8
         np1z=0.0_r8
         do js=1,2
            ppz(1:2,js)=ppz(1:2,js)+sp(1:2,i)*sp(js,j)
            nnz(3:4,js+2)=nnz(3:4,js+2)+sp(3:4,i)*sp(js+2,j)
         enddo
         np0z(1,3)=np0z(1,3)+sp(1,i)*sp(3,j)-sp(3,i)*sp(1,j)
         np0z(1,4)=np0z(1,4)+sp(1,i)*sp(4,j)-sp(3,i)*sp(2,j)
         np0z(2,3)=np0z(2,3)+sp(2,i)*sp(3,j)-sp(4,i)*sp(1,j)
         np0z(2,4)=np0z(2,4)+sp(2,i)*sp(4,j)-sp(4,i)*sp(2,j)
         np0z(3,1)=np0z(3,1)+sp(3,i)*sp(1,j)-sp(1,i)*sp(3,j)
         np0z(3,2)=np0z(3,2)+sp(3,i)*sp(2,j)-sp(1,i)*sp(4,j)
         np0z(4,1)=np0z(4,1)+sp(4,i)*sp(1,j)-sp(2,i)*sp(3,j)
         np0z(4,2)=np0z(4,2)+sp(4,i)*sp(2,j)-sp(2,i)*sp(4,j)
         np1z(1,3)=np1z(1,3)+sp(1,i)*sp(3,j)+sp(3,i)*sp(1,j)
         np1z(1,4)=np1z(1,4)+sp(1,i)*sp(4,j)+sp(3,i)*sp(2,j)
         np1z(2,3)=np1z(2,3)+sp(2,i)*sp(3,j)+sp(4,i)*sp(1,j)
         np1z(2,4)=np1z(2,4)+sp(2,i)*sp(4,j)+sp(4,i)*sp(2,j)
         np1z(3,1)=np1z(3,1)+sp(3,i)*sp(1,j)+sp(1,i)*sp(3,j)
         np1z(3,2)=np1z(3,2)+sp(3,i)*sp(2,j)+sp(1,i)*sp(4,j)
         np1z(4,1)=np1z(4,1)+sp(4,i)*sp(1,j)+sp(2,i)*sp(3,j)
         np1z(4,2)=np1z(4,2)+sp(4,i)*sp(2,j)+sp(2,i)*sp(4,j)
         np0(ij)=0.5_r8*sum(d2b(:,:,ij)*np0z(:,:))
         np1(ij)=0.5_r8*sum(d2b(:,:,ij)*np1z(:,:))
         pp(ij)=sum(d2b(:,:,ij)*ppz(:,:))
         nn(ij)=sum(d2b(:,:,ij)*nnz(:,:))
      enddo
   enddo
   end subroutine op2val

   subroutine calpot(cvs,v2,v3,v4,v5,v6)
   real(kind=r8) :: v2(:,:),v3(:,:),v4(:,:),v5(:,:,:,:),v6(:,:,:,:)
   complex(kind=r8) :: cvs(:)
   integer(kind=i4) :: it,is,js,ij,i,j
   cvs=czero
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         do it=1,3
            cvs(2)=cvs(2)+tau(it,it,ij)*v2(i,j)
            do is=1,3
               cvs(4)=cvs(4)+sigtau(is,is,it,it,ij)*v4(i,j)
            enddo
         enddo
         do is=1,3
            cvs(3)=cvs(3)+sigma(is,is,ij)*v3(i,j)
            do js=1,3
               cvs(5)=cvs(5)+sigma(is,js,ij)*v5(is,i,js,j)
               do it=1,3
                  cvs(6)=cvs(6)+sigtau(is,js,it,it,ij)*v6(is,i,js,j)
               enddo
            enddo
         enddo
      enddo
   enddo
   end subroutine calpot

   subroutine calpotcoul(vcoul,vcoulc,vcoulmat)
   complex(kind=r8) :: vcoul
   real(kind=r8) :: vcoulc,vcoulmat(:,:)
   integer(kind=i4) :: i,j,ij
   vcoul=0.0_r8
   vcoulc=0.0_r8
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         vcoulc=vcoulc+vcoulmat(i,j)
         vcoul=vcoul+(tau1(3,i)+tau1(3,j)+tau(3,3,ij))*vcoulmat(i,j)
      enddo
   enddo
   vcoulc=0.25_r8*vcoulc
   vcoul=0.25_r8*vcoul
   end subroutine calpotcoul

   subroutine getop1(sig1in,tau1in,sigtau1in)
   complex(kind=r8) :: sig1in(:,:),tau1in(:,:),sigtau1in(:,:,:)
   sig1in=sigma1
   tau1in=tau1
   sigtau1in=sigtau1
   end subroutine getop1

   subroutine getop2(sig2in,tau2in,sigtau2in,np0in,np1in,ppin,nnin)
   complex(Kind=r8), dimension(:,:,:) :: sig2in,tau2in
   complex(kind=r8), dimension(:,:,:,:,:) :: sigtau2in
   complex(kind=r8), dimension(:) :: np0in,np1in,ppin,nnin
   sig2in=sigma
   tau2in=tau
   sigtau2in=sigtau
   np0in=np0
   np1in=np1
   ppin=pp
   nnin=nn
   end subroutine getop2

   subroutine getoptens(x,tin,ttauin)
   real(kind=r8) :: x(:,:),dx(3)
   complex(kind=r8), dimension(:) :: tin,ttauin
   integer(kind=i4) :: ij,i,j,is,js,it
   tin=0.0_r8
   ttauin=0.0_r8
   ij=1
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx*eli)
         dx=dx/sqrt(sum(dx**2))
         do is=1,3
            do js=1,3
               tin(ij)=tin(ij)+sigma(is,js,ij)*3.0_r8*dx(is)*dx(js)
               if (is.eq.js) tin(ij)=tin(ij)-sigma(is,js,ij)
               do it=1,3
                  ttauin(ij)=ttauin(ij)+sigtau(is,js,it,it,ij)*3.0_r8*dx(is)*dx(js)
                  if (is.eq.js) ttauin(ij)=ttauin(ij)-sigtau(is,js,it,it,ij)
               enddo
            enddo
         enddo
         ij=ij+1
      enddo
   enddo
   end subroutine getoptens

   subroutine v3bval(x,spx,d2b,d3b,vc,v3a,v3tm,v3c,v3cxd,v3cdd,v3d1,v3d2,v3e,v3axdel,v3adeldel,v3asc,v3axdelsc,v3adeldelsc)
   use v3bpot
   real(kind=r8) :: x(:,:)
   complex(kind=r8), intent(in) :: d2b(4,4,npair),d3b(4,4,4,ntrip),spx(:,:,:)
   complex(kind=r8), intent(out) :: v3a,v3tm,v3c,v3cxd,v3cdd,v3d1,v3d2,v3asc,v3e
   complex(kind=r8) :: v3axdel,v3adeldel,v3axdelsc,v3adeldelsc
   real(kind=r8) :: vc,xpi(3,npart,3,npart),xxpi(3,npart,3,npart)
   real(kind=r8) :: xd(3,npart,3,npart),xdel(3,npart,3,npart)
   real(kind=r8) :: xpic(3,npart,3,npart)
   real(kind=r8) :: delta(npart,npart),ddelta(npart,npart)
   real(kind=r8) :: vfacxx,vfacxd,vfacdd
   real(kind=r8) :: g2s3b(3,npart,npart),vtm
   real(kind=r8) :: delsum(npart)
   complex(kind=r8) :: v2tmp1(4,4),v2tmp2(4,4),v2tmp3a(4,4),v2tmp3b(4,4),v2tmp4(4,4)
   complex(kind=r8) :: v2xdel(4,4),v2deldel(4,4)
   complex(kind=r8) :: v3tmp1(4,4,4),v3tmp2(4,4,4),v3tmp3(4,4,4)
   integer(kind=i4) :: i,j,k,ic,jc,kc,ic3,jc3,kc3,js,ks,ij,ijk,it
   real(kind=r8) :: a2p3b,a2s3b,avd,acfac,a2psc,a2pxdsc,a2pddsc,ave,dfac
   complex(kind=r8) :: spxfac(4)
   real(kind=r8) :: rscal(3) ! match the dimension of rscal in jastrowtabop
   rscal=1.0_r8
   vc=0.0_r8
   if (dov3.eq.0) then
      v3asc=0.0_r8
      v3a=0.0_r8
      v3axdelsc=0.0_r8
      v3axdel=0.0_r8
      v3adeldelsc=0.0_r8
      v3adeldel=0.0_r8
      v3tm=0.0_r8
      v3d1=0.0_r8
      v3d2=0.0_r8
      v3e=0.0_r8
      v3c=0.0_r8
      v3cxd=0.0_r8
      v3cdd=0.0_r8
      return
   endif
   call hstnipot(x,vc,xpi,g2s3b,delta,a2p3b,a2s3b,avd,ave,acfac,a2psc,a2pxdsc,a2pddsc,dfac,rscal)
   delsum=0.0_r8
   xd=0.0_r8
   xxpi=0.0_r8
   xxpi=reshape(matmul(reshape(xpi,(/3*npart,3*npart/)), &
      reshape(xpi,(/3*npart,3*npart/))),(/3,npart,3,npart/))
   xdel=0.0_r8
   ddelta=0.0_r8
   if (dov3.eq.2) then
      do i=1,npart-1
         do j=i+1,npart
            delsum(i)=delsum(i)+delta(i,j)
            delsum(j)=delsum(j)+delta(i,j)
            xd(1,i,1,j)=-delta(i,j)
            xd(2,i,2,j)=-delta(i,j)
            xd(3,i,3,j)=-delta(i,j)
            xd(1,j,1,i)=-delta(i,j)
            xd(2,j,2,i)=-delta(i,j)
            xd(3,j,3,i)=-delta(i,j)
         enddo
      enddo
      xdel=reshape(matmul(reshape(xpi,(/3*npart,3*npart/)), &
         reshape(xd,(/3*npart,3*npart/))),(/3,npart,3,npart/))
      xdel=xdel+reshape(matmul(reshape(xd,(/3*npart,3*npart/)), &
         reshape(xpi,(/3*npart,3*npart/))),(/3,npart,3,npart/))
      ddelta=matmul(delta,delta)
   endif
   v3a=czero
   v3axdel=czero
   v3adeldel=czero
   v3tm=czero
   v3d1=czero
   v3d2=czero
   v3e=czero
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         v2tmp1=0.0_r8
         v2tmp2=0.0_r8
         v2tmp3a=0.0_r8
         v2tmp3b=0.0_r8
         v2tmp4=0.0_r8
         v2xdel=0.0_r8
         v2deldel=0.0_r8
         do ic=1,3
            ic3=(ic+1)*3
            do jc=1,3
               jc3=(jc+1)*3
               vtm=sum(g2s3b(ic,:,i)*g2s3b(jc,:,j))
               do it=1,3
                  do js=1,4
                     v2tmp1(:,js)=v2tmp1(:,js)+xxpi(ic,i,jc,j) &
                        *spx(:,ic3+it,i)*spx(js,jc3+it,j)  ! anticommutator
                     if (dov3.eq.2) then
                        v2xdel(:,js)=v2xdel(:,js)+xdel(ic,i,jc,j) &
                           *spx(:,ic3+it,i)*spx(js,jc3+it,j)  ! Vd,Cc3
                        v2tmp2(:,js)=v2tmp2(:,js)+vtm &
                           *spx(:,ic3+it,i)*spx(js,jc3+it,j)  ! Tucson-Melbourne
!                       v2tmp3a(:,js)=v2tmp3a(:,js)+xpi(ic,i,jc,j)*( &  ! Vd,1 Kevin's version
!                          delsum(j)+delsum(i)-2.0_r8*delta(i,j)) &
!                          *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                        v2tmp3a(:,js)=v2tmp3a(:,js)-xdel(ic,i,jc,j) & ! Vd,1 Ingo/Joel's version
                                     *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                        if (jc.eq.ic) then
                           v2deldel(:,js)=v2deldel(:,js)+ddelta(i,j) &
                                *spx(:,ic3+it,i)*spx(js,jc3+it,j)  ! Ve,Cc3
!                          v2tmp3b(:,js)=v2tmp3b(:,js)+xd(ic,i,ic,j) &  ! Vd,2 Kevin's version
!                               *(delsum(j)+delsum(i)-2.0_r8*delta(i,j)) & 
!                               *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                           v2tmp3b(:,js)=v2tmp3b(:,js)-2.0_r8*ddelta(i,j) &  ! Vd,2 Ingo/Joel's version
                                *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                        endif
                     endif
                  enddo
               enddo
            enddo
            if (dov3.eq.2) then
               do js=1,4
                  v2tmp4(:,js)=v2tmp4(:,js)+ddelta(i,j)*spx(:,ic+3,i)*spx(js,ic+3,j) ! Ve
               enddo
            endif
         enddo
         v3a=v3a+sum(v2tmp1*d2b(:,:,ij))
         if (dov3.eq.2) then
            v3axdel=v3axdel+sum(v2xdel*d2b(:,:,ij))
            v3adeldel=v3adeldel+sum(v2deldel*d2b(:,:,ij))
            v3tm=v3tm+sum(v2tmp2*d2b(:,:,ij))
            v3d1=v3d1+sum(v2tmp3a*d2b(:,:,ij))
            v3d2=v3d2+sum(v2tmp3b*d2b(:,:,ij))
            v3e=v3e+sum(v2tmp4*d2b(:,:,ij))
         endif
      enddo
   enddo
   v3asc=4.0_r8*a2psc*a2p3b*v3a  !anticommutator (propagator)
   v3a=4.0_r8*a2p3b*v3a  !anticommutator
   v3axdelsc=4.0_r8*a2p3b*a2pxdsc*v3axdel  !Vd,cc3 (propagator)
   v3axdel=4.0_r8*a2p3b*v3axdel  !Vd,cc3
   v3adeldelsc=4.0_r8*a2p3b*a2pddsc*v3adeldel !Ve,cc3 (propagator)
   v3adeldel=4.0_r8*a2p3b*v3adeldel !Ve,cc3
   v3tm=a2s3b*v3tm  !Tucson-Melbourne
   v3d1=avd/dfac*v3d1  !Vd,1
   v3d2=avd/dfac*v3d2  !Vd,2
   v3e=ave/dfac**2*v3e  !Ve
   xpic=0.0_r8
   v3c=czero
   v3cxd=czero
   v3cdd=czero
   if (a2p3b.ne.0.0_r8.and.acfac.ne.0.0_r8) then
      do i=1,npart-1
         do j=i+1,npart
            xpic(:,i,:,j)=xpi(:,i,:,j)
            xpic(1,i,1,j)=xpic(1,i,1,j)-delta(i,j)
            xpic(2,i,2,j)=xpic(2,i,2,j)-delta(i,j)
            xpic(3,i,3,j)=xpic(3,i,3,j)-delta(i,j)
            xpic(:,j,:,i)=xpic(:,i,:,j)
         enddo
      enddo
      ijk=0
      do i=1,npart-2
         do j=i+1,npart-1
            do k=j+1,npart
               ijk=ijk+1
               if (dotrip(ijk)) then
                  v3tmp1=czero
                  v3tmp2=czero
                  v3tmp3=czero
                  do ic=1,3
                     ic3=(ic+1)*3
                     do jc=1,3
                        jc3=(jc+1)*3
                        do kc=1,3
                           kc3=(kc+1)*3
                           vfacxx=xpi(ic,i,levi(1,jc),j)*xpi(levi(2,jc),j,kc,k) &
                                 -xpi(ic,i,levi(2,jc),j)*xpi(levi(1,jc),j,kc,k) &
                                 +xpi(jc,j,levi(1,kc),k)*xpi(levi(2,kc),k,ic,i) &
                                 -xpi(jc,j,levi(2,kc),k)*xpi(levi(1,kc),k,ic,i) &
                                 +xpi(kc,k,levi(1,ic),i)*xpi(levi(2,ic),i,jc,j) &
                                 -xpi(kc,k,levi(2,ic),i)*xpi(levi(1,ic),i,jc,j)
                           if (dov3.eq.2) then
                              vfacxd=xpi(ic,i,levi(1,jc),j)*xd(levi(2,jc),j,kc,k) &
                                    -xpi(ic,i,levi(2,jc),j)*xd(levi(1,jc),j,kc,k) &
                                    +xpi(jc,j,levi(1,kc),k)*xd(levi(2,kc),k,ic,i) &
                                    -xpi(jc,j,levi(2,kc),k)*xd(levi(1,kc),k,ic,i) &
                                    +xpi(kc,k,levi(1,ic),i)*xd(levi(2,ic),i,jc,j) &
                                    -xpi(kc,k,levi(2,ic),i)*xd(levi(1,ic),i,jc,j) &
                                    +xd(ic,i,levi(1,jc),j)*xpi(levi(2,jc),j,kc,k) &
                                    -xd(ic,i,levi(2,jc),j)*xpi(levi(1,jc),j,kc,k) &
                                    +xd(jc,j,levi(1,kc),k)*xpi(levi(2,kc),k,ic,i) &
                                    -xd(jc,j,levi(2,kc),k)*xpi(levi(1,kc),k,ic,i) &
                                    +xd(kc,k,levi(1,ic),i)*xpi(levi(2,ic),i,jc,j) &
                                    -xd(kc,k,levi(2,ic),i)*xpi(levi(1,ic),i,jc,j)
                              vfacdd=xd(ic,i,levi(1,jc),j)*xd(levi(2,jc),j,kc,k) &
                                    -xd(ic,i,levi(2,jc),j)*xd(levi(1,jc),j,kc,k) &
                                    +xd(jc,j,levi(1,kc),k)*xd(levi(2,kc),k,ic,i) &
                                    -xd(jc,j,levi(2,kc),k)*xd(levi(1,kc),k,ic,i) &
                                    +xd(kc,k,levi(1,ic),i)*xd(levi(2,ic),i,jc,j) &
                                    -xd(kc,k,levi(2,ic),i)*xd(levi(1,ic),i,jc,j)
                           endif
                           do it=1,3
                              do ks=1,4
                                 do js=1,4
                                    spxfac=spx(:,ic3+it,i) &
                                        *(spx(js,jc3+levi(1,it),j)*spx(ks,kc3+levi(2,it),k) &
                                         -spx(js,jc3+levi(2,it),j)*spx(ks,kc3+levi(1,it),k))
                                    v3tmp1(:,js,ks)=v3tmp1(:,js,ks)+spxfac*vfacxx
                                    if (dov3.eq.2) then
                                       v3tmp2(:,js,ks)=v3tmp2(:,js,ks)+spxfac*vfacxd
                                       v3tmp3(:,js,ks)=v3tmp3(:,js,ks)+spxfac*vfacdd
                                    endif
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  v3c=v3c+sum(v3tmp1(:,:,:)*d3b(:,:,:,ijk))
                  v3cxd=v3cxd+sum(v3tmp2(:,:,:)*d3b(:,:,:,ijk))
                  v3cdd=v3cdd+sum(v3tmp3(:,:,:)*d3b(:,:,:,ijk))
               endif
            enddo
         enddo
      enddo
      v3c=4.0_r8*a2p3b*acfac*v3c
      v3cxd=4.0_r8*a2p3b*acfac*v3cxd
      v3cdd=4.0_r8*a2p3b*acfac*v3cdd
   endif
   end subroutine v3bval

   subroutine setxspx(w,rcut,acut,noprot)
   use stack
   use random
   type (walker) :: w
   real(kind=r8) :: dxij(3),dxjk(3),dxik(3),rn(1),r,rcut,acut,prob
   integer(kind=i4) :: ijk,i,j,k
   logical :: noprot
   ijk=0
   dotrip=.false.
   if (dov3.eq.0.or.noprot) then
      probinvijk=1.0_r8
      return
   endif
!  dotrip=.true.
!  probinvijk=1.0_r8
! return
   call setrn(w%irn)
   do i=1,npart-2
      do j=i+1,npart-1
         dxij=w%x(:,i)-w%x(:,j)
         dxij=dxij-el*nint(dxij*eli)
         do k=j+1,npart
            ijk=ijk+1
            dxjk=w%x(:,j)-w%x(:,k)
            dxik=w%x(:,i)-w%x(:,k)
            dxjk=dxjk-el*nint(dxjk*eli)
            dxik=dxik-el*nint(dxik*eli)
            r=sqrt(sum(dxij**2)+sum(dxjk**2)+sum(dxik**2))
            if (r.le.rcut) then
               prob=1.0_r8
            else
               prob=exp(-acut*(r-rcut))
            endif
            rn=randn(1)
            if (rn(1).lt.prob) then
               dotrip(ijk)=.true.
               probinvijk(ijk)=1.0_r8/prob
            endif
         enddo
      enddo
   enddo
   call savern(w%irn)
   end subroutine setxspx
end module correlator
