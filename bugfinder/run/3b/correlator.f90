module correlator
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   real(kind=r8), private, parameter :: tiny=1e-5_r8
   integer(kind=i4), private, save :: npart,npair,ntrip
   complex(kind=r8), private, save, allocatable :: sxmallz(:,:,:),sp(:,:) &
      ,spx(:,:,:),sxmall(:,:,:)
   complex(kind=r8), private, save, allocatable :: vsz(:,:,:),vcoulz(:,:,:)
   real(kind=r8), private, save, allocatable :: fstvec(:,:,:),fstval(:,:) &
      ,ft(:),vsvec(:,:,:),vsval(:,:),vstvec(:,:,:),vstval(:,:)
   logical, private, save, allocatable :: dofst(:),doft(:) &
      ,dovs(:)
   logical, private, save :: docoul
   integer(kind=i4), private, parameter :: levi(2,3) = &
      reshape((/2,3, 3,1, 1,2/),(/2,3/))
!   private :: opmult
   complex(kind=r8), private, save, allocatable :: tau(:,:,:),tauz(:,:,:,:,:) &
      ,sigma(:,:,:),sigmaz(:,:,:,:,:),sigtau(:,:,:,:,:),sigtauz(:,:,:,:,:,:,:)
   real(kind=r8), private, save, allocatable :: v2(:),v3(:),v4(:),v5(:,:,:) &
      ,v6(:,:,:)
   logical, private, save :: isdiag=.true. ! calculate only tau.tau operators, not other components
   logical, private, save :: opcalc=.true. ! store all the spin/isospin operators
   complex(kind=r8), private, save, allocatable :: sigma1(:,:),tau1(:,:),sigtau1(:,:,:)
contains
   subroutine initcormod(npartin,docoulin)
   integer(kind=i4) :: npartin
   logical :: docoulin
   docoul=docoulin
   npart=npartin
   npair=(npart*(npart-1))/2
   ntrip=(npart*(npart-1)*(npart-2))/6
   if (allocated(sxmallz)) then
      deallocate(sxmallz,sp,spx,sxmall,vstvec,vstval)
      deallocate(vsvec,vsval,fstvec,fstval,ft)
      deallocate(dofst,doft,dovs)
      deallocate(vsz,vcoulz,v2,v3,v4,v5,v6)
      if (opcalc) deallocate(tau,tauz,sigma,sigmaz,sigtau,sigtauz)
   endif
   allocate(sxmallz(npart,4,npart),sp(4,npart),spx(4,15,npart))
   allocate(sxmall(npart,15,npart),vstvec(3,3,npair),vstval(3,npair))
   allocate(vsvec(3,3,npair),vsval(3,npair),fstvec(3,3,npair),fstval(3,npair))
   allocate(ft(npair),dofst(npair),doft(npair),dovs(npair))
   allocate(vsz(4,4,npair),vcoulz(4,4,npair))
   if (opcalc) then
      allocate(tau(3,3,npair),tauz(4,4,3,3,npair))
      allocate(sigma(3,3,npair),sigmaz(4,4,3,3,npair))
      allocate(sigtau(3,3,3,3,npair),sigtauz(4,4,3,3,3,3,npair))
   endif
   allocate(v2(npair),v3(npair),v4(npair))
   allocate(v5(3,3,npair),v6(3,3,npair))
   allocate(sigma1(3,npart),tau1(3,npart),sigtau1(3,3,npart))
   end subroutine initcormod

   subroutine calfop(ftauin,fsigtauin,cut)
   use matrixmod
   real(kind=r8) :: ftauin(:,:),fsigtauin(:,:,:,:),cut
   integer(kind=i4) :: i,j,ij
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         dofst(ij)=maxval(abs(fsigtauin(:,i,:,j))).gt.cut
         doft(ij)=abs(ftauin(i,j)).gt.cut
         if (dofst(ij)) then
            fstvec(:,:,ij)=fsigtauin(:,i,:,j)
            call eigenrs(fstvec(:,:,ij),fstval(:,ij),3)
         else
            fstval(:,ij)=0.0_r8
         endif
         if (doft(ij)) then
            ft(ij)=ftauin(i,j)
         else
            ft(ij)=0.0_r8
         endif
      enddo
   enddo
   end subroutine calfop

   subroutine calvop(spin,v2in,v3in,v4in,v5in,v6in,vcoulin,cut)
   use matrixmod
   complex(kind=r8) :: spin(:,:)
   real(kind=r8) :: vtau(npair),vsig(3,3,npair),vsigtau(3,3,npair),cut
   real(kind=r8) :: v2in(:,:),v3in(:,:),v4in(:,:),v5in(:,:,:,:),v6in(:,:,:,:)
   real(kind=r8) :: vcoulin(:,:)
   integer(kind=i4) :: i,j,ij,ic,jc,js,it,jt
   spx=opmult(spin)
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         v2(ij)=v2in(i,j)
         v3(ij)=v3in(i,j)
         v4(ij)=v4in(i,j)
         v5(:,:,ij)=v5in(:,i,:,j)
         v6(:,:,ij)=v6in(:,i,:,j)
      enddo
   enddo
   vtau=v2
   vsig=0.0_r8
   vsigtau=0.0_r8
   do ic=1,3
      vsig(ic,ic,:)=v3
      vsigtau(ic,ic,:)=v4
   enddo
   vsig=vsig+v5
   vsigtau=vsigtau+v6
   vsz=czero
   vcoulz=czero
   if (opcalc) then
      tauz=czero
      sigmaz=czero
      sigtauz=czero
   endif
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         if (.not.opcalc) then
            dovs(ij)=maxval(abs(vsig(:,:,ij))).gt.cut.or. &
               maxval(abs(vsigtau(:,:,ij))).gt.cut.or.abs(vtau(ij)).gt.cut
         endif
         if (dovs(ij).or.opcalc) then
            do js=1,4
               do ic=1,3
                  do jc=1,3
                     if (opcalc) then
                        sigmaz(:,js,ic,jc,ij)=sigmaz(:,js,ic,jc,ij)+spx(:,ic,i)*spx(js,jc,j)
                     else
                        vsz(:,js,ij)=vsz(:,js,ij) &
                           +spx(:,ic,i)*spx(js,jc,j)*vsig(ic,jc,ij)
                     endif
                  enddo
               enddo
            enddo
            do js=1,4
               do ic=1,3
                  do jc=1,3
                     if (opcalc) then
                        do it=1,3
                           if (isdiag) then
                              sigtauz(:,js,ic,jc,it,it,ij)=sigtauz(:,js,ic,jc,it,it,ij) &
                                 +  spx(:,3*(ic-1)+it+6,i)*spx(js,3*(jc-1)+it+6,j)
                           else
                              do jt=1,3
                                 sigtauz(:,js,ic,jc,it,jt,ij)=sigtauz(:,js,ic,jc,it,jt,ij) &
                                    +  spx(:,3*(ic-1)+it+6,i)*spx(js,3*(jc-1)+jt+6,j)
                              enddo
                           endif
                        enddo
                     else
                        vsz(:,js,ij)=vsz(:,js,ij) &
                           +(spx(:,3*ic+4,i)*spx(js,3*jc+4,j) &
                            +spx(:,3*ic+5,i)*spx(js,3*jc+5,j) &
                            +spx(:,3*ic+6,i)*spx(js,3*jc+6,j))*vsigtau(ic,jc,ij)
                     endif
                  enddo
               enddo
            enddo
            do js=1,4
               if (opcalc) then
                  do it=1,3
                     if (isdiag) then
                        tauz(:,js,it,it,ij)=tauz(:,js,it,it,ij) &
                              +spx(:,it+3,i)*spx(js,it+3,j)
                     else
                        do jt=1,3
                           tauz(:,js,it,jt,ij)=tauz(:,js,it,jt,ij) &
                                 +spx(:,it+3,i)*spx(js,jt+3,j)
                        enddo
                     endif
                  enddo
               else
                  vsz(:,js,ij)=vsz(:,js,ij)+(spx(:,4,i)*spx(js,4,j) &
                      +spx(:,5,i)*spx(js,5,j)+spx(:,6,i)*spx(js,6,j))*vtau(ij)
               endif
            enddo
         endif
         if (docoul) then
            do js=1,2
               vcoulz(1:2,js,ij)=vcoulz(1:2,js,ij) &
                  +spin(1:2,i)*spin(js,j)*vcoulin(i,j)
            enddo
         endif
      enddo
   enddo
   end subroutine calvop

   subroutine cordet(detrat,sxmallzin,spin)
! calfop must be called first
   complex(kind=r8) :: detrat,sxmallzin(:,:,:),spin(:,:)
   complex(kind=r8) :: ctmp1,d1,d2,d3,d4
   integer(kind=i4) :: i,j,ij,is,it
   sp=spin
   sxmallz=sxmallzin
   spx=opmult(sp)
   do i=1,npart
      sxmall(:,:,i)=matmul(sxmallz(:,:,i),spx(:,:,i))
   enddo
   ij=0
   detrat=czero
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         if (dofst(ij)) then
!dir$ ivdep
            do is=1,3
               ctmp1=czero
!dir$ ivdep
               do it=1,3
                  d1=sum(sxmall(i,it+6:it+12:3,i)*fstvec(:,is,ij))
                  d2=sum(sxmall(j,it+6:it+12:3,j)*fstvec(:,is,ij))
                  d3=sum(sxmall(i,it+6:it+12:3,j)*fstvec(:,is,ij))
                  d4=sum(sxmall(j,it+6:it+12:3,i)*fstvec(:,is,ij))
                  ctmp1=ctmp1+d1*d2-d3*d4
               enddo
               detrat=detrat+ctmp1*fstval(is,ij)
            enddo
         endif
         if (doft(ij)) then
            ctmp1=czero
!dir$ ivdep
            do it=4,6
               d1=sxmall(i,it,i)
               d2=sxmall(j,it,j)
               d3=sxmall(i,it,j)
               d4=sxmall(j,it,i)
               ctmp1=ctmp1+d1*d2-d3*d4
            enddo
            detrat=detrat+ctmp1*ft(ij)
         endif
      enddo
   enddo
   detrat=detrat+cone
   end subroutine cordet

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

   subroutine v6val(cvs,cvcoul,d2b)
! calvop must be called first.
   complex(kind=r8), intent(out) :: cvs,cvcoul
   complex(kind=r8), intent(in) :: d2b(:,:,:)
   integer(kind=i4) :: i,j,ij
   cvs=czero
   cvcoul=czero
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         if (dovs(ij)) cvs=cvs+sum(d2b(:,:,ij)*vsz(:,:,ij))
         if (docoul) cvcoul=cvcoul+sum(d2b(:,:,ij)*vcoulz(:,:,ij))
      enddo
   enddo
   end subroutine v6val

   subroutine g3bval(d3b,sxz,fij)
   complex(kind=r8), intent(inout):: d3b(4,4,4,ntrip)
   complex(kind=r8), intent(in):: sxz(4,npart,npart),fij
   integer(kind=i4) :: i,j,k,js,ks,ijk
   ijk=0
   do i=1,npart-2
      do j=i+1,npart-1
         do k=j+1,npart
            ijk=ijk+1
            do ks=1,4
               do js=1,4
                  d3b(:,js,ks,ijk)= d3b(:,js,ks,ijk)+fij*( &
                      sxz(:,i,i)*(sxz(js,j,j)*sxz(ks,k,k) &
                     -sxz(js,j,k)*sxz(ks,k,j))+sxz(:,i,j) &
                     *(sxz(js,j,k)*sxz(ks,k,i)-sxz(js,j,i)*sxz(ks,k,k)) &
                     +sxz(:,i,k)*(sxz(js,j,i)*sxz(ks,k,j) &
                                -sxz(js,j,j)*sxz(ks,k,i)))
               enddo
            enddo
         enddo
      enddo
   enddo
   end subroutine g3bval

   subroutine v6tot(cvs,cvcoul)
! calvop, calfop and cordet must be called first.
   complex(kind=r8) :: cvs(:),cvcoul
   complex(kind=r8) :: di(npart,12),dj(npart,12),fij(12)
   complex(kind=r8) :: sxi(4,npart,12),sxj(4,npart,12),detrati(12)
   complex(kind=r8) :: detrat(12),sxz(4,npart,npart)
   complex(kind=r8) :: sinvijz(4,npart,npart)
   complex(kind=r8) :: sx15(4,15,npart,npart)
   complex(kind=r8) :: d1b(4,npart),d2b(4,4,npair),d3b(4,4,4,ntrip)
   integer(kind=i4) :: i,j,ij,iop,m,n,ktau,kval
   logical :: doops(12)
   real(kind=r8) :: vec(3,3)
   sxz=reshape(transpose(reshape(sxmallz,(/npart,4*npart/))),shape(sxz))
   d1b=czero
   d2b=czero
   call g1bval(d1b,sxz,cone)
   call g2bval(d2b,sxz,cone)
   d3b=czero
   call g3bval(d3b,sxz,cone)
   do i=1,npart
      sx15(:,:,:,i)=conjg(opmult(conjg(sxz(:,i,:))))
   enddo
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         fij=czero !if this remains zero, the rest is not calculated
         if ((.not.doft(ij)).and.(.not.dofst(ij))) cycle
         if (doft(ij)) fij(1:3)=ft(ij)
         if (dofst(ij)) then
            vec=fstvec(:,:,ij)
            fij(4:6)=fstval(1,ij)
            fij(7:9)=fstval(2,ij)
            fij(10:12)=fstval(3,ij)
         endif
         doops(1:3)=doft(ij)
         doops(4:12)=dofst(ij)
         if (doft(ij)) then
            do ktau=1,3
               sxi(:,:,ktau)=sx15(:,ktau+3,:,i)
               sxj(:,:,ktau)=sx15(:,ktau+3,:,j)
            enddo
            di(:,1:3)=sxi(1,:,1:3)*sp(1,i)+sxi(2,:,1:3)*sp(2,i) &
                +sxi(3,:,1:3)*sp(3,i)+sxi(4,:,1:3)*sp(4,i)
            dj(:,1:3)=sxj(1,:,1:3)*sp(1,j)+sxj(2,:,1:3)*sp(2,j) &
                +sxj(3,:,1:3)*sp(3,j)+sxj(4,:,1:3)*sp(4,j)
            detrat(1:3)=di(i,1:3)*dj(j,1:3)-di(j,1:3)*dj(i,1:3)
            detrati(1:3)=cone/detrat(1:3)
            fij(1:3)=detrat(1:3)*fij(1:3)
         endif
         if (dofst(ij)) then
            do concurrent(m=1:npart,ktau=1:3,kval=1:3)
               sxi(:,m,ktau+3*kval)=sx15(:,ktau+6,m,i)*vec(1,kval) &
                 +sx15(:,ktau+9,m,i)*vec(2,kval) &
                 +sx15(:,ktau+12,m,i)*vec(3,kval)
               sxj(:,m,ktau+3*kval)=sx15(:,ktau+6,m,j)*vec(1,kval) &
                 +sx15(:,ktau+9,m,j)*vec(2,kval) &
                 +sx15(:,ktau+12,m,j)*vec(3,kval)
               di(m,ktau+3*kval)=sum(sxi(:,m,ktau+3*kval)*sp(:,i))
               dj(m,ktau+3*kval)=sum(sxj(:,m,ktau+3*kval)*sp(:,j))
            enddo
            detrat(4:12)=di(i,4:12)*dj(j,4:12)-di(j,4:12)*dj(i,4:12)
            detrati(4:12)=cone/detrat(4:12)
            fij(4:12)=detrat(4:12)*fij(4:12)
         endif
         do iop=1,12
         if (.not.doops(iop)) cycle
!dir$ ivdep
            do m=1,npart
               sinvijz(:,:,m)= &
                  (di(i,iop)*dj(j,iop)-dj(i,iop)*di(j,iop))*sxz(:,:,m) &
                 -(di(i,iop)*dj(m,iop)-dj(i,iop)*di(m,iop))*sxz(:,:,j) &
                 +(di(j,iop)*dj(m,iop)-dj(j,iop)*di(m,iop))*sxz(:,:,i)
            enddo
            do m=1,npart
               sinvijz(:,i,m)= &
                  (di(i,iop)*dj(j,iop)-dj(i,iop)*di(j,iop))*sxi(:,m,iop) &
                 -(di(i,iop)*dj(m,iop)-dj(i,iop)*di(m,iop))*sxi(:,j,iop) &
                 +(di(j,iop)*dj(m,iop)-dj(j,iop)*di(m,iop))*sxi(:,i,iop)
               sinvijz(:,j,m)= &
                  (di(i,iop)*dj(j,iop)-dj(i,iop)*di(j,iop))*sxj(:,m,iop) &
                 -(di(i,iop)*dj(m,iop)-dj(i,iop)*di(m,iop))*sxj(:,j,iop) &
                 +(di(j,iop)*dj(m,iop)-dj(j,iop)*di(m,iop))*sxj(:,i,iop)
            enddo
            do n=1,npart
               sinvijz(:,n,i)=dj(j,iop)*sxz(:,n,i)-dj(i,iop)*sxz(:,n,j)
               sinvijz(:,n,j)=di(i,iop)*sxz(:,n,j)-di(j,iop)*sxz(:,n,i)
            enddo
            sinvijz(:,i,i)=dj(j,iop)*sxi(:,i,iop)-dj(i,iop)*sxi(:,j,iop)
            sinvijz(:,j,j)=di(i,iop)*sxj(:,j,iop)-di(j,iop)*sxj(:,i,iop)
            sinvijz(:,j,i)=dj(j,iop)*sxj(:,i,iop)-dj(i,iop)*sxj(:,j,iop)
            sinvijz(:,i,j)=di(i,iop)*sxi(:,j,iop)-di(j,iop)*sxi(:,i,iop)
            sinvijz(:,:,:)=detrati(iop)*sinvijz(:,:,:)
            call g1bval(d1b,sinvijz,fij(iop))
            call g2bval(d2b,sinvijz,fij(iop))
            call g3bval(d3b,sinvijz,fij(iop))
         enddo
      enddo
   enddo
   cvs=czero
   if (opcalc) then
      call opval(tau,sigma,sigtau,d2b,cvcoul)
      call calpot(cvs,v2,v3,v4,v5,v6)
   else
      call v6val(cvs(1),cvcoul,d2b)
   endif
   call op1val(d1b)
   end subroutine v6tot

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

   subroutine opval(ttau,ssigma,ssigtau,d2b,cvcoul)
! calvop must be called first.
   complex(kind=r8), intent(out) :: ttau(:,:,:),ssigma(:,:,:),ssigtau(:,:,:,:,:)
   complex(kind=r8), intent(in) :: d2b(:,:,:)
   complex(kind=r8) :: cvcoul
   integer(kind=i4) :: i,j,ic,jc,it,jt,ij
   ttau=czero
   ssigma=czero
   ssigtau=czero
   cvcoul=czero
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         do ic=1,3
            do jc=1,3
               ttau(ic,jc,ij)=ttau(ic,jc,ij)+sum(d2b(:,:,ij)*tauz(:,:,ic,jc,ij))
               ssigma(ic,jc,ij)=ssigma(ic,jc,ij)+sum(d2b(:,:,ij)*sigmaz(:,:,ic,jc,ij))
               do it=1,3
                  if (isdiag) then
                     ssigtau(ic,jc,it,it,ij)=ssigtau(ic,jc,it,it,ij) &
                              +sum(d2b(:,:,ij)*sigtauz(:,:,ic,jc,it,it,ij))
                  else
                     do jt=1,3
                        ssigtau(ic,jc,it,jt,ij)=ssigtau(ic,jc,it,jt,ij) &
                                 +sum(d2b(:,:,ij)*sigtauz(:,:,ic,jc,it,jt,ij))
                     enddo
                  endif
               enddo
            enddo
         enddo
         if (docoul) cvcoul=cvcoul+sum(d2b(:,:,ij)*vcoulz(:,:,ij))
      enddo
   enddo
   end subroutine opval

   subroutine calpot(cvs,v2,v3,v4,v5,v6)
   real(kind=r8) :: v2(:),v3(:),v4(:),v5(:,:,:),v6(:,:,:)
   complex(kind=r8) :: cvs(:)
   integer(kind=i4) :: it,is,js,ij
   cvs=czero
   do ij=1,npair
      do it=1,3
         cvs(2)=cvs(2)+tau(it,it,ij)*v2(ij)
         do is=1,3
            cvs(4)=cvs(4)+sigtau(is,is,it,it,ij)*v4(ij)
         enddo
      enddo
      do is=1,3
         cvs(3)=cvs(3)+sigma(is,is,ij)*v3(ij)
         do js=1,3
            cvs(5)=cvs(5)+sigma(is,js,ij)*v5(is,js,ij)
            do it=1,3
               cvs(6)=cvs(6)+sigtau(is,js,it,it,ij)*v6(is,js,ij)
            enddo
         enddo
      enddo
   enddo
   end subroutine calpot

   subroutine g1bval(d1b,sxz,fij)
   complex(kind=r8), intent(inout) :: d1b(:,:)
   complex(kind=r8), intent(in) :: sxz(:,:,:),fij
   integer(kind=i4) :: i
   do i=1,npart
      d1b(:,i)=d1b(:,i)+fij*sxz(:,i,i)
   enddo
   end subroutine g1bval

   subroutine op1val(d1b)
   complex(kind=r8), intent(in) :: d1b(:,:)
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

   subroutine getop1(sig1in,tau1in,sigtau1in)
   complex(kind=r8) :: sig1in(:,:),tau1in(:,:),sigtau1in(:,:,:)
   sig1in=sigma1
   tau1in=tau1
   sigtau1in=sigtau1
   end subroutine getop1

   subroutine getop2(sig2in,tau2in,sigtau2in)
   complex(Kind=r8), dimension(:,:,:) :: sig2in,tau2in
   complex(kind=r8), dimension(:,:,:,:,:) :: sigtau2in
   sig2in=sigma
   tau2in=tau
   sigtau2in=sigtau
   end subroutine getop2

end module correlator
