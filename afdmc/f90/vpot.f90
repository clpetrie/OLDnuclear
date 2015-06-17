module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   real(kind=r8), private, parameter :: tiny=1e-8_r8
   integer(kind=i4), private, save :: npart
contains
   subroutine initvpot(npartin)
   integer(kind=i4) :: npartin
   npart=npartin
   end subroutine initvpot
 
   subroutine v6potcalc(cvs,cvst,cvt,vcoulomb,sxmall,sxxmall,vasig,vasigtau,vatau,fasigtau,fatau,docoul,vcmat)
   use v6pot
   use matrixmod
   complex(kind=r8) :: cvs,cvst,cvt
   complex(kind=r8) :: sxmall(:,:,:),sxxmall(:,:,:,:)
   real(kind=r8) :: fasigtau(:,:,:,:),fatau(:,:)
   real(kind=r8) :: vasig(:,:,:,:),vasigtau(:,:,:,:),vatau(:,:)
   complex(kind=r8) :: d1,d2,d3,d4
   integer(kind=i4) :: ic1,i1,ic2,i2,i,j,ic,it,jc1,jc2,jc,k,l,kt,kc,lc1,lc,kc1
   complex(kind=r8) :: M3v(3,3),M2v(2,2),M4(4,4),M3(3,3),M2(2,2),dM4,dM3,dM2
   complex(kind=r8) :: Ma(2,2),Mb(2,2),Mc(2,2),Md(2,2),dMa,dMc,dm
   complex(kind=r8) :: vcoulomb(:),vcoul,vcoulc
   real(kind=r8) :: vcmat(:,:)
   logical :: docoul
   vcoul=czero
   vcoulc=czero
   if (docoul) then
      do i=1,npart-1
         do j=i+1,npart
            vcoulc=vcoulc+vcmat(i,j)
         enddo
      enddo
   endif
   cvs=czero
   cvst=czero
   cvt=czero
   do ic1=1,3
      do i1=1,npart-1
         d1=sxmall(i1,ic1,i1)
         do ic2=1,3
            do i2=i1+1,npart
               d2=sxmall(i2,ic2,i2)
               d3=sxmall(i1,ic2,i2)
               d4=sxmall(i2,ic1,i1)
               cvs=cvs+(d1*d2-d3*d4)*vasig(ic1,i1,ic2,i2)
            enddo
         enddo
      enddo
   enddo
   do it=1,3
      do ic1=1,3
         jc1=3*(ic1-1)+6+it
         do i1=1,npart-1
            d1=sxmall(i1,jc1,i1)
            do ic2=1,3
               jc2=3*(ic2-1)+6+it
               do i2=i1+1,npart
                  d2=sxmall(i2,jc2,i2)
                  d3=sxmall(i1,jc2,i2)
                  d4=sxmall(i2,jc1,i1)
                  cvst=cvst+(d1*d2-d3*d4)*vasigtau(ic1,i1,ic2,i2)
               enddo
            enddo
         enddo
      enddo
   enddo
   do it=4,6
      do i1=1,npart-1
         d1=sxmall(i1,it,i1)
         do i2=i1+1,npart
            d2=sxmall(i2,it,i2)
            d3=sxmall(i1,it,i2)
            d4=sxmall(i2,it,i1)
            cvt=cvt+(d1*d2-d3*d4)*vatau(i1,i2)
            if (docoul.and.it.eq.6) vcoul=vcoul+(d1*d2-d3*d4)*vcmat(i1,i2)
         enddo
      enddo
   enddo
! v_ij f_jk i < j != k < l
   do i=1,npart-1
     do j=i+1,npart
       do it=1,3
         do ic=1,3
           ic1=3*(ic-1)+6+it
           M4(1,1)=sxmall(i,ic1,i)     
           do jc=1,3
             if (abs(fasigtau(ic,i,jc,j)).lt.tiny) cycle
             jc1=3*(jc-1)+6+it
             M4(1,2)=sxmall(i,jc1,j)
             M4(2,2)=sxmall(j,jc1,j)
             M4(2,1)=sxmall(j,ic1,i)
             Ma=M4(1:2,1:2)
             call cmatinv2(Ma,dMa)
             do k=1,npart-1
               if(k.eq.i.or.k.eq.j) cycle
               M4(3,1)=sxmall(k,ic1,i)
               M4(3,2)=sxmall(k,jc1,j)
               do l=k+1,npart
                 if(l.eq.i.or.l.eq.j) cycle
                 M4(4,1)=sxmall(l,ic1,i)
                 M4(4,2)=sxmall(l,jc1,j)
                 Mc=M4(3:4,1:2)
                 Mc=matmul(Mc,Ma)
                 do kc=1,3
                   M4(1,3)=sxmall(i,kc,k)
                   M4(2,3)=sxmall(j,kc,k)
                   M4(3,3)=sxmall(k,kc,k)
                   M4(4,3)=sxmall(l,kc,k)
                   do lc=1,3
                     if (abs(vasig(kc,k,lc,l)).lt.tiny) cycle
                     M4(1,4)=sxmall(i,lc,l)
                     M4(2,4)=sxmall(j,lc,l)
                     M4(3,4)=sxmall(k,lc,l)
                     M4(4,4)=sxmall(l,lc,l)
                     if(dMa.ne.0.) then
                       Mb=M4(1:2,3:4)
                       Md=M4(3:4,3:4)
                       call cmatdet2(Md-matmul(Mc,Mb),dMc)
                       dm=dMa*dMc
                     else
                       call cmatdet4(M4,dM4)
                       dm=dM4
                     endif
                     cvs=cvs+dm*vasig(kc,k,lc,l)*fasigtau(ic,i,jc,j)
                   enddo
                 enddo
                 do kt=1,3
                   do kc=1,3
                     kc1=3*(kc-1)+6+kt
                     M4(1,3)=sxmall(i,kc1,k)
                     M4(2,3)=sxmall(j,kc1,k)
                     M4(3,3)=sxmall(k,kc1,k)
                     M4(4,3)=sxmall(l,kc1,k)
                     do lc=1,3
                       if (abs(vasigtau(kc,k,lc,l)).lt.tiny) cycle
                       lc1=3*(lc-1)+6+kt
                       M4(1,4)=sxmall(i,lc1,l)
                       M4(2,4)=sxmall(j,lc1,l)
                       M4(3,4)=sxmall(k,lc1,l)
                       M4(4,4)=sxmall(l,lc1,l)
                       if(dMa.ne.0.) then
                         Mb=M4(1:2,3:4)
                         Md=M4(3:4,3:4)
                         call cmatdet2(Md-matmul(Mc,Mb),dMc)
                         dm=dMa*dMc
                       else
                         call cmatdet4(M4,dM4)
                         dm=dM4
                       endif
                       cvst=cvst+dm*vasigtau(kc,k,lc,l)*fasigtau(ic,i,jc,j)
                     enddo
                   enddo
                 enddo
                 do kt=4,6
                   if (abs(vatau(k,l)).lt.tiny) cycle
                   M4(1,3)=sxmall(i,kt,k)
                   M4(2,3)=sxmall(j,kt,k)
                   M4(3,3)=sxmall(k,kt,k)
                   M4(4,3)=sxmall(l,kt,k)
                   M4(1,4)=sxmall(i,kt,l)
                   M4(2,4)=sxmall(j,kt,l)
                   M4(3,4)=sxmall(k,kt,l)
                   M4(4,4)=sxmall(l,kt,l)
                   if(dMa.ne.0.) then
                     Mb=M4(1:2,3:4)
                     Md=M4(3:4,3:4)
                     call cmatdet2(Md-matmul(Mc,Mb),dMc)
                     dm=dMa*dMc
                   else
                     call cmatdet4(M4,dM4)
                     dm=dM4
                   endif
                   cvt=cvt+dm*vatau(k,l)*fasigtau(ic,i,jc,j)
                   if (docoul.and.kt.eq.6) vcoul=vcoul+dm*vcmat(k,l)*fasigtau(ic,i,jc,j)
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo

! f_ij v_kj (j=l)
   do i=1,npart
    do it=1,3
     do ic=1,3
       ic1=3*(ic-1)+6+it
         M3v(1,1)=sxmall(i,ic1,i)
         do jc=1,3
           jc1=3*(jc-1)+6+it
           do j=1,npart
             if (abs(fasigtau(ic,i,jc,j)).lt.tiny) cycle
             if(j.eq.i) cycle
             M3v(2,1)=sxmall(j,ic1,i)
             do k=1,npart
               if(k.eq.i.or.k.eq.j) cycle
               M3v(3,1)=sxmall(k,ic1,i)
               do kc=1,3
                 M3v(1,3)=sxmall(i,kc,k)
                 M3v(2,3)=sxmall(j,kc,k)
                 M3v(3,3)=sxmall(k,kc,k)
                 do lc=1,3
                   if (abs(vasig(kc,k,lc,j)).lt.tiny) cycle
                   M3v(1,2)=sxxmall(i,jc1,lc,j)
                   M3v(2,2)=sxxmall(j,jc1,lc,j)
                   M3v(3,2)=sxxmall(k,jc1,lc,j)
                   call cmatdet3(M3v,dM3)
                   cvs=cvs+dM3*vasig(kc,k,lc,j)*fasigtau(ic,i,jc,j)
                 enddo
               enddo
               do kt=1,3
                 do kc=1,3
                   kc1=3*(kc-1)+6+kt
                   M3v(1,3)=sxmall(i,kc1,k)
                   M3v(2,3)=sxmall(j,kc1,k)
                   M3v(3,3)=sxmall(k,kc1,k)
                   do lc=1,3
                     if (abs(vasigtau(kc,k,lc,j)).lt.tiny) cycle
                     lc1=3*(lc-1)+6+kt
                     M3v(1,2)=sxxmall(i,jc1,lc1,j)
                     M3v(2,2)=sxxmall(j,jc1,lc1,j)
                     M3v(3,2)=sxxmall(k,jc1,lc1,j)
                     call cmatdet3(M3v,dM3)
                     cvst=cvst+dM3*vasigtau(kc,k,lc,j)*fasigtau(ic,i,jc,j)
                   enddo
                 enddo
               enddo
               do kt=4,6
                 if (abs(vatau(k,j)).lt.tiny) cycle
                 M3v(1,2)=sxxmall(i,jc1,kt,j)
                 M3v(1,3)=sxmall(i,kt,k)
                 M3v(2,2)=sxxmall(j,jc1,kt,j)
                 M3v(3,2)=sxxmall(k,jc1,kt,j)
                 M3v(2,3)=sxmall(j,kt,k)
                 M3v(3,3)=sxmall(k,kt,k)
                 call cmatdet3(M3v,dM3)
                 cvt=cvt+dM3*vatau(k,j)*fasigtau(ic,i,jc,j)
                 if (docoul.and.kt.eq.6) vcoul=vcoul+dM3*vcmat(k,j)*fasigtau(ic,i,jc,j)
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo

!f_ij v_ij (i=k,j=l)
   do it=1,3
     do ic=1,3
       ic1=3*(ic-1)+6+it
       do i=1,npart-1
         do jc=1,3
           jc1=3*(jc-1)+6+it
           do j=i+1,npart
             if (abs(fasigtau(ic,i,jc,j)).lt.tiny) cycle
             do kc=1,3
               M2v(1,1)=sxxmall(i,ic1,kc,i)
               M2v(2,1)=sxxmall(j,ic1,kc,i)
               do lc=1,3
                 if (abs(vasig(kc,i,lc,j)).lt.tiny) cycle
                 M2v(1,2)=sxxmall(i,jc1,lc,j)
                 M2v(2,2)=sxxmall(j,jc1,lc,j)
                 call cmatdet2(M2v,dM2)
                 cvs=cvs+dM2*vasig(kc,i,lc,j)*fasigtau(ic,i,jc,j)
               enddo
             enddo
             do kt=1,3
               do kc=1,3
                 kc1=3*(kc-1)+6+kt
                 M2v(1,1)=sxxmall(i,ic1,kc1,i)
                 M2v(2,1)=sxxmall(j,ic1,kc1,i)
                 do lc=1,3
                   if (abs(vasigtau(kc,i,lc,j)).lt.tiny) cycle
                   lc1=3*(lc-1)+6+kt
                   M2v(1,2)=sxxmall(i,jc1,lc1,j)
                   M2v(2,2)=sxxmall(j,jc1,lc1,j)
                   call cmatdet2(M2v,dM2)
                   cvst=cvst+dM2*vasigtau(kc,i,lc,j)*fasigtau(ic,i,jc,j)
                 enddo
               enddo
             enddo
             do kt=4,6
               if (abs(vatau(i,j)).lt.tiny) cycle
               M2v(1,1)=sxxmall(i,ic1,kt,i)
               M2v(2,1)=sxxmall(j,ic1,kt,i)
               M2v(1,2)=sxxmall(i,jc1,kt,j)
               M2v(2,2)=sxxmall(j,jc1,kt,j)
               call cmatdet2(M2v,dM2)
               cvt=cvt+dM2*vatau(i,j)*fasigtau(ic,i,jc,j)
               if (docoul.and.kt.eq.6) vcoul=vcoul+dM2*vcmat(i,j)*fasigtau(ic,i,jc,j)
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! tau correlation 
! v_ij f_kl
   do i=1,npart-1
     do it=4,6
       M4(1,1)=sxmall(i,it,i)     
       do j=i+1,npart
         if (abs(fatau(i,j)).lt.tiny) cycle
         M4(1,2)=sxmall(i,it,j)
         M4(2,2)=sxmall(j,it,j)
         M4(2,1)=sxmall(j,it,i)
         Ma=M4(1:2,1:2)
         call cmatinv2(Ma,dMa)
         do k=1,npart-1
           if(k.eq.i.or.k.eq.j) cycle
           M4(3,1)=sxmall(k,it,i)
           M4(3,2)=sxmall(k,it,j)
           do l=k+1,npart
             if(l.eq.i.or.l.eq.j) cycle
             M4(4,1)=sxmall(l,it,i)
             M4(4,2)=sxmall(l,it,j)
             Mc=M4(3:4,1:2)
             Mc=matmul(Mc,Ma)
             do kc=1,3  !sigma S
               M4(1,3)=sxmall(i,kc,k)
               M4(2,3)=sxmall(j,kc,k)
               M4(3,3)=sxmall(k,kc,k)
               M4(4,3)=sxmall(l,kc,k)
               do lc=1,3
                 if (abs(vasig(kc,k,lc,l)).lt.tiny) cycle
                 M4(1,4)=sxmall(i,lc,l)
                 M4(2,4)=sxmall(j,lc,l)
                 M4(3,4)=sxmall(k,lc,l)
                 M4(4,4)=sxmall(l,lc,l)
                 if(dMa.ne.0.) then
                   Mb=M4(1:2,3:4)
                   Md=M4(3:4,3:4)
                   call cmatdet2(Md-matmul(Mc,Mb),dMc)
                   cvs=cvs+dMa*dMc*vasig(kc,k,lc,l)*fatau(i,j)
                 else
                   call cmatdet4(M4,dM4)
                   cvs=cvs+dM4*vasig(kc,k,lc,l)*fatau(i,j)
                 endif
               enddo
             enddo
             do kt=1,3 !sigma-tau S-tau
               do kc=1,3
                 kc1=3*(kc-1)+6+kt
                 M4(1,3)=sxmall(i,kc1,k)
                 M4(2,3)=sxmall(j,kc1,k)
                 M4(3,3)=sxmall(k,kc1,k)
                 M4(4,3)=sxmall(l,kc1,k)
                 do lc=1,3
                   if (abs(vasigtau(kc,k,lc,l)).lt.tiny) cycle
                   lc1=3*(lc-1)+6+kt
                   M4(1,4)=sxmall(i,lc1,l)
                   M4(2,4)=sxmall(j,lc1,l)
                   M4(3,4)=sxmall(k,lc1,l)
                   M4(4,4)=sxmall(l,lc1,l)
                   if(dMa.ne.0.) then
                     Mb=M4(1:2,3:4)
                     Md=M4(3:4,3:4)
                     call cmatdet2(Md-matmul(Mc,Mb),dMc)
                     cvst=cvst+dMa*dMc*vasigtau(kc,k,lc,l)*fatau(i,j)
                   else
                     call cmatdet4(M4,dM4)
                     cvst=cvst+dM4*vasigtau(kc,k,lc,l)*fatau(i,j)
                   endif
                 enddo
               enddo
             enddo
             do kt=4,6 !tau
               if (abs(vatau(k,l)).lt.tiny) cycle
               M4(1,3)=sxmall(i,kt,k)
               M4(2,3)=sxmall(j,kt,k)
               M4(3,3)=sxmall(k,kt,k)
               M4(4,3)=sxmall(l,kt,k)
               M4(1,4)=sxmall(i,kt,l)
               M4(2,4)=sxmall(j,kt,l)
               M4(3,4)=sxmall(k,kt,l)
               M4(4,4)=sxmall(l,kt,l)
               if(dMa.ne.0.) then
                 Mb=M4(1:2,3:4)
                 Md=M4(3:4,3:4)
                 call cmatdet2(Md-matmul(Mc,Mb),dMc)
                 dm=dMa*dMc
               else
                 call cmatdet4(M4,dM4)
                 dm=dM4
               endif 
               cvt=cvt+dm*vatau(k,l)*fatau(i,j)
               if (docoul.and.kt.eq.6) vcoul=vcoul+dm*vcmat(k,l)*fatau(i,j)
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! f_ij v_kj (j=l)
   do i=1,npart
     do it=4,6
       M3v(1,1)=sxmall(i,it,i)
       do j=1,npart
         if(j.eq.i) cycle
         if (abs(fatau(i,j)).lt.tiny) cycle
         M3v(2,1)=sxmall(j,it,i)
         do k=1,npart
           if(k.eq.i.or.k.eq.j) cycle
           M3v(3,1)=sxmall(k,it,i)
           do kc=1,3 !sigma S
             M3v(1,3)=sxmall(i,kc,k)
             M3v(2,3)=sxmall(j,kc,k)
             M3v(3,3)=sxmall(k,kc,k)
             do lc=1,3
               if (abs(vasig(kc,k,lc,j)).lt.tiny) cycle
               M3v(1,2)=sxxmall(i,it,lc,j)
               M3v(2,2)=sxxmall(j,it,lc,j)
               M3v(3,2)=sxxmall(k,it,lc,j)
               call cmatdet3(M3v,dM3)
               cvs=cvs+dM3*vasig(kc,k,lc,j)*fatau(i,j)
             enddo
           enddo
           do kt=1,3 !sigma-tau S-tau
             do kc=1,3
               kc1=3*(kc-1)+6+kt
               M3v(1,3)=sxmall(i,kc1,k)
               M3v(2,3)=sxmall(j,kc1,k)
               M3v(3,3)=sxmall(k,kc1,k)
               do lc=1,3
                 if (abs(vasigtau(kc,k,lc,j)).lt.tiny) cycle
                 lc1=3*(lc-1)+6+kt
                 M3v(1,2)=sxxmall(i,it,lc1,j)
                 M3v(2,2)=sxxmall(j,it,lc1,j)
                 M3v(3,2)=sxxmall(k,it,lc1,j)
                 call cmatdet3(M3v,dM3)
                 cvst=cvst+dM3*vasigtau(kc,k,lc,j)*fatau(i,j)
               enddo 
             enddo
           enddo
           do kt=4,6 !tau
             if (abs(vatau(k,j)).lt.tiny) cycle
             M3v(1,2)=sxxmall(i,it,kt,j)
             M3v(1,3)=sxmall(i,kt,k)
             M3v(2,2)=sxxmall(j,it,kt,j)
             M3v(3,2)=sxxmall(k,it,kt,j)
             M3v(2,3)=sxmall(j,kt,k)
             M3v(3,3)=sxmall(k,kt,k)
             call cmatdet3(M3v,dM3)
             cvt=cvt+dM3*vatau(k,j)*fatau(i,j)
             if (docoul.and.kt.eq.6) vcoul=vcoul+dM3*vcmat(k,j)*fatau(i,j)
           enddo
         enddo
       enddo
     enddo
   enddo
!f_ij v_ij
   do it=4,6
     do i=1,npart-1
       do j=i+1,npart
         if (abs(fatau(i,j)).lt.tiny) cycle
         do kc=1,3 !sigma S
           M2v(1,1)=sxxmall(i,it,kc,i)
           M2v(2,1)=sxxmall(j,it,kc,i)
           do lc=1,3
             if (abs(vasig(kc,i,lc,j)).lt.tiny) cycle
             M2v(1,2)=sxxmall(i,it,lc,j)
             M2v(2,2)=sxxmall(j,it,lc,j)
             call cmatdet2(M2v,dM2)
             cvs=cvs+dM2*vasig(kc,i,lc,j)*fatau(i,j)
           enddo
         enddo
         do kt=1,3 !sigma-tau S-tau
           do kc=1,3
             kc1=3*(kc-1)+6+kt
             M2v(1,1)=sxxmall(i,it,kc1,i)
             M2v(2,1)=sxxmall(j,it,kc1,i)
             do lc=1,3
               if (abs(vasigtau(kc,i,lc,j)).lt.tiny) cycle
               lc1=3*(lc-1)+6+kt
               M2v(1,2)=sxxmall(i,it,lc1,j)
               M2v(2,2)=sxxmall(j,it,lc1,j)
               call cmatdet2(M2v,dM2)
               cvst=cvst+dM2*vasigtau(kc,i,lc,j)*fatau(i,j)
             enddo
           enddo
         enddo 
         do kt=4,6 !tau
           if (abs(vatau(i,j)).lt.tiny) cycle
           M2v(1,1)=sxxmall(i,it,kt,i)
           M2v(2,1)=sxxmall(j,it,kt,i)
           M2v(1,2)=sxxmall(i,it,kt,j)
           M2v(2,2)=sxxmall(j,it,kt,j)
           call cmatdet2(M2v,dM2)
           cvt=cvt+dM2*vatau(i,j)*fatau(i,j)
           if (docoul.and.kt.eq.6) vcoul=vcoul+dM2*vcmat(i,j)*fatau(i,j)
         enddo
       enddo
     enddo
   enddo
   vcoulomb(1)=vcoulc
   vcoulomb(2)=vcoul
   contains

   subroutine cmatdet2(a,det)
   complex (kind=r8) :: a(:,:),det
   det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   end subroutine cmatdet2

   subroutine cmatdet3(a,det)
   complex (kind=r8) :: a(:,:),det
   det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-&
       a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+&
       a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
   end subroutine cmatdet3

   subroutine cmatdet4(a,det)
   complex (kind=r8) :: a(:,:),det
   det=(a(1,1)*a(2,2)-a(2,1)*a(1,2))*(a(3,3)*a(4,4)-a(4,3)*a(3,4)) &
      -(a(1,1)*a(2,3)-a(2,1)*a(1,3))*(a(3,2)*a(4,4)-a(4,2)*a(3,4)) &
      +(a(1,1)*a(2,4)-a(2,1)*a(1,4))*(a(3,2)*a(4,3)-a(4,2)*a(3,3)) &
      +(a(1,2)*a(2,3)-a(2,2)*a(1,3))*(a(3,1)*a(4,4)-a(4,1)*a(3,4)) &
      -(a(1,2)*a(2,4)-a(2,2)*a(1,4))*(a(3,1)*a(4,3)-a(4,1)*a(3,3)) &
      +(a(1,3)*a(2,4)-a(2,3)*a(1,4))*(a(3,1)*a(4,2)-a(4,1)*a(3,2))
   end subroutine cmatdet4

   subroutine cmatinv2(a,det)
   complex (kind=r8) :: a(:,:),det,temp,deti
   det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   deti=1.0_r8/det
   temp=a(1,1)
   a(1,1)=a(2,2)*deti
   a(2,2)=temp*deti
   a(1,2)=-a(1,2)*deti
   a(2,1)=-a(2,1)*deti
   end subroutine cmatinv2
   end subroutine v6potcalc

   subroutine opcalc(sigma,tau,sigmatau,sxmall,sxxmall,fasigtau,fatau)
   use matrixmod
   complex(kind=r8) :: sigma(:,:),tau(:,:),sigmatau(:,:,:)
   complex(kind=r8) :: sxmall(:,:,:),sxxmall(:,:,:,:)
   real(kind=r8) :: fasigtau(:,:,:,:),fatau(:,:)
   complex(kind=r8) :: d1,d2,d3,d4
   integer(kind=i4) :: ic1,i1,ic2,i2,i,j,ic,it,jc1,jc2,jc,k,l,kt,kc,lc1,lc,kc1,kc2
   complex(kind=r8) :: M3v(3,3),M2v(2,2),M4(4,4),M3(3,3),M2(2,2),dM4,dM3,dM2
   complex(kind=r8) :: Ma(2,2),Mb(2,2),Mc(2,2),Md(2,2),dMa,dMc,dm
write (50,*) 'called opcalc'
call flush(50)
   sigma=czero
   tau=czero
   sigmatau=czero
   do i=1,npart
      do ic=1,3
         sigma(ic,i)=sxmall(i,ic,i)
         tau(ic,i)=sxmall(i,ic+3,i)
         do jc=1,3
            kc=jc+3*(ic-1)+6
            sigmatau(ic,jc,i)=sxmall(i,kc,i)
         enddo
      enddo
   enddo
! O_k f_ij i < j != k
   do it=1,3
     do ic=1,3
       ic1=3*(ic-1)+6+it
       do i=1,npart-1
         M3(1,1)=sxmall(i,ic1,i)
         do jc=1,3
           jc1=3*(jc-1)+6+it
           do j=i+1,npart
             M3(1,2)=sxmall(i,jc1,j)
             M3(2,1)=sxmall(j,ic1,i)
             M3(2,2)=sxmall(j,jc1,j)
             do k=1,npart
               if(k.eq.i.or.k.eq.j) cycle
               M3(3,1)=sxmall(k,ic1,i)
               M3(3,2)=sxmall(k,jc1,j)
               do kc=1,3
                 M3(1,3)=sxmall(i,kc,k)
                 M3(2,3)=sxmall(j,kc,k)
                 M3(3,3)=sxmall(k,kc,k)
                 call cmatdet3(M3,dM3)
                 sigma(kc,k)=sigma(kc,k)+dM3*fasigtau(ic,i,jc,j)
                 M3(1,3)=sxmall(i,kc+3,k)
                 M3(2,3)=sxmall(j,kc+3,k)
                 M3(3,3)=sxmall(k,kc+3,k)
                 call cmatdet3(M3,dM3)
                 tau(kc,k)=tau(kc,k)+dM3*fasigtau(ic,i,jc,j)
                 do kc1=1,3
                   kc2=kc1+3*(kc-1)+6
                   M3(1,3)=sxmall(i,kc2,k)
                   M3(2,3)=sxmall(j,kc2,k)
                   M3(3,3)=sxmall(k,kc2,k)
                   call cmatdet3(M3,dM3)
                   sigmatau(kc,kc1,k)=sigmatau(kc,kc1,k)+dM3*fasigtau(ic,i,jc,j)
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! O_j f_ij i!= (j=k)
   do it=1,3
     do ic=1,3
       ic1=3*(ic-1)+6+it
       do i=1,npart
         M2(1,1)=sxmall(i,ic1,i)
         do jc=1,3
           jc1=3*(jc-1)+6+it
           do j=1,npart
             if(i.eq.j) cycle
             M2(2,1)=sxmall(j,ic1,i)
             do kc=1,3
               M2(1,2)=sxxmall(i,jc1,kc,j)
               M2(2,2)=sxxmall(j,jc1,kc,j)
               call cmatdet2(M2,dM2)
               sigma(kc,j)=sigma(kc,j)+dM2*fasigtau(ic,i,jc,j)
               M2(1,2)=sxxmall(i,jc1,kc+3,j)
               M2(2,2)=sxxmall(j,jc1,kc+3,j)
               call cmatdet2(M2,dM2)
               tau(kc,j)=tau(kc,j)+dM2*fasigtau(ic,i,jc,j)
               do kc1=1,3
                 kc2=kc1+3*(kc-1)+6
                 M2(1,2)=sxxmall(i,jc1,kc2,j)
                 M2(2,2)=sxxmall(j,jc1,kc2,j)
                 call cmatdet2(M2,dM2)
                 sigmatau(kc,kc1,j)=sigmatau(kc,kc1,j)+dM2*fasigtau(ic,i,jc,j)
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! O_k f_ij i < j != k
   do it=4,6
     do i=1,npart-1
       M3(1,1)=sxmall(i,it,i)
       do j=i+1,npart
         M3(1,2)=sxmall(i,it,j)
         M3(2,1)=sxmall(j,it,i)
         M3(2,2)=sxmall(j,it,j)
         do k=1,npart
           if(k.eq.i.or.k.eq.j) cycle
           M3(3,1)=sxmall(k,it,i)
           M3(3,2)=sxmall(k,it,j)
           do kc=1,3
             M3(1,3)=sxmall(i,kc,k)
             M3(2,3)=sxmall(j,kc,k)
             M3(3,3)=sxmall(k,kc,k)
             call cmatdet3(M3,dM3)
             sigma(kc,k)=sigma(kc,k)+dM3*fatau(i,j)
             M3(1,3)=sxmall(i,kc+3,k)
             M3(2,3)=sxmall(j,kc+3,k)
             M3(3,3)=sxmall(k,kc+3,k)
             call cmatdet3(M3,dM3)
             tau(kc,k)=tau(kc,k)+dM3*fatau(i,j)
             do kc1=1,3
               kc2=kc1+3*(kc-1)+6
               M3(1,3)=sxmall(i,kc2,k)
               M3(2,3)=sxmall(j,kc2,k)
               M3(3,3)=sxmall(k,kc2,k)
               call cmatdet3(M3,dM3)
               sigmatau(kc,kc1,k)=sigmatau(kc,kc1,k)+dM3*fatau(i,j)
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! O_j f_ij i!= (j=k)
   do it=4,6
     do i=1,npart
       M2(1,1)=sxmall(i,it,i)
       do j=1,npart
         if(i.eq.j) cycle
         M2(2,1)=sxmall(j,it,i)
         do kc=1,3
           M2(1,2)=sxxmall(i,it,kc,j)
           M2(2,2)=sxxmall(j,it,kc,j)
           call cmatdet2(M2,dM2)
           sigma(kc,j)=sigma(kc,j)+dM2*fatau(i,j)
           M2(1,2)=sxxmall(i,it,kc+3,j)
           M2(2,2)=sxxmall(j,it,kc+3,j)
           call cmatdet2(M2,dM2)
           tau(kc,j)=tau(kc,j)+dM2*fatau(i,j)
           do kc1=1,3
             kc2=kc1+3*(kc-1)+6
             M2(1,2)=sxxmall(i,it,kc2,j)
             M2(2,2)=sxxmall(j,it,kc2,j)
             call cmatdet2(M2,dM2)
             sigmatau(kc,kc1,j)=sigmatau(kc,kc1,j)+dM2*fatau(i,j)
           enddo
         enddo
       enddo
     enddo
   enddo
   contains

   subroutine cmatdet2(a,det)
   complex (kind=r8) :: a(:,:),det
   det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   end subroutine cmatdet2

   subroutine cmatdet3(a,det)
   complex (kind=r8) :: a(:,:),det
   det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-&
       a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+&
       a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
   end subroutine cmatdet3

   end subroutine opcalc

   subroutine op2calc(cvs,cvst,cvt,sxmall,sxxmall,fasigtau,fatau)
   use matrixmod
   complex(kind=r8) :: cvs,cvst,cvt
   complex(kind=r8) :: sxmall(:,:,:),sxxmall(:,:,:,:)
   real(kind=r8) :: fasigtau(:,:,:,:),fatau(:,:)
   complex(kind=r8) :: d1,d2,d3,d4
   integer(kind=i4) :: ic1,i1,ic2,i2,i,j,ic,it,jc1,jc2,jc,k,l,kt,kc,lc1,lc,kc1
   complex(kind=r8) :: M3v(3,3),M2v(2,2),M4(4,4),M3(3,3),M2(2,2),dM4,dM3,dM2
   complex(kind=r8) :: Ma(2,2),Mb(2,2),Mc(2,2),Md(2,2),dMa,dMc,dm
   complex(kind=r8) :: sigma(3,3,npart,npart)
   complex(kind=r8) :: sigmatau(9,9,npart,npart)
   complex(kind=r8) :: tau(3,3,npart,npart)
   sigma=czero
   sigmatau=czero
   tau=czero
   do ic1=1,3
      do i1=1,npart-1
         d1=sxmall(i1,ic1,i1)
         do ic2=1,3
            do i2=i1+1,npart
               d2=sxmall(i2,ic2,i2)
               d3=sxmall(i1,ic2,i2)
               d4=sxmall(i2,ic1,i1)
               sigma(ic1,ic2,i1,i2)=sigma(ic1,ic2,i1,i2)+d1*d2-d3*d4
            enddo
         enddo
      enddo
   enddo
   do it=1,3
      do ic1=1,3
         jc1=3*(ic1-1)+6+it
         do i1=1,npart-1
            d1=sxmall(i1,jc1,i1)
            do ic2=1,3
               jc2=3*(ic2-1)+6+it
               do i2=i1+1,npart
                  d2=sxmall(i2,jc2,i2)
                  d3=sxmall(i1,jc2,i2)
                  d4=sxmall(i2,jc1,i1)
                  sigmatau(jc1,jc2,i1,i2)=sigmatau(jc1,jc2,i1,i2)+d1*d2-d3*d4
               enddo
            enddo
         enddo
      enddo
   enddo
! note, only the diagonal part of tau.tau is calculated!!!
   do it=4,6
      do i1=1,npart-1
         d1=sxmall(i1,it,i1)
         do i2=i1+1,npart
            d2=sxmall(i2,it,i2)
            d3=sxmall(i1,it,i2)
            d4=sxmall(i2,it,i1)
            tau(it-3,it-3,i1,i2)=tau(it-3,it-3,i1,i2)+d1*d2-d3*d4
         enddo
      enddo
   enddo
! v_ij f_jk i < j != k < l
   do i=1,npart-1
     do j=i+1,npart
       do it=1,3
         do ic=1,3
           ic1=3*(ic-1)+6+it
           M4(1,1)=sxmall(i,ic1,i)     
           do jc=1,3
             jc1=3*(jc-1)+6+it
             M4(1,2)=sxmall(i,jc1,j)
             M4(2,2)=sxmall(j,jc1,j)
             M4(2,1)=sxmall(j,ic1,i)
             Ma=M4(1:2,1:2)
             call cmatinv2(Ma,dMa)
             do k=1,npart-1
               if(k.eq.i.or.k.eq.j) cycle
               M4(3,1)=sxmall(k,ic1,i)
               M4(3,2)=sxmall(k,jc1,j)
               do l=k+1,npart
                 if(l.eq.i.or.l.eq.j) cycle
                 M4(4,1)=sxmall(l,ic1,i)
                 M4(4,2)=sxmall(l,jc1,j)
                 Mc=M4(3:4,1:2)
                 Mc=matmul(Mc,Ma)
                 do kc=1,3
                   M4(1,3)=sxmall(i,kc,k)
                   M4(2,3)=sxmall(j,kc,k)
                   M4(3,3)=sxmall(k,kc,k)
                   M4(4,3)=sxmall(l,kc,k)
                   do lc=1,3
                     M4(1,4)=sxmall(i,lc,l)
                     M4(2,4)=sxmall(j,lc,l)
                     M4(3,4)=sxmall(k,lc,l)
                     M4(4,4)=sxmall(l,lc,l)
                     if(dMa.ne.0.) then
                       Mb=M4(1:2,3:4)
                       Md=M4(3:4,3:4)
                       call cmatdet2(Md-matmul(Mc,Mb),dMc)
                       dm=dMa*dMc
                     else
                       call cmatdet4(M4,dM4)
                       dm=dM4
                     endif
                     sigma(kc,lc,k,l)=sigma(kc,lc,k,l)+dm*fasigtau(ic,i,jc,j)
                   enddo
                 enddo
                 do kt=1,3
                   do kc=1,3
                     kc1=3*(kc-1)+6+kt
                     M4(1,3)=sxmall(i,kc1,k)
                     M4(2,3)=sxmall(j,kc1,k)
                     M4(3,3)=sxmall(k,kc1,k)
                     M4(4,3)=sxmall(l,kc1,k)
                     do lc=1,3
                       lc1=3*(lc-1)+6+kt
                       M4(1,4)=sxmall(i,lc1,l)
                       M4(2,4)=sxmall(j,lc1,l)
                       M4(3,4)=sxmall(k,lc1,l)
                       M4(4,4)=sxmall(l,lc1,l)
                       if(dMa.ne.0.) then
                         Mb=M4(1:2,3:4)
                         Md=M4(3:4,3:4)
                         call cmatdet2(Md-matmul(Mc,Mb),dMc)
                         dm=dMa*dMc
                       else
                         call cmatdet4(M4,dM4)
                         dm=dM4
                       endif
                       sigmatau(kc,lc,k,l)=sigmatau(kc,lc,k,l)+dm*fasigtau(ic,i,jc,j)
                     enddo
                   enddo
                 enddo
                 do kt=4,6
                   M4(1,3)=sxmall(i,kt,k)
                   M4(2,3)=sxmall(j,kt,k)
                   M4(3,3)=sxmall(k,kt,k)
                   M4(4,3)=sxmall(l,kt,k)
                   M4(1,4)=sxmall(i,kt,l)
                   M4(2,4)=sxmall(j,kt,l)
                   M4(3,4)=sxmall(k,kt,l)
                   M4(4,4)=sxmall(l,kt,l)
                   if(dMa.ne.0.) then
                     Mb=M4(1:2,3:4)
                     Md=M4(3:4,3:4)
                     call cmatdet2(Md-matmul(Mc,Mb),dMc)
                     tau(kt-3,kt-3,k,l)=tau(kt-3,kt-3,k,l)+dMa*dMc*fasigtau(ic,i,jc,j)
                   else
                     call cmatdet4(M4,dM4)
                     tau(kt-3,kt-3,k,l)=tau(kt-3,kt-3,k,l)+dM4*fasigtau(ic,i,jc,j)
                   endif
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! f_ij v_kj (j=l)
   do i=1,npart
    do it=1,3
     do ic=1,3
       ic1=3*(ic-1)+6+it
         M3v(1,1)=sxmall(i,ic1,i)
         do jc=1,3
           jc1=3*(jc-1)+6+it
           do j=1,npart
             if(j.eq.i) cycle
             M3v(2,1)=sxmall(j,ic1,i)
             do k=1,npart
               if(k.eq.i.or.k.eq.j) cycle
               M3v(3,1)=sxmall(k,ic1,i)
               do kc=1,3
                 M3v(1,3)=sxmall(i,kc,k)
                 M3v(2,3)=sxmall(j,kc,k)
                 M3v(3,3)=sxmall(k,kc,k)
                 do lc=1,3
                   M3v(1,2)=sxxmall(i,jc1,lc,j)
                   M3v(2,2)=sxxmall(j,jc1,lc,j)
                   M3v(3,2)=sxxmall(k,jc1,lc,j)
                   call cmatdet3(M3v,dM3)
                   sigma(kc,lc,k,l)=sigma(kc,lc,k,l)+dm3*fasigtau(ic,i,jc,j)
                 enddo
               enddo
               do kt=1,3
                 do kc=1,3
                   kc1=3*(kc-1)+6+kt
                   M3v(1,3)=sxmall(i,kc1,k)
                   M3v(2,3)=sxmall(j,kc1,k)
                   M3v(3,3)=sxmall(k,kc1,k)
                   do lc=1,3
                     lc1=3*(lc-1)+6+kt
                     M3v(1,2)=sxxmall(i,jc1,lc1,j)
                     M3v(2,2)=sxxmall(j,jc1,lc1,j)
                     M3v(3,2)=sxxmall(k,jc1,lc1,j)
                     call cmatdet3(M3v,dM3)
                     sigmatau(kc,lc,k,j)=sigmatau(kc,lc,k,j)+dm3*fasigtau(ic,i,jc,j)
                   enddo
                 enddo
               enddo
               do kt=4,6
                 M3v(1,2)=sxxmall(i,jc1,kt,j)
                 M3v(1,3)=sxmall(i,kt,k)
                 M3v(2,2)=sxxmall(j,jc1,kt,j)
                 M3v(3,2)=sxxmall(k,jc1,kt,j)
                 M3v(2,3)=sxmall(j,kt,k)
                 M3v(3,3)=sxmall(k,kt,k)
                 call cmatdet3(M3v,dM3)
                 tau(kt-3,kt-3,k,j)=tau(kt-3,kt-3,k,j)+dM3*fasigtau(ic,i,jc,j)
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
!f_ij v_ij (i=k,j=l)
   do it=1,3
     do ic=1,3
       ic1=3*(ic-1)+6+it
       do i=1,npart-1
         do jc=1,3
           jc1=3*(jc-1)+6+it
           do j=i+1,npart
             do kc=1,3
               M2v(1,1)=sxxmall(i,ic1,kc,i)
               M2v(2,1)=sxxmall(j,ic1,kc,i)
               do lc=1,3
                 M2v(1,2)=sxxmall(i,jc1,lc,j)
                 M2v(2,2)=sxxmall(j,jc1,lc,j)
                 call cmatdet2(M2v,dM2)
                 sigma(kc,lc,i,j)=sigma(kc,lc,i,j)+dm2*fasigtau(ic,i,jc,j)
               enddo
             enddo
             do kt=1,3
               do kc=1,3
                 kc1=3*(kc-1)+6+kt
                 M2v(1,1)=sxxmall(i,ic1,kc1,i)
                 M2v(2,1)=sxxmall(j,ic1,kc1,i)
                 do lc=1,3
                   lc1=3*(lc-1)+6+kt
                   M2v(1,2)=sxxmall(i,jc1,lc1,j)
                   M2v(2,2)=sxxmall(j,jc1,lc1,j)
                   call cmatdet2(M2v,dM2)
                   sigmatau(kc,lc,i,j)=sigmatau(kc,lc,i,j)+dm2*fasigtau(ic,i,jc,j)
                 enddo
               enddo
             enddo
             do kt=4,6
               M2v(1,1)=sxxmall(i,ic1,kt,i)
               M2v(2,1)=sxxmall(j,ic1,kt,i)
               M2v(1,2)=sxxmall(i,jc1,kt,j)
               M2v(2,2)=sxxmall(j,jc1,kt,j)
               call cmatdet2(M2v,dM2)
               tau(kt-3,kt-3,i,j)=tau(kt-3,kt-3,i,j)+dM2*fasigtau(ic,i,jc,j)
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! tau correlation 
! v_ij f_kl
   do i=1,npart-1
     do it=4,6
       M4(1,1)=sxmall(i,it,i)     
       do j=i+1,npart
         M4(1,2)=sxmall(i,it,j)
         M4(2,2)=sxmall(j,it,j)
         M4(2,1)=sxmall(j,it,i)
         Ma=M4(1:2,1:2)
         call cmatinv2(Ma,dMa)
         do k=1,npart-1
           if(k.eq.i.or.k.eq.j) cycle
           M4(3,1)=sxmall(k,it,i)
           M4(3,2)=sxmall(k,it,j)
           do l=k+1,npart
             if(l.eq.i.or.l.eq.j) cycle
             M4(4,1)=sxmall(l,it,i)
             M4(4,2)=sxmall(l,it,j)
             Mc=M4(3:4,1:2)
             Mc=matmul(Mc,Ma)
             do kc=1,3  !sigma S
               M4(1,3)=sxmall(i,kc,k)
               M4(2,3)=sxmall(j,kc,k)
               M4(3,3)=sxmall(k,kc,k)
               M4(4,3)=sxmall(l,kc,k)
               do lc=1,3
                 M4(1,4)=sxmall(i,lc,l)
                 M4(2,4)=sxmall(j,lc,l)
                 M4(3,4)=sxmall(k,lc,l)
                 M4(4,4)=sxmall(l,lc,l)
                 if(dMa.ne.0.) then
                   Mb=M4(1:2,3:4)
                   Md=M4(3:4,3:4)
                   call cmatdet2(Md-matmul(Mc,Mb),dMc)
                   sigma(kc,lc,k,l)=sigma(kc,lc,k,l)+dma*dmc*fatau(i,j)
                 else
                   call cmatdet4(M4,dM4)
                   sigma(kc,lc,k,l)=sigma(kc,lc,k,l)+dm4*fatau(i,j)
                 endif
               enddo
             enddo
             do kt=1,3 !sigma-tau S-tau
               do kc=1,3
                 kc1=3*(kc-1)+6+kt
                 M4(1,3)=sxmall(i,kc1,k)
                 M4(2,3)=sxmall(j,kc1,k)
                 M4(3,3)=sxmall(k,kc1,k)
                 M4(4,3)=sxmall(l,kc1,k)
                 do lc=1,3
                   lc1=3*(lc-1)+6+kt
                   M4(1,4)=sxmall(i,lc1,l)
                   M4(2,4)=sxmall(j,lc1,l)
                   M4(3,4)=sxmall(k,lc1,l)
                   M4(4,4)=sxmall(l,lc1,l)
                   if(dMa.ne.0.) then
                     Mb=M4(1:2,3:4)
                     Md=M4(3:4,3:4)
                     call cmatdet2(Md-matmul(Mc,Mb),dMc)
                     sigmatau(kc,lc,k,l)=sigmatau(kc,lc,k,l)+dma*dmc*fatau(i,j)
                   else
                     call cmatdet4(M4,dM4)
                     sigmatau(kc,lc,k,l)=sigmatau(kc,lc,k,l)+dm4*fatau(i,j)
                   endif
                 enddo
               enddo
             enddo
             do kt=4,6 !tau
               M4(1,3)=sxmall(i,kt,k)
               M4(2,3)=sxmall(j,kt,k)
               M4(3,3)=sxmall(k,kt,k)
               M4(4,3)=sxmall(l,kt,k)
               M4(1,4)=sxmall(i,kt,l)
               M4(2,4)=sxmall(j,kt,l)
               M4(3,4)=sxmall(k,kt,l)
               M4(4,4)=sxmall(l,kt,l)
               if(dMa.ne.0.) then
                 Mb=M4(1:2,3:4)
                 Md=M4(3:4,3:4)
                 call cmatdet2(Md-matmul(Mc,Mb),dMc)
                 tau(kt-3,kt-3,k,l)=tau(kt-3,kt-3,k,l)+dma*dmc*fatau(i,j)
               else
                 call cmatdet4(M4,dM4)
                 tau(kt-3,kt-3,k,l)=tau(kt-3,kt-3,k,l)+dm4*fatau(i,j)
! TO HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!
               endif 
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo
! f_ij v_kj (j=l)
   do i=1,npart
     do it=4,6
       M3v(1,1)=sxmall(i,it,i)
       do j=1,npart
         if(j.eq.i) cycle
         M3v(2,1)=sxmall(j,it,i)
         do k=1,npart
           if(k.eq.i.or.k.eq.j) cycle
           M3v(3,1)=sxmall(k,it,i)
           do kc=1,3 !sigma S
             M3v(1,3)=sxmall(i,kc,k)
             M3v(2,3)=sxmall(j,kc,k)
             M3v(3,3)=sxmall(k,kc,k)
             do lc=1,3
               M3v(1,2)=sxxmall(i,it,lc,j)
               M3v(2,2)=sxxmall(j,it,lc,j)
               M3v(3,2)=sxxmall(k,it,lc,j)
               call cmatdet3(M3v,dM3)
!              cvs=cvs+dM3*vasig(kc,k,lc,j)*fatau(i,j)
             enddo
           enddo
           do kt=1,3 !sigma-tau S-tau
             do kc=1,3
               kc1=3*(kc-1)+6+kt
               M3v(1,3)=sxmall(i,kc1,k)
               M3v(2,3)=sxmall(j,kc1,k)
               M3v(3,3)=sxmall(k,kc1,k)
               do lc=1,3
                 lc1=3*(lc-1)+6+kt
                 M3v(1,2)=sxxmall(i,it,lc1,j)
                 M3v(2,2)=sxxmall(j,it,lc1,j)
                 M3v(3,2)=sxxmall(k,it,lc1,j)
                 call cmatdet3(M3v,dM3)
!                cvst=cvst+dM3*vasigtau(kc,k,lc,j)*fatau(i,j)
               enddo 
             enddo
           enddo
           do kt=4,6 !tau
             M3v(1,2)=sxxmall(i,it,kt,j)
             M3v(1,3)=sxmall(i,kt,k)
             M3v(2,2)=sxxmall(j,it,kt,j)
             M3v(3,2)=sxxmall(k,it,kt,j)
             M3v(2,3)=sxmall(j,kt,k)
             M3v(3,3)=sxmall(k,kt,k)
             call cmatdet3(M3v,dM3)
!            cvt=cvt+dM3*vatau(k,j)*fatau(i,j)
           enddo
         enddo
       enddo
     enddo
   enddo
!f_ij v_ij
   do it=4,6
     do i=1,npart-1
       do j=i+1,npart
         do kc=1,3 !sigma S
           M2v(1,1)=sxxmall(i,it,kc,i)
           M2v(2,1)=sxxmall(j,it,kc,i)
           do lc=1,3
             M2v(1,2)=sxxmall(i,it,lc,j)
             M2v(2,2)=sxxmall(j,it,lc,j)
             call cmatdet2(M2v,dM2)
!            cvs=cvs+dM2*vasig(kc,i,lc,j)*fatau(i,j)
           enddo
         enddo
         do kt=1,3 !sigma-tau S-tau
           do kc=1,3
             kc1=3*(kc-1)+6+kt
             M2v(1,1)=sxxmall(i,it,kc1,i)
             M2v(2,1)=sxxmall(j,it,kc1,i)
             do lc=1,3
               lc1=3*(lc-1)+6+kt
               M2v(1,2)=sxxmall(i,it,lc1,j)
               M2v(2,2)=sxxmall(j,it,lc1,j)
               call cmatdet2(M2v,dM2)
!              cvst=cvst+dM2*vasigtau(kc,i,lc,j)*fatau(i,j)
             enddo
           enddo
         enddo 
         do kt=4,6 !tau
           M2v(1,1)=sxxmall(i,it,kt,i)
           M2v(2,1)=sxxmall(j,it,kt,i)
           M2v(1,2)=sxxmall(i,it,kt,j)
           M2v(2,2)=sxxmall(j,it,kt,j)
           call cmatdet2(M2v,dM2)
!          cvt=cvt+dM2*vatau(i,j)*fatau(i,j)
         enddo
       enddo
     enddo
   enddo
   contains

   subroutine cmatdet2(a,det)
   complex (kind=r8) :: a(:,:),det
   det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   end subroutine cmatdet2

   subroutine cmatdet3(a,det)
   complex (kind=r8) :: a(:,:),det
   det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-&
       a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+&
       a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
   end subroutine cmatdet3

   subroutine cmatdet4(a,det)
   complex (kind=r8) :: a(:,:),det
   det=(a(1,1)*a(2,2)-a(2,1)*a(1,2))*(a(3,3)*a(4,4)-a(4,3)*a(3,4)) &
      -(a(1,1)*a(2,3)-a(2,1)*a(1,3))*(a(3,2)*a(4,4)-a(4,2)*a(3,4)) &
      +(a(1,1)*a(2,4)-a(2,1)*a(1,4))*(a(3,2)*a(4,3)-a(4,2)*a(3,3)) &
      +(a(1,2)*a(2,3)-a(2,2)*a(1,3))*(a(3,1)*a(4,4)-a(4,1)*a(3,4)) &
      -(a(1,2)*a(2,4)-a(2,2)*a(1,4))*(a(3,1)*a(4,3)-a(4,1)*a(3,3)) &
      +(a(1,3)*a(2,4)-a(2,3)*a(1,4))*(a(3,1)*a(4,2)-a(4,1)*a(3,2))
   end subroutine cmatdet4

   subroutine cmatinv2(a,det)
   complex (kind=r8) :: a(:,:),det,temp,deti
   det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   deti=1.0_r8/det
   temp=a(1,1)
   a(1,1)=a(2,2)*deti
   a(2,2)=temp*deti
   a(1,2)=-a(1,2)*deti
   a(2,1)=-a(2,1)*deti
   end subroutine cmatinv2
   end subroutine op2calc


end module potential
