   program drive3b
   use random
   use matrixmod
   use correlator
   use v6pot
   use v3bpot
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), parameter :: czero=(0.0_r8,0.0_r8)
   integer(kind=i8) :: irn
   integer(kind=i4) :: irnd,i,j,ic,is1,is2,is3,js1,js2,js3,npart
   integer(kind=i4) :: is,js,it,jt
   real(kind=r8) :: x(3,3),small,vfact(8)
   complex(kind=r8) :: ph(3,4,3),spsmall(4,3),sp(4,3),smati(3,3),det
   complex(kind=r8) :: sp0(4,3),v3a,v3b,v3c,v3d
   complex(kind=r8) :: sxmallz(3,4,3),spx(4,15,3),sxmall(3,15,3)
   complex(kind=r8) :: sxz(4,3,3)
   complex(kind=r8) :: d2b(4,4,3),d3b(4,4,4,1)
   real(kind=r8) :: c,a2p,a2s,a3,ar
   integer(kind=i4) :: v3bflag,iunit,ic1,jc1,icharge,ncut
   logical :: dov3b2
   real(kind=r8) :: sigmatau1(3,3,3,3),sigmatau2(3,3,3,3)
   real(kind=r8) :: delta(3,3),delta2(3,3)
   real(kind=r8) :: vc
   real(kind=r8) :: xpi(3,3,3,3),xxpi(3,3,3,3)
   real(kind=r8) :: a2p3b,a3p3b
   complex(kind=r8) :: v32ppwa,v32psw,v32ppwc,v33paa,v33pss,ss3pi
   complex(kind=r8) :: v3e,v3d2,v3ec4
   integer(kind=i4), parameter :: isospin(3,3,0:3) = &
       reshape((/3,3,3, -1,-1,-1, -1,-1,-1, &
         1,3,3, 3,1,3, 3,3,1, &
         1,1,3, 1,3,1, 3,1,1, &
         1,1,1, -1,-1,-1, -1,-1,-1/),(/3,3,4/))
   integer(kind=i4), parameter :: ncharge(0:3)=(/1,3,3,1/)
   open(unit=10,file='anticomm.nnn.out',status='unknown')
   open(unit=11,file='anticomm.nnp.out',status='unknown')
   open(unit=12,file='anticomm.npp.out',status='unknown')
   open(unit=13,file='anticomm.ppp.out',status='unknown')
   open(unit=14,file='comm.nnn.out',status='unknown')
   open(unit=15,file='comm.nnp.out',status='unknown')
   open(unit=16,file='comm.npp.out',status='unknown')
   open(unit=17,file='comm.ppp.out',status='unknown')
   open(unit=18,file='tma.nnn.out',status='unknown')
   open(unit=19,file='tma.nnp.out',status='unknown')
   open(unit=20,file='tma.npp.out',status='unknown')
   open(unit=21,file='tma.ppp.out',status='unknown')
   open(unit=22,file='ve.nnn.out',status='unknown')
   open(unit=23,file='ve.nnp.out',status='unknown')
   open(unit=24,file='ve.npp.out',status='unknown')
   open(unit=25,file='ve.ppp.out',status='unknown')
   open(unit=26,file='vd2.nnn.out',status='unknown')
   open(unit=27,file='vd2.nnp.out',status='unknown')
   open(unit=28,file='vd2.npp.out',status='unknown')
   open(unit=29,file='vd2.ppp.out',status='unknown')
   open(unit=30,file='vcc4ve.nnn.out',status='unknown')
   open(unit=31,file='vcc4ve.nnp.out',status='unknown')
   open(unit=32,file='vcc4ve.npp.out',status='unknown')
   open(unit=33,file='vcc4ve.ppp.out',status='unknown')
   do i=10,33
      rewind i
      do j=1,13
         write (i,*)
      enddo
   enddo
   npart=3
   vfact=0.0_r8
   call v6potinit(npart,0,1.0e10_r8,1,1,vfact)
   call initcormod(npart,.false.)
   read (*,*) irn
   write (*,*) 'random seed ',irn
   call setrn(irn)
   read (*,*) small
   read (*,*) irnd
   if (irnd.eq.1) then
      x=reshape(2.0_r8*(randn(9)-0.5_r8),shape(x))
   else
      read (5,'(3g25.18)') x
   endif
   write (*,'(''particle positions'')') 
   write (*,'(3g25.18)') x
   call v3bpotinit(npart,1.0e10_r8,'uix')
   spsmall=small*reshape(cmplx(randn(12)-0.5_r8,randn(12)-0.5_r8),shape(sp))
!write (*,*) sp
   do icharge=0,3
      do it=1,ncharge(icharge)
         do jt=1,ncharge(icharge)
            do is=0,7
               do js=0,7
                  is1=isospin(1,it,icharge)-iand(is,1)+1
                  is2=isospin(2,it,icharge)-iand(shiftr(is,1),1)+1
                  is3=isospin(3,it,icharge)-iand(shiftr(is,2),1)+1
                  js1=isospin(1,jt,icharge)-iand(js,1)+1
                  js2=isospin(2,jt,icharge)-iand(shiftr(js,1),1)+1
                  js3=isospin(3,jt,icharge)-iand(shiftr(js,2),1)+1
                     ph=czero
                     ph(1,is1,1)=cone
                     ph(2,is2,2)=cone
                     ph(3,is3,3)=cone
                     sp0=czero
                     sp0(js1,1)=cone
                     sp0(js2,2)=cone
                     sp0(js3,3)=cone
                     sp=spsmall
                     sp(js1,1)=cone
                     sp(js2,2)=cone
                     sp(js3,3)=cone
                     do i=1,3
                        smati(:,i)=matmul(ph(:,:,i),sp(:,i))
                     enddo
                     call matinv(smati,det,3)
                     sxmallz=reshape(matmul(smati,reshape(ph(:,:,:) &
                        ,(/npart,4*npart/))),shape(sxmallz))
                     sxz=reshape(transpose(reshape(sxmallz,(/npart,4*npart/))) &
                       ,shape(sxz))
                     d2b=czero
                     call g2bval(d2b,sxz,cone)
                     d3b=czero
                     call g3bval(d3b,sxz,cone)
!write (*,'(2f10.5)')  d3b
!stop
                     spx=opmult(sp0)
                     call setxspx(x,spx)
                     call v3bval(d2b,d3b,vc,v3a,v3b,v3c,v3d)
                     write(10+icharge,'(4i2,2f15.8)') it,jt,is,js,det*v3a
                     write(14+icharge,'(4i2,2f15.8)') it,jt,is,js,det*v3c
                     write(18+icharge,'(4i2,2f15.8)') it,jt,is,js,det*v3b
!                     write(22+icharge,'(4i2,2f15.8)') it,jt,is,js,det*v3e
!                     write(26+icharge,'(4i2,2f15.8)') it,jt,is,js,det*v3d2
!                     write(30+icharge,'(4i2,2f15.8)') it,jt,is,js,det*v3ec4
!                     write (14+icharge,'(6i2,2g20.10)') is1,is2,is3,js1,js2,js3,v32ppwc*det
               enddo
            enddo
         enddo
      enddo
   enddo
   end program drive3b
