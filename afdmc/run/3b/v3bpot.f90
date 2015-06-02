module v3bpot
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   integer(kind=i4), private, parameter :: levi(2,3) = &
      reshape((/2,3, 3,1, 1,2/),(/2,3/))
   real(kind=r8), private, parameter :: amu=0.69953054_r8
   character(len=3), private, save :: v3b
   integer(kind=i4), private, save :: npart,ncut,ncutt,npair,ntrip,iv3bflag
   real(kind=r8), private, save :: el,eli,pi
   real(kind=r8), private, save :: c3b,a2p3ba,a2p3bc,a2s3b,a3p3b,ar3b,dfac,dnorm
   real(kind=r8), private, save :: a1pcon
   real(kind=r8), private, save, allocatable :: x(:,:)
   complex(kind=r8), private, save, allocatable :: spx(:,:,:)

contains
   subroutine v3bpotinit(npartin,elin,v3bin)
   use mympi
   integer(kind=i4), intent(in) :: npartin
   real(kind=r8), intent(in) :: elin
   character(len=3), intent(in) :: v3bin
   npart=npartin
   npair=(npart*(npart-1))/2
   ntrip=(npair*(npart-2))/3
   el=elin
   eli=1.0_r8/el
   v3b=v3bin
   pi=4.0_r8*atan(1.0_r8)
   if(v3b.eq.'no3') then
      ncut=1
      ncutt=2
      c3b=0.0_r8
      a2p3ba=0.0_r8
      a2p3bc=0.0_r8
      a2s3b=0.0_r8
      a3p3b=0.0_r8
      ar3b=0.0_r8
      a1pcon=0.0_r8
      dfac=0.0_r8
      iv3bflag=0
   else if(v3b.eq.'ch1') then
      ncut=4
      ncutt=1
      c3b=(1.0_r8/1.2_r8)**ncut
      a2p3ba=-0.0293_r8
      a2p3bc=-0.0293_r8*0.25_r8
      a2s3b=1.0_r8
      a3p3b=0.0_r8
      ar3b=0.0048_r8
      a1pcon=1.0_r8
      dfac=4.0_r8*pi/amu**3
      iv3bflag=1
   else if(v3b.eq.'uix') then
      ncut=2
      ncutt=2
      c3b=2.1_r8
      a2p3ba=-0.0293_r8
      a2p3bc=-0.0293_r8*0.25_r8
      a2s3b=0.0_r8
      a3p3b=0.0_r8
      ar3b=0.00480_r8
      a1pcon=0.0_r8
      dfac=0.0_r8
      iv3bflag=1
   else if(v3b.eq.'il1') then
      ncut=2
      ncutt=2
      c3b=2.1_r8
      a2p3ba=-0.0385_r8
      a2p3bc=-0.0385_r8*0.25_r8
      a2s3b=0.1_r8
      a3p3b=0.0026_r8
      ar3b=0.00705_r8
      a1pcon=0.0_r8
      dfac=0.0_r8
      iv3bflag=2
   else if(v3b.eq.'il2') then
      ncut=2
      ncutt=2
      c3b=2.1_r8
      a2p3ba=-0.0370_r8
      a2p3bc=-0.0370_r8*0.25_r8
      a2s3b=-1.0_r8
      a3p3b=0.0026_r8
      ar3b=0.00705_r8
      a1pcon=0.0_r8
      dfac=0.0_r8
      iv3bflag=2
   else if(v3b.eq.'il3') then
      ncut=2
      ncutt=2
      c3b=1.5_r8
      a2p3ba=-0.070_r8
      a2p3bc=-0.070_r8*0.25_r8
      a2s3b=-1.0_r8
      a3p3b=0.0065_r8
      ar3b=0.032_r8
      a1pcon=0.0_r8
      dfac=0.0_r8
      iv3bflag=2
   else if(v3b.eq.'il4')then
      ncut=2
      ncutt=2
      c3b=2.1_r8
      a2p3ba=-0.0280_r8
      a2p3bc=-0.0280_r8*0.25_r8
      a2s3b=-1.0_r8
      a3p3b=0.0021_r8
      ar3b=0.0039_r8
      a1pcon=0.0_r8
      dfac=0.0_r8
      iv3bflag=2
   else if(v3b.eq.'rd3') then
      if (myrank().eq.0) then
         read (5,*) ncut
         read (5,*) ncutt
         read (5,*) c3b
         read (5,*) a2p3ba
         read (5,*) a2p3bc
         read (5,*) a2s3b
         read (5,*) ar3b
         read (5,*) a1pcon
         read (5,*) dfac
         read (5,*) iv3bflag
      endif
      call bcast(ncut)
      call bcast(ncutt)
      call bcast(c3b)
      call bcast(a2p3ba)
      call bcast(a2p3bc)
      call bcast(a2s3b)
      call bcast(ar3b)
      call bcast(a1pcon)
      call bcast(dfac)
      call bcast(iv3bflag)
   else
      if (myrank().eq.0) then
         write (6,'(1x,''serious error pot3bpar'')')
         write (6,'(1x,''bad three body potential name  '',a3)') v3b
         write (6,'(1x,''allowed names: '')')
         write (6,'(16x,a3)') 'no3'
         write (6,'(16x,a3)') 'ch1'
         write (6,'(16x,a3)') 'uix'
         write (6,'(16x,a3)') 'il1'
         write (6,'(16x,a3)') 'il2'
         write (6,'(16x,a3)') 'il3'
         write (6,'(16x,a3)') 'il4'
         write (6,'(16x,a21)') 'rd3 to read in values'
         write (6,'(1x,''stopping the calculation '')')
         call abort
      else
         call barrier !wait for 0 to print
      endif
   endif
   if (myrank().eq.0.and.iv3bflag.ne.0) then
      write (6,'(''three-body potential name ='',t35,a3)') v3b
      write (6,'(''anticommutator strength ='',t35,f10.5)') a2p3ba
      write (6,'(''commutator strength ='',t35,f10.5)') a2p3bc
      write (6,'(''Tucson-Melbourne strength ='',t35,f10.5)') a2s3b
      write (6,'(''Central strength ='',t35,f10.5)') ar3b
      write (6,'(''1-pion contact strength ='',t35,f10.5)') a1pcon
      write (6,'(''Pion propagator Delta strength ='',t35,f10.5)') dfac
      write (6,'(''cut off exponent power ='',t35,i10)') ncut
      write (6,'(''tensor cut off power ='',t35,i10)') ncutt
      write (6,'(''cut off exponent strength ='',t35,f10.5)') c3b
   endif
   dnorm=ncut/(4.0_r8*pi*gamma(3.0_r8/ncut))
   dnorm=c3b**(3.0_r8/ncut)*dnorm
   if (allocated(x)) then
      deallocate(x,spx)
   endif
   allocate(x(3,npart),spx(4,15,npart))
   end subroutine v3bpotinit

   subroutine setxspx(xin,spxin)
   real(kind=r8), intent(in) :: xin(3,npart)
   complex(kind=r8), intent(in) :: spxin(4,15,npart)
   x=xin
   spx=spxin
   end subroutine setxspx

   subroutine v3bval(d2b,d3b,vc,v3a,v3tm,v3c,v3d)
   complex(kind=r8), intent(in) :: d2b(4,4,npair),d3b(4,4,4,ntrip)
   complex(kind=r8), intent(out) :: v3a,v3tm,v3c,v3d
   real(kind=r8) :: vc,xpi(3,npart,3,npart),xxpi(3,npart,3,npart),cut
   real(kind=r8) :: delta(npart,npart),delta2(npart,npart),dx(3),r,vfac
   real(kind=r8) :: gr3b(npart),g2s3b(3,npart,npart),ex,exc,y2,t2,z2,vtm
   real(kind=r8) :: delsum(npart)
   complex(kind=r8) :: v2tmp1(4,4),v2tmp2(4,4),v2tmp3(4,4)
   complex(kind=r8) :: v3tmp1(4,4,4),v3tmp2(4,4,4)
   integer(kind=i4) :: i,j,k,ic,jc,kc,ic3,jc3,kc3,js,ks,ij,ijk,it
   v3a=czero
   v3tm=czero
   v3c=czero
   v3d=czero
   vc=0.0_r8
   if (iv3bflag.eq.0) return
   gr3b=0.0_r8
   g2s3b=0.0_r8
   xpi=0.0_r8
   delta=0.0_r8
   delsum=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dx=x(:,i)-x(:,j)
         dx=dx-el*nint(dx*eli)
         r=sqrt(sum(dx**2))
         dx=dx/r
         ex=exp(-amu*r)
         exc=exp(-c3b*r**ncut)
         cut=1.0_r8-exc
         y2=cut*ex/(amu*r)
         t2=(1.0_r8+3.0_r8/(r*amu)+3.0_r8/(r*amu)**2)*y2*cut**(ncutt-1)
         z2=amu*r*(y2-t2)/3.0_r8
         delta(i,j)=dnorm*exc
         delta(j,i)=delta(i,j)
         delsum(i)=delsum(i)+delta(i,j)
         delsum(j)=delsum(j)+delta(i,j)
         xpi(:,i,1,j)=3.0_r8*t2*dx(1)*dx
         xpi(:,i,2,j)=3.0_r8*t2*dx(2)*dx
         xpi(:,i,3,j)=3.0_r8*t2*dx(3)*dx
         xpi(1,i,1,j)=xpi(1,i,1,j)+(y2-t2)-dfac*delta(i,j)
         xpi(2,i,2,j)=xpi(2,i,2,j)+(y2-t2)-dfac*delta(i,j)
         xpi(3,i,3,j)=xpi(3,i,3,j)+(y2-t2)-dfac*delta(i,j)
         xpi(:,j,:,i)=xpi(:,i,:,j)
         g2s3b(:,j,i)=dx*z2
         g2s3b(:,i,j)=-g2s3b(:,j,i)
         gr3b(i)=gr3b(i)+t2**2
         gr3b(j)=gr3b(j)+t2**2
         vc=vc-ar3b*t2**4
      enddo
   enddo
   vc=vc+0.5_r8*ar3b*sum(gr3b**2)
   delta2=matmul(delta,delta)
   xxpi=reshape(matmul(reshape(xpi,(/3*npart,3*npart/)), &
      reshape(xpi,(/3*npart,3*npart/))),(/3,npart,3,npart/))
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         v2tmp1=0.0_r8
         v2tmp2=0.0_r8
         v2tmp3=0.0_r8
         do ic=1,3
            ic3=(ic+1)*3
            do jc=1,3
               jc3=(jc+1)*3
               vtm=sum(g2s3b(ic,:,i)*g2s3b(jc,:,j))
               do it=1,3
                  do js=1,4
                     v2tmp1(:,js)=v2tmp1(:,js)+xxpi(ic,i,jc,j) &
                        *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                     v2tmp2(:,js)=v2tmp2(:,js)+vtm &
                        *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                     v2tmp3(:,js)=v2tmp3(:,js)+xpi(ic,i,jc,j)*( &
                        delsum(j)+delsum(i)-2.0_r8*delta(i,j)) &
                        *spx(:,ic3+it,i)*spx(js,jc3+it,j)
                  enddo
               enddo
            enddo
         enddo
         v3a=v3a+sum(v2tmp1*d2b(:,:,ij))
         v3tm=v3tm+sum(v2tmp2*d2b(:,:,ij))
         v3d=v3d+sum(v2tmp3*d2b(:,:,ij))
      enddo
   enddo
   v3a=4.0_r8*a2p3ba*v3a !anticommutator
   v3tm=a2s3b*v3tm      !Tucson-Melbourne
   v3d=a1pcon*v3d       !1-pion + contact
   ijk=0
   do i=1,npart-2
      do j=i+1,npart-1
         do k=j+1,npart
            ijk=ijk+1
            v3tmp1=czero
            do ic=1,3
               ic3=(ic+1)*3
               do jc=1,3
                  jc3=(jc+1)*3
                  do kc=1,3
                     kc3=(kc+1)*3
                     vfac=xpi(ic,i,levi(1,jc),j)*xpi(levi(2,jc),j,kc,k) &
                        -xpi(ic,i,levi(2,jc),j)*xpi(levi(1,jc),j,kc,k) &
                        +xpi(jc,j,levi(1,kc),k)*xpi(levi(2,kc),k,ic,i) &
                        -xpi(jc,j,levi(2,kc),k)*xpi(levi(1,kc),k,ic,i) &
                        +xpi(kc,k,levi(1,ic),i)*xpi(levi(2,ic),i,jc,j) &
                        -xpi(kc,k,levi(2,ic),i)*xpi(levi(1,ic),i,jc,j)
                     do it=1,3
                        do ks=1,4
                           do js=1,4
                              v3tmp1(:,js,ks)=v3tmp1(:,js,ks)+spx(:,ic3+it,i) &
                                  *(spx(js,jc3+levi(1,it),j) &
                                   *spx(ks,kc3+levi(2,it),k) &
                                   -spx(js,jc3+levi(2,it),j) &
                                   *spx(ks,kc3+levi(1,it),k))*vfac
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            v3c=v3c+sum(v3tmp1(:,:,:)*d3b(:,:,:,ijk))
         enddo
      enddo
   enddo
   v3c=a2p3bc*4.0_r8*v3c
   end subroutine

   function iv3bflg()
   integer(kind=i4) :: iv3bflg
   iv3bflg=iv3bflag
   end function iv3bflg
   
end module v3bpot
