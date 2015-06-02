module v3bpot
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   integer(kind=i4), private, parameter :: levi(2,3) = &
      reshape((/2,3, 3,1, 1,2/),(/2,3/))
   real(kind=r8), private, parameter :: hbarc=197.327_r8 ! MeV*fm
   real(kind=r8), private, parameter :: mpi=138.03 ! MeV
   real(kind=r8), private, parameter :: fpi=92.4_r8 ! MeV
   real(kind=r8), private, parameter :: ga=1.267_r8
   character(len=3), private, save :: v3b
   integer(kind=i4), private, save :: npart,ncut,ncutt
   real(kind=r8), private, save :: el,eli,pi
   real(kind=r8), private, save :: c3b,amu,a2p3b,a2s3b,avd,ave,ar3b,dfac,dnorm
   real(kind=r8), private, save :: acfac,a2pscal,a2pxdscal,a2pddscal
   integer(kind=i4), private, save :: dov3
contains
   subroutine v3bpotinit(npartin,elin,v3bin,dov3in,a2pscalin,a2pxdscalin,a2pddscalin)
   use mympi
   integer(kind=i4) :: npartin,dov3in
   real(kind=r8), intent(in) :: elin,a2pscalin,a2pxdscalin,a2pddscalin
   real(kind=r8) :: lamchi,c1,c3,c4,cd,ce,r0
   character(len=3), intent(in) :: v3bin
   npart=npartin
   el=elin
   if (el.gt.0.0_r8) then
      eli=1.0_r8/el
   else
      el=0.0_r8
      eli=0.0_r8
   endif
   v3b=v3bin
   pi=4.0_r8*atan(1.0_r8)
   if(v3b.eq.'no3') then
      ncut=1
      ncutt=1
      c3b=0.0_r8
      amu=0.0_r8
      a2p3b=0.0_r8
      a2s3b=0.0_r8
      avd=0.0_r8
      ave=0.0_r8
      ar3b=0.0_r8
      acfac=0.0_r8
      dfac=0.0_r8
      dnorm=0.0_r8
      dov3=0
   else if(v3b.eq.'uix') then
      ncut=2
      ncutt=2
      c3b=2.1_r8
      amu=0.69953054_r8
      a2p3b=-0.0293_r8
      a2s3b=0.0_r8
      avd=0.0_r8
      ave=0.0_r8
      ar3b=0.00480_r8
      acfac=0.25_r8
      dfac=1.0_r8
      dnorm=0.0_r8
      dov3=1
   else if(v3b.eq.'eft') then
      ncut=4
      ncutt=1
      r0=1.0_r8 ! cutoff
      c3b=(1.0_r8/r0)**ncut
      amu=0.69949880_r8
      lamchi=700.0_r8 ! MeV
      c1=-0.00081_r8 ! MeV**-1
      c3=-0.0034_r8 ! MeV**-1
      c4=0.0034_r8 ! MeV**-1
      cd=-1.0_r8 ! something
      ce=-0.1_r8 ! something
 cd=0.7_r8
 ce=-0.63_r8
      a2p3b=0.5_r8*(ga/fpi**2)**2*(1.0_r8/(4.0_r8*pi))**2*mpi**6/9.0_r8*c3/4.0_r8
      a2s3b=(ga/(2.0_r8*fpi))**2*(1.0_r8/(4.0_r8*pi))**2*4.0_r8*mpi**6/fpi**2*c1
      avd=mpi**3/(12.0_r8*pi)*ga/(8.0_r8*fpi**2)*1.0_r8/(fpi**2*lamchi)*cd*hbarc**3
      ave=ce/(fpi**4*lamchi)*hbarc**6
      ar3b=0.0_r8
      acfac=-c4/(2.0_r8*c3)
      dfac=4.0_r8*pi/mpi**3*hbarc**3
      dnorm=ncut/(4.0_r8*pi*gamma(3.0_r8/ncut)*r0**3)
      dov3=2
   else if(v3b.eq.'rd3') then
      if (myrank().eq.0) then
         read (5,*) ncut
         read (5,*) ncutt
         read (5,*) c3b
         read (5,*) a2p3b
         read (5,*) a2s3b
         read (5,*) ar3b
         read (5,*) dfac
      endif
      call bcast(ncut)
      call bcast(ncutt)
      call bcast(c3b)
      call bcast(a2p3b)
      call bcast(a2s3b)
      call bcast(ar3b)
      call bcast(dfac)
   else
      if (myrank().eq.0) then
         write (6,'(1x,''serious error pot3bpar'')')
         write (6,'(1x,''bad three body potential name  '',a3)') v3b
         write (6,'(1x,''allowed names: '')')
         write (6,'(16x,a3)') 'no3'
         write (6,'(16x,a3)') 'uix'
         write (6,'(16x,a3)') 'eft'
         write (6,'(16x,a21)') 'rd3 to read in values'
         write (6,'(1x,''stopping the calculation '')')
         call abort
      else
         call barrier !wait for 0 to print
      endif
   endif
   a2pscal=a2pscalin
   a2pxdscal=a2pxdscalin
   a2pddscal=a2pddscalin
   if (myrank().eq.0) then
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''Do V3'',t40,l10)') dov3.gt.0
      write (6,'(''Three body force'',t47,a3)') v3b
      write (6,'(''three body a2piPWa  ='',t40,f20.12)') a2p3b
      write (6,'(''three body a2piPwc ='',t40,f20.12)') a2p3b*acfac
      write (6,'(''three body a2piSw ='',t40,f20.12)') a2s3b
      write (6,'(''three body Vd,C,c3 ='',t40,f20.12)') 4.0_r8*a2p3b*dfac !dfac included in the delta function 
      write (6,'(''three body Ve,C,c3 ='',t40,f20.12)') 4.0_r8*a2p3b*dfac**2 !dfac**2 included in the delta function
      write (6,'(''three body A_R ='',t40,f20.15)') ar3b
      write (6,'(''three body Vd1 ='',t40,f20.12)') avd
      write (6,'(''three body Vd2 ='',t40,f20.12)') 2.0_r8*dfac*avd
      write (6,'(''three body Ve ='',t40,f20.12)') ave
      write (6,'(''propagator: multiply the antic by this ='',t40,f20.15)') a2pscal
      write (6,'(''propagator: multiply the antic-xd by this ='',t40,f20.15)') a2pxdscal
      write (6,'(''propagator: multiply the antic-dd by this ='',t40,f20.15)') a2pddscal
      write (6,'(''propagator three body a2piPW ='',t40,f10.5)') a2p3b*a2pscal
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
   endif
   dov3in=dov3
   end subroutine v3bpotinit

   subroutine hstnipot(x,vc,xpi,g2s3b,delta,a2pin,a2sin,avdin,avein,acfacin,a2pscalin,a2pxdscalin,a2pddscalin,dfacin,rscal)
   real(kind=r8) :: x(:,:)
   real(kind=r8) :: vc,xpi(3,npart,3,npart),xxpi(3,npart,3,npart),cut
   real(kind=r8) :: delta(npart,npart),dx(3),r
   real(kind=r8) :: gr3b(npart),g2s3b(3,npart,npart),exc,y2,t2,z2
   integer(kind=i4) :: i,j
   real(kind=r8) :: a2pin,a2sin,avdin,acfacin,a2pscalin,a2pxdscalin,a2pddscalin,rscal,avein,dfacin
   vc=0.0_r8
   gr3b=0.0_r8
   g2s3b=0.0_r8
   xpi=0.0_r8
   xxpi=0.0_r8
   delta=0.0_r8
   a2pin=a2p3b
   a2sin=a2s3b
   avdin=avd !/(4.0_r8*pi/mpi**3*hbarc**3)
   avein=ave
   acfacin=acfac
   a2pscalin=a2pscal
   a2pxdscalin=a2pxdscal
   a2pddscalin=a2pddscal
   dfacin=dfac
   if (dov3.eq.0) return
   do i=1,npart-1
      do j=i+1,npart
         dx=x(:,i)-x(:,j)
         dx=rscal*dx
         dx=dx-el*nint(dx*eli)
         r=sqrt(sum(dx**2))
         dx=dx/r
         exc=exp(-c3b*r**ncut)
         cut=1.0_r8-exc
         y2=cut*exp(-amu*r)/(amu*r)
         t2=(1.0_r8+3.0_r8/(r*amu)+3.0_r8/(r*amu)**2)*y2*cut**(ncutt-1)
         z2=amu*r*(y2-t2)/3.0_r8
         delta(i,j)=dfac*dnorm*exc
         delta(j,i)=delta(i,j)
         xpi(:,i,1,j)=3.0_r8*t2*dx(1)*dx
         xpi(:,i,2,j)=3.0_r8*t2*dx(2)*dx
         xpi(:,i,3,j)=3.0_r8*t2*dx(3)*dx
         xpi(1,i,1,j)=xpi(1,i,1,j)+(y2-t2)
         xpi(2,i,2,j)=xpi(2,i,2,j)+(y2-t2)
         xpi(3,i,3,j)=xpi(3,i,3,j)+(y2-t2)
         xpi(:,j,:,i)=xpi(:,i,:,j)
         g2s3b(:,j,i)=dx*z2
         g2s3b(:,i,j)=-g2s3b(:,j,i)
         gr3b(i)=gr3b(i)+t2**2
         gr3b(j)=gr3b(j)+t2**2
         vc=vc-ar3b*t2**4
      enddo
   enddo
   vc=vc+0.5_r8*ar3b*sum(gr3b**2)
   end subroutine hstnipot

   subroutine hstnimat(x,vc,ast,asttm,astvd1,astvd2,astvdc3,astvec3,atau,a2pscalin,a2pxdscalin,a2pddscalin,rscal)
   real(kind=r8) :: x(3,npart),rscal
   real(kind=r8), dimension (3,npart,3,npart) :: ast,asttm,astvd1,astvd2,astvdc3,astvec3
   real(kind=r8) :: atau(npart,npart)
   real(kind=r8) :: vc,xpi(3,npart,3,npart),xxpi(3,npart,3,npart)
   real(kind=r8) :: xd(3,npart,3,npart),xdel(3,npart,3,npart),ddelta(npart,npart)
   real(kind=r8) :: delta(npart,npart)
   real(kind=r8) :: g2s3b(3,npart,npart)
   real(kind=r8) :: delsum(npart)
   real(kind=r8) :: a2pscalin,a2pxdscalin,a2pddscalin
   real(kind=r8) :: dummy
   integer(kind=i4) :: i,j,k,jc,kc,ic
   if (dov3.eq.0) then
      ast=0.0_r8
      asttm=0.0_r8
      astvd1=0.0_r8
      astvd2=0.0_r8
      astvdc3=0.0_r8
      astvec3=0.0_r8
      atau=0.0_r8
      return
   endif
   call hstnipot(x,vc,xpi,g2s3b,delta,dummy,dummy,dummy,dummy,dummy,a2pscalin,a2pxdscalin,a2pddscalin,dummy,rscal)
   delsum=0.0_r8
   xd=0.0_r8
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
   xxpi=0.0_r8
   xxpi=reshape(matmul(reshape(xpi,(/3*npart,3*npart/)), &
      reshape(xpi,(/3*npart,3*npart/))),(/3,npart,3,npart/))
   xdel=0.0_r8
   ddelta=0.0_r8
   if (dov3.eq.2) then
      xdel=reshape(matmul(reshape(xpi,(/3*npart,3*npart/)), &
         reshape(xd,(/3*npart,3*npart/))),(/3,npart,3,npart/))
      xdel=xdel+reshape(matmul(reshape(xd,(/3*npart,3*npart/)), &
         reshape(xpi,(/3*npart,3*npart/))),(/3,npart,3,npart/))
      ddelta=matmul(delta,delta)
   endif
   ast=0.0_r8
   asttm=0.0_r8
   astvd1=0.0_r8
   astvd2=0.0_r8
   astvdc3=0.0_r8
   astvec3=0.0_r8
   atau=0.0_r8
   do j=2,npart
      do k=1,j-1
         do jc=1,3
            do kc=1,3
               ast(kc,k,jc,j)=ast(kc,k,jc,j)+4.0_r8*a2p3b*xxpi(kc,k,jc,j)  ! anticommutator
               if (dov3.eq.2) then
                  astvdc3(kc,k,jc,j)=astvdc3(kc,k,jc,j)+4.0_r8*a2p3b*xdel(kc,k,jc,j)  ! Vd,Cc3
                  asttm(kc,k,jc,j)=asttm(kc,k,jc,j)+a2s3b*sum(g2s3b(jc,:,j)*g2s3b(kc,:,k))  ! Tucson-Melbourne
                  astvd1(kc,k,jc,j)=astvd1(kc,k,jc,j)+avd/dfac*xpi(kc,k,jc,j)*(delsum(k)+delsum(j)-2.0_r8*delta(k,j)) ! Vd,1 - Kevin
!                 astvd1(kc,k,jc,j)=astvd1(kc,k,jc,j)+avd/dfac*xdel(kc,k,jc,j) ! Vd,1 - Ingo/Joel
                  if (kc.eq.jc) then
                     astvec3(kc,k,jc,j)=astvec3(kc,k,jc,j)+4.0_r8*a2p3b*ddelta(k,j) ! Ve,Cc3
                     astvd2(kc,k,jc,j)=astvd2(kc,k,jc,j)-avd/dfac*delta(k,j)*(delsum(k)+delsum(j)-2.0_r8*delta(k,j)) ! Vd,2 - Kevin
!                    astvd2(kc,k,jc,j)=astvd2(kc,k,jc,j)-avd/dfac*ddelta(k,j) ! Vd,2 - Ingo
                  endif
               endif
            enddo
         enddo
         if (dov3.eq.2) atau(k,j)=atau(k,j)+ave/dfac**2*sum(delta(k,:)*delta(j,:)) ! Ve
      enddo
   enddo
   do i=2,npart
      do ic=1,3
         ast(ic,i,:,1:i-1)=ast(:,1:i-1,ic,i)
         asttm(ic,i,:,1:i-1)=asttm(:,1:i-1,ic,i)
         astvd1(ic,i,:,1:i-1)=astvd1(:,1:i-1,ic,i)
         astvd2(ic,i,:,1:i-1)=astvd2(:,1:i-1,ic,i)
         astvdc3(ic,i,:,1:i-1)=astvdc3(:,1:i-1,ic,i)
         astvec3(ic,i,:,1:i-1)=astvec3(:,1:i-1,ic,i)
      enddo
      atau(i,1:i-1)=atau(1:i-1,i)
   enddo
   end subroutine hstnimat
end module v3bpot
