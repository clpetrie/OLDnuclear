module phimod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci = (0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone = (1.0_r8,0.0_r8)

   integer(kind=i4), private, save :: jmax,npart,ntab,nrad,ndet
   real(kind=r8), private, save, allocatable :: orb(:,:),dorb(:,:),d2orb(:,:)
   integer(kind=i4), private, save, allocatable :: iorb(:,:),jpval(:,:)
   integer(kind=i4), private, save, allocatable :: mpval(:,:),lval(:,:)
   integer(kind=i4), private, save, allocatable :: isoval(:,:),lorb(:)
   real(kind=r8), private, save :: rmax
   real(kind=r8), private, save, allocatable :: sscale(:)

contains
   subroutine setphi(ntabin,jmaxin,nradin,npartin,ndetin,jpvalin,mpvalin, &
      lvalin,isovalin,iorbin,lorbin,orbin,rmaxin,fscal)
   integer(kind=i4) :: ntabin,jmaxin,nradin,npartin,ndetin
   integer(kind=i4) :: jpvalin(:,:),mpvalin(:,:),lvalin(:,:),isovalin(:,:)
   integer(kind=i4) :: iorbin(:,:),lorbin(:),i
   real(kind=r8) :: orbin(:,:),rmaxin,dr,fscal(:)
   real(kind=r8), allocatable :: ddr(:)
   ntab=ntabin
   jmax=jmaxin
   nrad=nradin
   npart=npartin
   ndet=ndetin
   rmax=rmaxin
   allocate(sscale(nrad),ddr(nrad))
   sscale=ntab/rmax
   ddr=1.0_r8/sscale
   ddr=fscal*ddr
   sscale=1.0_r8/ddr
   if (allocated(orb)) deallocate(orb,d2orb,jpval,mpval,lval,isoval,iorb)
   nrad=maxval(iorbin)
   allocate(orb(ntab,nrad),d2orb(ntab,nrad))
   allocate(jpval(npart,ndet),mpval(npart,ndet),lval(npart,ndet))
   allocate(isoval(npart,ndet),iorb(npart,ndet),lorb(nrad))
   orb=orbin
   do i=1,nrad
      dr=ddr(i)
!     anorm=0.0_r8
!     do j=1,ntab
!        r=j*dr
!        anorm=anorm+(r*orb(j,i))**2
!     enddo
!     orb(:,i)=orb(:,i)/sqrt(anorm*dr)
      orb(ntab,i)=orb(ntab-1,i)
      d2orb(1,i)=(11*orb(5,i)-56*orb(4,i)+114*orb(3,i)-104*orb(2,i) &
         +35*orb(1,i))/(12*dr*dr)
      d2orb(2,i)=(11*orb(6,i)-56*orb(5,i)+114*orb(4,i)-104*orb(3,i) &
         +35*orb(2,i))/(12*dr*dr)
      d2orb(3,i)=(11*orb(7,i)-56*orb(6,i)+114*orb(5,i)-104*orb(4,i) &
         +35*orb(3,i))/(12*dr*dr)
      d2orb(4:ntab-3,i)=(2*orb(7:ntab,i)-27*orb(6:ntab-1,i) &
         +270*orb(5:ntab-2,i)-490*orb(4:ntab-3,i)+270*orb(3:ntab-4,i) &
         -27*orb(2:ntab-5,i)+2*orb(1:ntab-6,i))/(180*dr*dr)
      d2orb(ntab-2,i)=(35*orb(ntab-2,i)-104*orb(ntab-3,i)+114*orb(ntab-4,i) &
         -56*orb(ntab-5,i)+11*orb(ntab-6,i))/(dr*dr)
      d2orb(ntab-1,i)=(35*orb(ntab-1,i)-104*orb(ntab-2,i)+114*orb(ntab-3,i) &
         -56*orb(ntab-4,i)+11*orb(ntab-5,i))/(dr*dr)
      d2orb(ntab,i)=(35*orb(ntab,i)-104*orb(ntab-1,i)+114*orb(ntab-2,i) &
         -56*orb(ntab-3,i)+11*orb(ntab-4,i))/(dr*dr)
   enddo
   jpval=jpvalin
   mpval=mpvalin
   lval=lvalin
   isoval=isovalin
   iorb=iorbin
   lorb=lorbin
!  do i=1,nrad
!     do j=1,ntab
!        r=j*ddr(i)
!        write(97,'(i5,3f25.15)') i,r,orb(j,i),d2orb(j,i)
!     enddo
!     write(97,*) ' '
!     write(97,*) ' '
!  enddo
   end subroutine setphi

   subroutine getphi(x,ph,dph)
   real(kind=r8) :: x(3)
   complex(kind=r8) :: ph(:,:,:),dph(:,:,:,:)
   complex(kind=r8) :: chi(2,2*jmax*(jmax+1)),dchi(3,2,2*jmax*(jmax+1))
   real(kind=r8) :: r,ri,dx(3),rx,drx,a,b,c,d,yi,yii,ddyi,ddyii,ril
   real(kind=r8) :: rad(nrad),radp(nrad)
   real(kind=r8) :: scale
   integer(kind=i4) :: i,j,l,ix,io,jm,it
   ph=czero
   dph=czero
   r=sqrt(sum(x**2))
   if (r.gt.rmax) return
   ri=1.0_r8/r
   dx=x*ri
   call chical(x,jmax,chi,dchi)
   do i=1,nrad
      scale=sscale(i)
      rx=r*scale
      ix=nint(rx)
      ix=max(1,min(ix,ntab-1))
      drx=rx-ix
      b=drx
      a=1.0_r8-drx
      c=a*(a**2-1.0_r8)/(6.0_r8*scale**2)
      d=b*(b**2-1.0_r8)/(6.0_r8*scale**2)
      yi=orb(ix,i)
      yii=orb(ix+1,i)
      ddyi=d2orb(ix,i)
      ddyii=d2orb(ix+1,i)
      rad(i)=a*yi+b*yii+c*ddyi+d*ddyii
      radp(i)=(yii-yi)*scale+((3*b**2-1)*ddyii-(3*a**2-1)*ddyi)/(6*scale)
      l=lorb(i)
      ril=ri**l
      radp(i)=ril*(radp(i)-l*ri*rad(i))
      rad(i)=ril*rad(i)
   enddo
   do j=1,ndet
      do i=1,npart
         l=jpval(i,j)+lval(i,j)
         jm=2*jpval(i,j)*l+mpval(i,j)
         it=isoval(i,j)
         io=iorb(i,j)
         ph(i,it,j)=rad(io)*chi(1,jm)
         ph(i,it+1,j)=rad(io)*chi(2,jm)
         dph(i,:,it,j)=dx(:)*radp(io)*chi(1,jm)+rad(io)*dchi(:,1,jm)
         dph(i,:,it+1,j)=dx(:)*radp(io)*chi(2,jm)+rad(io)*dchi(:,2,jm)
      enddo
   enddo
   end subroutine getphi

   subroutine chical(x,jpmax,chi,dchi)
   use ylmmod
   integer(kind=i4) :: jpmax
   real(kind=r8) :: x(3)
   complex(kind=r8) :: chi(2,2*jpmax*(jpmax+1)),dchi(3,2,2*jpmax*(jpmax+1))
   complex(kind=r8) :: ylm(0:(jpmax+1)*(jpmax+1)-1)
   complex(kind=r8) :: dylm(3,0:(jpmax+1)*(jpmax+1)-1)
   integer(kind=i4) :: jp,mp,k,lm
   real(kind=r8) :: cmup,cmdown,cpup,cpdown,ajp,emp
!
! Indices for |J,M, L=J-1/2,S=1/2> are given by j'=J+1/2 = 1,2,3,...,
! and m'=M+J+1 = 1,2,3,..., and the value is indexed as
! 2j'(j'-1)+m' = (2J+1)*(J-1/2)+M+J+1.
! Indices for |J,M, L=J+1/2,S=1/2> are 2j'^2+m'=2(J+1/2)^2+M+J+1
! The Clebsch-Gordon Coefficients are:
! <L=J-1/2,M-1/2 up|J,M, L=J-1/2,S=1/2> = sqrt((L+M+1/2)/(2L+1))
!                                       = sqrt((m'-1)/(2j'-1))
! <L=J-1/2,M+1/2 down |J,M, L=J-1/2,S=1/2> = sqrt((L-M+1/2)/(2L+1))
!                                          = sqrt((2j'-m')/(2j'-1))
! <L=J+1/2,M-1/2 up|J,M, L=J+1/2,S=1/2> = -sqrt((L-M+1/2)/(2L+1))
!                                       = -sqrt((2j'-m'+1)/(2j'+1))
! <L=J+1/2,M+1/2 down |J,M, L=J+1/2,S=1/2> = sqrt((L+M+1/2)/(2L+1))
!                                          = sqrt(m'/(2j'+1))
!
! jpmax = maximum value of j'=J+1/2.
!
   call ylmcal(x,jpmax,ylm,dylm)
!
! use complex conjugate to be consistent with wavenuclei psi_T(R,S)^*
!
   ylm=conjg(ylm)
   dylm=conjg(dylm)
   k=0
   do jp=1,jpmax
      ajp=jp
!
! L=J-1/2, M=-J special case then loop
!
      k=k+1
      lm=(jp-1)**2 ! l=jp-1, lm=l*(l+1)+m=(jp-1)*jp-(jp-1)=jp^2-2jp+1
      chi(1,k)=czero
      chi(2,k)=ylm(lm)
      dchi(:,1,k)=czero
      dchi(:,2,k)=dylm(:,lm)
      do mp=2,2*jp-1
         k=k+1
         emp=mp
         cmup=sqrt((emp-1.0_r8)/(2.0_r8*ajp-1.0_r8))
         cmdown=sqrt((2.0_r8*ajp-emp)/(2.0_r8*ajp-1.0_r8))
         lm=jp*(jp-2)+mp-1 ! l*(l+1)+m l=jp-1, m=M-1/2=mp-jp-1
         chi(1,k)=cmup*ylm(lm)
         chi(2,k)=cmdown*ylm(lm+1)
         dchi(:,1,k)=cmup*dylm(:,lm)
         dchi(:,2,k)=cmdown*dylm(:,lm+1)
      enddo
!
! special case L=J-1/2, M=J
!
      k=k+1
      lm=jp**2-1
      chi(1,k)=ylm(lm)
      chi(2,k)=czero
      dchi(:,1,k)=dylm(:,lm)
      dchi(:,2,k)=czero
!
! L=J+1/2
!
      do mp=1,2*jp
         emp=mp
         cpup=-sqrt((2.0_r8*ajp-emp+1.0_r8)/(2.0_r8*ajp+1.0_r8))
         cpdown=sqrt(emp/(2.0_r8*ajp+1.0_r8))
         k=k+1
         lm=jp*jp+mp-1  ! l*(l+1)+m l=jp, m=M-1/2 = mp-jp-1
         chi(1,k)=cpup*ylm(lm)
         chi(2,k)=cpdown*ylm(lm+1)
         dchi(:,1,k)=cpup*dylm(:,lm)
         dchi(:,2,k)=cpdown*dylm(:,lm+1)
      enddo
   enddo
   end subroutine chical

end module phimod
