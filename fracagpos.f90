!    This code generates the position coordinates x_i, y_i, z_i, i=1,npart,
!    for the npart spheres in a fractal aggregate.   The fractal scaling
!    is given by
!
!       npart = k0 (Rg/a)^df
!
!    where k0: structure factor, df: fractal dimension,, a: sphere radius,
!    Rg: radius of gyration:
!
!    Rg^2 = (1/npart) sum r_i^2
!
!    with r_i: distance from ith sphere to center of mass of cluster.
!
!    the sphere radii are taken to be unity in the code.  The positions
!    are therefore given in units of sphere radii, and two touching spheres
!    will be two units apart.
!
!    The code uses a pseudo-random algorithm to mimic cluster-cluster
!    aggregation, subject to the constraint that the fractal scaling is
!    identically satisfied for the cluster.
!
!    input parameters:  all should be obvious except nsamp:  this is the
!    initial number of spheres from which the npart--sized custer will be
!    generated.   nsamp should be >= npart.   I've used nsamp around
!    1 - 2 * npart.
!
!    iccmod=0: the code uses a cluster-cluster algorithm.
!    iccmod=1: the original DLA algorithm:  one sphere at a time
!    is added to the cluster.    Not realistic of typical soot dynamics.


      implicit none
      integer :: npart,nsamp,iccmod,iseed,i,more
      real(4) :: rk,df,rkf1
      real(4), allocatable :: xp(:,:)
      character :: fout*30

10    write(*,'('' npart, nsamp, k0, df, iccmod, output file:'',$)')
      read(*,*) npart, nsamp, rk,df,iccmod, fout
      !write(*,'('' output file:'',$)')
      !read(*,'(a)') fout

      if(fout.ne.' ') open(1,file=fout)
      allocate(xp(3,nsamp))
      iseed=0
      rkf1=(1./rk)**(1./df)/2.d0

      call ppclus(npart,nsamp,rkf1,df,iccmod,xp,iseed)

      if(fout.ne.' ') then
         do i=1,npart
            write(1,'(4f10.4)') 1.d0,xp(1,i),xp(2,i),xp(3,i)
         enddo
         close(1)
      endif

      end
!
!  random number generator from numerical recipes, Press et al.
!

      function ran3(idum)
      implicit none
      integer :: idum
      integer :: mbig,mseed,mz
      real(4) ::  fac, ran3
      parameter(mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
      integer :: i,ii,k
      integer :: mj,mk
      integer, save:: iff,inext,inextp,ma(55)
      data iff/0/
      if(idum.lt.0.or.iff.eq.0) then
         iff=1
         mj=mseed-iabs(idum)
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
         enddo
         do k=1,4
            do i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
            enddo
         enddo
         inext=0
         inextp=31
         idum=1
      endif
 20   inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      if(ran3.gt.1.e-9.and.ran3.lt.0.99999999) return
      goto 20
      end function ran3

      subroutine ppclus(nptot,nsamp,rkf1,df,iccmod,xp,iseed0)
      implicit none
      integer :: nptot,nsamp,iccmod,iseed0,iarray(3),iadd(nsamp),npc(nsamp), &
                 ihr,imin,isec,ihsec,nc,i,j,iaddic,iaddjc,npic,npjc,k,nk,ic, &
                 iseed,jc,jp,irun
      real(4) :: ran,ran3,xp(3,nsamp),xpc(3,nsamp),xpt(3,nsamp),rkf1,df
!
! the intrinsic time function is used to seed the random generator
!
! insert the corresponding function call for your compiler
!
      if(iseed0.le.0) then
!         call gettim(ihr,imin,isec,ihsec)
!         iseed=360000*ihr+6000*imin+1000*isec+ihsec
        call itime(iarray)
        iseed=iarray(1)+iarray(2)+iarray(3)
      else
         iseed=iseed0
      endif

      ran=ran3(iseed)

      do irun=1,1000
         do i=1,nsamp
            xpt(:,i)=0.d0
            iadd(i)=i
            npc(i)=1
         enddo
         nc=nsamp
         do while (nc.gt.1)
            i=1
            j=i
            do while(i.eq.j)
               if(iccmod.eq.0) then
                  i=int(real(nc)*ran3(1))+1
               else
                  i=1
               endif
               j=int(real(nc)*ran3(1))+1
            enddo
            ic=min(i,j)
            jc=max(i,j)
            iaddic=iadd(ic)
            iaddjc=iadd(jc)
            npic=npc(ic)
            npjc=npc(jc)

            do i=1,npjc
               xpc(:,i)=xpt(:,i+iaddjc-1)
            enddo

            do k=jc-1,ic+1,-1
               do nk=npc(k),1,-1
                  i=nk+iadd(k)-1
                  j=i+npjc
                  xpt(:,j)=xpt(:,i)
               enddo
               iadd(k)=iadd(k)+npjc
            enddo

            call combine(npic,xpt(1,iaddic),npjc,xpc,rkf1,df)
            npc(ic)=npic+npjc

            if(npc(ic).eq.nptot) then
               do i=1,nptot
                  j=i+iaddic-1
                  xp(:,i)=xpt(:,j)
               enddo
               return
            endif

            do j=jc,nc-1
               jp=j+1
               iadd(j)=iadd(jp)
               npc(j)=npc(jp)
            enddo
            nc=nc-1
         enddo
      enddo
      write(*,'('' process did not work.  Try again'')')

      return
      end subroutine ppclus

      subroutine combine(np1,xp,npc,xp2,rkf1,df)
      implicit none

      integer :: np1,npc,nits,it,np3,i,m,ic1,ic2
      real(4) :: xp(3,*),xp2(3,*),ran3,xc(3),xpc(3,np1+npc),x1(3),x2(3), &
            rc1(3),rc2(3),pi,ct,phi,st,rkf1,df,ctc,stc,pc,rgc,rg1,rg3,c, &
            r12,r1,r2,beta,alpha,gamma,tc,betac,r12c,fatan

      nits=10000
      pi=4.*atan(1.d0)

      if(np1.eq.1.and.npc.eq.1) then
         ct=1.-2.*ran3(1)
         phi=2.*pi*ran3(1)
         st=sqrt((1.-ct)*(1.+ct))
         xp(1,1)=st*cos(phi)
         xp(2,1)=st*sin(phi)
         xp(3,1)=ct
         xp(:,2)=-xp(:,1)
         return
      endif

      if(npc.eq.1) then
         call addone(np1,xp,rkf1,df)
         return
      endif

      rgc=dble(npc)**(1./df)*2.*rkf1
      np3=np1+npc
      rg1=dble(np1)**(1./df)*2.*rkf1
      rg3=dble(np3)**(1./df)*2.*rkf1
      c=sqrt(dble(np3)*(np3*rg3*rg3-np1*rg1*rg1-npc*rgc*rgc) &
         /dble(npc*(np3-npc)))

      do it=1,nits
         ctc=1.d0-2.*ran3(1)
         stc=sqrt((1.d0-ctc)*(1.d0+ctc))
         pc=2.*pi*ran3(1)
         xc(1)=c*stc*cos(pc)
         xc(2)=c*stc*sin(pc)
         xc(3)=c*ctc
         do i=1,npc
            xpc(:,i)=xp2(:,i)+xc(:)
         enddo

         call finmin(0,np1,xp,npc,xpc,ic1,ic2,r12)

         if(r12.gt.4.d0) cycle

         r1=0.
         r2=0.
         x1(:)=xp(:,ic1)-xc(:)
         x2(:)=xp2(:,ic2)
         r1=dot_product(x1,x1)
         r2=dot_product(x2,x2)
         r1=sqrt(r1)
         r2=sqrt(r2)

         if(abs(r1-r2).gt.2.d0) cycle

         call ctos(x2,rc2)
         beta=acos(rc2(2))
         alpha=rc2(3)
         call rotate(0,alpha,beta,0.,x1)
         gamma=fatan(x1(2),x1(1))

         call rotate(0,0.,0.,gamma,x1)
         call ctos(x1,rc1)
         ctc=(rc1(1)*rc1(1)+rc2(1)*rc2(1)-4.d0)/(2.*rc1(1)*rc2(1))
         if(abs(ctc).gt.1.d0) cycle
         tc=acos(ctc)
         betac=tc-acos(rc1(2))

         do i=1,npc
            call rotate(0,alpha,beta,gamma,xp2(1,i))
            call rotate(0,0.,betac,0.,xp2(1,i))
            call rotate(1,alpha,beta,gamma,xp2(1,i))
            do m=1,3
               xpc(m,i)=xp2(m,i)+xc(m)
            enddo
         enddo

         call finmin(0,np1,xp,npc,xpc,ic1,ic2,r12c)

         if(r12c.ge.2.d0) exit

      enddo

      if(r12c.lt.2.d0) then
         write(*,'('' clusters did not combine'')')
      endif

      do i=1,npc
         xp(:,i+np1)=xpc(:,i)
      enddo

      xc=0.d0
      do i=1,np3
         xc(:)=xc(:)+xp(:,i)
      enddo
      xc=xc/dble(np3)
      do i=1,np3
         xp(:,i)=xp(:,i)-xc(:)
      enddo

      end subroutine combine

      function fatan(y,x)
      implicit none
      real(4) :: x,y,fatan
      if(x.eq.0..and.y.eq.0.) then
         fatan=0.
      else
         fatan=atan2(y,x)
      endif
      return
      end function fatan

      subroutine ctos(x,r)
      implicit none
      real(4) :: x(3),r(3),fatan
      r(1)=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      if(r(1).eq.0.) then
         r(2)=1.d0
         r(3)=0.d0
      else
        r(2)=x(3)/r(1)
        r(3)=fatan(x(2),x(1))
      endif
      return
      end subroutine ctos

      subroutine rotate(idir,alpha,beta,gamma,x)
      implicit none
      integer :: idir
      real(4) x(3),xt(3),alpha,beta,gamma,sa,ca,sb,cb,sg,cg

      sa=sin(alpha)
      ca=cos(alpha)
      sb=sin(beta)
      cb=cos(beta)
      sg=sin(gamma)
      cg=cos(gamma)

      if(idir.eq.0) then
         xt(1)=(ca*cb*cg-sa*sg)*x(1)+(cb*cg*sa+ca*sg)*x(2) &
               -cg*sb*x(3)
         xt(2)=(-cg*sa-ca*cb*sg)*x(1)+(ca*cg-cb*sa*sg)*x(2) &
               +sb*sg*x(3)
         xt(3)=ca*sb*x(1)+sa*sb*x(2)+cb*x(3)
      else
         xt(1)=(ca*cb*cg-sa*sg)*x(1)-(cb*sg*ca+sa*cg)*x(2) &
               +ca*sb*x(3)
         xt(2)=(sg*ca+sa*cb*cg)*x(1)+(ca*cg-cb*sa*sg)*x(2) &
               +sb*sa*x(3)
         xt(3)=-cg*sb*x(1)+sg*sb*x(2)+cb*x(3)
      endif

      x=xt
      end subroutine rotate

      subroutine finmin(isame,np1,xp,np2,xpc,ic1,ic2,r12)
      implicit none
      integer :: isame,np1,np2,ic1,ic2,i,j,m
      real(4) :: xp(3,*),xpc(3,*),rmin,rij,r12,x12

      rmin=10000.
      do i=1,np1
         do j=1,np2
            if(isame.eq.0.or.i.ne.j) then
               rij=0.
               do m=1,3
                  x12=xp(m,i)-xpc(m,j)
                  rij=rij+x12*x12
               enddo
               if(rij.lt.rmin) then
                  rmin=rij
                  ic1=i
                  ic2=j
               endif
            endif
         enddo
      enddo
      r12=sqrt(rmin)
      end subroutine finmin

      subroutine addone(nptot,xp,rkf1,df)
      implicit none
      integer :: nptot,ijp(nptot),itmax,np3,nj,i,j,it,icon,m,k,iseed,ij,n
      real(4) :: ran3,xp(3,*),rp(nptot),rkf1,df,rgn2,rmax,ri,rgn,rg3,rn2,rn, &
                 rg32,rj2,rj,rjn,rt,ctj,stj,phij,sphij,cphij,ctp,stp,ran1,phi,zpp, &
                 xpp,ypp,z,x,y,xi,yi,zi,ri2,x0,y0,z0

      itmax=20000

      rgn2=0.
      rmax=0.
      do i=1,nptot
         ri=0.
         do m=1,3
            ri=ri+xp(m,i)*xp(m,i)
         enddo
         rgn2=rgn2+ri
         rmax=max(rmax,ri)
      enddo
      rgn2=rgn2/dble(nptot)
      rgn=sqrt(rgn2)
      rmax=sqrt(ri)

      np3=nptot+1
      rg3=dble(np3)**(1./df)*2.*rkf1
      rn2=np3*(np3/real(nptot)*rg3*rg3-rgn2)
      rn=sqrt(rn2)
!
!  check if rn is too big
!
      if(rn-rmax.gt.2.) then
         rn=rmax+1.8
         rn2=rn*rn
         rg32=(rn2/real(np3)+rgn2)*(nptot)/real(np3)
         rg3=sqrt(rg32)
      endif
!
!  find particles that intersect with rn
!
   20 i=1
      do j=1,nptot
         rj2=0.
         do k=1,3
            rj2=rj2+xp(k,j)*xp(k,j)
         enddo
         rj=sqrt(rj2)
         rjn=abs(rj-rn)
         if(rjn.le.2.) then
            ijp(i)=j
            rp(i)=rj
            i=i+1
         endif
      enddo
      nj=i-1
      do i=1,nj
         ran1=ran3(iseed)
         j=nj*ran1+1
         it=ijp(i)
         ijp(i)=ijp(j)
         ijp(j)=it
         rt=rp(i)
         rp(i)=rp(j)
         rp(j)=rt
      enddo

      outerloop : do ij=1,nj
         j=ijp(ij)
         rj=rp(ij)
         rj2=rj*rj
         if(rj+rn.lt.2.) cycle
         ctj=xp(3,j)/rj
         stj=sqrt((1.-ctj)*(1.+ctj))
         phij=atan2(xp(2,j),xp(1,j))
         sphij=sin(phij)
         cphij=cos(phij)
         ctp=(rn2+rj2-4.)/2./rn/rj
         stp=sqrt((1.-ctp)*(1.+ctp))
!
!  randomly find attachment point
!
         it=0
         do while(it.lt.itmax)
            it=it+1
            ran1=ran3(iseed)
            phi=2.*3.141592654*ran1
            zpp=rn*ctp
            xpp=rn*stp*cos(phi)
            ypp=rn*stp*sin(phi)
            z=zpp*ctj-xpp*stj
            x=(zpp*stj+xpp*ctj)*cphij-ypp*sphij
            y=(zpp*stj+xpp*ctj)*sphij+ypp*cphij
            icon=0
            do i=1,nptot
               xi=x-xp(1,i)
               yi=y-xp(2,i)
               zi=z-xp(3,i)
               ri2=xi*xi+yi*yi+zi*zi
               if(ri2.lt.3.999) icon=1
            enddo
            if(icon.eq.0) exit outerloop
         enddo
      enddo outerloop

      if(ij.gt.nj.and.icon.ne.0) then
         rn=rn+0.01
         rn2=rn*rn
         rg32=(rn2/real(np3)+rgn2)*(nptot)/real(np3)
         rg3=sqrt(rg32)
         write(*,'(''+particle '',i3,'' does not fit.'', &
                 &'' new rn, rgn:'',2f8.2)') n,rn,rgn
         goto 20
      endif
!
!  shift coordinates to new origin
!
      xp(1,np3)=x
      xp(2,np3)=y
      xp(3,np3)=z
      x0=0.
      y0=0.
      z0=0.
      do i=1,np3
         x0=x0+xp(1,i)
         y0=y0+xp(2,i)
         z0=z0+xp(3,i)
      enddo
      x0=x0/real(np3)
      y0=y0/real(np3)
      z0=z0/real(np3)
      do i=1,np3
         xp(1,i)=xp(1,i)-x0
         xp(2,i)=xp(2,i)-y0
         xp(3,i)=xp(3,i)-z0
      enddo
      end subroutine addone
