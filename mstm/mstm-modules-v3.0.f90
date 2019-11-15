!
!  numerical constants
!
!
!  last revised: 15 January 2011
!
      module numconstants
      implicit none
      integer :: print_intermediate_results
      integer, allocatable :: monen(:)
      integer, private :: nmax=0
      real(8) :: pi
      real(8), allocatable :: bcof(:,:),fnr(:),vwh_coef(:,:,:,:)
      real(8), allocatable :: vcc_const(:,:,:),fnm1_const(:,:),fn_const(:,:),fnp1_const(:,:)
      data pi/3.141592653589793/

      contains

         subroutine init(notd)
         implicit none
         integer :: notd,l,n,ierr,nbc,m,mm1,mp1,np1,nm1,nn1
         real(8) :: fnorm1,fnorm2
!
!  bcof(n,l)=((n+l)!/(n!l!))^(1/2)
!
         if(notd.le.nmax) return
         nmax=max(nmax,notd)
         nbc=6*notd+6
         if(allocated(fnr)) deallocate(monen,fnr,bcof)
         allocate (monen(0:2*notd),bcof(0:nbc,0:nbc),fnr(0:2*nbc),stat=ierr)
!         write(*,'('' nmax, bcof status:'',2i5)') nmax,ierr
         do n=0,2*notd
            monen(n)=(-1)**n
         enddo
         fnr(0)=0.d0
         do n=1,2*nbc
            fnr(n)=dsqrt(dble(n))
         enddo
         bcof(0,0)=1.d0
         do n=0,nbc-1
            do l=n+1,nbc
               bcof(n,l)=fnr(n+l)*bcof(n,l-1)/fnr(l)
               bcof(l,n)=bcof(n,l)
            enddo
            bcof(n+1,n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n,n)/fnr(n+1)/fnr(n+1)
         enddo
         if(allocated(vwh_coef)) deallocate(vwh_coef)
         allocate(vwh_coef(-notd:notd,1:notd,-1:1,-1:1))
!
!  constants used for calculation of svwf functions.
!
         do n=1,notd
            nn1=n*(n+1)
            np1=n+1
            nm1=n-1
            fnorm1=-.5d0/fnr(n+n+1)/fnr(n)/fnr(n+1)
            fnorm2=-.5d0*fnr(n+n+1)/fnr(n)/fnr(n+1)
            m=-n
            mp1=m+1
            mm1=m-1
            vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
            vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
            vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
            vwh_coef(m,n,-1,-1)=0.d0
            vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
            vwh_coef(m,n, 0,-1)=0.d0
            vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
            vwh_coef(m,n,-1, 0)=-0.d0
            vwh_coef(m,n, 0, 0)=-fnorm2*m
            do m=-n+1,-1
               mp1=m+1
               mm1=m-1
               vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
               vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
               vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
               vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
               vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
               vwh_coef(m,n, 0,-1)=fnorm1*np1*fnr(n+m)*fnr(n-m)
               vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
               vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
               vwh_coef(m,n, 0, 0)=-fnorm2*m
            enddo
            do m=0,n-1
               mp1=m+1
               mm1=m-1
               vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
               vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
               vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
               vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
               vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
               vwh_coef(m,n, 0,-1)=fnorm1*np1*fnr(n+m)*fnr(n-m)
               vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
               vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
               vwh_coef(m,n, 0, 0)=-fnorm2*m
            enddo
            m=n
            mp1=m+1
            mm1=m-1
            vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
            vwh_coef(m,n, 1,-1)=0.d0
            vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
            vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
            vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
            vwh_coef(m,n, 0,-1)=0.d0
            vwh_coef(m,n, 1, 0)=-0.d0
            vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
            vwh_coef(m,n, 0, 0)=-fnorm2*m
         enddo
         end subroutine init

      end module numconstants
!
!  special function for the multiple sphere problem
!
      module specialfuncs
      implicit none
      contains

         subroutine timewrite(iunit,char1,time)
         use intrinsics
         implicit none
         integer :: iunit
         real(8) :: time,time2
         character(*) :: char1
         if(time.gt.3600.d0) then
            time2=time/3600.d0
            write(iunit,'(a,f9.3,'' hours'')') char1,time2
         elseif(time.gt.60.d0) then
            time2=time/60.d0
            write(iunit,'(a,f9.2,'' min'')') char1,time2
         else
            write(iunit,'(a,f9.2,'' sec'')') char1,time
         endif
         call flush(iunit)
         end subroutine timewrite
!
!  ricatti-bessel function psi(n), real argument
!
         subroutine ricbessel(n,ds,eps,nmax,psi)
         implicit none
         integer :: n,nmax,ns,i
         real(8) :: ds,dns,sn,psi(0:n),psit,ds2,sum,eps,err
         if(int(ds).lt.n) then
            ns=nint(ds+4.*(ds**.3333d0)+17)
            ns=max(n+10,ns)
            dns=0.d0
            do i=ns-1,n,-1
               sn=dble(i+1)/ds
               dns=sn-1.d0/(dns+sn)
            enddo
            psi(n)=dns
            psi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
            do i=n-2,1,-1
               sn=dble(i+1)/ds
               psi(i)=sn-1.d0/(psi(i+1)+sn)
            enddo
            psit=dsin(ds)
            psi(0)=psit
            ds2=ds*ds
            sum=psit*psit/ds2
            do i=1,n
               psit=psit/(dble(i)/ds+psi(i))
               sum=sum+dble(i+i+1)*psit*psit/ds2
               err=dabs(1.d0-sum)
               psi(i)=psit
               if(err.lt.eps) then
                  nmax=i
                  return
               endif
            enddo
            nmax=n
         else
            psi(0)=dsin(ds)
            psi(1)=psi(0)/ds-dcos(ds)
            do i=1,n-1
               sn=dble(i+i+1)/ds
               psi(i+1)=sn*psi(i)-psi(i-1)
            enddo
            nmax=n
         endif
         end subroutine ricbessel
!
!  ricatti-hankel function xi(n), real argument
!
!
!  last revised: 15 January 2011
!
         subroutine richankel(n,ds,xi)
         implicit none
         integer :: n,i,ns
         real(8) :: ds,dns,sn,chi0,chi1,chi2,psi,psi0,psi1
         complex(8) :: xi(0:n)
         if(int(ds).lt.n) then
            ns=nint(ds+4.*(ds**.3333)+17)
            ns=max(n+10,ns)
            dns=0.d0
            do i=ns-1,n,-1
               sn=dble(i+1)/ds
               dns=sn-1.d0/(dns+sn)
            enddo
            xi(n)=dns
            xi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
            do i=n-2,1,-1
               sn=dble(i+1)/ds
               xi(i)=sn-1.d0/(xi(i+1)+sn)
            enddo
            chi0=-dcos(ds)
            psi=dsin(ds)
            chi1=chi0/ds-psi
            xi(0)=dcmplx(psi,chi0)
            do i=1,n
               chi2=dble(i+i+1)/ds*chi1-chi0
               psi=psi/(dble(i)/ds+xi(i))
               xi(i)=dcmplx(psi,chi1)
               chi0=chi1
               chi1=chi2
            enddo
            return
         else
            chi0=-dcos(ds)
            psi0=dsin(ds)
            chi1=chi0/ds-psi0
            psi1=psi0/ds+chi0
            xi(0)=dcmplx(psi0,chi0)
            xi(1)=dcmplx(psi1,chi1)
            do i=1,n-1
               sn=dble(i+i+1)/ds
               xi(i+1)=sn*xi(i)-xi(i-1)
            enddo
            return
         endif
         end subroutine richankel
!
!  ricatti-bessel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!
         subroutine cricbessel(n,ds,psi)
         implicit none
         integer :: n,i
         complex(8) :: ds,psi(0:n),chi(0:n)
         call cspherebessel(n,ds,psi,chi)
         do i=0,n
            psi(i)=psi(i)*ds
         enddo
         return
         end subroutine cricbessel
!
!  ricatti-hankel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!  March 2013
!  The condition abs(xi(i))/abs(psi(0)) << 1.d-6
!  implies an argument with large imag part, and use of xi = psi + i chi will have
!  round off problems.   Upwards recurrence is used in this case.
!
         subroutine crichankel(n,ds,xi)
         implicit none
         integer :: n,i
         complex(8) :: ds,psi(0:n),chi(0:n),xi(0:n),ci,&
                       psi0
         data ci/(0.d0,1.d0)/
         xi(0)=-ci*cdexp(ci*ds)
         psi0=cdsin(ds)
         if(cdabs(xi(0))/cdabs(psi0).lt.1.d-6) then
            xi(1)=-cdexp(ci*ds)*(ci+ds)/ds
            do i=2,n
               xi(i)=dble(i+i+1)/ds*xi(i-1)-xi(i-2)
            enddo
         else
            call cspherebessel(n,ds,psi,chi)
            do i=1,n
               xi(i)=(psi(i)+ci*chi(i))*ds
            enddo
         endif
         end subroutine crichankel
!
!     ==========================================================
!     Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!              for a complex argument
!     Input :  z --- Complex argument
!              n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!     Output:  CSJ(n) --- jn(z)
!              CSY(n) --- yn(z)
!              NM --- Highest order computed
!     Routines called:
!              MSTA1 and MSTA2 for computing the starting
!              point for backward recurrence
!     ==========================================================
!
!    obtained from, and copywrited by, Jian-Ming Jin
!    http://jin.ece.uiuc.edu/
!
!
!  last revised: 15 January 2011
!
         subroutine cspherebessel(n,z,csj,csy)
         implicit none
         integer :: n,nm,k,m
         real(8) :: a0
         complex(8) :: z,csj(0:n),csy(0:n),csa,csb,cs,cf0,cf1,cf
         a0=cdabs(z)
         nm=n
         if (a0.lt.1.0d-60) then
            csj=(0.d0,0.d0)
            csy=(-1.d300,0.d0)
            csy(0)=(1.d0,0.d0)
            return
         endif
         csj=(0.d0,0.d0)
         csj(0)=cdsin(z)/z
         csj(1)=(csj(0)-cdcos(z))/z
         if (n.ge.2) then
            csa=csj(0)
            csb=csj(1)
            m=msta1(a0,200)
            if (m.lt.n) then
               nm=m
            else
               m=msta2(a0,n,15)
            endif
            cf0=0.0d0
            cf1=1.0d0-100
            do k=m,0,-1
               cf=(2.0d0*k+3.0d0)*cf1/z-cf0
               if (k.le.nm) csj(k)=cf
               cf0=cf1
               cf1=cf
            enddo
            if (cdabs(csa).gt.cdabs(csb)) cs=csa/cf
            if (cdabs(csa).le.cdabs(csb)) cs=csb/cf0
            do k=0,min(nm,n)
               csj(k)=cs*csj(k)
            enddo
         endif
         csy=(1.d200,0.d0)
         csy(0)=-cdcos(z)/z
         csy(1)=(csy(0)-cdsin(z))/z
         do k=2,min(nm,n)
            if (cdabs(csj(k-1)).gt.cdabs(csj(k-2))) then
               csy(k)=(csj(k)*csy(k-1)-1.0d0/(z*z))/csj(k-1)
            else
               csy(k)=(csj(k)*csy(k-2)-(2.0d0*k-1.0d0)/z**3)/csj(k-2)
            endif
         enddo
         end subroutine cspherebessel
!
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that the magnitude of
!              Jn(x) at that point is about 10^(-MP)
!     Input :  x     --- Argument of Jn(x)
!              MP    --- Value of magnitude
!     Output:  MSTA1 --- Starting point
!     ===================================================
!
!
!  last revised: 15 January 2011
!
         integer function msta1(x,mp)
         implicit none
         integer :: mp,n0,n1,it,nn
         real(8) :: x, a0,f1,f,f0
         a0=dabs(x)
         n0=int(1.1*a0)+1
         f0=envj(n0,a0)-mp
         n1=n0+5
         f1=envj(n1,a0)-mp
         do it=1,20
            nn=n1-(n1-n0)/(1.0d0-f0/f1)
            f=envj(nn,a0)-mp
            if(abs(nn-n1).lt.1) exit
            n0=n1
            f0=f1
            n1=nn
            f1=f
         enddo
         msta1=nn
         end function msta1
!
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that all Jn(x) has MP
!              significant digits
!     Input :  x  --- Argument of Jn(x)
!              n  --- Order of Jn(x)
!              MP --- Significant digit
!     Output:  MSTA2 --- Starting point
!     ===================================================
!
!
!  last revised: 15 January 2011
!
         integer function msta2(x,n,mp)
         implicit none
         integer :: n,mp,n0,n1,it,nn
         real(8) :: x,a0,hmp,ejn,obj,f0,f1,f
         a0=dabs(x)
         hmp=0.5d0*dble(mp)
         ejn=envj(n,a0)
         if (ejn.le.hmp) then
            obj=mp
            n0=int(1.1*a0)
         else
            obj=hmp+ejn
            n0=n
         endif
         f0=envj(n0,a0)-obj
         n1=n0+5
         f1=envj(n1,a0)-obj
         do it=1,20
            nn=n1-(n1-n0)/(1.0d0-f0/f1)
            f=envj(nn,a0)-obj
            if (abs(nn-n1).lt.1) exit
            n0=n1
            f0=f1
            n1=nn
            f1=f
         enddo
         msta2=nn+10
         end function msta2

         real(8) function envj(n,x)
         implicit none
         integer :: n
         real(8) :: x
         n=max(1,abs(n))
         envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
         end function envj
!
!    vector coupling coefficients vc(w) = C(m,n|k,l|m+k,w), w = |n-l|,... n+l
!    uses downwards and upwards recurrence
!
!
!  last revised: 15 January 2011
!
         subroutine vcfunc(m,n,k,l,vcn)
         use numconstants
         implicit none
         integer :: m,n,k,l,wmax,wmin,w,mk
         real(8) :: vcn(0:n+l),t1,t2,t3,vcmax,vctest,rat
         vcn=0.d0
         wmax=n+l
         wmin=max(abs(n-l),abs(m+k))
         vcn(wmax)=bcof(n+m,l+k)*bcof(n-m,l-k)/bcof(n+n,l+l)
         if(wmin.eq.wmax) return
         vcn(wmax-1)=vcn(wmax)*(l*m-k*n)*fnr(2*(l+n)-1)/fnr(l)/fnr(n)&
        &  /fnr(n+l+m+k)/fnr(n+l-m-k)
         if(wmin.eq.wmax-1) return
         mk=m+k
         vcmax=abs(vcn(wmax))+abs(vcn(wmax-1))
!
!  a downwards recurrence is used initially
!
         do w=wmax,wmin+2,-1
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)&
        &     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))&
        &    /dble(2*w*(w-1))
            t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1)&
        &     *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3)&
        &     *fnr(2*w-1))
            vcn(w-2)=(t2*vcn(w-1)-vcn(w)/t1)/t3
            if(mod(wmax-w,2).eq.1) then
               vctest=abs(vcn(w-2))+abs(vcn(w-1))
               vcmax=max(vcmax,vctest)
               rat=vctest/vcmax
!
!  if/when the coefficients start to decrease in magnitude, an upwards recurrence takes over
!
               if(rat.lt.0.01d0) exit
            endif
         enddo
         if(w-2.gt.wmin) then
            wmax=w-3
            call vcfuncuprec(m,n,k,l,wmax,vcn)
         endif
         end subroutine vcfunc
!
!  upwards VC coefficient recurrence
!
!
!  last revised: 15 January 2011
!
         subroutine vcfuncuprec(m,n,k,l,wmax,vcn)
         use numconstants
         implicit none
         integer :: m,n,k,l,wmax,w,mk,nl,m1,n1,l1,k1,w1,w2
         real(8) :: vcn(0:n+l),t1,t2,t3,vc1
         mk=abs(m+k)
         nl=abs(n-l)
         if(nl.ge.mk) then
            w=nl
            if(n.ge.l) then
               m1=m
               n1=n
               l1=l
               k1=k
            else
               m1=k
               n1=l
               k1=m
               l1=n
            endif
            vc1=(-1)**(k1+l1)*bcof(l1+k1,w-m1-k1) &
               *bcof(l1-k1,w+m1+k1)/bcof(l1+l1,w+w+1)
         else
            w=mk
            if(m+k.ge.0) then
               vc1=(-1)**(n+m)*bcof(n-l+w,l-k)*bcof(l-n+w,n-m) &
                  /bcof(w+w+1,n+l-w)
            else
               vc1=(-1)**(l+k)*bcof(n-l+w,l+k)*bcof(l-n+w,n+m) &
                 /bcof(w+w+1,n+l-w)
            endif
         endif
         w1=w
         vcn(w)=vc1
         w=w1+1
         mk=m+k
         w2=min(wmax,n+l)
         if(w2.gt.w1) then
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk) &
              *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            if(w1.eq.0) then
               t2=.5*dble(m-k)
            else
               t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1)) &
                 /dble(2*w*(w-1))
            endif
            vcn(w)=t1*t2*vcn(w1)
         endif
         do w=w1+2,w2
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk) &
              *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1)) &
             /dble(2*w*(w-1))
            t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1) &
              *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3) &
              *fnr(2*w-1))
            vcn(w)=t1*(t2*vcn(w-1)-t3*vcn(w-2))
         enddo
         end subroutine vcfuncuprec
!
!  Normalized associated legendre functions
!
!
!  last revised: 15 January 2011
!
         subroutine normalizedlegendre(cbe,mmax,nmax,dc)
         use numconstants
         implicit none
         integer :: nmax,mmax,m,n,im
         real(8) :: dc(-mmax:mmax,0:nmax),cbe,sbe
         sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
         dc=0.d0
         do m=0,mmax
            dc(m,m)=(-1)**m*(0.5d0*sbe)**m*bcof(m,m)
            if(m.eq.nmax) exit
            dc(m,m+1)=fnr(m+m+1)*cbe*dc(m,m)
            do n=m+1,nmax-1
               dc(m,n+1)=(-fnr(n-m)*fnr(n+m)*dc(m,n-1)+dble(n+n+1)*cbe*dc(m,n)) &
                         /(fnr(n+1-m)*fnr(n+1+m))
            enddo
         enddo
         do m=1,mmax
            im=(-1)**m
            do n=m,nmax
               dc(-m,n)=im*dc(m,n)
            enddo
         enddo
         end subroutine normalizedlegendre
!
!  Generalized spherical functions
!
!  dc(m,n*(n+1)+k)=(-1)^(m + k)((n - k)!(n + k)!/(n - m)!/(n + m)!)^(1/2)
!  ((1 + x)/2)^((m + k)/2)((1 - x)/2)^((k - m)/2)JacobiP[n - k, k - m, k + m, x]
!
!  for |m| <= kmax, n=0,1,...nmax, |k| <= n
!
!
!  last revised: 15 January 2011
!
         subroutine rotcoef(cbe,kmax,nmax,dc)
         use numconstants
         implicit none
         integer :: kmax,nmax,k,m,in,n,knmax,nn1,kn,im,m1
         real(8) :: cbe,sbe,dc(-kmax:kmax,0:nmax*(nmax+2)),cbe2,sbe2,dk0(-nmax-1:nmax+1),&
                    dk01(-nmax-1:nmax+1),sben,dkt,fmn,dkm0,dkm1,dkn1
         sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         in=1
         dk0(0)=1.d0
         sben=1.d0
         dc(0,0)=1.d0
         dk01(0)=0.d0
         do n=1,nmax
            knmax=min(n,kmax)
            nn1=n*(n+1)
            in=-in
            sben=sben*sbe/2.d0
            dk0(n)=in*sben*bcof(n,n)
            dk0(-n)=in*dk0(n)
            dk01(n)=0.d0
            dk01(-n)=0.d0
            dc(0,nn1+n)=dk0(n)
            dc(0,nn1-n)=dk0(-n)
            do k=-n+1,n-1
               kn=nn1+k
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)&
                     /(fnr(n+k)*fnr(n-k))
               dc(0,kn)=dk0(k)
            enddo
            im=1
            do m=1,knmax
               im=-im
               fmn=1.d0/fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.d0
               do k=-n,n
                  kn=nn1+k
                  dkm1=dkm0
                  dkm0=dc(m1,kn)
                  if(k.eq.n) then
                     dkn1=0.d0
                  else
                     dkn1=dc(m1,kn+1)
                  endif
                  dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                          -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1  &
                          -dble(k)*sbe*dc(m1,kn))*fmn
                  dc(-m,nn1-k)=dc(m,kn)*(-1)**(k)*im
               enddo
            enddo
         enddo
         end subroutine rotcoef

!
!  tau are the vector spherical harmonic functions, normalized
!
!
!  last revised: 15 January 2011
!
         subroutine taufunc(cb,nmax,tau)
         use numconstants
         implicit none
         integer :: nmax,n,m,nn1,mn
         real(8) :: drot(-1:1,0:nmax*(nmax+2)),tau(0:nmax+1,nmax,2),cb,fnm
         call rotcoef(cb,1,nmax,drot)
         do n=1,nmax
            nn1=n*(n+1)
            fnm=sqrt(dble(n+n+1)/2.d0)/4.d0
            do m=-n,-1
               mn=nn1+m
               tau(n+1,-m,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(n+1,-m,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
            do m=0,n
               mn=nn1+m
               tau(m,n,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(m,n,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
         enddo
         end subroutine taufunc
!
! vector spherical harmonic function
! november 2011
! april 2012: lr formulation
!

         subroutine pifunc(cb,ephi,nmax,ndim,pivec)
         use numconstants
         implicit none
         integer :: nmax,n,m,nn1,ndim
         real(8) :: drot(-1:1,0:nmax*(nmax+2)),cb
         complex(8) :: pivec(0:ndim+1,ndim,2),ephi,ephim(-nmax:nmax),cin
         call rotcoef(cb,1,nmax,drot)
         ephim(0)=1.d0
         do m=1,nmax
            ephim(m)=ephi*ephim(m-1)
            ephim(-m)=dconjg(ephim(m))
         enddo
         do n=1,nmax
            cin=(0.d0,-1.d0)**(n)*fnr(n+n+1)
            nn1=n*(n+1)
            pivec(n+1,n:1:-1,1)=cin*drot(1,nn1-n:nn1-1)*ephim(-n:-1)
            pivec(n+1,n:1:-1,2)=cin*drot(-1,nn1-n:nn1-1)*ephim(-n:-1)
            pivec(0:n,n,1)=cin*drot(1,nn1:nn1+n)*ephim(0:n)
            pivec(0:n,n,2)=cin*drot(-1,nn1:nn1+n)*ephim(0:n)
         enddo
         end subroutine pifunc
!
!  regular vswf expansion coefficients for a plane wave.
!  alpha, beta: incident azimuth and polar angles.
!
!
!  last revised: 15 January 2011
!  april 2012: lr formulation
!
         subroutine planewavecoef(alpha,beta,nodr,pmnp0)
         use numconstants
         implicit none
         integer :: nodr,m,n,p,sp
         real(8) :: alpha,beta,cb,sb,ca,sa
         real(8), allocatable :: tau(:,:,:),taulr(:,:,:)
         complex(8) :: ealpha,ci,cin
         complex(8), allocatable :: ealpham(:)
         complex(8) :: pmnp0(0:nodr+1,nodr,2,2)
         data ci/(0.d0,1.d0)/
         call init(nodr)
         allocate(ealpham(-nodr:nodr))
         allocate(tau(0:nodr+1,nodr,2),taulr(0:nodr+1,nodr,2))
         cb=cos(beta)
         sb=sqrt((1.d0-cb)*(1.d0+cb))
         ca=cos(alpha)
         sa=sin(alpha)
         ealpha=dcmplx(ca,sa)
         call taufunc(cb,nodr,tau)
         taulr(:,:,1)=(tau(:,:,1)+tau(:,:,2))*.5d0
         taulr(:,:,2)=(tau(:,:,1)-tau(:,:,2))*.5d0
         call ephicoef(ealpha,nodr,ealpham)
         do n=1,nodr
            cin=4.d0*ci**(n+1)
            do p=1,2
               sp=-(-1)**p
               do m=-n,-1
                  pmnp0(n+1,-m,p,1)=-cin*taulr(n+1,-m,p)*ealpham(-m)
                  pmnp0(n+1,-m,p,2)=sp*ci*cin*taulr(n+1,-m,p)*ealpham(-m)
               enddo
               do m=0,n
                  pmnp0(m,n,p,1)=-cin*taulr(m,n,p)*ealpham(-m)
                  pmnp0(m,n,p,2)=sp*ci*cin*taulr(m,n,p)*ealpham(-m)
               enddo
            enddo
         enddo
         deallocate(ealpham,tau)
         end subroutine planewavecoef
!
!  regular vswf expansion coefficients for a gaussian beam, localized approximation.
!  cbeam = 1/(k omega)
!
!
!  last revised: 15 January 2011
!
         subroutine gaussianbeamcoef(alpha,beta,cbeam,nodr,pmnp0)
         use numconstants
         implicit none
         integer :: nodr,m,n,p,k
         real(8) :: alpha,beta,cbeam,gbn
         complex(8) :: pmnp0(0:nodr+1,nodr,2,2)
         call planewavecoef(alpha,beta,nodr,pmnp0)
         do n=1,nodr
            gbn=dexp(-((dble(n)+.5d0)*cbeam)**2.)
            do p=1,2
               do k=1,2
                  do m=-n,-1
                     pmnp0(n+1,-m,p,k)=pmnp0(n+1,-m,p,k)*gbn
                  enddo
                  do m=0,n
                     pmnp0(m,n,p,k)=pmnp0(m,n,p,k)*gbn
                  enddo
               enddo
            enddo
         enddo
         end subroutine gaussianbeamcoef
!
!  plane wave expansion coefficients at sphere origins.  uses a phase shift.
!
!
!  last revised: 15 January 2011
!
         subroutine sphereplanewavecoef(nsphere,neqns,nodr,nodrmax,alpha,beta, &
                    rpos,hostsphere,numberfieldexp,rimedium,pmnp)
         implicit none
         integer :: m,n,p,nsphere,i,l,nodr(nsphere),nodrmax,neqns,k, &
                    hostsphere(nsphere),numberfieldexp(nsphere),j
         real(8) :: alpha,beta,cb,sb,ca,sa,rpos(3,nsphere)
         complex(8) :: ci,phasefac,pmnp(neqns*2),rimedium(2),rib
         complex(8) :: pmnp0(0:nodrmax+1,nodrmax,2,2)
         data ci/(0.d0,1.d0)/
         rib=2./(1.d0/rimedium(1)+1.d0/rimedium(2))
         call planewavecoef(alpha,beta,nodrmax,pmnp0)
         cb=cos(beta)
         sb=sqrt((1.d0-cb)*(1.d0+cb))
         ca=cos(alpha)
         sa=sin(alpha)
         l=0
         do i=1,nsphere
            phasefac=cdexp(ci*rib*((ca*rpos(1,i)+sa*rpos(2,i))*sb+rpos(3,i)*cb))
            do j=1,numberfieldexp(i)
               do k=1,2
                  do p=1,2
                     do n=1,nodr(i)
                        do m=0,nodr(i)+1
                           l=l+1
                           if(hostsphere(i).eq.0.and.j.eq.1) then
                              pmnp(l)=phasefac*pmnp0(m,n,p,k)
                           else
                              pmnp(l)=0.d0
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         end subroutine sphereplanewavecoef
!
!  axial translation coefficients calculated by the diamond recurrence formula
!  new: 10 october 2011
!  april 2012: lr formulation
!  may 2012: new ordering scheme:
!  input:  itype : 1 or 3 (regular, outgoing)
!     r: axial translation distance (positive)
!     ri: rank 2 complex array: L and R refractive indices of medium
!     nmax, lmax: largest row and column orders.
!     ndim: dimension of ac
!  output:
!     ac:  rank 1 complex array, dimension ndim, containing the matrix elements.
!  storage scheme:   for each degree m, with ordering m=0, -1, 1, -2, 2, ..min(nmax,lmax),
!  the elements for degree m are stored
!
         subroutine axialtrancoefrecurrence(itype,r,ri,nmax,lmax,ndim,ac)
         use numconstants
         implicit none
         integer :: itype,nmax,lmax,n,l,m,p,nlmin, &
                    wmin,wmax,ml,m1,np1,nm1,lm1,lp1,sp
         integer :: iadd,nlmax,iadd0,iadd1,ndim
         integer :: ma,blockdim
         integer, save :: nlmax0
         real(8) :: r
         complex(8) :: ri(2),ci,z(2),xi(0:nmax+lmax,2)
         complex(8) :: ac(ndim),act(nmax,lmax,2),actt(2,2)
         data ci,nlmax0/(0.d0,1.d0),0/
         nlmax=max(nmax,lmax)
         nlmin=min(nmax,lmax)
         if(nlmax.gt.nlmax0) then
            nlmax0=nlmax
            call axialtrancoefinit(nlmax)
         endif

         if(r.eq.0.d0) then
            ac=(0.d0,0.d0)
            if(itype.ne.1) return
            iadd0=0
            do ma=0,nlmin
               m1=max(1,ma)
               do m=-ma,ma,2*m1
                  blockdim=(nmax-m1+1)*(lmax-m1+1)*2
                  iadd1=iadd0+blockdim
                  do l=m1,nlmax
                     act(l,l,1)=1.d0
                     act(l,l,2)=1.d0
                  enddo
                  ac(iadd0+1:iadd1)=reshape(act(m1:nmax,m1:lmax,1:2),(/blockdim/))
                  iadd0=iadd1
               enddo
            enddo
            return
         endif
         z=r*ri
         do p=1,2
            if(itype.eq.1) then
               call cricbessel(nmax+lmax,z(p),xi(0:,p))
            else
               call crichankel(nmax+lmax,z(p),xi(0:,p))
            endif
            xi(0:,p)=xi(0:,p)/z(p)
         enddo
         lm1=lmax-1

         iadd0=0
         do ma=0,nlmin
            m1=max(1,ma)
            lp1=m1+1
            do m=-ma,ma,2*m1
               blockdim=2*(nmax-m1+1)*(lmax-m1+1)
               iadd1=iadd0+blockdim
               n=m1
               do l=m1,lmax
                  wmin=abs(n-l)
                  wmax=n+l
                  iadd=iadd+1
                  ml=l*(l+1)+m
                  do p=1,2
                     actt(1,p)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2,p))
                     actt(2,p)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2,p))
                  enddo
                  act(n,l,1)=actt(1,1)+actt(2,1)
                  act(n,l,2)=actt(1,2)-actt(2,2)
               enddo
               l=lmax
               ml=l*(l+1)+m
               do n=m1+1,nmax
                  wmin=abs(n-l)
                  wmax=n+l
                  do p=1,2
                     actt(1,p)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2,p))
                     actt(2,p)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2,p))
                  enddo
                  act(n,l,1)=actt(1,1)+actt(2,1)
                  act(n,l,2)=actt(1,2)-actt(2,2)
               enddo

               if(m1.lt.nlmin) then
                  do n=m1,nmax-1
                     np1=n+1
                     nm1=n-1
                     do p=1,2
                        sp=-(-1)**p
                        act(np1,m1:lmax-1,p)= &
                          - act(n,m1+1:lmax,p)*fnp1_const(m,m1:lm1) &
                          + sp*(fn_const(m,m1:lm1)-fn_const(m,n))*ci*act(n,m1:lm1,p)
                        act(np1,m1+1:lm1,p)=act(np1,m1+1:lm1,p) &
                          + act(n,m1:lmax-2,p)*fnm1_const(m,lp1:lm1)
                        if(n.gt.m1) then
                           act(np1,m1:lm1,p)=act(np1,m1:lm1,p) &
                             + act(nm1,m1:lm1,p)*fnm1_const(m,n)
                        endif
                        act(np1,m1:lm1,p)=act(np1,m1:lm1,p)/fnp1_const(m,n)
                     enddo
                  enddo
               endif
               ac(iadd0+1:iadd1)=reshape(act(m1:nmax,m1:lmax,1:2),(/blockdim/))
               iadd0=iadd1
            enddo
         enddo
         end subroutine axialtrancoefrecurrence

!
!  constants for translation coefficient calculation
!
         subroutine axialtrancoefinit(nmax)
         use numconstants
         implicit none
         integer :: nmax,m,n,l,w,n21,ml,ll1,wmin,wmax,nlmin,lp1,lm1
         real(8) :: c1,c2,vc1(0:2*nmax),vc2(0:2*nmax)
         complex(8) :: ci,inlw
         data ci/(0.d0,1.d0)/
         if(allocated(vcc_const)) deallocate(vcc_const,fnm1_const,fn_const,fnp1_const)
         allocate(vcc_const(nmax,nmax*(nmax+2),0:2*nmax),fnm1_const(-nmax:nmax,nmax), &
                  fn_const(-nmax:nmax,nmax),fnp1_const(-nmax:nmax,nmax))
         do n=1,nmax
            n21=n+n+1
            do l=1,nmax
               c1=fnr(n21)*fnr(l+l+1)
               ll1=l*(l+1)
               call vcfunc(-1,n,1,l,vc2)
               wmin=abs(n-l)
               wmax=n+l
               nlmin=min(l,n)
               do m=-nlmin,nlmin
                  ml=ll1+m
                  c2=-c1*(-1)**m
                  call vcfunc(-m,n,m,l,vc1)
                  do w=wmin,wmax
                     inlw=ci**(n-l+w)
                     vcc_const(n,ml,w)=c2*vc1(w)*vc2(w)*(dble(inlw)+dimag(inlw))
                  enddo
               enddo
            enddo
         enddo
         fnm1_const=0.
         fn_const=0.
         fnp1_const=0.
         do m=-nmax,nmax
            do l=max(1,abs(m)),nmax
               lp1=l+1
               lm1=l-1
               fnm1_const(m,l)=fnr(lm1)*fnr(lp1)*fnr(l-m)*fnr(l+m)/fnr(lm1+l)/fnr(l+lp1)/dble(l)
               fn_const(m,l)=dble(m)/dble(l)/dble(lp1)
               fnp1_const(m,l)=fnr(l)*fnr(l+2)*fnr(lp1-m)*fnr(lp1+m)/fnr(l+lp1)/fnr(l+l+3)/dble(lp1)
            enddo
         enddo
         end subroutine axialtrancoefinit
!
!  test to determine convergence of regular vswf addition theorem for max. order lmax
!  and translation distance r w/ refractive index ri.
!
!
!  last revised: 15 January 2011
!
         subroutine tranordertest(r,ri,lmax,eps,nmax)
         use numconstants
         implicit none
         integer :: nmax,lmax,n,l,m,w,n21,wmin,wmax
         integer, parameter :: nlim=200
         real(8) :: r,alnw,sum,eps
         real(8) :: vc1(0:nlim+lmax)
         complex(8) :: ri,ci,z,a,b,c
         complex(8) :: xi(0:nlim+lmax)
         data ci/(0.d0,1.d0)/
         if(r.eq.0.d0) then
            nmax=lmax
            return
         endif
         z=r*ri
         sum=0.d0
         do n=1,nlim
            call init(n+lmax)
            call cricbessel(n+lmax,z,xi)
            do l=0,n+lmax
               xi(l)=xi(l)/z*ci**l
            enddo
            n21=n+n+1
            l=lmax
            c=fnr(n21)*fnr(l+l+1)*ci**(n-l)
            call vcfunc(-1,n,1,l,vc1)
            wmin=abs(n-l)
            wmax=n+l
            m=1
            a=0.
            b=0.
            do w=wmin,wmax
               alnw=vc1(w)*vc1(w)
               if(mod(n+l+w,2).eq.0) then
                  a=a+alnw*xi(w)
               else
                  b=b+alnw*xi(w)
               endif
            enddo
            a=c*a
            b=c*b
            sum=sum+a*conjg(a)+b*conjg(b)
            if(abs(1.d0-sum).lt.eps) exit
         enddo
         nmax=min(n,nlim)
         nmax=max(nmax,lmax)
         end subroutine tranordertest
!
!  address for axial translation coefficient
!
!
!  last revised: 15 January 2011
!

         integer function atcadd(m,n,ntot)
         implicit none
         integer :: m,n,ntot
         atcadd=n-ntot+(max(1,m)*(1+2*ntot-max(1,m)))/2+ntot*min(1,m)
         end function atcadd

         integer function atcdim(ntot,ltot)
         implicit none
         integer :: ntot,ltot,nmin,nmax
         nmin=min(ntot,ltot)
         nmax=max(ntot,ltot)
         atcdim=2*(nmin*(1- nmin*nmin + 3*nmax*(2 + nmin)))/3
         end function atcdim
!
! the offset (integer) for the ntot X ltot translation matrix for degree m
!
         integer function moffset(m,ntot,ltot)
         implicit none
         integer :: m,ntot,ltot
         if(m.eq.0) then
            moffset=0
         elseif(m.lt.0) then
            moffset=2*(-((1+m)*(2+m)*(3+2*m +3*ntot)) &
                   -3*ltot*(2+ntot+m*(3+m+2*ntot)))/3
         else
            moffset=2*(-3*ltot*(-1 + m)**2 +6*ltot*m*ntot &
                   +(-1+m)*(m*(-4+2*m-3*ntot)+3*(1 + ntot)))/3
         endif
         end function moffset
!
!   gentrancoef: calculates the vwh translation coefficients for
!   a general translation from one origin to another
!
!   input: itype: integer, =1, regular, =3, outgoing type harmonics
!          xptran: real, dim 3 vector: x,y,z components of translation, in units
!                   of 1/k
!          ri: complex, refractive index of medium
!          nrow0,nrow1,ncol0,ncol1: integer, starting and stopping row and column order
!          iaddrow0,iaddcol0: address offset for row and column order (see below)
!   output: ac(p,mn,kl): complex translation matrix.  calculated for mode p=1,2 (A or B type),
!           order n=nrow0,nrow1, degree m=-n,n
!           order l=ncol0,ncol1, degree k=-n,n
!           address is given by
!           mn=m+n*(n+1)-(nrow0-1)*(nrow0+1)+iaddrow0
!           kl=k+l*(l+1)-(ncol0-1)*(ncol0+1)+iaddcol0
!           that is, if iaddrow0=0 the address is mn=1 for n=nrow0 and m=-n.
!
!
!  last revised: 15 January 2011
!  april 2012: lr formulation
!
         subroutine gentrancoef(itype,xptran,ri,nrow0,nrow1,ncol0,ncol1, &
                               iaddrow0,iaddcol0,ac)
         use numconstants
         implicit none
         integer :: itype,nrow0,nrow1,ncol0,ncol1,iaddrow0,iaddcol0,kmax
         integer :: ntot,nblkr0,nblkr1,nblkc0,nblkc1
         integer :: v,vw,w,wmax,wmin,n,l,m,k,p,nn1,ll1,mn,kl,m1m
         real(8) :: vc1(0:nrow1+ncol1),vc2(0:nrow1+ncol1),&
                    xptran(3),r,ct,ct0
         real(8) :: drot(0:0,0:(nrow1+ncol1)*(nrow1+ncol1+2))
         complex(8) :: ri(2),ci,ephi,ac(2,nrow1*(nrow1+2)-(nrow0-1)*(nrow0+1)-iaddrow0,&
                       ncol1*(ncol1+2)-(ncol0-1)*(ncol0+1)-iaddcol0),&
                       z(2),c,a,b,atc(2,2)
         complex(8) :: ephim(-(nrow1+ncol1):nrow1+ncol1),jnc(0:nrow1+ncol1,2)
         data ci/(0.d0,1.d0)/
         call cartosphere(xptran,r,ct,ephi)
         ntot=nrow1+ncol1
         nblkr0=(nrow0-1)*(nrow0+1)
         nblkr1=nrow1*(nrow1+2)
         nblkc0=(ncol0-1)*(ncol0+1)
         nblkc1=ncol1*(ncol1+2)
         if(r.eq.0.d0) then
            do n=nblkr0+1,nblkr1
               mn=n-nblkr0+iaddrow0
               do l=nblkc0+1,nblkc1
                  kl=l-nblkc0+iaddcol0
                  do p=1,2
                     ac(p,mn,kl)=0.d0
                  enddo
               enddo
               if(n.gt.nblkc0.and.n.le.nblkc1.and.itype.eq.1) then
                  ac(1,mn,n-nblkc0+iaddcol0)=1.d0
                  ac(2,mn,n-nblkc0+iaddcol0)=1.d0
               endif
            enddo
            return
         endif
         kmax=0
         ct0=ct
         call rotcoef(ct0,kmax,ntot,drot)
         call ephicoef(ephi,ntot,ephim)
         z=ri*r
         do p=1,2
            if(itype.eq.1) then
               call cricbessel(ntot,z(p),jnc(0,p))
            else
               call crichankel(ntot,z(p),jnc(0,p))
            endif
            do n=0,ntot
               c=ci**n
               jnc(n,p)=c*jnc(n,p)/z(p)
            enddo
         enddo
         do l=ncol0,ncol1
            ll1=l*(l+1)
            do n=nrow0,nrow1
               nn1=n*(n+1)
               wmax=n+l
               call vcfunc(-1,n,1,l,vc2)
               c=-ci**(n-l)*fnr(n+n+1)*fnr(l+l+1)
               do k=-l,l
                  kl=ll1+k-nblkc0+iaddcol0
                  do m=-n,n
                     m1m=(-1)**m
                     mn=nn1+m-nblkr0+iaddrow0
                     v=k-m
                     call vcfunc(-m,n,k,l,vc1)
                     a=0.
                     b=0.
                     wmin=max(abs(v),abs(n-l))
                     do p=1,2
                        do w=wmax,wmin,-1
                           vw=w*(w+1)+v
                           if(mod(wmax-w,2).eq.0) then
                              a=a+vc1(w)*vc2(w)*jnc(w,p)*drot(0,vw)
                           else
                              b=b+vc1(w)*vc2(w)*jnc(w,p)*drot(0,vw)
                           endif
                        enddo
                        atc(p,1)=a
                        atc(p,2)=b
                     enddo
                     ac(1,mn,kl)=(atc(1,1)+atc(1,2))*c*m1m*ephim(v)
                     ac(2,mn,kl)=(atc(2,1)-atc(2,2))*c*m1m*ephim(v)
                  enddo
               enddo
            enddo
         enddo
         return
         end subroutine gentrancoef
!
! cartosphere takes the cartesian point (x,y,z) = xp(1), xp(2), xp(3)
! and converts to polar form: r: radius, ct: cos(theta), ep = exp(i phi)
!
!
!  last revised: 15 January 2011
!
         subroutine cartosphere(xp,r,ct,ep)
         implicit none
         real(8) :: xp(3),r,ct
         complex(8) :: ep
         r=xp(1)*xp(1)+xp(2)*xp(2)+xp(3)*xp(3)
         if(r.eq.0.d0) then
            ct=1.d0
            ep=(1.d0,0.d0)
            return
         endif
         r=sqrt(r)
         ct=xp(3)/r
         if(xp(1).eq.0.d0.and.xp(2).eq.0.d0) then
            ep=(1.d0,0.d0)
         else
            ep=dcmplx(xp(1),xp(2))/sqrt(xp(1)*xp(1)+xp(2)*xp(2))
         endif
         return
         end subroutine cartosphere
!
! euler rotation of a point (x,y,z) = xp(1), xp(2), xp(3)
! November 2012
!
         subroutine eulerrotation(xp,eulerangf,dir,xprot)
         implicit none
         integer :: dir
         real(8) :: xp(3),eulerangf(3),eulerang(3),cang(3),sang(3), &
                    mat1(3,3),mat2(3,3),mat3(3,3),xprot(3),xpt(3)
         xpt=xp
         if(dir.eq.1) then
            eulerang=eulerangf
         else
            eulerang(1:3)=-eulerangf(3:1:-1)
         endif
         cang=cos(eulerang)
         sang=sin(eulerang)
         mat1(1,:) = (/cang(1),sang(1),0.d0/)
         mat1(2,:) = (/-sang(1),cang(1),0.d0/)
         mat1(3,:) = (/0.d0,0.d0,1.d0/)
         mat2(1,:) = (/cang(2),0.d0,-sang(2)/)
         mat2(2,:) = (/0.d0,1.d0,0.d0/)
         mat2(3,:) = (/sang(2),0.d0,cang(2)/)
         mat3(1,:) = (/cang(3),sang(3),0.d0/)
         mat3(2,:) = (/-sang(3),cang(3),0.d0/)
         mat3(3,:) = (/0.d0,0.d0,1.d0/)
         xpt=matmul(mat1,xpt)
         xpt=matmul(mat2,xpt)
         xpt=matmul(mat3,xpt)
         xprot=xpt
         end subroutine eulerrotation
!
! ephicoef returns the complex array epm(m) = exp(i m phi) for
! m=-nodr,nodr.   ep =exp(i phi), and epm is dimensioned epm(-nd:nd)
!
!
!  last revised: 15 January 2011
!
         subroutine ephicoef(ep,nodr,epm)
         implicit none
         integer :: nodr,m
         complex(8) :: ep,epm(-nodr:nodr)
         epm(0)=(1.d0,0.d0)
         do m=1,nodr
            epm(m)=ep*epm(m-1)
            epm(-m)=dconjg(epm(m))
         enddo
         return
         end subroutine ephicoef
!
!  test to determine max order of vswf expansion of a plane wave at distance r
!
!
!  last revised: 15 January 2011
!
         subroutine planewavetruncationorder(r,rimedium,eps,nodr)
         implicit none
         integer :: nodr,n1,n
         real(8) :: r,eps,err
         complex(8), allocatable :: jn(:)
         complex(8) :: sum, ci,eir,rimedium(2),rri,rib
         data ci/(0.d0,1.d0)/
         rib=2.d0/(1.d0/rimedium(1)+1.d0/rimedium(2))
         n1=max(10,int(3.*r+1))
         allocate(jn(0:n1))
         rri=r*rib
         call cricbessel(n1,rri,jn)
         jn(0:n1)=jn(0:n1)/rri
         eir=cdexp(-ci*rri)
         sum=jn(0)*eir
         do n=1,n1
            sum=sum+ci**n*dble(n+n+1)*jn(n)*eir
            err=cdabs(1.d0-sum)
            if(err.lt.eps) then
               nodr=n
               deallocate(jn)
               return
            endif
         enddo
         nodr=n1
         deallocate(jn)
         end subroutine planewavetruncationorder
!
!  calculates the cartesian components of the vswf at position rpos, in ref. index ri.
!
!
!  original: 15 January 2011
!  revised: 23 February 2011: multiplied by root 2
!  april 2012: lr formulation
!
         subroutine vwhcalc(rpos,ri,nodr,itype,vwh)
         use numconstants
         implicit none
         integer :: nodr,itype,n,nodrp1,nodrm1,nn1,np1,nm1,p,sp
         integer, save :: nodrmax
         real(8) ::  rpos(3),r,ct
         real(8) pmn(0:0,0:(nodr+1)*(nodr+3))
         complex(8) :: ci,vwh(3,2,1:*),ri(2),ephi,a(2)
         complex(8)  :: a1vec(-nodr:nodr), &
                       b1vec(-nodr:nodr),z1vec(-nodr:nodr),a2vec(-nodr:nodr), &
                       b2vec(-nodr:nodr),z2vec(-nodr:nodr)
         complex(8) :: umn(-nodr-2:nodr+2,0:nodr+1,2), hn(0:nodr+1,2), ephim(-nodr-1:nodr+1)
         data ci,nodrmax/(0.d0,1.d0),0/
         if(nodr.gt.nodrmax) then
            nodrmax=nodr
            call init(nodr+2)
         endif
         call cartosphere(rpos,r,ct,ephi)
         if(r.le.1.d-4) then
            vwh(:,:,1:nodr*(nodr+2))=(0.d0,0.d0)
            if(itype.eq.3) return
            do p=1,2
               vwh(1,p,1)=.5d0*fnr(2)/fnr(3)
               vwh(2,p,1)=-.5d0*ci*fnr(2)/fnr(3)
               vwh(3,p,2)=1.d0*fnr(2)/fnr(6)
               vwh(1,p,3)=-.5d0*fnr(2)/fnr(3)
               vwh(2,p,3)=-.5d0*ci*fnr(2)/fnr(3)
            enddo
            return
         endif
         nodrp1=nodr+1
         nodrm1=nodr-1
!
! this is now a vector operation w/ l/r form
!
         a=ri*r
         do p=1,2
            if(itype.eq.1) then
               call cricbessel(nodrp1,a(p),hn(0,p))
            else
               call crichankel(nodrp1,a(p),hn(0,p))
            endif
            hn(0:nodrp1,p)=hn(0:nodrp1,p)/a(p)
         enddo
         call rotcoef(ct,0,nodrp1,pmn)
         call ephicoef(ephi,nodrp1,ephim)
         umn=0.d0
         do p=1,2
            umn(0,0,p)=hn(0,p)*fnr(2)
            do n=1,nodrp1
               nn1=n*(n+1)
               umn(-n:n,n,p)=fnr(2)*pmn(0,nn1-n:nn1+n)*ephim(-n:n)*hn(n,p)
               umn(-n-1,n,p)=0.d0
               umn(n+1,n,p)=0.d0
            enddo
         enddo
         do p=1,2
            sp=-(-1)**p
            do n=1,nodr
               nn1=n*(n+1)
               np1=n+1
               nm1=n-1
               a1vec(-n:n)=vwh_coef(-n:n,n,1,1)*umn(-nm1:np1,np1,p) &
                  +vwh_coef(-n:n,n,1,-1)*umn(-nm1:np1,nm1,p)
               b1vec(-n:n)=vwh_coef(-n:n,n,-1,1)*umn(-np1:nm1,np1,p) &
                  +vwh_coef(-n:n,n,-1,-1)*umn(-np1:nm1,nm1,p)
               z1vec(-n:n)=vwh_coef(-n:n,n,0,1)*umn(-n:n,np1,p) &
                  +vwh_coef(-n:n,n,0,-1)*umn(-n:n,nm1,p)
               a2vec(-n:n)=vwh_coef(-n:n,n,1,0)*umn(-nm1:np1,n,p)
               b2vec(-n:n)=vwh_coef(-n:n,n,-1,0)*umn(-np1:nm1,n,p)
               z2vec(-n:n)=vwh_coef(-n:n,n,0,0)*umn(-n:n,n,p)
               vwh(1,p,nn1-n:nn1+n)=-0.5d0*(a1vec(-n:n)+b1vec(-n:n)) &
                      -sp*0.5d0*ci*(a2vec(-n:n)+b2vec(-n:n))
               vwh(2,p,nn1-n:nn1+n)=-0.5d0*ci*(-a1vec(-n:n)+b1vec(-n:n)) &
                      -sp*0.5d0*(a2vec(-n:n)-b2vec(-n:n))
               vwh(3,p,nn1-n:nn1+n)=-z1vec(-n:n) &
                      -sp*ci*z2vec(-n:n)
            enddo
         enddo
         end subroutine vwhcalc
!
!  svwf calculation for an axial translation
!
!
!  original: 15 January 2011
!  revised: 23 February 2011: multiplied by root 2
!  april 2012: l/r formulation
!
         subroutine vwhaxialcalc(rpos,ri,nodr,itype,vwh)
         use numconstants
         implicit none
         integer :: nodr,itype,m,n,p,nodrp1,nodrm1,np1,nm1,mp1,mm1,sp,k
         integer, save :: nodrmax
         real(8) ::  rpos(3),r,ct
         real(8) pmn(-2:2,0:nodr+1)
         complex(8) :: ci,vwh(3,2,2,1:nodr),ri(2),ephi,a(2),a1,b1,z1,&
                       a2,b2,z2,umn(-2:2,0:nodr+1,2), hn(0:nodr+1,2), &
                       ephim(-2:2)
         data ci,nodrmax/(0.d0,1.d0),0/
         if(nodr.gt.nodrmax) then
            nodrmax=nodr
            call init(nodr+2)
         endif
         call cartosphere(rpos,r,ct,ephi)
         if(r.le.1.d-4) then
            vwh(:,:,:,1:nodr)=(0.d0,0.d0)
            if(itype.eq.3) return
            vwh(1,1,1,1)=.5d0*fnr(2)/fnr(3)
            vwh(2,1,1,1)=-.5d0*ci*fnr(2)/fnr(3)
            vwh(1,1,2,1)=-.5d0*fnr(2)/fnr(3)
            vwh(2,1,2,1)=-.5d0*ci*fnr(2)/fnr(3)
            return
         endif
         nodrp1=nodr+1
         nodrm1=nodr-1
         a=ri*r
         do p=1,2
            if(itype.eq.1) then
               call cricbessel(nodrp1,a(p),hn(0,p))
            else
               call crichankel(nodrp1,a(p),hn(0,p))
            endif
            hn(0:nodrp1,p)=hn(0:nodrp1,p)/a(p)
         enddo
         call normalizedlegendre(ct,2,nodrp1,pmn)
         call ephicoef(ephi,2,ephim)
         umn(-2:2,0:nodrp1,1:2)=0.d0
         umn(0,0,1:2)=hn(0,1:2)*fnr(2)
         do n=1,nodrp1
            p=min(n,2)
            do m=-p,p
               umn(m,n,1:2)=fnr(2)*pmn(m,n)*ephim(m)*hn(n,1:2)
            enddo
         enddo
         vwh(:,:,:,1:nodr)=0.d0
         do p=1,2
            sp=-(-1)**p
            do n=1,nodr
               np1=n+1
               nm1=n-1
               do m=-1,1,2
                  k=(m+1)/2+1
                  mp1=m+1
                  mm1=m-1
                  a1=vwh_coef(m,n,1,1)*umn(mp1,np1,p) &
                     +vwh_coef(m,n,1,-1)*umn(mp1,nm1,p)
                  b1=vwh_coef(m,n,-1,1)*umn(mm1,np1,p) &
                     +vwh_coef(m,n,-1,-1)*umn(mm1,nm1,p)
                  z1=vwh_coef(m,n,0,1)*umn(m,np1,p) &
                     +vwh_coef(m,n,0,-1)*umn(m,nm1,p)
                  a2=vwh_coef(m,n,1,0)*umn(mp1,n,p)
                  b2=vwh_coef(m,n,-1,0)*umn(mm1,n,p)
                  z2=vwh_coef(m,n,0,0)*umn(m,n,p)
                  vwh(1,p,k,n)=-0.5d0*(a1+b1)-sp*0.5d0*ci*(a2+b2)
                  vwh(2,p,k,n)=-0.5d0*ci*(-a1+b1)-sp*0.5d0*(a2-b2)
                  vwh(3,p,k,n)=-z1-sp*ci*z2
               enddo
            enddo
         enddo
         return
         end subroutine vwhaxialcalc
!
! inverse of a 2 X 2 complex matrix.
! March 2013
!
         subroutine twobytwoinverse(mat,imat)
         implicit none
         integer :: s,t,ss,st
         complex(8) :: mat(2,2),imat(2,2),tmat(2,2),det
         tmat=mat
         det=mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)
         do s=1,2
            ss=(-1)**s
            do t=1,2
               st=(-1)**t
               imat(s,t)=ss*st*tmat(3-t,3-s)/det
            enddo
         enddo
         end subroutine twobytwoinverse
!
! move between unequal size matrices.
! March 2013
!
         subroutine transfer(nin,nout,cin,cout)
         implicit none
         integer :: nin,nout,nmin
         complex(8) :: cin(0:nin+1,nin,2),cout(0:nout+1,nout,2)
         cout=0.d0
         nmin=min(nin,nout)
         cout(0:nmin+1,1:nmin,1:2)=cin(0:nmin+1,1:nmin,1:2)
         end subroutine transfer

      end module specialfuncs
!
!  module mpidata
!
!
!  last revised: 15 January 2011
!
      module mpidata
      implicit none
      integer :: group_comm,root_group_comm,base_rank,group_rank,root_group_rank, &
                 base_group,number_groups,proc_per_group,number_proc
      integer, allocatable :: mpi_sphere_index(:), mpi_sphere_number(:)

      contains
!
! allocates the processors into groups
!
!  last revised: 15 January 2011: original
!  20 April 2011: fixedorran=0 now looks for 2 groups.
!  10 october 2011: option for not storing matrices.  If fixorran=0, 2 groups, else
!     nproc groups
!  november 2011: near and far field translation differentiation
!  february 2013: completely rewritten.    now it sets a new group containing
!  min(nsphere*(nsphere-1)/2, numprocs) processors: this is the optimum size
!
         subroutine mpisetup(nsphere,newnumprocs,mpicomm)
         use mpidefs
         implicit none
         integer :: numprocs,rank,nsphere,mpicomm,newnumprocs,&
                    newgroup,i
         integer, allocatable :: grouplist(:)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='group',mpi_group=base_group)
         newnumprocs=nsphere*(nsphere-1)/2
         newnumprocs=max(newnumprocs,1)
         if(newnumprocs.ge.numprocs) then
            mpicomm=mpi_comm_world
            newnumprocs=numprocs
         else
            allocate(grouplist(0:newnumprocs-1))
            grouplist=(/(i,i=0,newnumprocs-1)/)
            call mstm_mpi(mpi_command='incl',mpi_group=base_group, &
            mpi_size=newnumprocs,mpi_new_group_list=grouplist, &
            mpi_new_group=newgroup)
            call mstm_mpi(mpi_command='create',mpi_group=newgroup, &
            mpi_new_comm=mpicomm)
         endif
         number_proc=numprocs
         proc_per_group=numprocs
         group_comm=mpi_comm_world
         root_group_comm=mpi_comm_world
         group_rank=rank
         root_group_rank=rank
         number_groups=1
         proc_per_group=numprocs
         end subroutine mpisetup

      end module mpidata

!
! module spheredata: used to 1) input sphere data, 2) dimension sphere data
! arrays, and 3) provide common access to the data in other subroutines.
!
!
!  last revised: 15 January 2011
!
!  30 March 2011: added optical activity
!  March 2013: a few more options
!
      module spheredata
      use specialfuncs
      use mpidata
      use numconstants
      implicit none
      integer, private :: numberspheres,numberiterations,fixedorrandom,numbertheta, &
               calcnf,nfplane,calctmatrix,runprintunit,calcamn,maxmemperproc, &
               trackiterations,nfoutdata,normalizesm,storetranmat,niterstep, &
               infinitemedium,excitedsphere,appendfile,numprocs,rank,&
               appendafile,appendnffile,smnumberprocessors,writespheredata,&
               incidentortargetframe,numbercommentlines,azimuthaverage
      integer, allocatable, private :: hostsphere(:)
      real(8), private :: lengthscalefactor,realriscalefactor,imriscalefactor,epsmie, &
                 epstran,epssoln,phideg,thetamindeg,thetamaxdeg,alphadeg, &
                 betadeg,epstcon,nfplanepos,nfplanevert(2,2),deltax,gammadeg,epspw, &
                 cgaussbeam,gaussbeamfocus(3),realchiralfactor,imchiralfactor,nfdistance, &
                 realrimedium,imagrimedium,eulerdeg(3),realchiralfactormedium, &
                 imagchiralfactormedium,phimindeg,phimaxdeg,deltathetadeg
      character(60), private :: positionfile,outputfile,nfoutputfile,tmatrixfile,printfile, &
                                amnfile,scatteringanglefile
      character(60), private, allocatable :: tmfile(:)
      character(128), private, allocatable :: commentline(:)
      real(8), private :: xspmax,xvsp
      real(8), private, allocatable :: rpos(:,:),xsp(:),scatteringanglearray(:,:)
      complex(8), private, allocatable :: ri(:,:)
      complex(8), private :: rimedium(2),chiralfactormedium
      logical, private, allocatable :: tmonfile(:)
      data numberiterations,fixedorrandom,numbertheta/2000,0,181/
      data calcamn,trackiterations,niterstep/1,1,20/
      data lengthscalefactor,realriscalefactor,imriscalefactor,epsmie, &
                 epstran,epssoln,phideg,thetamindeg,thetamaxdeg,alphadeg, &
                 betadeg,epstcon/1.d0,1.d0,1.d0,1.d-4,1.d-6,1.d-10,0.d0,0.d0, &
                 180.d0,0.d0,0.d0,1.d-6/
      data phimindeg,phimaxdeg,deltathetadeg/0.d0,0.d0,1.d0/
      data incidentortargetframe/0/
      data realchiralfactor,imchiralfactor/0.d0,0.d0/
      data normalizesm,storetranmat,nfdistance/0,0,1.0d8/
      data azimuthaverage/0/
      data maxmemperproc/1500/
      data cgaussbeam/0.d0/
      data gaussbeamfocus/0.d0,0.d0,0.d0/
      data calcnf,calctmatrix,nfoutdata/0,1,2/
      data nfplane,nfplanepos/1,0.d0/
      data nfplanevert/10.d0,10.d0,0.d0,0.d0/
      data runprintunit/6/
      data positionfile,outputfile,tmatrixfile,printfile/' ', &
           'mstm_default_out.dat','tm_default.dat',' '/
      data nfoutputfile/'nf_default.dat'/
      data scatteringanglefile/' '/
      data amnfile/' '/
      data realrimedium,imagrimedium/1.d0,0.d0/
      data realchiralfactormedium,imagchiralfactormedium/0.d0,0.d0/
      data infinitemedium,excitedsphere/0,0/
      data appendfile,appendafile,appendnffile/0,0,0/
      data eulerdeg/3*0.d0/
      data smnumberprocessors/10/
      data writespheredata,numbercommentlines/1,0/

      contains
!
! reads property and position information for the spheres.
! now includes t matrix file options
! february 2013
!
         subroutine readpositions(iunit,nsphere)
         implicit none
         integer :: iunit,i,nsphere,numrec,itemp
         real(8) :: rireal,riimag,betareal,betaimag,sdat(8)
         complex(8) :: ribulk,beta
         character*128 :: parmid

         if(allocated(xsp)) deallocate(xsp,rpos,ri,hostsphere,tmfile,tmonfile)
         allocate(xsp(0:nsphere),rpos(3,0:nsphere), &
                  ri(2,0:nsphere),hostsphere(nsphere),tmfile(nsphere), &
                  tmonfile(nsphere))
         rpos(:,0)=0.d0
         tmfile=' '
         tmonfile=.false.
         do i=1,nsphere
            sdat=1.d0
            read(iunit,'(a)',end=20) parmid
            call numberinstring(parmid,numrec)
            if(numrec.eq.4) then
               read(parmid,*) sdat(1:4)
            elseif(numrec.eq.6) then
               read(parmid,*) sdat(1:6)
            elseif(numrec.eq.8) then
               read(parmid,*) sdat(1:8)
            elseif(numrec.eq.5) then
               read(parmid,*) sdat(1:4),tmfile(i)
               tmonfile(i)=.true.
               sdat(5:8)=0.d0
            else
               exit
            endif
            if(tmonfile(i)) then
               open(2,file=tmfile(i))
               read(2,*) itemp,itemp,itemp
               read(2,*) itemp,xsp(i)
               close(2)
            else
               xsp(i)=sdat(1)*lengthscalefactor
            endif
            rpos(1:3,i)=sdat(2:4)*lengthscalefactor
            rireal=sdat(5)*realriscalefactor
            riimag=sdat(6)*imriscalefactor
            betareal=sdat(7)*realchiralfactor
            betaimag=sdat(8)*imchiralfactor
            ribulk=dcmplx(rireal,riimag)
            beta=dcmplx(betareal,betaimag)
            if(beta.eq.(0.d0,0.d0)) then
               ri(1,i)=ribulk
               ri(2,i)=ribulk
            else
               ri(1,i)=ribulk/(1.d0-beta*ribulk)
               ri(2,i)=ribulk/(1.d0+beta*ribulk)
            endif
         enddo
20       nsphere=min(nsphere,i-1)
         end subroutine readpositions

         subroutine setfileposition(iunit,string,num)
         implicit none
         integer :: iunit,num,i
         character :: string*35,frec*60
         i=0
         do while(i.lt.num)
            read(iunit,'(a)',end=20) frec
            frec=frec(:index(frec,' '))
            if(frec.eq.string) i=i+1
         enddo
         return
20       string='error'
         end subroutine setfileposition
!
! number of data points in a string
! march 2013
!
         subroutine numberinstring(parmid,numrec)
         implicit none
         integer :: i,l,numrec
         character*1 :: a
         character*128 :: parmid
         numrec=0
         l=len_trim(parmid)+1
         i=1
         do
            a=parmid(i:i)
            if(a.ne.' '.and.a.ne.','.and.ichar(a).ne.9) then
!
! start of a number
!
               numrec=numrec+1
!
! look for the delimeter
!
               do
                  i=i+1
                  if(i.ge.l) exit
                  a=parmid(i:i)
                  if(a.eq.' '.or.a.eq.','.or.ichar(a).eq.9) exit
               enddo
            endif
            i=i+1
            if(i.ge.l) exit
         enddo
         end subroutine numberinstring
!
!  Find the number of data points in input unit iunit, and reposition the unit to the
!  point after record containing parmid
!
!
!  last revised: 15 January 2011
!
         subroutine numberinrecord(iunit,parmid,numrec)
         implicit none
         integer :: numrec,iunit
         character*1 :: a
         character*35 :: parmid
         character*10 :: rec
         numrec=0
         do
            read(iunit,"(a)",advance="no",err=100,eor=100) a
            if(a.ne.' '.and.a.ne.','.and.ichar(a).ne.9) then
!
! start of a number
!
               numrec=numrec+1
!
! look for the delimeter
!
               do
                  read(iunit,"(a)",advance="no",err=100,eor=100) a
                  if(a.eq.' '.or.a.eq.','.or.ichar(a).eq.9) exit
               enddo
            endif
         enddo
100      if(parmid.eq.'rewind') then
            rewind(iunit)
         elseif(parmid.eq.'backspace') then
            backspace(iunit)
         else
            backspace(iunit)
            backspace(iunit)
            backspace(iunit)
            do
               read(iunit,'(a10)') rec
               if(rec.eq.parmid(1:10)) exit
            enddo
         endif
         end subroutine numberinrecord
!
!  inputdata:  reads parameters from inputfile
!              reads sphere data from position file
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fix output file initialization.
!  30 March 2011: added optical activity
!  November 2012: rewrote to allow repeated runs.    MPI is now initiallized
!     prior to calling inputdata, followed by all ranks sequentially calling
!     inputdata.
!  february 2013: tmatrix file options.
!  april 2013: nested multiple run option.
!              selection of coordinate frame for scattering matrix.
!              azimuth averaged scattering matrix.
!
         subroutine inputdata(inputfile,printdata,run_num,more_runs)
         use mpidefs
         implicit none
         integer :: imax,i,j,iunit,printdata,runnum,k, &
                    nextsphere,rank,numprocs,isphere,numberphi
         integer, save :: nsphere,nummultivars(3),posfilenumber,numbertheta0, &
                          nummultiloops,nestedmultiplerun
         integer, optional :: run_num
         logical :: moreruns,newoutfile,newafile,newnffile,posfilepres, &
                    multirun,loopdone,repeatread, &
                    newpositions,deltathetaspec,numthetaspec, &
                    autonfplanevert,switchloop,nestedloop
         logical, optional :: more_runs
         real(8) :: rmax,rtoi,rposmean(3),rij,xij(3),rijmax,euler(3),ftemp,&
                    thetatemp,phitemp,rposmax(3),rposmin(3),nfpad(2)
         real(8), save :: gaussbeamfocus0(3)
         real(8), save :: multivar(3,4),multivar1(3,4), &
                  multivar2(3,4),multistep(3,4)
         real(8), save, allocatable :: rpos0(:,:)
         complex(8) :: ribulk,beta
         character*35 :: parmid,tempparmid
         character*35, save :: multiparmid(3),endparmid,varmultiparmid
         character*30 :: inputfile
         data multirun,posfilepres/.false.,.false./
         data deltathetaspec,numthetaspec,autonfplanevert/.false.,.false.,.true./
         data switchloop,nestedloop/.true.,.true./
         data posfilenumber/0/
         data nummultiloops,nestedmultiplerun/0,1/
         data gaussbeamfocus0/3*0.d0/
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
!
!  run_num is the run number for the given job, and run_num=1 is the initial run.
!  If new_run is encountered in the input file the subroutine will set more_runs=.true.
!
         if(present(run_num)) then
            runnum=run_num
         else
            runnum=1
         endif
         newoutfile=.false.
         newafile=.false.
         newnffile=.false.
         moreruns=.false.
         repeatread=.true.
         newpositions=.false.
!
!  cycle through parameter input operations
!
         do while(repeatread)
            repeatread=.false.
            if(multirun) then
               open(1,file='multirun.inp',status='replace')
               do i=1,nummultiloops
                  write(1,'(a)') multiparmid(i)
                  write(1,*) multivar(i,1:nummultivars(i))
               enddo
!
! this is a specific case for fixed orientation, internal field calculation, in situations
! where the solution is unaltered from the first run
!
               if((varmultiparmid.eq.'near_field_plane_position' &
                .or.varmultiparmid.eq.'polarization_angle_deg' &
                .or.varmultiparmid.eq.'near_field_plane_vertices') &
                .and.runnum.gt.1.and.(.not.switchloop)) then
                  fixedorrandom=3
               else
                  fixedorrandom=0
               endif
               switchloop=.false.
               do i=nummultiloops,1,-1
                  varmultiparmid=multiparmid(i)
                  multivar(i,1:nummultivars(i))=multivar(i,1:nummultivars(i))+multistep(i,1:nummultivars(i))
                  loopdone=.true.
                  do j=1,nummultivars(i)
                     if(multistep(i,j).ne.0.d0) then
                        if(abs(multivar(i,j)-multivar1(i,j)).le.abs(multivar2(i,j)-multivar1(i,j))) then
                           loopdone=.false.
                        endif
                     endif
                  enddo
                  if(loopdone) then
                     multivar(i,1:nummultivars(i))=multivar1(i,1:nummultivars(i))
                     switchloop=.true.
                  elseif(nestedloop) then
                     exit
                  endif
               enddo
               if(loopdone) then
                  write(1,'(''end_of_options'')')
                  nummultiloops=0
                  multirun=.false.
               else
                  write(1,'(''new_run'')')
               endif
               close(1)
               open(1,file='multirun.inp')
            else
               open(1,file=inputfile)
!
!  find the run_num-1 occurence of new_run, and position the input file to
!  the next line.
!
               if(runnum.gt.1) then
                  parmid='new_run'
                  nsphere=numberspheres
                  call setfileposition(1,parmid,runnum-1)
                  if(parmid.eq.'error') then
                     if(rank.eq.0) then
                        write(runprintunit,'('' error finding data for run '',&
                        & i5)') runnum
                     endif
                     stop
                  endif
               endif
            endif
!
! the main loop of reading input parameters
!
            do
               read(1,'(a)',end=10) parmid
               parmid=parmid(:index(parmid,' '))
               if(parmid.eq.'number_spheres') then
                  read(1,*) ftemp
                  numberspheres=nint(ftemp)
                  nsphere=numberspheres
                  cycle
               endif
               if(parmid.eq.'sphere_position_file') then
                  read(1,'(a)') positionfile
                  positionfile=positionfile(:index(positionfile,' '))
                  if(positionfile.eq.'at_bottom') positionfile=' '
                  if(positionfile.ne.' ') then
                     posfilepres=.true.
                  endif
                  cycle
               endif
               if(parmid.eq.'output_file') then
                  read(1,'(a)') outputfile
                  outputfile=outputfile(:index(outputfile,' '))
                  if(outputfile.eq.' ') then
                     outputfile='mstm_default_out.dat'
                  endif
                  newoutfile=.true.
                  cycle
               endif
               if(parmid.eq.'run_print_file') then
                  read(1,'(a)') printfile
                  printfile=printfile(:index(printfile,' '))
                  if(printdata.eq.1) then
                     if((printfile.eq.' '.or.printfile.eq.'console')) then
                        printfile=' '
                        runprintunit=6
                     else
                        runprintunit=4
                        if(rank.eq.0) then
                           open(runprintunit,file=printfile)
                        endif
                     endif
                  else
                     runprintunit=6
                  endif
                  cycle
               endif
               if(parmid.eq.'append_output_file') then
                  read(1,*) appendfile
                  cycle
               endif
               if(parmid.eq.'write_sphere_data') then
                  read(1,*) writespheredata
                  cycle
               endif
               if(parmid.eq.'length_scale_factor') then
                  read(1,*) lengthscalefactor
                  cycle
               endif
               if(parmid.eq.'real_ref_index_scale_factor') then
                  read(1,*) realriscalefactor
                  cycle
               endif
               if(parmid.eq.'imag_ref_index_scale_factor') then
                  read(1,*) imriscalefactor
                  cycle
               endif
               if(parmid.eq.'medium_real_ref_index') then
                  read(1,*) realrimedium
                  cycle
               endif
               if(parmid.eq.'medium_imag_ref_index') then
                  read(1,*) imagrimedium
                  cycle
               endif
               if(parmid.eq.'medium_real_chiral_factor') then
                  read(1,*) realchiralfactormedium
                  cycle
               endif
               if(parmid.eq.'medium_imag_chiral_factor') then
                  read(1,*) imagchiralfactormedium
                  cycle
               endif
               if(parmid.eq.'real_chiral_factor') then
                  read(1,*) realchiralfactor
                  cycle
               endif
               if(parmid.eq.'imag_chiral_factor') then
                  read(1,*) imchiralfactor
                  cycle
               endif
               if(parmid.eq.'target_euler_angles_deg') then
                  read(1,*) eulerdeg
                  cycle
               endif
               if(parmid.eq.'mie_epsilon') then
                  read(1,*) epsmie
                  cycle
               endif
               if(parmid.eq.'translation_epsilon') then
                  read(1,*) epstran
                  cycle
               endif
               if(parmid.eq.'solution_epsilon') then
                  read(1,*) epssoln
                  cycle
               endif
               if(parmid.eq.'max_number_iterations') then
                  read(1,*) numberiterations
                  cycle
               endif
               if(parmid.eq.'max_memory_per_processor') then
                  read(1,*) maxmemperproc
                  cycle
               endif
               if(parmid.eq.'store_translation_matrix') then
                  read(1,*) storetranmat
                  cycle
               endif
               if(parmid.eq.'near_field_distance') then
                  read(1,*) nfdistance
                  cycle
               endif
               if(parmid.eq.'iterations_per_correction') then
                  read(1,*) niterstep
                  cycle
               endif
               if(parmid.eq.'sm_number_processors') then
                  read(1,*) smnumberprocessors
                  cycle
               endif
               if(parmid.eq.'fixed_or_random_orientation') then
                  read(1,*) fixedorrandom
                  cycle
               endif
               if(parmid.eq.'incident_or_target_frame') then
                  read(1,*) incidentortargetframe
                  cycle
               endif
               if(parmid.eq.'scattering_angle_file') then
                  read(1,'(a)') scatteringanglefile
                  scatteringanglefile=scatteringanglefile(:index(scatteringanglefile,' '))
                  cycle
               endif
               if(parmid.eq.'scattering_plane_angle_deg') then
                  read(1,*) phideg
                  phimindeg=phideg
                  phimaxdeg=phideg
                  cycle
               endif
               if(parmid.eq.'min_scattering_angle_deg') then
                  read(1,*) thetamindeg
                  cycle
               endif
               if(parmid.eq.'max_scattering_angle_deg') then
                  read(1,*) thetamaxdeg
                  cycle
               endif
               if(parmid.eq.'min_scattering_plane_angle_deg') then
                  read(1,*) phimindeg
                  cycle
               endif
               if(parmid.eq.'max_scattering_plane_angle_deg') then
                  read(1,*) phimaxdeg
                  cycle
               endif
               if(parmid.eq.'number_scattering_angles') then
                  read(1,*) numbertheta0
                  numbertheta=numbertheta0
                  deltathetaspec=.false.
                  cycle
               endif
               if(parmid.eq.'delta_scattering_angle_deg') then
                  read(1,*) deltathetadeg
                  deltathetaspec=.true.
                  cycle
               endif
               if(parmid.eq.'normalize_scattering_matrix') then
                  read(1,*) normalizesm
                  cycle
               endif
               if(parmid.eq.'azimuth_average_scattering_matrix') then
                  read(1,*) azimuthaverage
                  cycle
               endif
               if(parmid.eq.'incident_azimuth_angle_deg') then
                  read(1,*) alphadeg
                  cycle
               endif
               if(parmid.eq.'incident_polar_angle_deg') then
                  read(1,*) betadeg
                  cycle
               endif
               if(parmid.eq.'calculate_scattering_coefficients') then
                  read(1,*) calcamn
                  cycle
               endif
               if(parmid.eq.'scattering_coefficient_file') then
                  read(1,'(a)') amnfile
                  amnfile=amnfile(:index(amnfile,' '))
                  newafile=.true.
                  cycle
               endif
               if(parmid.eq.'track_iterations') then
                  read(1,*) trackiterations
                  cycle
               endif
               if(parmid.eq.'calculate_near_field') then
                  read(1,*) calcnf
                  cycle
               endif
               if(parmid.eq.'near_field_plane_coord') then
                  read(1,*) nfplane
                  cycle
               endif
               if(parmid.eq.'near_field_plane_position') then
                  read(1,*) nfplanepos
                  cycle
               endif
               if(parmid.eq.'near_field_plane_vertices') then
                  call numberinrecord(1,parmid,i)
                  if(i.eq.2) then
                     read(1,*) nfpad
                     autonfplanevert=.true.
                  else
                     read(1,*) nfplanevert
                     autonfplanevert=.false.
                  endif
                  cycle
               endif
               if(parmid.eq.'spacial_step_size') then
                  read(1,*) deltax
                  cycle
               endif
               if(parmid.eq.'polarization_angle_deg') then
                  read(1,*) gammadeg
                  cycle
               endif
               if(parmid.eq.'near_field_output_file') then
                  read(1,'(a)') nfoutputfile
                  if(nfoutputfile.eq.' ') then
                     nfoutputfile='mstm_nf_default_out.dat'
                  else
                     nfoutputfile=nfoutputfile(:index(nfoutputfile,' '))
                  endif
                  newnffile=.true.
                  cycle
               endif
               if(parmid.eq.'near_field_output_data') then
                  read(1,*) nfoutdata
                  cycle
               endif
               if(parmid.eq.'plane_wave_epsilon') then
                  read(1,*) epspw
                  cycle
               endif
               if(parmid.eq.'gaussian_beam_constant') then
                  read(1,*) cgaussbeam
                  cycle
               endif
               if(parmid.eq.'gaussian_beam_focal_point') then
                  read(1,*) gaussbeamfocus0(1),gaussbeamfocus0(2), &
                            gaussbeamfocus0(3)
                  cycle
               endif
               if(parmid.eq.'t_matrix_convergence_epsilon') then
                  read(1,*) epstcon
                  cycle
               endif
               if(parmid.eq.'calculate_t_matrix') then
                  read(1,*) calctmatrix
                  cycle
               endif
               if(parmid.eq.'infinite_medium') then
                  read(1,*) infinitemedium
                  cycle
               endif
               if(parmid.eq.'excited_sphere') then
                  read(1,*) excitedsphere
                  cycle
               endif
               if(parmid.eq.'sm_number_processors') then
                  read(1,*) smnumberprocessors
                  cycle
               endif
               if(parmid.eq.'t_matrix_file') then
                  read(1,'(a)') tmatrixfile
                  if(tmatrixfile.eq.' ') then
                     tmatrixfile='tmatrix-temp.dat'
                  else
                     tmatrixfile=tmatrixfile(:index(tmatrixfile,' '))
                  endif
                  cycle
               endif
               if(parmid.eq.'sphere_sizes_and_positions') then
                  if(posfilepres) then
                     if(rank.eq.0) then
                        write(runprintunit,'('' error: cannot specify both '', &
                        &'' sphere_sizes_and_positions and sphere_position_file'')')
                     endif
                     stop
                  endif
                  posfilepres=.false.
                  do isphere=1,nsphere
                     read(1,'(a)',end=10) parmid
                  enddo
                  posfilenumber=posfilenumber+1
                  cycle
               endif
               if(parmid.eq.'new_run') then
                  moreruns=.true.
                  if(multirun) then
                     endparmid=parmid
                  endif
                  exit
               endif
               if(parmid.eq.'end_of_options') then
                  moreruns=.false.
                  endparmid=parmid
                  exit
               endif
               if(parmid.eq.'multiple_run') then
                  nummultiloops=nummultiloops+1
                  multirun=.true.
                  moreruns=.true.
                  read(1,'(a)') tempparmid
                  tempparmid=tempparmid(:index(tempparmid,' '))
                  varmultiparmid=tempparmid
                  switchloop=.true.
                  multiparmid(nummultiloops)=tempparmid
                  call numberinrecord(1,tempparmid,nummultivars(nummultiloops))
                  nummultivars(nummultiloops)=nummultivars(nummultiloops)/3
                  i=nummultivars(nummultiloops)
                  read(1,*) multivar1(nummultiloops,1:i), &
                     multivar2(nummultiloops,1:i),multistep(nummultiloops,1:i)
                  multivar(nummultiloops,1:i)=multivar1(nummultiloops,1:i)
                  repeatread=.true.
                  cycle
               endif
               if(parmid.eq.'nested_multiple_run') then
                  read(1,*) nestedmultiplerun
                  if(nestedmultiplerun.eq.1) then
                     nestedloop=.true.
                  else
                     nestedloop=.false.
                  endif
                  cycle
               endif
               if(parmid.eq.'begin_comment') then
                  numbercommentlines=0
                  do while(parmid.ne.'end_comment')
                     numbercommentlines=numbercommentlines+1
                     read(1,'(a)',end=10) parmid
                     parmid=parmid(:index(parmid,' '))
                  enddo
                  numbercommentlines=numbercommentlines-1
                  if(numbercommentlines.gt.0.and.rank.eq.0) then
                     if(allocated(commentline)) deallocate(commentline)
                     allocate(commentline(numbercommentlines))
                     do i=1,numbercommentlines+1
                        backspace(1)
                     enddo
                     do i=1,numbercommentlines
                        read(1,'(a)') commentline(i)
                     enddo
                     read(1,'(a)',end=10) parmid
                  endif
                  cycle
               endif
               write(runprintunit,'('' warning: unknown parameter ID:'',a35)') parmid
            enddo
!
!  end of parameter input options.   Input of sphere data follows
!
10          close(1)
            if(present(more_runs)) then
               more_runs=moreruns
            endif
!
!  reading of parameters is repeated only if multirun was set: the repeat
!  initializes the multirun variable
!
         enddo
!
!  automatically append output file for multiple runs
!
         if((.not.newoutfile).and.runnum.gt.1) then
            appendfile=1
         endif
         if((.not.newafile).and.runnum.gt.1) then
            appendafile=1
         else
            appendafile=0
         endif
         if((.not.newnffile).and.runnum.gt.1) then
            appendnffile=1
         else
            appendnffile=0
         endif
!
! read sphere properties from position file
!
         iunit=1
         nsphere=numberspheres
!
! this loop runs once if numberspheres = number of listings in position file,
! and twice if it does not
!
         do
            if(posfilepres) then
               open(iunit,file=positionfile)
            else
               open(iunit,file=inputfile)
               parmid='sphere_sizes_and_positions'
               call setfileposition(iunit,parmid,posfilenumber)
            endif
!
! readpositions reads the sphere data from either input or sphere file.  The
! positions are store in module variable rpos
!
            call readpositions(iunit,nsphere)
            if(nsphere.eq.numberspheres) then
               close(iunit)
               exit
            endif
            numberspheres=nsphere
            if(rank.eq.0) then
               write(runprintunit,'('' insufficient points in position file. '',&
               &'' nsphere changed to '',i6)') numberspheres
            endif
            close(iunit)
         enddo
         newpositions=.true.
!
!  record the original sphere positions
!
         if(newpositions) then
            if(allocated(rpos0)) deallocate(rpos0)
            allocate(rpos0(3,0:nsphere))
            rpos0=rpos
         else
            rpos=rpos0
         endif
!
!  find the host spheres for the set
!
         call findhostspheres(nsphere,xsp(1:nsphere),rpos(:,1:nsphere),hostsphere)
!
! check for overlapping spheres, and find maximum translation
!
         rijmax=0.
         do i=1,nsphere
            do j=i+1,nsphere
               xij=rpos(:,i)-rpos(:,j)
               rij=sqrt(dot_product(xij,xij))
               if((rij.lt.xsp(i)+xsp(j)-0.0001d0).and.(rij.gt.abs(xsp(i)-xsp(j))+0.0001d0)) then
                  if(rank.eq.0) then
                     write(runprintunit,'('' warning: spheres '',i4,'' and '',&
                        & i4 '' overlap. '',&
                        & '' scaled distance:'' f8.4)') i,j,rij/(xsp(i)+xsp(j))
                  endif
               endif
               if(hostsphere(i).eq.hostsphere(j)) then
                  rijmax=max(rijmax,rij)
               endif
            enddo
         enddo
!
!  find the mean position of the external spheres and the volume equiv. radius
!
         rposmean=0.d0
         nextsphere=0
         xvsp=0.d0
         do i=1,nsphere
            if(hostsphere(i).eq.0) then
               rposmean=rposmean+rpos(:,i)
               nextsphere=nextsphere+1
               xvsp=xvsp+xsp(i)**3.d0
            endif
         enddo
         rposmean=rposmean/dble(nextsphere)
         xvsp=xvsp**(1.d0/3.d0)
!
!  the target origin is defined as the GB focal point.
!  shift the positions accordingly, and find the circumscribing radius
!
         gaussbeamfocus=gaussbeamfocus0*lengthscalefactor
         rmax=0.d0
         do i=1,nsphere
            rtoi=sqrt(dot_product(rpos(1:3,i)-rposmean(1:3), &
                   rpos(1:3,i)-rposmean(1:3)))+xsp(i)
            rpos(1:3,i)=rpos(1:3,i)-gaussbeamfocus(1:3)
            if(rtoi.gt.rmax) then
               rmax=rtoi
               imax=i
            endif
         enddo
         xspmax=rmax
!
!  xsp(0) is the circumscribing sphere size parameter
!
         xsp(0)=xspmax
         rpos(:,0)=rposmean(:)-gaussbeamfocus(:)
!
!  medium optical properties
!
         ribulk=dcmplx(realrimedium,imagrimedium)
         beta=dcmplx(realchiralfactormedium,imagchiralfactormedium)
         if(beta.eq.(0.d0,0.d0)) then
            ri(1,0)=ribulk
            ri(2,0)=ribulk
         else
            ri(1,0)=ribulk/(1.d0-beta*ribulk)
            ri(2,0)=ribulk/(1.d0+beta*ribulk)
         endif
         rimedium=ri(:,0)
!
!  coordinate rotation option
!
         if(eulerdeg(1).ne.0.d0.or.eulerdeg(2).ne.0.d0.or.eulerdeg(3).ne.0.d0) then
            euler=4.d0*datan(1.d0)/180.d0*eulerdeg
            do i=1,nsphere
               call eulerrotation(rpos(:,i),euler,1,rpos(:,i))
            enddo
         endif
!
! scattering angle array determination: new in march 2013
!
         if(scatteringanglefile.ne.' ') then
!
! if file is present, check number of points
!
            open(1,file=scatteringanglefile)
            i=0
            j=0
            do while(j.eq.0)
               i=i+1
               read(1,*,iostat=j) ftemp,ftemp
            enddo
            close(1)
            numbertheta=min(i-1,numbertheta0)
            if(rank.eq.0) then
               write(runprintunit,'('' number of angles in file:'',i10)') numbertheta
               call flush(runprintunit)
            endif
         else
!
! determine the number of points for automatic selection
!
            if(azimuthaverage.eq.1.or.fixedorrandom.eq.1) then
               phimindeg=0.d0
               phimaxdeg=0.d0
            endif
            if(deltathetaspec) then
               numbertheta0=nint((thetamaxdeg-thetamindeg)/deltathetadeg)+1
            else
               deltathetadeg=(thetamaxdeg-thetamindeg)/max(1,numbertheta0-1)
            endif
            numbertheta=0
            do i=1,numbertheta0
               thetatemp=thetamindeg &
                 +(i-1)*(thetamaxdeg-thetamindeg)/max(1.d0,dble(numbertheta0-1))
               numberphi=nint((phimaxdeg-phimindeg)*sin(thetatemp*pi/180.d0)/deltathetadeg)+1
               numbertheta=numbertheta+numberphi
            enddo
         endif
         if(allocated(scatteringanglearray)) deallocate(scatteringanglearray)
         allocate(scatteringanglearray(2,numbertheta))
         if(scatteringanglefile.ne.' ') then
            open(1,file=scatteringanglefile)
            do i=1,numbertheta
               read(1,*) scatteringanglearray(1:2,i)
            enddo
            close(1)
         else
            k=0
            do i=1,numbertheta0
               thetatemp=thetamindeg &
                 +(i-1)*(thetamaxdeg-thetamindeg)/max(1.d0,dble(numbertheta0-1))
               numberphi=nint((phimaxdeg-phimindeg)*sin(thetatemp*pi/180.d0)/deltathetadeg)+1
               do j=1,numberphi
                  phitemp=phimindeg+(phimaxdeg-phimindeg)*dble((j-1))/dble(max(1,numberphi-1))
                  k=k+1
                  scatteringanglearray(1:2,k) = (/thetatemp,phitemp/)
               enddo
            enddo
         endif
!
!  automatic near field plane specification
!
         if(autonfplanevert) then
            rposmax=-1.d-8
            rposmin=1.d-8
            do i=1,nsphere
               rposmax=max(rposmax,rpos0(:,i)+xsp(i))
               rposmin=min(rposmin,rpos0(:,i)-xsp(i))
            enddo
            if(nfplane.eq.1) then
               i=2
               j=3
            elseif(nfplane.eq.2) then
               i=3
               j=1
            else
               i=1
               j=2
            endif
            nfplanevert(1:2,1)=(/rposmin(i)-nfpad(1),rposmin(j)-nfpad(2)/)
            nfplanevert(1:2,2)=(/rposmax(i)+nfpad(1),rposmax(j)+nfpad(2)/)
         endif
!
!  write run data to run file and output file
!
         if(printdata.eq.1.and.rank.eq.0) then
            call writerundata(runprintunit,input_file=inputfile, &
                 run_num=runnum)
            call flush(runprintunit)
            if(appendfile.eq.0) then
               open(1,file=outputfile,status='replace',action='write')
            else
               open(1,file=outputfile,position='append')
            endif
            call writerundata(1,input_file=inputfile,run_num=runnum)
            close(1)
            numbercommentlines=0
         endif
         end subroutine inputdata
!
!  findhostspheres finds the host sphere of each sphere in the set.   host=0 for
!  external sphere.
!
!  december 2011
!  march 2013: something changed
!
         subroutine findhostspheres(nsphere_t,xsp_t,rpos_t,hostsphere_t)
         implicit none
         integer :: nsphere_t,hostsphere_t(nsphere_t),i,j
         real(8) :: xsp_t(nsphere_t),rpos_t(3,nsphere_t),xij(3),rij,xspmin

         hostsphere_t=0
         do i=1,nsphere_t
            xspmin=1.d6
            do j=1,nsphere_t
               if(xsp_t(j).gt.xsp_t(i)) then
                  xij(:)=rpos_t(:,j)-rpos_t(:,i)
                  rij=sqrt(dot_product(xij,xij))
                  if(rij.le.xsp_t(j)-xsp_t(i).and.xsp_t(j).lt.xspmin) then
                     hostsphere_t(i)=j
                     xspmin=xsp_t(j)
                  endif
               endif
            enddo
         enddo
         end subroutine findhostspheres
!
!  writes run data to output unit iunit
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  march 2013: more bells and whistles.
!
         subroutine writerundata(iunit,input_file,run_num)
         use mpidefs
         implicit none
         integer :: iunit,runnum,i,numprocs
         integer, optional :: run_num
         character*30, optional :: input_file
         character*1 :: lf
         if(iunit.ne.1) then
            lf = ' '
         else
            lf = '/'
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
         if(present(run_num)) then
            runnum=run_num
         else
            runnum=1
         endif
         if(present(input_file).and.runnum.eq.1) then
            write(iunit, '('' input file is '' '//lf//',a)') input_file
         endif
         if(runnum.gt.1) then
            write(iunit,'('' *****************************************************'')')
         endif
         if(numbercommentlines.gt.0) then
            write(iunit,'('' input file comments follow:'')')
            do i=1,numbercommentlines
               write(iunit,'(a)') trim(commentline(i))
            enddo
         endif
         if(runnum.eq.1) then
            write(iunit,'('' number of processors:'' '//lf//',i4)') numprocs
         endif
         write(iunit,'('' input parameters for run number '' '//lf//',i5)') runnum
         write(iunit,'('' number of spheres, volume size parameter:'' '//lf//',i5,e13.5)') &
                       numberspheres,xvsp
         write(iunit,'('' position file:'' '//lf//',a)') positionfile
         write(iunit,'('' output file:'' '//lf//',a)') outputfile
         write(iunit,'('' length, ref. indx. scale factors:'' '//lf//',3f8.3)') lengthscalefactor, &
                      realriscalefactor,imriscalefactor
         write(iunit,'('' chiral factors:'' '//lf//',2e13.5)')  &
                      realchiralfactor,imchiralfactor
         if(incidentortargetframe.eq.0) then
            write(iunit,'('' scattering matrix directions based on incident frame'')')
         else
            write(iunit,'('' scattering matrix directions based on target frame'')')
         endif
         if(scatteringanglefile.eq.' ') then
            write(iunit,'('' automatic scattering angle calculation'')')
            write(iunit,'('' thetamin, thetamax:'' '//lf//',2f9.1)') &
                        thetamindeg,thetamaxdeg
            if(fixedorrandom.eq.0.and.azimuthaverage.eq.0) then
               write(iunit,'('' phimin, phimax:'' '//lf//',2f9.1)') &
                        phimindeg,phimaxdeg
            endif
         else
            write(iunit,'('' scattering angle file:'' '//lf//',a)') scatteringanglefile
         endif
         if(fixedorrandom.eq.0.and.azimuthaverage.eq.1) then
            write(iunit,'('' scattering matrix is averaged over azimuthal angle'')')
         endif
         write(iunit,'('' number scattering angles:'' '//lf//',i9)') numbertheta
         write(iunit,'('' epsmie, epssoln, max number iterations:'' '//lf//',2e12.4,i5)') epsmie, &
                      epssoln, numberiterations
         write(iunit,'('' medium refractive index:'' '//lf//',2f8.4)') ri(1,0)
         write(iunit,'('' target euler rotation angles:'' '//lf//',3f10.4)') eulerdeg
         write(iunit,'('' far field kr, iterations/correction:'' '//lf//',e12.4,i5)') &
                 nfdistance,niterstep
         if(cgaussbeam.ne.0.d0) then
            write(iunit,'('' gaussian incident beam: 1/width:'' '//lf//',f9.4,)') cgaussbeam
            write(iunit,'('' beam focal point:'' '//lf//',3f9.3,)') gaussbeamfocus
         else
            write(iunit,'('' plane wave incidence'')')
         endif
         if(fixedorrandom.ne.1) then
            write(iunit,'('' fixed orientation calculations'')')
            write(iunit,'('' incident azimuth and polar angles:'' '//lf//',2f10.3)') &
                    alphadeg,betadeg
            if(scatteringanglefile.ne.' ') then
               write(iunit,'('' scattering plane angle:'' '//lf//',f10.3)') &
                    phideg
            endif
            write(iunit,'('' common expansion epsilon:'' '//lf//',e12.4)') epstran
            if(fixedorrandom.eq.0) then
               if(calcamn.eq.0) then
                  write(iunit,'('' scattering coefficients read from file '' '//lf//',a)') amnfile
               else
                  if(amnfile.ne.' ') then
                     write(iunit,'('' scattering coefficients calculated, stored in file '' '//lf//',a)')&
                         amnfile
                  else
                     write(iunit,'('' scattering coefficients calculated '')')
                  endif
               endif
            else
               write(iunit,'('' scattering coefficients same as previous run'')')
            endif
            if(calcnf.eq.1) then
               if(appendnffile.eq.0) then
                  write(iunit,'('' near field calculated, stored in file '' '//lf//',a)') nfoutputfile
               else
                  write(iunit,'('' near field calculated, appended to file '' '//lf//',a)') nfoutputfile
               endif
               write(iunit,'('' near field data output option: '' '//lf//',i4)') nfoutdata
               write(iunit,'('' near field plane, position: '' '//lf//', i4,f9.3)') nfplane, nfplanepos
               write(iunit,'('' near field plane vertices: '' '//lf//',4f9.3)') nfplanevert
               write(iunit,'('' spacial step size:'' '//lf//',f9.4)') deltax
               write(iunit,'('' polarization angle, deg.:'' '//lf//',f9.2)') gammadeg
               write(iunit,'('' plane wave epsilon:'' '//lf//',e13.5)') epspw
            endif
         else
            write(iunit,'('' random orientation calculations'')')
            if(calctmatrix.eq.0) then
               write(iunit,'('' t matrix read from file '' '//lf//',a)') tmatrixfile
            elseif(calctmatrix.eq.1) then
               write(iunit,'('' t matrix calculated, stored in file '' '//lf//',a)') tmatrixfile
               write(iunit,'('' t matrix convergence epsilon:'' '//lf//',e12.4)') epstcon
               if(excitedsphere.eq.0) then
                  write(iunit,'('' all spheres excited (total T matrix) '')')
               else
                  write(iunit,'('' T(i) calculated for sphere:'' '//lf//',i5)') excitedsphere
               endif
            else
               write(iunit,'('' t matrix calculated from end of file '' '//lf//',a)') tmatrixfile
               write(iunit,'('' t matrix convergence epsilon:'' '//lf//',e12.4)') epstcon
               if(excitedsphere.eq.0) then
                  write(iunit,'('' all spheres excited (total T matrix) '')')
               else
                  write(iunit,'('' T(i) calculated for sphere:'' '//lf//',i5)') excitedsphere
               endif
            endif
            if(infinitemedium.eq.1) then
               write(iunit,'('' infinite medium correction applied '')')
            endif
         endif
         end subroutine writerundata
!
!  getspheredata: retrieves sphere data
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  february 2013: particle T matrix file options.
!
         subroutine getspheredata(number_spheres, sphere_size_parameters, sphere_positions, &
              sphere_refractive_indices, volume_size_parameter,host_spheres, &
              tmatrix_on_file,tmatrix_file)
         implicit none
         integer, optional :: number_spheres,host_spheres(numberspheres)
         real(8), optional :: sphere_size_parameters(numberspheres), &
                              sphere_positions(3,numberspheres), volume_size_parameter
         complex(8), optional :: sphere_refractive_indices(2,0:numberspheres)
         character(60), optional :: tmatrix_file(numberspheres)
         logical, optional :: tmatrix_on_file(numberspheres)
         if (present(number_spheres)) number_spheres=numberspheres
         if (present(sphere_size_parameters)) sphere_size_parameters(1:numberspheres)=xsp(1:numberspheres)
         if (present(sphere_positions)) sphere_positions(:,1:numberspheres)=rpos(:,1:numberspheres)
         if (present(sphere_refractive_indices)) &
               sphere_refractive_indices(:,0:numberspheres)=ri(:,0:numberspheres)
         if (present(volume_size_parameter)) volume_size_parameter=xvsp
         if (present(host_spheres)) host_spheres(1:numberspheres)=hostsphere(1:numberspheres)
         if (present(tmatrix_file)) tmatrix_file(1:numberspheres)=tmfile(1:numberspheres)
         if (present(tmatrix_on_file)) tmatrix_on_file(1:numberspheres)=tmonfile(1:numberspheres)
         end subroutine getspheredata

         subroutine getspheredataone(sphere,sphere_size_parameter, sphere_position, &
              sphere_refractive_index,host_sphere,tmatrix_file,tmatrix_on_file)
         implicit none
         integer :: sphere
         integer, optional :: host_sphere
         real(8), optional :: sphere_size_parameter,sphere_position(3)
         complex(8), optional :: sphere_refractive_index(2)
         character(60), optional :: tmatrix_file
         logical, optional :: tmatrix_on_file
         if (present(sphere_size_parameter)) sphere_size_parameter=xsp(sphere)
         if (present(sphere_position)) sphere_position(:)=rpos(:,sphere)
         if (present(sphere_refractive_index)) &
               sphere_refractive_index(:)=ri(:,sphere)
         if (present(host_sphere)) host_sphere=hostsphere(sphere)
         if (present(tmatrix_file)) tmatrix_file=tmfile(sphere)
         if (present(tmatrix_on_file)) tmatrix_on_file=tmonfile(sphere)
         end subroutine getspheredataone
!
!  setspheredata: sets sphere data
!
         subroutine setspheredata(number_spheres, sphere_size_parameters, sphere_positions, &
              sphere_refractive_indices, volume_size_parameter,host_spheres)
         implicit none
         integer, optional :: number_spheres,host_spheres(*)
         real(8), optional :: sphere_size_parameters(*), &
                              sphere_positions(3,*), volume_size_parameter
         complex(8), optional :: sphere_refractive_indices(2,0:*)
         if (present(number_spheres)) then
            numberspheres=number_spheres
            if(allocated(xsp)) deallocate(xsp,rpos,ri,hostsphere)
            allocate(xsp(0:numberspheres),rpos(3,0:numberspheres),ri(2,0:numberspheres), &
                     hostsphere(numberspheres))
         endif
         if (present(sphere_size_parameters))    xsp(1:numberspheres)       =sphere_size_parameters(1:numberspheres)
         if (present(sphere_positions))          rpos(:,1:numberspheres)    =sphere_positions(:,1:numberspheres)
         if (present(sphere_refractive_indices)) ri(:,0:numberspheres)      =sphere_refractive_indices(:,0:numberspheres)
         if (present(sphere_size_parameters))    xsp(1:numberspheres)       =sphere_size_parameters(1:numberspheres)
         if (present(host_spheres))              hostsphere(1:numberspheres)=host_spheres(1:numberspheres)
         if (present(volume_size_parameter))     xvsp                       =volume_size_parameter
         end subroutine setspheredata
!
!  getrunparameters: retrieves run parameters read from input file
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  march 2013: more options
!
         subroutine getrunparameters(number_spheres,sphere_position_file,output_file, &
                       length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       max_number_iterations,fixed_or_random_orientation,scattering_plane_angle_deg, &
                       min_scattering_angle_deg,max_scattering_angle_deg,number_scattering_angles, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,calculate_near_field, &
                       near_field_plane_coord,near_field_plane_position,near_field_plane_vertices, &
                       spacial_step_size,polarization_angle_deg,near_field_output_file, &
                       plane_wave_epsilon,t_matrix_convergence_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point,calculate_t_matrix,t_matrix_file,run_print_file, &
                       run_print_unit,calculate_scattering_coefficients,scattering_coefficient_file, &
                       max_memory_per_processor,track_iterations,near_field_output_data, &
                       real_chiral_factor,imag_chiral_factor,normalize_scattering_matrix, &
                       store_translation_matrix,near_field_distance,iterations_per_correction,&
                       medium_ref_index,infinite_medium,excited_sphere,append_a_file,append_nf_file, &
                       target_euler_angles_deg,sm_number_processors,write_sphere_data, &
                       scattering_angle_array,incident_or_target_frame,azimuth_average_scattering_matrix)
         implicit none
         integer, optional :: number_spheres,max_number_iterations,fixed_or_random_orientation, &
                              number_scattering_angles,calculate_near_field,near_field_plane_coord, &
                              calculate_t_matrix,run_print_unit,calculate_scattering_coefficients, &
                              max_memory_per_processor,track_iterations,near_field_output_data, &
                              normalize_scattering_matrix,store_translation_matrix, &
                              iterations_per_correction,infinite_medium,excited_sphere, &
                              append_a_file,append_nf_file,sm_number_processors, &
                              write_sphere_data,incident_or_target_frame,azimuth_average_scattering_matrix
         real(8), optional :: length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       scattering_plane_angle_deg, &
                       min_scattering_angle_deg,max_scattering_angle_deg, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,t_matrix_convergence_epsilon, &
                       near_field_plane_position,near_field_plane_vertices(2,2),spacial_step_size, &
                       polarization_angle_deg,plane_wave_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point(3),real_chiral_factor,imag_chiral_factor, &
                       near_field_distance,target_euler_angles_deg(3), &
                       scattering_angle_array(2,*)
         complex(8), optional :: medium_ref_index
         character*60, optional :: sphere_position_file,output_file,near_field_output_file, &
                                   t_matrix_file,run_print_file,scattering_coefficient_file
         if(present(number_spheres))                    number_spheres                    =numberspheres
         if(present(sphere_position_file))              sphere_position_file              =positionfile
         if(present(output_file))                       output_file                       =outputfile
         if(present(length_scale_factor))               length_scale_factor               =lengthscalefactor
         if(present(real_ref_index_scale_factor))       real_ref_index_scale_factor       =realriscalefactor
         if(present(imag_ref_index_scale_factor))       imag_ref_index_scale_factor       =imriscalefactor
         if(present(mie_epsilon))                       mie_epsilon                       =epsmie
         if(present(translation_epsilon))               translation_epsilon               =epstran
         if(present(solution_epsilon))                  solution_epsilon                  =epssoln
         if(present(max_number_iterations))             max_number_iterations             =numberiterations
         if(present(track_iterations))                  track_iterations                  =trackiterations
         if(present(max_memory_per_processor))          max_memory_per_processor          =maxmemperproc
         if(present(fixed_or_random_orientation))       fixed_or_random_orientation       =fixedorrandom
         if(present(incident_or_target_frame))          incident_or_target_frame          =incidentortargetframe
         if(present(scattering_plane_angle_deg))        scattering_plane_angle_deg        =phideg
         if(present(min_scattering_angle_deg))          min_scattering_angle_deg          =thetamindeg
         if(present(max_scattering_angle_deg))          max_scattering_angle_deg          =thetamaxdeg
         if(present(number_scattering_angles))          number_scattering_angles          =numbertheta
         if(present(scattering_angle_array))            scattering_angle_array(1:2,1:numbertheta) &
                    =scatteringanglearray(1:2,1:numbertheta)
         if(present(normalize_scattering_matrix))       normalize_scattering_matrix       =normalizesm
         if(present(incident_azimuth_angle_deg))        incident_azimuth_angle_deg        =alphadeg
         if(present(incident_polar_angle_deg))          incident_polar_angle_deg          =betadeg
         if(present(t_matrix_convergence_epsilon))      t_matrix_convergence_epsilon      =epstcon
         if(present(calculate_near_field))              calculate_near_field              =calcnf
         if(present(near_field_plane_coord))            near_field_plane_coord            =nfplane
         if(present(near_field_plane_position))         near_field_plane_position         =nfplanepos
         if(present(near_field_plane_vertices))         near_field_plane_vertices         =nfplanevert
         if(present(spacial_step_size))                 spacial_step_size                 =deltax
         if(present(polarization_angle_deg))            polarization_angle_deg            =gammadeg
         if(present(near_field_output_file))            near_field_output_file            =nfoutputfile
         if(present(near_field_output_data))            near_field_output_data            =nfoutdata
         if(present(plane_wave_epsilon))                plane_wave_epsilon                =epspw
         if(present(gaussian_beam_constant))            gaussian_beam_constant            =cgaussbeam
         if(present(gaussian_beam_focal_point))         gaussian_beam_focal_point         =gaussbeamfocus
         if(present(t_matrix_file))                     t_matrix_file                     =tmatrixfile
         if(present(calculate_t_matrix))                calculate_t_matrix                =calctmatrix
         if(present(run_print_file))                    run_print_file                    =printfile
         if(present(run_print_unit))                    run_print_unit                    =runprintunit
         if(present(calculate_scattering_coefficients)) calculate_scattering_coefficients =calcamn
         if(present(scattering_coefficient_file))       scattering_coefficient_file       =amnfile
         if(present(real_chiral_factor))                real_chiral_factor                =realchiralfactor
         if(present(imag_chiral_factor))                imag_chiral_factor                =imchiralfactor
         if(present(store_translation_matrix))          store_translation_matrix          =storetranmat
         if(present(near_field_distance))               near_field_distance               =nfdistance
         if(present(iterations_per_correction))         iterations_per_correction         =niterstep
         if(present(medium_ref_index))                  medium_ref_index                  =rimedium(1)
         if(present(infinite_medium))                   infinite_medium                   =infinitemedium
         if(present(excited_sphere))                    excited_sphere                    =excitedsphere
         if(present(append_a_file))                     append_a_file                     =appendafile
         if(present(append_nf_file))                    append_nf_file                    =appendnffile
         if(present(target_euler_angles_deg))           target_euler_angles_deg           =eulerdeg
         if(present(sm_number_processors))              sm_number_processors              =smnumberprocessors
         if(present(write_sphere_data))                 write_sphere_data                 =writespheredata
         if(present(azimuth_average_scattering_matrix)) azimuth_average_scattering_matrix =azimuthaverage
         end subroutine getrunparameters
!
!  set run parameters: set run parameters
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine setrunparameters(number_spheres,sphere_position_file,output_file, &
                       length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       max_number_iterations,fixed_or_random_orientation,scattering_plane_angle_deg, &
                       min_scattering_angle_deg,max_scattering_angle_deg,number_scattering_angles, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,calculate_near_field, &
                       near_field_plane_coord,near_field_plane_position,near_field_plane_vertices, &
                       spacial_step_size,polarization_angle_deg,near_field_output_file, &
                       plane_wave_epsilon,t_matrix_convergence_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point,calculate_t_matrix,t_matrix_file,run_print_file, &
                       run_print_unit,calculate_scattering_coefficients,scattering_coefficient_file, &
                       max_memory_per_processor,track_iterations,near_field_output_data, &
                       real_chiral_factor,imag_chiral_factor,store_translation_matrix, &
                       near_field_distance,iterations_per_correction,medium_ref_index)
         implicit none
         integer, optional :: number_spheres,max_number_iterations,fixed_or_random_orientation, &
                              number_scattering_angles,calculate_near_field,near_field_plane_coord, &
                              calculate_t_matrix,run_print_unit,calculate_scattering_coefficients, &
                              max_memory_per_processor,track_iterations,near_field_output_data, &
                              store_translation_matrix,iterations_per_correction
         real(8), optional :: length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       scattering_plane_angle_deg,near_field_distance,&
                       min_scattering_angle_deg,max_scattering_angle_deg, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,t_matrix_convergence_epsilon, &
                       near_field_plane_position,near_field_plane_vertices(2,2),spacial_step_size, &
                       polarization_angle_deg,plane_wave_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point(3),real_chiral_factor,imag_chiral_factor
         complex(8), optional :: medium_ref_index
         character*30, optional :: sphere_position_file,output_file,near_field_output_file, &
                                   t_matrix_file,run_print_file,scattering_coefficient_file

         if(present(number_spheres))                     numberspheres        =number_spheres
         if(present(sphere_position_file))               positionfile         =sphere_position_file
         if(present(output_file))                        outputfile           =output_file
         if(present(length_scale_factor))                lengthscalefactor    =length_scale_factor
         if(present(real_ref_index_scale_factor))        realriscalefactor    =real_ref_index_scale_factor
         if(present(imag_ref_index_scale_factor))        imriscalefactor      =imag_ref_index_scale_factor
         if(present(mie_epsilon))                        epsmie               =mie_epsilon
         if(present(translation_epsilon))                epstran              =translation_epsilon
         if(present(solution_epsilon))                   epssoln              =solution_epsilon
         if(present(max_number_iterations))              numberiterations     =max_number_iterations
         if(present(track_iterations))                   trackiterations      =track_iterations
         if(present(max_memory_per_processor))           maxmemperproc        =max_memory_per_processor
         if(present(fixed_or_random_orientation))        fixedorrandom        =fixed_or_random_orientation
         if(present(scattering_plane_angle_deg))         phideg               =scattering_plane_angle_deg
         if(present(min_scattering_angle_deg))           thetamindeg          =min_scattering_angle_deg
         if(present(max_scattering_angle_deg))           thetamaxdeg          =max_scattering_angle_deg
         if(present(number_scattering_angles))           numbertheta          =number_scattering_angles
         if(present(incident_azimuth_angle_deg))         alphadeg             =incident_azimuth_angle_deg
         if(present(incident_polar_angle_deg))           betadeg              =incident_polar_angle_deg
         if(present(t_matrix_convergence_epsilon))       epstcon              =t_matrix_convergence_epsilon
         if(present(calculate_near_field))               calcnf               =calculate_near_field
         if(present(near_field_plane_coord))             nfplane              =near_field_plane_coord
         if(present(near_field_plane_position))          nfplanepos           =near_field_plane_position
         if(present(near_field_plane_vertices))          nfplanevert          =near_field_plane_vertices
         if(present(spacial_step_size))                  deltax               =spacial_step_size
         if(present(polarization_angle_deg))             gammadeg             =polarization_angle_deg
         if(present(near_field_output_file))             nfoutputfile         =near_field_output_file
         if(present(near_field_output_data))             nfoutdata            =near_field_output_data
         if(present(plane_wave_epsilon))                 epspw                =plane_wave_epsilon
         if(present(gaussian_beam_constant))             cgaussbeam           =gaussian_beam_constant
         if(present(gaussian_beam_focal_point))          gaussbeamfocus       =gaussian_beam_focal_point
         if(present(t_matrix_file))                      tmatrixfile          =t_matrix_file
         if(present(calculate_t_matrix))                 calctmatrix          =calculate_t_matrix
         if(present(run_print_file))                     printfile            =run_print_file
         if(present(run_print_unit))                     runprintunit         =run_print_unit
         if(present(calculate_scattering_coefficients))  calcamn              =calculate_scattering_coefficients
         if(present(scattering_coefficient_file))        amnfile              =scattering_coefficient_file
         if(present(real_chiral_factor))                 realchiralfactor     =real_chiral_factor
         if(present(imag_chiral_factor))                 imchiralfactor       =imag_chiral_factor
         if(present(store_translation_matrix))           storetranmat         =store_translation_matrix
         if(present(near_field_distance))                nfdistance           =near_field_distance
         if(present(iterations_per_correction))          niterstep            =iterations_per_correction
         if(present(medium_ref_index)) then
            rimedium=medium_ref_index
            ri(1,0)=rimedium(1)
            ri(2,0)=rimedium(2)
         endif
         end subroutine setrunparameters

      end module spheredata
!
! module miecoefdata: used to 1) calculate single sphere mie coefficient values,
! 2) store values in an allocated array, 3) provide common access to values, and
! 4) perform multiplication of coefficient values with vectors containing VWH scattering
! coefficients.
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
      module miecoefdata
      implicit none
      integer, private :: number_eqns,max_mie_order,number_spheres
      integer, allocatable, private :: mie_order(:),mie_offset(:),mie_block(:), &
               mie_block_offset(:),mie_number_field_expansions(:),tm_offset(:)
      logical, allocatable, private :: is_optically_active(:),tm_on_file(:)
      real(8), allocatable, private :: qext_mie(:),qabs_mie(:)
      complex(8), allocatable, private :: an_mie(:),cn_mie(:),un_mie(:),vn_mie(:), &
               dn_mie(:),an_inv_mie(:),stored_tm(:)
      interface getmiedata
            module procedure getmiedataall, getmiedataone
      end interface getmiedata

      contains
!
!  calculation of the max order of sphere expansions and storage of mie coefficients
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  april 2012: all spheres are assumed OA: l/r formulation
!  february 2013: tmatrix file option.
!
         subroutine miecoefcalc(nsphere,xsp,ri,hostsphere,numberfieldexp,qeps, &
         tmatrix_file)
         use mpidefs
         use spheredata
         implicit none
         integer :: i,nodrn,nsphere,ntermstot,nblktot,nterms,hostsphere(nsphere), &
                    n1,n2,numberfieldexp(nsphere),rank,sizetmstore,n, &
                    j,tint(1),iunit
         logical :: newtm
         real(8) :: qext,qsca,qeps,xsp(nsphere),qabs
         complex(8) :: ri(2,0:nsphere),rihost(2)
         complex(8), allocatable :: anp(:,:,:),cnp(:,:,:),unp(:,:,:), &
                     vnp(:,:,:),dnp(:,:,:),anpinv(:,:,:)
         character(60) :: tmtemp
         character(60), optional :: tmatrix_file(nsphere)
         if(allocated(mie_order)) deallocate(mie_order,mie_offset,mie_block, &
                      mie_block_offset,qext_mie,qabs_mie,is_optically_active, &
                      mie_number_field_expansions,tm_on_file,tm_offset)
         allocate(mie_order(nsphere),mie_offset(nsphere+1),mie_block(nsphere), &
                      mie_block_offset(nsphere+1),qext_mie(nsphere), &
                      qabs_mie(nsphere),is_optically_active(nsphere), &
                      mie_number_field_expansions(nsphere),tm_on_file(nsphere), &
                      tm_offset(nsphere))
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call getrunparameters(run_print_unit=iunit)
         ntermstot=0
         nblktot=0
         max_mie_order=0
         mie_number_field_expansions=numberfieldexp
         number_spheres=nsphere
         tm_on_file=.false.
         sizetmstore=0
!
!  executed if tmatrix file present.   These lines
!    read and set orders, offsets for spheres
!    sets up storage space
!    identifies repeated files and compresses accordingly
!
         if(present(tmatrix_file)) then
            do i=1,nsphere
               tmtemp=tmatrix_file(i)
               tmtemp=tmtemp(:index(tmtemp,' '))
               if(tmtemp.ne.' ') then
                  if(numberfieldexp(i).eq.2) then
                     if(rank.eq.0) then
                        write(iunit,'('' error: sphere '',i5,'' with TM file option '', &
                      &''cannot have inclusions'')') i
                     endif
                     stop
                  endif
                  tm_on_file(i)=.true.
                  newtm=.true.
                  do j=1,i-1
                     if(tm_on_file(j).and. &
                        (tmatrix_file(j).eq.tmatrix_file(i))) then
                        newtm=.false.
                        exit
                     endif
                  enddo
                  if(newtm) then
                     if(rank.eq.0) then
                        open(3,file=tmtemp)
                        read(3,*) n,nodrn,n
                        close(3)
                     endif
                     tint(1)=nodrn
                     call mstm_mpi(mpi_command='bcast',mpi_number=1, &
                          mpi_send_buf_i=tint(1),mpi_rank=0)
                     nodrn=tint(1)
                     tm_offset(i)=sizetmstore
                     sizetmstore=sizetmstore+(2*nodrn*(nodrn+2))*(2*nodrn*(nodrn+2))
                  else
                     tm_offset(i)=tm_offset(j)
                     nodrn=mie_order(j)
                  endif
                  numberfieldexp(i)=1
                  mie_order(i)=nodrn
                  mie_block(i)=nodrn*(nodrn+2)*2
                  mie_number_field_expansions(i)=1
               endif
            enddo
            if(sizetmstore.gt.0) then
               if(allocated(stored_tm)) deallocate(stored_tm)
               allocate(stored_tm(sizetmstore))
            endif
         endif
!
! march 2013: calculates orders, and
!             forces host to have at least same order as constituents.
!
         do i=1,nsphere
            if(.not.tm_on_file(i)) then
               rihost=ri(:,hostsphere(i))
               call mieoa(xsp(i),ri(1,i),nodrn,qeps,qext,qsca,qabs, &
                       ri_medium=rihost)
               mie_order(i)=nodrn
            endif
         enddo
         do i=1,nsphere
            j=hostsphere(i)
            do while(j.ne.0)
               mie_order(j)=max(mie_order(j),mie_order(i))
               j=hostsphere(j)
            enddo
         enddo
!
!  calculate the order limits and efficiencies
!
         do i=1,nsphere
            rihost=ri(:,hostsphere(i))
            if(ri(1,i).eq.ri(2,i)) then
               is_optically_active(i)=.false.
            else
               is_optically_active(i)=.true.
            endif
!
!  this step reads the t matrix if present and if needed
!
            if(tm_on_file(i)) then
               newtm=.true.
               do j=1,i-1
                  if(tm_on_file(j).and. &
                     (tmatrix_file(j).eq.tmatrix_file(i))) then
                     newtm=.false.
                     exit
                  endif
               enddo
               if(newtm) then
                  call readtmatrix(tmatrix_file(i),i,nodrn,qext,qabs)
                  qext_mie(i)=qext*2./xsp(i)/xsp(i)
                  qabs_mie(i)=qabs*2./xsp(i)/xsp(i)
               else
                  qext_mie(i)=qext_mie(j)
                  qabs_mie(i)=qabs_mie(j)
               endif
            else
               nodrn=mie_order(i)
               call mieoa(xsp(i),ri(1,i),nodrn,0.d0,qext,qsca,qabs, &
                       ri_medium=rihost)
               nterms=4*nodrn
               mie_offset(i)=ntermstot
               ntermstot=ntermstot+nterms
               qext_mie(i)=qext
               qabs_mie(i)=qabs
               mie_block_offset(i)=nblktot
               mie_order(i)=nodrn
               mie_block(i)=nodrn*(nodrn+2)*2
            endif
            nblktot=nblktot+mie_block(i)*numberfieldexp(i)
            max_mie_order=max(max_mie_order,mie_order(i))
         enddo
         mie_offset(nsphere+1)=ntermstot
         mie_block_offset(nsphere+1)=nblktot
         number_eqns=nblktot
!
! calculate the mie coefficients, and store in memory
!
         if(allocated(an_mie)) deallocate(an_mie,cn_mie,un_mie,vn_mie,dn_mie, &
                               an_inv_mie)
         allocate(an_mie(ntermstot),cn_mie(ntermstot),un_mie(ntermstot), &
                  vn_mie(ntermstot),dn_mie(ntermstot),an_inv_mie(ntermstot))
         do i=1,nsphere
            if(.not.tm_on_file(i)) then
               nodrn=mie_order(i)
               rihost=ri(:,hostsphere(i))
               allocate(anp(2,2,nodrn),cnp(2,2,nodrn),unp(2,2,nodrn), &
                  vnp(2,2,nodrn),dnp(2,2,nodrn),anpinv(2,2,nodrn))
               call mieoa(xsp(i),ri(1,i),nodrn,0.d0,qext,qsca,qabs, &
                    anp_mie=anp,cnp_mie=cnp,dnp_mie=dnp, &
                    unp_mie=unp,vnp_mie=vnp,anp_inv_mie=anpinv, &
                    ri_medium=rihost)
               nterms=4*nodrn
               n1=mie_offset(i)+1
               n2=mie_offset(i)+nterms
               an_mie(n1:n2)=reshape(anp(1:2,1:2,1:nodrn),(/nterms/))
               cn_mie(n1:n2)=reshape(cnp(1:2,1:2,1:nodrn),(/nterms/))
               dn_mie(n1:n2)=reshape(dnp(1:2,1:2,1:nodrn),(/nterms/))
               un_mie(n1:n2)=reshape(unp(1:2,1:2,1:nodrn),(/nterms/))
               vn_mie(n1:n2)=reshape(vnp(1:2,1:2,1:nodrn),(/nterms/))
               an_inv_mie(n1:n2)=reshape(anpinv(1:2,1:2,1:nodrn),(/nterms/))
               deallocate(anp,cnp,unp,vnp,dnp,anpinv)
            endif
         enddo
         end subroutine miecoefcalc
!
! reads and stores a particle T matrix, for subsequent use in multmiecoeffmult
! the file structure of the T matrix is identical to that written during a T
! matrix solution.   It is TE-TM based.
! February 2013.
!
         subroutine readtmatrix(tmfile,sphere,nodrt,qext,qsca)
         use mpidefs
         implicit none
         integer :: nodrt,nblkt,m,n,p,k,l,q,kl,klm,sizetm,nstop,i, &
                    mn,mnm,lt,kt,qt,nt,mt,ikm,nspheret,ma,na,la,ka,sphere, &
                    rank,numprocs,tint(1)
         real(8) :: qext,qsca,qabs,fc(4)
         complex(8) :: tm22(2,2)
         complex(8), allocatable :: tc(:,:,:,:),tct(:,:,:,:,:,:)
         character :: tmfile*60
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
         if(rank.eq.0) then
            open(3,file=tmfile)
            read(3,*) m,nodrt,n
            read(3,*) nspheret,fc(1)
         endif
         if(numprocs.gt.1) then
            tint(1)=nodrt
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_i=tint, &
              mpi_number=1,mpi_rank=0)
            nodrt=tint(1)
         endif
         nblkt=nodrt*(nodrt+2)
         sizetm=4*nblkt*nblkt
         if(rank.eq.0) then
            allocate(tc(2,nblkt,2,nblkt), &
               tct(0:nodrt+1,nodrt,2,0:nodrt+1,nodrt,2))
            tc=(0.,0.)
            do l=1,nodrt
               nstop=l
               do k=-l,l
                  kl=l*(l+1)+k
                  klm=l*(l+1)-k
                  do q=1,2
                     read(3,*) lt,kt,qt
                     do n=1,nstop
                        do m=-n,n
                           mn=n*(n+1)+m
                           mnm=n*(n+1)-m
                           read(3,*) nt,mt,fc
                           tc(1,mn,q,kl)=cmplx(fc(1),fc(2))
                           tc(2,mn,q,kl)=cmplx(fc(3),fc(4))
                           if(n.lt.l) then
                              ikm=(-1)**(m+k)
                              do p=1,2
                                 tc(q,klm,p,mnm)=tc(p,mn,q,kl)*ikm
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               enddo
               do i=1,nspheret
                  read(3,*) nt,qext,qabs,qsca
               enddo
            enddo
            close(3)
            qext=0.d0
            qsca=0.d0
            do n=1,nodrt
               do m=-n,n
                  mn=n*(n+1)+m
                  if(m.ge.0) then
                     ma=m
                     na=n
                  else
                     ma=n+1
                     na=-m
                  endif
                  do p=1,2
                     qext=qext-tc(p,mn,p,mn)
                  enddo
                  do l=1,nodrt
                     do k=-l,l
                        kl=l*(l+1)+k
                        do p=1,2
                           do q=1,2
                              qsca=qsca+tc(p,mn,q,kl)*conjg(tc(p,mn,q,kl))
                           enddo
                        enddo
                        if(k.ge.0) then
                           ka=k
                           la=l
                        else
                           ka=l+1
                           la=-k
                        endif
                        tm22(1:2,1:2)=tc(1:2,mn,1:2,kl)
                        call lrtomodetran(tm22,tm22)
                        tct(ma,na,1:2,ka,la,1:2)=tm22
                     enddo
                  enddo
               enddo
            enddo
            qabs=qext-qsca
            stored_tm(tm_offset(sphere)+1:tm_offset(sphere)+sizetm) &
              = reshape(tct,(/sizetm/))
         endif
         if(numprocs.gt.1) then
            fc(1)=qext
            fc(2)=qsca
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=fc(1), &
               mpi_number=1,mpi_rank=0)
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=fc(2), &
               mpi_number=1,mpi_rank=0)
            qext=fc(1)
            qsca=fc(2)
            call mstm_mpi(mpi_command='bcast', &
               mpi_send_buf_dc=stored_tm(tm_offset(sphere)+1:tm_offset(sphere)+sizetm), &
               mpi_number=sizetm,mpi_rank=0)
         endif
         end subroutine readtmatrix
!
! transformation between lr and te tm basis
!
         subroutine lrtomodetran(at,am)
         implicit none
         complex(8) :: a(2,2),am(2,2),at(2,2)
         a=at
         am(1,1)=(a(1,1) + a(1,2) + a(2,1) + a(2,2))/2.
         am(1,2)=(a(1,1) - a(1,2) + a(2,1) - a(2,2))/2.
         am(2,1)=(a(1,1) + a(1,2) - a(2,1) - a(2,2))/2.
         am(2,2)=(a(1,1) - a(1,2) - a(2,1) + a(2,2))/2.
         end subroutine lrtomodetran
!
! optically active lorenz/mie coefficients
! original 30 March 2011
! April 2012: generalized LR formulation, generalized mie coefficients
!
         subroutine mieoa(x,ri,nodr0,qeps,qext,qsca,qabs,anp_mie,dnp_mie, &
           unp_mie,vnp_mie,cnp_mie,ri_medium,anp_inv_mie)
         use specialfuncs
         implicit none
         integer :: nstop,n,i,p,q,nodr0,s,t,ss,st
         real(8) :: x,qeps,qext,qsca,fn1,err,qextt,qabs
         real(8), allocatable :: qext1(:)
         complex(8) :: ci,ri(2),xri(2,2),rii(2,2),ribulk(2),psip(2,2), &
                     xip(2,2),gmatinv(2,2),bmatinv(2,2),gmat(2,2),bmat(2,2),&
                     amat(2,2),dmat(2,2),umat(2,2),vmat(2,2),amatinv(2,2)
         complex(8), optional :: anp_mie(2,2,*),dnp_mie(2,2,*),ri_medium(2), &
                                 unp_mie(2,2,*),vnp_mie(2,2,*),cnp_mie(2,2,*), &
                                 anp_inv_mie(2,2,*)
         complex(8), allocatable :: psi(:,:,:),xi(:,:,:)
         data ci/(0.d0,1.d0)/

         if(present(ri_medium)) then
            rii(:,1)=ri_medium
         else
            rii(:,1)=(/(1.d0,0.d0),(1.d0,0.d0)/)
         endif
         rii(:,2)=ri
         ribulk(:)=2.d0/(1.d0/rii(1,:)+1.d0/rii(2,:))
         xri=x*rii
         if(qeps.gt.0.) then
            nstop=nint(x+4.*x**(1./3.))+5.
         elseif(qeps.lt.0.) then
            nstop=ceiling(-qeps)
            nodr0=nstop
         else
            nstop=nodr0
         endif
         allocate(psi(0:nstop+1,2,2),xi(0:nstop+1,2,2),qext1(nstop))

         do i=1,2
            do p=1,2
               call cricbessel(nstop+1,xri(p,i),psi(0,p,i))
               call crichankel(nstop+1,xri(p,i),xi(0,p,i))
            enddo
         enddo
         qabs=0.d0
         qsca=0.d0
         qext=0.d0
         do n=1,nstop
            psip(:,:) = psi(n-1,:,:)-dble(n)*psi(n,:,:)/xri(:,:)
            xip(:,:) = xi(n-1,:,:)-dble(n)*xi(n,:,:)/xri(:,:)
            do s=1,2
               ss=(-1)**s
               do t=1,2
                  st=(-1)**t
                  gmat(s,t)=(ss*ribulk(1)+st*ribulk(2))/(rii(s,2)*rii(t,1)) &
                    *(psip(s,2)*xi(n,t,1)-ss*st*psi(n,s,2)*xip(t,1))
                  amat(s,t)=(ss*ribulk(1)+st*ribulk(2))/(rii(s,2)*rii(t,1)) &
                    *(psip(s,2)*psi(n,t,1)-ss*st*psi(n,s,2)*psip(t,1))
!                  umat(s,t)=(ss*ribulk(2)+st*ribulk(2))/(rii(s,2)*rii(t,2)) &
!                    *(xip(s,2)*psi(n,t,2)-ss*st*xi(n,s,2)*psip(t,2))
                  umat(s,t)=(ss*ribulk(2)+st*ribulk(2))/(rii(s,2)*rii(t,2))*ci
                  bmat(s,t)=(ss*ribulk(2)+st*ribulk(1))/(rii(s,1)*rii(t,2)) &
                    *(xip(s,1)*psi(n,t,2)-ss*st*xi(n,s,1)*psip(t,2))
!                  dmat(s,t)=(ss*ribulk(1)+st*ribulk(1))/(rii(s,1)*rii(t,1)) &
!                    *(psip(s,1)*xi(n,t,1)-ss*st*psi(n,s,1)*xip(t,1))
                  dmat(s,t)=-(ss*ribulk(1)+st*ribulk(1))/(rii(s,1)*rii(t,1))*ci
                  vmat(s,t)=(ss*ribulk(2)+st*ribulk(1))/(rii(s,1)*rii(t,2)) &
                    *(xip(s,1)*xi(n,t,2)-ss*st*xi(n,s,1)*xip(t,2))
               enddo
            enddo

            call twobytwoinverse(gmat,gmatinv)
            call twobytwoinverse(bmat,bmatinv)

            amat=-matmul(gmatinv,amat)
            umat=-matmul(gmatinv,umat)
            dmat=-matmul(bmatinv,dmat)
            vmat=-matmul(bmatinv,vmat)

            if(present(anp_mie)) then
               anp_mie(:,:,n)=amat(:,:)
            endif
            if(present(dnp_mie)) then
               dnp_mie(:,:,n)=dmat(:,:)
            endif
            if(present(unp_mie)) then
               unp_mie(:,:,n)=umat(:,:)
            endif
            if(present(vnp_mie)) then
               vnp_mie(:,:,n)=vmat(:,:)
            endif
            if(present(cnp_mie)) then
               call twobytwoinverse(amat,amatinv)
               cnp_mie(:,:,n)=matmul(dmat,amatinv)
            endif
            if(present(anp_inv_mie)) then
               call twobytwoinverse(amat,amatinv)
               anp_inv_mie(:,:,n)=amatinv(:,:)
            endif
            qext1(n)=0.d0
            fn1=n+n+1
            do p=1,2
               do q=1,2
                  qsca=qsca+fn1*cdabs(amat(p,q))*cdabs(amat(p,q))
               enddo
               qext1(n)=qext1(n)-fn1*dble(amat(p,p))
            enddo
            qext=qext+qext1(n)
         enddo

         if(qeps.gt.0.d0) then
            qextt=qext
            qext=0.
            do n=1,nstop
               qext=qext+qext1(n)
               err=abs(1.d0-qext/qextt)
               if(err.lt.qeps.or.n.eq.nstop) exit
            enddo
            nodr0=n
         endif
         qsca=2./x/x*qsca
         qext=2./x/x*qext
         nstop=min(n,nstop)
         return
         end subroutine mieoa
!
!  retrieve the array of mie data
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine getmiedataall(sphere_order, sphere_block, &
                     sphere_order_offset, sphere_block_offset, sphere_qext, &
                     sphere_qabs, sphere_mie_coefficients, sphere_int_mie_coefficients, &
                     number_equations, max_order,optically_active, &
                     number_field_expansions)
         use spheredata
         implicit none
         integer, optional :: sphere_order(:), sphere_block(:), sphere_order_offset(:), &
                     sphere_block_offset(:),number_equations, max_order, &
                     number_field_expansions(:)
         integer :: nsphere
         logical, optional :: optically_active(:)
         real(8), optional :: sphere_qext(:), sphere_qabs(:)
         complex(8), optional :: sphere_mie_coefficients(:), &
                       sphere_int_mie_coefficients(:)
         call getspheredata(number_spheres=nsphere)
         if(present(sphere_order)) sphere_order=mie_order
         if(present(number_field_expansions)) &
             number_field_expansions=mie_number_field_expansions
         if(present(sphere_block)) sphere_block=mie_block
         if(present(sphere_order_offset)) sphere_order_offset=mie_offset
         if(present(sphere_block_offset)) sphere_block_offset=mie_block_offset
         if(present(sphere_qext)) sphere_qext=qext_mie
         if(present(sphere_qabs)) sphere_qabs=qabs_mie
         if(present(number_equations)) number_equations=number_eqns
         if(present(max_order)) max_order=max_mie_order
         if(present(optically_active)) optically_active=is_optically_active
         if(present(sphere_mie_coefficients)) sphere_mie_coefficients=an_mie
         if(present(sphere_int_mie_coefficients)) sphere_int_mie_coefficients=cn_mie
         end subroutine getmiedataall
!
!  retrieve mie data for a single sphere
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine getmiedataone(which_sphere, sphere_order, sphere_block, &
                     sphere_order_offset, sphere_block_offset, sphere_qext, &
                     sphere_qabs, sphere_mie_coefficients, sphere_int_mie_coefficients, &
                     number_equations, max_order, optically_active)
         use spheredata
         implicit none
         integer, optional :: sphere_order, sphere_block, sphere_order_offset, &
                     sphere_block_offset, number_equations, max_order
         integer :: which_sphere,n1,n2,nodrn
         integer :: i
         logical, optional :: optically_active
         real(8), optional :: sphere_qext, sphere_qabs
         complex(8), optional :: sphere_mie_coefficients(:), sphere_int_mie_coefficients(:)
         i=which_sphere
         if(present(sphere_order)) sphere_order=mie_order(i)
         if(present(sphere_block)) sphere_block=mie_block(i)
         if(present(sphere_order_offset)) sphere_order_offset=mie_offset(i)
         if(present(sphere_block_offset)) sphere_block_offset=mie_block_offset(i)
         if(present(sphere_qext)) sphere_qext=qext_mie(i)
         if(present(sphere_qabs)) sphere_qabs=qabs_mie(i)
         if(present(number_equations)) number_equations=number_eqns
         if(present(max_order)) max_order=max_mie_order
         if(present(optically_active)) optically_active=is_optically_active(i)
         if(present(sphere_mie_coefficients)) then
            nodrn=mie_order(i)
            n1=mie_offset(i)+1
            n2=mie_offset(i+1)
            sphere_mie_coefficients(1:n2-n1+1)=an_mie(n1:n2)
         endif
         if(present(sphere_int_mie_coefficients)) then
            nodrn=mie_order(i)
            n1=mie_offset(i)+1
            n2=mie_offset(i+1)
            sphere_int_mie_coefficients(1:n2-n1+1)=cn_mie(n1:n2)
         endif
         end subroutine getmiedataone
!
!
!  multiplies coefficients for sphere i by appropriate lm coefficient.
!  lr, oa model
!  april 2012
!
         subroutine onemiecoeffmult(i,nodr,cx,cy,mie_coefficient)
         use spheredata
         implicit none
         integer :: i,n,p,nodr,n1,n2,nterms
         complex(8) :: cx(0:nodr+1,nodr,2),cy(0:nodr+1,nodr,2)
         complex(8), allocatable :: an1(:,:,:)
         character*1, optional :: mie_coefficient
         character*1 :: miecoefficient

         if(present(mie_coefficient)) then
            miecoefficient=mie_coefficient
         else
            miecoefficient='a'
         endif
         nterms=4*nodr
         n1=mie_offset(i)+1
         n2=mie_offset(i)+nterms
         allocate(an1(2,2,nodr))
         if(miecoefficient.eq.'a') then
            an1=reshape(an_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'c') then
            an1=reshape(cn_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'d') then
            an1=reshape(dn_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'u') then
            an1=reshape(un_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'v') then
            an1=reshape(vn_mie(n1:n2),(/2,2,nodr/))
         endif
         do n=1,nodr
            do p=1,2
               cy(n+1,n:1:-1,p)=an1(p,1,n)*cx(n+1,n:1:-1,1) &
                 +an1(p,2,n)*cx(n+1,n:1:-1,2)
               cy(0:n,n,p)=an1(p,1,n)*cx(0:n,n,1) &
                 +an1(p,2,n)*cx(0:n,n,2)
            enddo
         enddo
         deallocate(an1)
         end subroutine onemiecoeffmult
!
! generalized mie coefficient mult:
!  (a,f) = (generalized mie matrix)*(g,b)
! idir not = 1 does the transpose.
! aout is written over in this one.
! february 2013: tmatrix file option
!
         subroutine multmiecoeffmult(neqns,nrhs,idir,ain,aout,rhs_list)
         implicit none
         integer :: neqns,i,n,p,q,nodri,nblki,noffi,n1,n2,b11,b12,b21,b22,idir, &
                    nterms,j,nrhs
         complex(8) :: ain(neqns*nrhs),aout(neqns*nrhs)
         complex(8), allocatable :: gin_t(:,:,:,:),aout_t(:,:,:,:), &
                     an1(:,:,:),dn1(:,:,:),un1(:,:,:),vn1(:,:,:), &
                     bin_t(:,:,:,:),fout_t(:,:,:,:),tmtemp(:,:), &
                     ain_t(:)
         logical, optional :: rhs_list(nrhs)
         logical :: rhslist(nrhs)
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif

         noffi=0
         do i=1,number_spheres
            nodri=mie_order(i)
            nblki=2*nodri*(nodri+2)
            b11=noffi+1
            if(tm_on_file(i)) then
               allocate(tmtemp(nblki,nblki),ain_t(nblki*nrhs))
               ain_t=ain(noffi+1:noffi+nblki*nrhs)
               b21=1
               n1=tm_offset(i)+1
               n2=tm_offset(i)+nblki*nblki
               tmtemp=reshape(stored_tm(n1:n2),(/nblki,nblki/))
               do j=1,nrhs
                  b12=b11+nblki-1
                  b22=b21+nblki-1
                  if(rhslist(j).and.idir.eq.1) then
                     aout(b11:b12)=matmul(tmtemp,ain_t(b21:b22))
                  elseif(rhslist(j).and.idir.ne.1) then
                     aout(b11:b12)=matmul(ain_t(b21:b22),tmtemp)
                  endif
                  b11=b11+nblki
                  b21=b21+nblki
               enddo
               deallocate(tmtemp,ain_t)
            else
               allocate(gin_t(0:nodri+1,nodri,2,nrhs),aout_t(0:nodri+1,nodri,2,nrhs), &
                        an1(2,2,nodri))
               b12=b11+nblki*nrhs-1
               gin_t=reshape(ain(b11:b12),(/nodri+2,nodri,2,nrhs/))
               nterms=4*nodri
               n1=mie_offset(i)+1
               n2=mie_offset(i)+nterms
               an1=reshape(an_mie(n1:n2),(/2,2,nodri/))
               if(mie_number_field_expansions(i).eq.1) then
                  aout_t=0.d0
                  do j=1,nrhs
                     if(rhslist(j).and.idir.eq.1) then
                        do n=1,nodri
                           do p=1,2
                              do q=1,2
                                 aout_t(n+1,n:1:-1,p,j)=aout_t(n+1,n:1:-1,p,j) &
                                     +an1(p,q,n)*gin_t(n+1,n:1:-1,q,j)
                                 aout_t(0:n,n,p,j)=aout_t(0:n,n,p,j)&
                                     +an1(p,q,n)*gin_t(0:n,n,q,j)
                              enddo
                           enddo
                        enddo
                     elseif(rhslist(j).and.idir.ne.1) then
                        do n=1,nodri
                           do p=1,2
                              do q=1,2
                                 aout_t(n+1,n:1:-1,p,j)=aout_t(n+1,n:1:-1,p,j) &
                                     +an1(q,p,n)*gin_t(n+1,n:1:-1,q,j)
                                 aout_t(0:n,n,p,j)=aout_t(0:n,n,p,j)&
                                     +an1(q,p,n)*gin_t(0:n,n,q,j)
                              enddo
                           enddo
                        enddo
                     endif
                  enddo
                  aout(b11:b12)&
                     =reshape(aout_t(0:nodri+1,1:nodri,1:2,1:nrhs),(/nblki*nrhs/))
               else
                  b21=noffi+nblki*nrhs+1
                  b22=b21+nblki*nrhs-1
                  allocate(bin_t(0:nodri+1,nodri,2,1:nrhs),&
                      fout_t(0:nodri+1,nodri,2,1:nrhs), &
                      dn1(2,2,nodri),un1(2,2,nodri),vn1(2,2,nodri))
                  bin_t=reshape(ain(b21:b22),(/nodri+2,nodri,2,nrhs/))
                  dn1=reshape(dn_mie(n1:n2),(/2,2,nodri/))
                  un1=reshape(un_mie(n1:n2),(/2,2,nodri/))
                  vn1=reshape(vn_mie(n1:n2),(/2,2,nodri/))
                  aout_t=0.d0
                  fout_t=0.d0
                  do j=1,nrhs
                     if(rhslist(j).and.idir.eq.1) then
                        do n=1,nodri
                           do p=1,2
                              do q=1,2
                                 aout_t(n+1,n:1:-1,p,j)=aout_t(n+1,n:1:-1,p,j) &
                                   + an1(p,q,n)*gin_t(n+1,n:1:-1,q,j) &
                                   + un1(p,q,n)*bin_t(n+1,n:1:-1,q,j)
                                 aout_t(0:n,n,p,j)=aout_t(0:n,n,p,j) &
                                   + an1(p,q,n)*gin_t(0:n,n,q,j) &
                                   + un1(p,q,n)*bin_t(0:n,n,q,j)
                                 fout_t(n+1,n:1:-1,p,j)=fout_t(n+1,n:1:-1,p,j) &
                                   + dn1(p,q,n)*gin_t(n+1,n:1:-1,q,j) &
                                   + vn1(p,q,n)*bin_t(n+1,n:1:-1,q,j)
                                 fout_t(0:n,n,p,j)=fout_t(0:n,n,p,j) &
                                   + dn1(p,q,n)*gin_t(0:n,n,q,j) &
                                   + vn1(p,q,n)*bin_t(0:n,n,q,j)
                              enddo
                           enddo
                        enddo
                     elseif(rhslist(j).and.idir.ne.1) then
                        do n=1,nodri
                           do p=1,2
                              do q=1,2
                                 aout_t(n+1,n:1:-1,p,j)=aout_t(n+1,n:1:-1,p,j) &
                                   + an1(q,p,n)*gin_t(n+1,n:1:-1,q,j) &
                                   + dn1(q,p,n)*bin_t(n+1,n:1:-1,q,j)
                                 aout_t(0:n,n,p,j)=aout_t(0:n,n,p,j) &
                                   + an1(q,p,n)*gin_t(0:n,n,q,j) &
                                   + dn1(q,p,n)*bin_t(0:n,n,q,j)
                                 fout_t(n+1,n:1:-1,p,j)=fout_t(n+1,n:1:-1,p,j) &
                                   + un1(q,p,n)*gin_t(n+1,n:1:-1,q,j) &
                                   + vn1(q,p,n)*bin_t(n+1,n:1:-1,q,j)
                                 fout_t(0:n,n,p,j)=fout_t(0:n,n,p,j) &
                                   + un1(q,p,n)*gin_t(0:n,n,q,j) &
                                   + vn1(q,p,n)*bin_t(0:n,n,q,j)

                              enddo
                           enddo
                        enddo
                     endif
                  enddo
                  aout(b11:b12)&
                     =reshape(aout_t(0:nodri+1,1:nodri,1:2,1:nrhs),(/nblki*nrhs/))
                  aout(b21:b22)&
                     =reshape(fout_t(0:nodri+1,1:nodri,1:2,1:nrhs),(/nblki*nrhs/))
                  deallocate(bin_t,fout_t,un1,vn1,dn1)
               endif
               deallocate(gin_t,aout_t,an1)
            endif
            noffi=noffi+nblki*nrhs*mie_number_field_expansions(i)
         enddo
         end subroutine multmiecoeffmult
!
! vector product for each rhs element of coefficient arrary.
! february 2013
!
         subroutine dotproduct(neqns,nrhs,gnp,err,rhs_list)
         implicit none
         integer :: neqns,nrhs,i,nodri,nblki,b11,b12, &
                    j,k
         real(8) :: err(nrhs)
         complex(8) :: gnp(neqns*nrhs)
         logical :: rhslist(nrhs)
         logical, optional :: rhs_list(nrhs)
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         b11=1
         err=0.d0
         do i=1,number_spheres
            nodri=mie_order(i)
            nblki=2*nodri*(nodri+2)
            do k=1,mie_number_field_expansions(i)
               do j=1,nrhs
                  if(rhslist(j)) then
                     b12=b11+nblki-1
                     err(j)=err(j)+sum((gnp(b11:b12)*conjg(gnp(b11:b12))))
                  endif
                  b11=b11+nblki
               enddo
            enddo
         enddo
         end subroutine dotproduct

      end module miecoefdata
!
!  module translation contains subroutines for VSWF translation and rotation
!
!
!  last revised: 15 January 2011
!  february 2013: generalized sphere configurations
!  march 2013: new mpi group convention
!
      module translation
      implicit none
      integer, private :: number_spheres,stored_max_order,store_tran_mat, &
                          number_host_spheres,number_interior_spheres, &
                          number_a_eqns,number_eqns,max_level
      integer, allocatable, private :: host_sphere(:),sphere_order(:), &
                 sphere_block(:), sphere_offset(:), interior_list(:), &
                 interior_number(:), interior_list_index(:),translation_order(:), &
                 number_field_expansions(:),number_sets(:),sphere_level(:), &
                 reordered_sphere_list(:),number_in_set(:,:),sphere_list_index(:,:)
      logical :: fftran_present
      logical, allocatable, private :: external_pair(:,:),ff_translation(:,:)
      real(8), private :: near_field_distance
      real(8), allocatable, private :: sphere_position(:,:)
      complex(8), allocatable, private :: ref_index(:,:)
      real(8), allocatable, private :: rotmat_store(:)
      complex(8), allocatable, private :: tranmat_store(:), ephimat_store(:)

      contains
!
! calculates lists for identifying host and interior spheres
!
! december 2011
!
         subroutine hostconfiguration(nsphere,hostsphere,numberfieldexp,maximum_level)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: nsphere,i,j,hostsphere(nsphere),k,numberfieldexp(nsphere),n,nset
         integer, optional :: maximum_level
         if(allocated(host_sphere)) deallocate(interior_list,interior_number, &
           interior_list_index,host_sphere,external_pair,number_field_expansions)
         allocate(host_sphere(nsphere),interior_list(nsphere),interior_number(nsphere), &
              interior_list_index(nsphere),external_pair(nsphere,nsphere), &
              number_field_expansions(nsphere))
         host_sphere=hostsphere
         k=1
         number_host_spheres=0
         do i=1,nsphere
            interior_number(i)=0
            interior_list_index(i)=k
            do j=1,nsphere
               if(host_sphere(j).eq.i) then
                  interior_number(i)=interior_number(i)+1
                  interior_list(k)=j
                  k=k+1
               endif
               if(host_sphere(i).eq.host_sphere(j).and.i.ne.j) then
                  external_pair(i,j)=.true.
               else
                  external_pair(i,j)=.false.
               endif
            enddo
            if(interior_number(i).gt.0) then
               number_host_spheres=number_host_spheres+1
            endif
         enddo
         do i=1,nsphere
            if(interior_number(i).eq.0) then
               number_field_expansions(i)=1
            else
               number_field_expansions(i)=2
            endif
         enddo
         numberfieldexp=number_field_expansions

         if(allocated(sphere_level)) deallocate(sphere_level)
         allocate(sphere_level(nsphere))
         max_level=0
         do i=1,nsphere
            j=i
            sphere_level(i)=0
            do while(host_sphere(j).ne.0)
               j=host_sphere(j)
               sphere_level(i)=sphere_level(i)+1
            enddo
            max_level=max(max_level,sphere_level(i))
         enddo
         if(present(maximum_level)) then
            maximum_level=max_level
         endif
         if(allocated(number_sets)) deallocate(number_sets,reordered_sphere_list, &
                  sphere_list_index,number_in_set)
         allocate(number_sets(0:max_level),reordered_sphere_list(nsphere), &
                  sphere_list_index(0:max_level,nsphere+1),number_in_set(0:max_level,nsphere))
         number_sets(0)=1
         k=1
         sphere_list_index(0,1)=1
         do i=1,nsphere
            if(hostsphere(i).eq.0) then
               reordered_sphere_list(k)=i
               k=k+1
            endif
         enddo
         sphere_list_index(0,2)=k
         number_in_set(0,1)=k-1
         do n=1,max_level
            number_sets(n)=0
            nset=1
            sphere_list_index(n,nset)=k
            do i=1,nsphere
               if((sphere_level(i).eq.n-1).and.(interior_number(i).ne.0)) then
                  do j=1,nsphere
                     if(hostsphere(j).eq.i) then
                        reordered_sphere_list(k)=j
                        k=k+1
                     endif
                  enddo
                  sphere_list_index(n,nset+1)=k
                  number_in_set(n,nset)=sphere_list_index(n,nset+1)-sphere_list_index(n,nset)
                  nset=nset+1
                  number_sets(n)=number_sets(n)+1
               endif
            enddo
         enddo
         end subroutine hostconfiguration
!
!  sets up the stored translation matrices and sets other constants used for translation.
!
!  last revised: 15 January 2011
!  november 2011: added near and far field translation
!  April 2012: all spheres are now assumed OA: l/r formulation
!  march 2013: mpicomm argument added
!
         subroutine rottranmtrxsetup(nsphere,nodr,rpos,ri,istore,&
                    nfdistance,ntran,runprintunit,fftranpresent,mpicomm)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: nsphere,nodr(nsphere),i,j,nodrmax, &
                    rank,runprintunit,isendok,tag,numprocs, &
                    nsend,istore,ntran(nsphere),rsize,tsize, &
                    esize,task,proc,noi,noj,nmin,nmax,mpicomm
         logical :: fftranpresent
         real(8) :: rpos(3,nsphere),xij(3),r,memused(1),memusedmax(1),memusedmin(1), &
                    nfdistance
         complex(8) :: ri(2,0:nsphere)
         data isendok,tag/0,1/
         if(mpicomm.eq.mpi_comm_null) return
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nodrmax=maxval(nodr)
         call init(nodrmax)
         number_spheres=nsphere
         store_tran_mat=istore
         near_field_distance=nfdistance
         if(allocated(sphere_position)) deallocate(sphere_position)
         if(allocated(ref_index)) deallocate(ref_index, &
            sphere_order,sphere_block,sphere_offset,translation_order,ff_translation)
         allocate(sphere_position(3,nsphere),ref_index(2,0:nsphere), &
                  sphere_order(nsphere),sphere_block(nsphere),sphere_offset(nsphere), &
                  translation_order(nsphere),ff_translation(nsphere,nsphere))
         sphere_position(:,1:nsphere)=rpos(:,1:nsphere)
         ref_index(:,0:nsphere)=ri(:,0:nsphere)
         translation_order=ntran
         ff_translation=.false.
         sphere_order(1)=nodr(1)
         sphere_block(1)=2*(nodr(1)*(nodr(1)+2))
         sphere_offset(1)=0
         do i=2,nsphere
            sphere_order(i)=nodr(i)
            sphere_block(i)=2*(nodr(i)*(nodr(i)+2))
            sphere_offset(i)=sphere_offset(i-1)  &
               +sphere_block(i-1)*number_field_expansions(i-1)
         enddo
         number_a_eqns=sum(sphere_block)
         number_eqns=sum(sphere_block*number_field_expansions)

         fftranpresent=.false.
         task=0
         rsize=0
         tsize=0
         esize=0
         do i=1,nsphere-1
            do j=i+1,nsphere
               if(host_sphere(j).eq.host_sphere(i)) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     noi=sphere_order(i)
                     noj=sphere_order(j)
                     xij(:)=sphere_position(:,i)-sphere_position(:,j)
                     r=sqrt(dot_product(xij,xij))
                     if(near_field_distance.lt.0.) then
                        nfdistance=(.5*(noi+noj))**2.
                     else
                        nfdistance=near_field_distance
                     endif
                     if(r.le.nfdistance) then
                        nmin=min(noi,noj)
                        nmax=max(noi,noj)
                        rsize=rsize+(2*nmin+1)*(1+nmax*(nmax+2))
                        tsize=tsize+atcdim(noi,noj)
                        esize=esize+2*nmax+1
                     else
                        fftranpresent=.true.
                        ff_translation(i,j)=.true.
                        ff_translation(j,i)=.true.
                     endif
                  endif
               endif
            enddo
         enddo

         fftran_present=fftranpresent
         if(store_tran_mat.eq.0) return

         stored_max_order=nodrmax
         memused(1)=dble(8*rsize+16*(tsize+esize))*1.d-6
         nsend=1
         call mstm_mpi(mpi_command='reduce',mpi_send_buf_dp=memused,mpi_recv_buf_dp=memusedmax,&
                     mpi_number=1,mpi_rank=0,mpi_operation=mstm_mpi_max,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='reduce',mpi_send_buf_dp=memused,mpi_recv_buf_dp=memusedmin,&
                     mpi_number=1,mpi_rank=0,mpi_operation=mstm_mpi_min,mpi_comm=mpicomm)
         if(rank.eq.0) then
            write(runprintunit,'('' maximum translation matrix storage:'',f9.4,'' MB'')') memusedmax
            write(runprintunit,'('' minimum translation matrix storage:'',f9.4,'' MB'')') memusedmin
            call flush(runprintunit)
         endif
!
!  calculate the matrices and store in memory
!
         if(allocated(rotmat_store)) deallocate(rotmat_store,tranmat_store,ephimat_store)
         allocate(rotmat_store(rsize))
         allocate(tranmat_store(tsize))
         allocate(ephimat_store(esize))
         end subroutine rottranmtrxsetup
!
!  clear the stored translation matrices
!
!  last revised: 15 January 2011
!
         subroutine rottranmtrxclear()
         implicit none
         if(allocated(rotmat_store)) deallocate(rotmat_store,tranmat_store,ephimat_store)
         if(allocated(sphere_position)) deallocate(sphere_position)
         end subroutine rottranmtrxclear

! The general sphere interaction driver
! LR formulation
! february 2013: number of rhs now a required argument list, and mpi comm is an option
!
         subroutine sphereinteraction(neqns,nrhs,ain,aout,ain_ct,aout_ct,initial_run, &
                    rhs_list,mpi_comm)
         use miecoefdata
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: neqns,rank,numprocs,nsphere,nrhs,npi1,npi2,nsend, &
                    mpicomm
         integer, optional :: mpi_comm
         logical :: twoop,initrun,rhslist(nrhs)
         logical, optional :: initial_run,rhs_list(nrhs)
         complex(8), optional :: ain_ct(neqns*nrhs),aout_ct(neqns*nrhs)
         complex(8)  :: ain_ct_t(neqns*nrhs),aout_t(neqns*nrhs), &
                        ain_t(neqns*nrhs),aout_ct_t(neqns*nrhs), &
                        ain(neqns*nrhs),aout(neqns*nrhs), &
                        aout_t2(neqns*nrhs),aout_ct_t2(neqns*nrhs)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         if(present(ain_ct)) then
            twoop=.true.
         else
            twoop=.false.
         endif
         if(present(initial_run)) then
            initrun=initial_run
         else
            initrun=.false.
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
!
! sphere-to-sphere H (i-j) interaction
!
         ain_t=ain
         aout=0.d0
         aout_t=0.
         if(twoop) then
            ain_ct_t=conjg(ain_ct)
            aout_ct=0.
            aout_ct_t=0.d0
            call multmiecoeffmult(neqns,nrhs,-1,ain_ct_t,ain_ct_t, &
                 rhs_list=rhslist)
         endif

         if(twoop) then
            call external_to_external_expansion(neqns,nrhs,ain_t,aout_t, &
                 ain_ct_t,aout_ct_t,initial_run=initrun, &
                 rhs_list=rhslist,mpi_comm=mpicomm)
         else
            call external_to_external_expansion(neqns,nrhs,ain_t,aout_t, &
                 initial_run=initrun,rhs_list=rhslist,mpi_comm=mpicomm)
         endif

         if(number_host_spheres.gt.0) then
            aout_t2=0.
            if(twoop) then
               aout_ct_t2=0.
               call external_to_internal_expansion(neqns,nrhs,ain_t,aout_t2, &
                    ain_ct_t,aout_ct_t2,rhs_list=rhslist,mpi_comm=mpicomm)
               aout_t=aout_t+aout_t2
               aout_ct_t=aout_ct_t+aout_ct_t2
            else
               call external_to_internal_expansion(neqns,nrhs,ain_t,aout_t2, &
                    rhs_list=rhslist,mpi_comm=mpicomm)
               aout_t=aout_t+aout_t2
            endif
         endif

         call multmiecoeffmult(neqns,nrhs,1,aout_t,aout_t,rhs_list=rhslist)
         aout=aout_t
         if(twoop) then
            aout_ct=conjg(aout_ct_t)
         endif
         if(numprocs.gt.1) then
            npi1=1
            nsend=neqns*nrhs
            npi2=nsend
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=aout(npi1:npi2), &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
            if(twoop) then
               call mstm_mpi(mpi_command='allreduce', &
                    mpi_recv_buf_dc=aout_ct(npi1:npi2), &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
            endif
         endif
         end subroutine sphereinteraction
!
! outgoing translation operation:  a(i) = H(i-j) a(j).
! February 2013: number of rhs is a required argument.   mpi comm option added.
! this does not perform an allgather on output arrays.  that operation will be needed
! to use the results
!
         subroutine external_to_external_expansion(neqns,nrhs,ain,gout,ain_ct,gout_ct, &
                    far_field_option,store_matrix_option,initial_run,rhs_list, &
                    mpi_comm)
         use miecoefdata
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: neqns,rank,numprocs,nsphere,nrhs,mpicomm, &
                    i,j,npi1,npi2,npj1,npj2,noj,noi,task,proc,noff, &
                    toffset,roffset,eoffset,nmin,nmax,tdim,rdim,edim
         integer, optional :: mpi_comm
         integer, allocatable :: noffi(:)
         logical :: twoop,ffopt,firstrun,smopt,rhslist(nrhs)
         logical, optional :: far_field_option,store_matrix_option,initial_run, &
                              rhs_list(nrhs)
         real(8) :: xij(3)
         real(8), allocatable :: rotmat(:)
         complex(8), optional :: ain_ct(neqns*nrhs),gout_ct(neqns*nrhs)
         complex(8)  :: ain(neqns*nrhs),gout(neqns*nrhs),rimedium(2)
         complex(8), allocatable :: tranmat(:),ephimat(:)
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         allocate(noffi(nsphere))
         gout=0.
         if(present(far_field_option)) then
            ffopt=far_field_option
         else
            ffopt=.true.
         endif
         if(present(store_matrix_option)) then
            smopt=store_matrix_option
         else
            smopt=.true.
         endif
         if(present(initial_run)) then
            firstrun=initial_run
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         if(present(ain_ct)) then
            twoop=.true.
            gout_ct=0.
         else
            twoop=.false.
         endif
!
!  compute offsets for scattering coefficients
!
         noff=0
         do i=1,nsphere
            noffi(i)=noff
            noff=noff+sphere_block(i) &
                 *nrhs*number_field_expansions(i)
         enddo
         roffset=0
         toffset=0
         eoffset=0
         task=0
         do i=1,nsphere-1
            do j=i+1,nsphere
               if(host_sphere(j).eq.host_sphere(i)) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     rimedium=ref_index(:,host_sphere(i))
                     noi=sphere_order(i)
                     npi1=noffi(i)+1
                     npi2=noffi(i)+sphere_block(i)*nrhs
                     noj=sphere_order(j)
                     npj1=noffi(j)+1
                     npj2=noffi(j)+sphere_block(j)*nrhs
                     xij(:)=sphere_position(:,i)-sphere_position(:,j)
                     if(ff_translation(i,j).and.ffopt) then
                        if(twoop) then
                           call rottranfarfield(nodrj=noj,nodri=noi,translation_vector=xij, &
                             refractive_index=rimedium,cxj=ain(npj1:npj2), &
                             cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),&
                             cyj=gout(npj1:npj2),cxj_t=ain_ct(npj1:npj2),&
                             cyi_t=gout_ct(npi1:npi2), cxi_t=ain_ct(npi1:npi2),&
                             cyj_t=gout_ct(npj1:npj2), &
                             number_rhs=nrhs,rhs_list=rhslist)
                        else
                           call rottranfarfield(nodrj=noj,nodri=noi,translation_vector=xij, &
                             refractive_index=rimedium,cxj=ain(npj1:npj2), &
                             cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),&
                             cyj=gout(npj1:npj2), &
                             number_rhs=nrhs,rhs_list=rhslist)
                        endif
                     else
                        if((store_tran_mat.eq.1).and.smopt) then
                           nmin=min(noi,noj)
                           nmax=max(noi,noj)
                           rdim=(2*nmin+1)*(1+nmax*(nmax+2))
                           tdim=atcdim(noi,noj)
                           edim=2*nmax+1
                           if(firstrun) then
                              allocate(rotmat(rdim),tranmat(tdim),ephimat(edim))
                              if(twoop) then
                                 call rottran(nodrj=noj,nodri=noi,translation_vector=xij, &
                                   refractive_index=rimedium,vswf_type=3,cxj=ain(npj1:npj2), &
                                   cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),&
                                   cyj=gout(npj1:npj2), &
                                   cxj_t=ain_ct(npj1:npj2),cyi_t=gout_ct(npi1:npi2), &
                                   cxi_t=ain_ct(npi1:npi2),cyj_t=gout_ct(npj1:npj2), &
                                   stored_rotation_matrix=rotmat, &
                                   stored_translation_matrix=tranmat, &
                                   stored_ephi_matrix=ephimat,number_rhs=nrhs,rhs_list=rhslist)
                              else
                                 call rottran(nodrj=noj,nodri=noi,translation_vector=xij, &
                                   refractive_index=rimedium,vswf_type=3,cxj=ain(npj1:npj2), &
                                   cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),&
                                   cyj=gout(npj1:npj2), &
                                   stored_rotation_matrix=rotmat,stored_translation_matrix=tranmat, &
                                   stored_ephi_matrix=ephimat,number_rhs=nrhs,rhs_list=rhslist)
                              endif
                              rotmat_store(roffset+1:roffset+rdim)=rotmat(1:rdim)
                              tranmat_store(toffset+1:toffset+tdim)=tranmat(1:tdim)
                              ephimat_store(eoffset+1:eoffset+edim)=ephimat(1:edim)
                              deallocate(rotmat,tranmat,ephimat)
                           else
                              if(twoop) then
                                 call rottran(nodrj=noj,nodri=noi,cxj=ain(npj1:npj2), &
                                   cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),&
                                   cyj=gout(npj1:npj2), &
                                   cxj_t=ain_ct(npj1:npj2),cyi_t=gout_ct(npi1:npi2), &
                                   cxi_t=ain_ct(npi1:npi2),cyj_t=gout_ct(npj1:npj2), &
                                   stored_rotation_matrix=rotmat_store(roffset+1:roffset+rdim), &
                                   stored_translation_matrix=tranmat_store(toffset+1:toffset+tdim), &
                                   stored_ephi_matrix=ephimat_store(eoffset+1:eoffset+edim),&
                                   number_rhs=nrhs,rhs_list=rhslist)
                              else
                                 call rottran(nodrj=noj,nodri=noi,cxj=ain(npj1:npj2), &
                                   cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),&
                                   cyj=gout(npj1:npj2), &
                                   stored_rotation_matrix=rotmat_store(roffset+1:roffset+rdim), &
                                   stored_translation_matrix=tranmat_store(toffset+1:toffset+tdim), &
                                   stored_ephi_matrix=ephimat_store(eoffset+1:eoffset+edim),&
                                   number_rhs=nrhs,rhs_list=rhslist)
                              endif
                           endif
                           roffset=roffset+rdim
                           toffset=toffset+tdim
                           eoffset=eoffset+edim
                        else
                           if(twoop) then
                              call rottran(nodrj=noj,nodri=noi,translation_vector=xij, &
                                refractive_index=rimedium,vswf_type=3,cxj=ain(npj1:npj2), &
                                cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),cyj=gout(npj1:npj2), &
                                cxj_t=ain_ct(npj1:npj2),cyi_t=gout_ct(npi1:npi2), &
                                cxi_t=ain_ct(npi1:npi2),cyj_t=gout_ct(npj1:npj2), &
                                number_rhs=nrhs,rhs_list=rhslist)
                           else
                              call rottran(nodrj=noj,nodri=noi,translation_vector=xij, &
                                refractive_index=rimedium,vswf_type=3,cxj=ain(npj1:npj2), &
                                cyi=gout(npi1:npi2),cxi=ain(npi1:npi2),cyj=gout(npj1:npj2), &
                                number_rhs=nrhs,rhs_list=rhslist)
                           endif
                        endif
                     endif
                  endif
               endif
            enddo
         enddo

         if(store_tran_mat.eq.1.and.firstrun) then
            firstrun=.false.
         endif

         end subroutine external_to_external_expansion
!
!  calculation of bmnp(i) = J(i-j) amnp(j) for i internal, host j=i, and
!  gmnp(i) = J(i-j) f(j), for host i = j.    This is the regular translation operation.   g and b
!  are returned ordered as (a,f) in the output array.
!  february 2013: number of rhs option
!  the routine does not perform an mpi reduce.
!
         subroutine external_to_internal_expansion(neqns,nrhs,ain,bout,ain_ct,&
                    bout_ct,rhs_list,mpi_comm)
         use miecoefdata
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: neqns,rank,numprocs,nsphere,nrhs,noff,mpicomm, &
                    i,j,task,proc,extsurf,intsurf,np1ext1,np1ext2,&
                    np2ext1,np2ext2,np1int1,np1int2,np2int1,np2int2,noext,noint
         integer, optional :: mpi_comm
         integer, allocatable :: noffi(:)
         logical :: twoop,rhslist(nrhs)
         logical, optional :: rhs_list(nrhs)
         real(8) :: xij(3)
         complex(8), optional :: ain_ct(neqns*nrhs),bout_ct(neqns*nrhs)
         complex(8)  :: ain(neqns*nrhs),bout(neqns*nrhs),rimedium(2)
         integer :: count
         data count/0/
         count=count+1
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         allocate(noffi(nsphere))
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         bout=0.d0
         if(present(ain_ct)) then
            twoop=.true.
            bout_ct=0.
         else
            twoop=.false.
         endif
         task=0
         noff=0
         do i=1,nsphere
            noffi(i)=noff
            noff=noff+sphere_block(i) &
                 *nrhs*number_field_expansions(i)
         enddo

         do i=1,nsphere-1
            do j=i+1,nsphere
               if(host_sphere(j).eq.i.or.host_sphere(i).eq.j) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     if(host_sphere(j).eq.i) then
                        extsurf=j
                        intsurf=i
                     else
                        extsurf=i
                        intsurf=j
                     endif
                     noext=sphere_order(extsurf)
                     noint=sphere_order(intsurf)
                     xij(:)=sphere_position(:,intsurf)-sphere_position(:,extsurf)
                     rimedium=ref_index(:,intsurf)
                     np1ext1=noffi(extsurf)+1
                     np1ext2=np1ext1-1+sphere_block(extsurf)*nrhs
                     np1int1=noffi(intsurf)+1+sphere_block(intsurf)*nrhs
                     np1int2=np1int1-1+sphere_block(intsurf)*nrhs
                     np2int1=noffi(intsurf)+1+sphere_block(intsurf)*nrhs
                     np2int2=np2int1-1+sphere_block(intsurf)*nrhs
                     np2ext1=noffi(extsurf)+1
                     np2ext2=np2ext1-1+sphere_block(extsurf)*nrhs
                     if(twoop) then
                        call rottran(nodrj=noext,nodri=noint,translation_vector=xij, &
                          refractive_index=rimedium,vswf_type=1,cxj=ain(np1ext1:np1ext2), &
                          cyi=bout(np2int1:np2int2),cxi=ain(np1int1:np1int2), &
                          cyj=bout(np2ext1:np2ext2),cxj_t=ain_ct(np1ext1:np1ext2), &
                          cyi_t=bout_ct(np2int1:np2int2),cxi_t=ain_ct(np1int1:np1int2), &
                          cyj_t=bout_ct(np2ext1:np2ext2),number_rhs=nrhs,rhs_list=rhslist)
                     else
                        call rottran(nodrj=noext,nodri=noint,translation_vector=xij, &
                          refractive_index=rimedium,vswf_type=1,cxj=ain(np1ext1:np1ext2), &
                          cyi=bout(np2int1:np2int2),cxi=ain(np1int1:np1int2), &
                          cyj=bout(np2ext1:np2ext2),number_rhs=nrhs,rhs_list=rhslist)
                     endif
                  endif
               endif
            enddo
         enddo
         end subroutine external_to_internal_expansion
!
! sign flipped for odd degrees
!
         subroutine m1_to_the_n(nodr,cx)
         implicit none
         integer :: nodr,n
         complex(8) :: cx(0:nodr+1,1:nodr,1:2)
         do n=1,nodr,2
            cx(0:n,n,1:2)=-cx(0:n,n,1:2)
            cx(n+1,1:n,1:2)=-cx(n+1,1:n,1:2)
         enddo
         end subroutine m1_to_the_n
!
!  far field formula for outgoing VSWF translation
!  October 2011
!  April 2011: lr formulation.   summation over mode is eliminated.   Definition of pi
!   is now based on lr model
!  february 2013: number rhs option
!
         subroutine rottranfarfield(nodrj,nodri,translation_vector,refractive_index, &
              cxj,cyi,cxi,cyj,cxj_t,cyi_t,cxi_t,cyj_t,number_rhs,rhs_list)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrj,nodri,nmax,nrhs,i
         integer, optional :: number_rhs
         logical :: twoop,tranop
         logical, optional :: rhs_list(*)
         logical, allocatable :: rhslist(:)
         real(8) :: translation_vector(3),r,ct
         complex(8) :: cxj(0:nodrj+1,nodrj,2,*),cyi(0:nodri+1,nodri,2,*),ri(2), &
                       ephi,phasefunc(2), sumx(2), &
                       pivec(0:max(nodri,nodrj)+1,max(nodri,nodrj),2), &
                       ctemp(0:max(nodri,nodrj)+1,max(nodri,nodrj),2)
         complex(8), optional :: refractive_index(2)
         complex(8), optional :: cxi(0:nodri+1,nodri,2,*),cyj(0:nodrj+1,nodrj,2,*), &
                     cxj_t(0:nodrj+1,nodrj,2,*),cyi_t(0:nodri+1,nodri,2,*), &
                     cxi_t(0:nodri+1,nodri,2,*),cyj_t(0:nodrj+1,nodrj,2,*)

         if(present(cxi)) then
            twoop=.true.
         else
            twoop=.false.
         endif
         if(present(cxj_t)) then
            tranop=.true.
         else
            tranop=.false.
         endif
         if(present(refractive_index)) then
            ri=refractive_index
         else
            ri=(/(1.d0,0.d0),(1.d0,0.d0)/)
         endif
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         if(allocated(rhslist)) deallocate(rhslist)
         allocate(rhslist(nrhs))
         if(present(rhs_list)) then
            rhslist=rhs_list(1:nrhs)
         else
            rhslist=.true.
         endif

         call cartosphere(translation_vector,r,ct,ephi)
         nmax=max(nodrj,nodri)
         call pifunc(ct,ephi,nmax,nmax,pivec)
         phasefunc(1)=cdexp((0.d0,1.d0)*ri(1)*r)/((0.d0,1.d0)*ri(1)*r)
         phasefunc(2)=cdexp((0.d0,1.d0)*ri(2)*r)/((0.d0,1.d0)*ri(2)*r)
!
! normal
!
         do i=1,nrhs
            if(rhslist(i)) then
               sumx(1)=sum(pivec(0:nodrj+1,1:nodrj,1)*cxj(0:nodrj+1,1:nodrj,1,i))
               sumx(2)=sum(pivec(0:nodrj+1,1:nodrj,2)*cxj(0:nodrj+1,1:nodrj,2,i))
               sumx=sumx*phasefunc
               cyi(0:nodri+1,1:nodri,1,i)=cyi(0:nodri+1,1:nodri,1,i) &
                  + sumx(1)*conjg(pivec(0:nodri+1,1:nodri,1))
               cyi(0:nodri+1,1:nodri,2,i)=cyi(0:nodri+1,1:nodri,2,i) &
                  + sumx(2)*conjg(pivec(0:nodri+1,1:nodri,2))
            endif
         enddo
! double

         if(twoop) then
            do i=1,nrhs
               if(rhslist(i)) then
                  ctemp(0:nodri+1,1:nodri,1:2)=cxi(0:nodri+1,1:nodri,1:2,i)
                  call m1_to_the_n(nodri,ctemp)
                  sumx(1)=sum(pivec(0:nodri+1,1:nodri,2)*ctemp(0:nodri+1,1:nodri,1))
                  sumx(2)=sum(pivec(0:nodri+1,1:nodri,1)*ctemp(0:nodri+1,1:nodri,2))
                  sumx=sumx*phasefunc
                  ctemp(0:nodrj+1,1:nodrj,1)=sumx(1)*conjg(pivec(0:nodrj+1,1:nodrj,2))
                  ctemp(0:nodrj+1,1:nodrj,2)=sumx(2)*conjg(pivec(0:nodrj+1,1:nodrj,1))
                  call m1_to_the_n(nodrj,ctemp)
                  cyj(0:nodrj+1,1:nodrj,1:2,i)=cyj(0:nodrj+1,1:nodrj,1:2,i) &
                     + ctemp(0:nodrj+1,1:nodrj,1:2)
               endif
            enddo
         endif
         if(tranop) then
            do i=1,nrhs
               if(rhslist(i)) then
                  ctemp(0:nodrj+1,1:nodrj,1:2)=cxj_t(0:nodrj+1,1:nodrj,1:2,i)
                  call m1_to_the_n(nodrj,ctemp)
                  sumx(1)=sum(conjg(pivec(0:nodrj+1,1:nodrj,2))*ctemp(0:nodrj+1,1:nodrj,1))
                  sumx(2)=sum(conjg(pivec(0:nodrj+1,1:nodrj,1))*ctemp(0:nodrj+1,1:nodrj,2))
                  sumx=sumx*phasefunc
                  ctemp(0:nodri+1,1:nodri,1)=sumx(1)*pivec(0:nodri+1,1:nodri,2)
                  ctemp(0:nodri+1,1:nodri,2)=sumx(2)*pivec(0:nodri+1,1:nodri,1)
                  call m1_to_the_n(nodri,ctemp)
                  cyi_t(0:nodri+1,1:nodri,1:2,i)=cyi_t(0:nodri+1,1:nodri,1:2,i) &
                     + ctemp(0:nodri+1,1:nodri,1:2)
               endif
            enddo
         endif
         if(twoop.and.tranop) then
            do i=1,nrhs
               if(rhslist(i)) then
                  sumx(1)=sum(conjg(pivec(0:nodri+1,1:nodri,1))*cxi_t(0:nodri+1,1:nodri,1,i))
                  sumx(2)=sum(conjg(pivec(0:nodri+1,1:nodri,2))*cxi_t(0:nodri+1,1:nodri,2,i))
                  sumx=sumx*phasefunc
                  cyj_t(0:nodrj+1,1:nodrj,1,i)=cyj_t(0:nodrj+1,1:nodrj,1,i) &
                     + sumx(1)*pivec(0:nodrj+1,1:nodrj,1)
                  cyj_t(0:nodrj+1,1:nodrj,2,i)=cyj_t(0:nodrj+1,1:nodrj,2,i) &
                     + sumx(2)*pivec(0:nodrj+1,1:nodrj,2)
               endif
            enddo
         endif
         end subroutine rottranfarfield
!
! correction term for hybrid bcgm solution: difference between exact and
! ff translation field
! november 2011
! april 2012: lr formulation (ri now vector length 2)
! february 2013: number rhs option
!
         subroutine farfieldtranslationerror(neqns,ain,gout,number_rhs,rhs_list, &
                    mpi_comm)
         use miecoefdata
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: neqns,rank,numprocs,nsphere,nrhs,noff,ntot,mpicomm, &
                    nsend,i,j,npi1,npi2,npj1,npj2,noj,noi,task,proc
         integer, optional :: number_rhs,mpi_comm
         integer, allocatable :: noffi(:)
         logical, allocatable :: rhslist(:)
         logical, optional :: rhs_list(*)
         real(8) :: xij(3),rij,nfdist
         complex(8)  :: ain(*),gout(*),rimedium(2)
         complex(8), allocatable :: goutff(:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         ntot=neqns*nrhs
         allocate(rhslist(nrhs),noffi(nsphere))
         if(present(rhs_list)) then
            rhslist=rhs_list(1:nrhs)
         else
            rhslist=.true.
         endif
         noff=0
         do i=1,nsphere
            noffi(i)=noff
            noff=noff+sphere_block(i) &
                 *nrhs*number_field_expansions(i)
         enddo
         allocate(goutff(ntot))
         gout(1:ntot)=0.
         goutff(1:ntot)=0.
         task=0
         do i=1,nsphere-1
            do j=i+1,nsphere
               if(host_sphere(j).eq.host_sphere(i)) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     rimedium=ref_index(:,host_sphere(i))
                     noi=sphere_order(i)
                     npi1=noffi(i)+1
                     npi2=noffi(i)+sphere_block(i)*nrhs
                     noj=sphere_order(j)
                     npj1=noffi(j)+1
                     npj2=noffi(j)+sphere_block(j)*nrhs
                     xij(:)=sphere_position(:,i)-sphere_position(:,j)
                     rij=sqrt(dot_product(xij,xij))
                     if(near_field_distance.lt.0.) then
                        nfdist=(.5*(noi+noj))**2.
                     else
                        nfdist=near_field_distance
                     endif
                     if(rij.gt.nfdist) then
                        call rottran(nodrj=noj,nodri=noi,translation_vector=xij, &
                          refractive_index=rimedium,vswf_type=3,cxj=ain(npj1:npj2), &
                          cyi=gout(npi1:npi2),cxi=ain(npi1:npi2), &
                          cyj=gout(npj1:npj2), &
                          number_rhs=nrhs,rhs_list=rhslist)
                        call rottranfarfield(nodrj=noj,nodri=noi,translation_vector=xij, &
                          refractive_index=rimedium,cxj=ain(npj1:npj2), &
                          cyi=goutff(npi1:npi2),cxi=ain(npi1:npi2), &
                          cyj=goutff(npj1:npj2),&
                          number_rhs=nrhs,rhs_list=rhslist)
                     endif
                  endif
               endif
            enddo
         enddo
         gout(1:ntot)=gout(1:ntot)-goutff(1:ntot)

         if(numprocs.gt.1) then
            npi1=1
            nsend=neqns*nrhs
            npi2=nsend
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=gout(npi1:npi2), &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
         deallocate(goutff)
         end subroutine farfieldtranslationerror
!
!  the vectorized rotation-translation-rotation operation
!
!
!  last revised: 15 January 2011
!  april 2012: lr formulation: this involves removing the mode coupling for translation
!  may 2012: completely rewritten, combines previous rottran matrix operations.
!  february 2013: number rhs option:   performs translations for each rhs column of
!  input vectors.
!
         subroutine rottran(nodrj,nodri,translation_vector,refractive_index, &
              vswf_type,stored_rotation_matrix,stored_translation_matrix, &
              stored_ephi_matrix,cxj,cyi,cxi,cyj,cxj_t,cyi_t,cxi_t,cyj_t, &
              number_rhs,rhs_list)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrj,nodri,nmin,nmax,ndim,m,n,nn1,nn2,m1,p,n1,offset, &
                    blocksize,itype,nrhs,i
         integer, optional :: vswf_type,number_rhs
         logical :: twoop,tranop
         logical, allocatable :: rhslist(:)
         logical, optional :: rhs_list(*)
         real(8) :: r,ct
         real(8), allocatable :: rotmat(:,:)
         real(8), optional :: stored_rotation_matrix(*),translation_vector(3)
         complex(8) :: cxj(0:nodrj+1,nodrj,2,*),cyi(0:nodri+1,nodri,2,*),ri(2), &
                       atc(nodri,nodrj,2),ephi
         complex(8), optional :: stored_ephi_matrix(*), &
                                 stored_translation_matrix(*), &
                                 refractive_index(2)
         complex(8), allocatable :: cxjt(:,:,:,:),cxit(:,:,:,:),cxjt_t(:,:,:,:), &
                                    cxit_t(:,:,:,:),ephimat(:),tranmat(:)
         complex(8), optional :: cxi(0:nodri+1,nodri,2,*),cyj(0:nodrj+1,nodrj,2,*), &
                     cxj_t(0:nodrj+1,nodrj,2,*),cyi_t(0:nodri+1,nodri,2,*), &
                     cxi_t(0:nodri+1,nodri,2,*),cyj_t(0:nodrj+1,nodrj,2,*)

         if(present(cxi)) then
            twoop=.true.
         else
            twoop=.false.
         endif
         if(present(cxj_t)) then
            tranop=.true.
         else
            tranop=.false.
         endif
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         if(allocated(rhslist)) deallocate(rhslist)
         allocate(rhslist(nrhs))
         if(present(rhs_list)) then
            rhslist=rhs_list(1:nrhs)
         else
            rhslist=.true.
         endif
         if(present(vswf_type)) then
            itype=vswf_type
         else
            itype=1
         endif

         if(present(translation_vector)) then
            if(dot_product(translation_vector,translation_vector).lt.0.001 &
               .and.itype.eq.1) then
               do i=1,nrhs
                  if(rhslist(i)) then
                     call transfer(nodrj,nodri,cxj(0,1,1,i),cyi(0,1,1,i))
                     if(tranop) call transfer(nodrj,nodri,cxj_t(0,1,1,i),cyi_t(0,1,1,i))
                     if(twoop) call transfer(nodri,nodrj,cxi(0,1,1,i),cyj(0,1,1,i))
                     if(tranop.and.twoop) call transfer(nodri,nodrj,cxi_t(0,1,1,i), &
                        cyj_t(0,1,1,i))
                  endif
               enddo
               return
            endif
         endif

         nmin=min(nodri,nodrj)
         nmax=max(nodri,nodrj)
         allocate(rotmat(-nmin:nmin,0:nmax*(nmax+2)))
         allocate(ephimat(-nmax:nmax))
         ndim=atcdim(nodri,nodrj)
         allocate(tranmat(1:ndim))
         if(present(translation_vector)) then
            call cartosphere(translation_vector,r,ct,ephi)
            if(present(refractive_index)) then
               ri=refractive_index
            else
               ri=(/(1.d0,0.d0),(1.d0,0.d0)/)
            endif
            call rotcoef(ct,nmin,nmax,rotmat)
            call axialtrancoefrecurrence(itype,r,ri,nodri,nodrj,ndim,tranmat)
            call ephicoef(ephi,nmax,ephimat)
            if(present(stored_rotation_matrix)) then
               n=(2*nmin+1)*(1+nmax*(nmax+2))
               stored_rotation_matrix(1:n) &
                  =reshape(rotmat(-nmin:nmin,0:nmax*(nmax+2)),(/n/))
            endif
            if(present(stored_translation_matrix)) then
               stored_translation_matrix(1:ndim)=tranmat(1:ndim)
            endif
            if(present(stored_ephi_matrix)) then
               stored_ephi_matrix(1:2*nmax+1)=ephimat(-nmax:nmax)
            endif
         else
            n=(2*nmin+1)*(1+nmax*(nmax+2))
            rotmat(-nmin:nmin,0:nmax*(nmax+2)) &
              =reshape(stored_rotation_matrix(1:n),(/2*nmin+1,1+nmax*(nmax+2)/))
            tranmat(1:ndim)=stored_translation_matrix(1:ndim)
            ephimat(-nmax:nmax)=stored_ephi_matrix(1:2*nmax+1)
         endif
!
! tranfer into temporary arrays and shift, as needed
!
!  regular case
         allocate(cxjt(-nmax:nmax,nmax,2,nrhs))
         do i=1,nrhs
            if(rhslist(i)) then
               cxjt(0,1:nodrj,1:2,i)=cxj(0,1:nodrj,1:2,i)
               do m=1,nodrj
                  cxjt(m,m:nodrj,1:2,i)=cxj(m,m:nodrj,1:2,i)*ephimat(m)
                  cxjt(-m,m:nodrj,1:2,i)=cxj(m+1:nodrj+1,m,1:2,i)*ephimat(-m)
               enddo
            endif
         enddo
! double case
         if(twoop) then
            allocate(cxit(-nmax:nmax,nmax,2,nrhs))
            do i=1,nrhs
               if(rhslist(i)) then
                  cxit(0,1:nodri,1:2,i)=cxi(0,1:nodri,1:2,i)
                  do m=1,nodri
                     cxit(m,m:nodri,1:2,i)=monen(m) &
                         *cxi(m+1:nodri+1,m,1:2,i)*ephimat(-m)
                     cxit(-m,m:nodri,1:2,i)=monen(m) &
                         *cxi(m,m:nodri,1:2,i)*ephimat(m)
                  enddo
               endif
            enddo
         endif
! transpose case
         if(tranop) then
            allocate(cxjt_t(-nmax:nmax,nmax,2,nrhs))
            do i=1,nrhs
               if(rhslist(i)) then
                  cxjt_t(0,1:nodrj,1:2,i)=cxj_t(0,1:nodrj,1:2,i)
                  do m=1,nodrj
                     cxjt_t(m,m:nodrj,1:2,i)=monen(m) &
                         *cxj_t(m+1:nodrj+1,m,1:2,i)*ephimat(m)
                     cxjt_t(-m,m:nodrj,1:2,i)=monen(m) &
                         *cxj_t(m,m:nodrj,1:2,i)*ephimat(-m)
                  enddo
               endif
            enddo
         endif
! double transpose case
         if(twoop.and.tranop) then
            allocate(cxit_t(-nmax:nmax,nmax,2,nrhs))
            do i=1,nrhs
               if(rhslist(i)) then
                  cxit_t(0,1:nodri,1:2,i)=cxi_t(0,1:nodri,1:2,i)
                  do m=1,nodri
                     cxit_t(m,m:nodri,1:2,i)=cxi_t(m,m:nodri,1:2,i)*ephimat(-m)
                     cxit_t(-m,m:nodri,1:2,i)=cxi_t(m+1:nodri+1,m,1:2,i)*ephimat(m)
                  enddo
               endif
            enddo
         endif
!
!  rotation to origin of target
!
         do n=1,nodrj
            nn1=n*(n+1)-n
            nn2=nn1+(2*n+1)-1
            n1=min(n,nodri)
            do i=1,nrhs
               if(rhslist(i)) then
                  cxjt(-n1:n1,n,1:2,i)=matmul(rotmat(-n1:n1,nn1:nn2),cxjt(-n:n,n,1:2,i))
               endif
            enddo
            if(tranop) then
               do i=1,nrhs
                  if(rhslist(i)) then
                     cxjt_t(-n1:n1,n,1:2,i) &
                     = matmul(rotmat(-n1:n1,nn1:nn2),cxjt_t(-n:n,n,1:2,i))
                  endif
               enddo
            endif
         enddo
         if(twoop) then
            do n=1,nodri
               nn1=n*(n+1)-n
               nn2=nn1+(2*n+1)-1
               n1=min(n,nodrj)
               do i=1,nrhs
                  if(rhslist(i)) then
                     cxit(-n1:n1,n,1:2,i) &
                     =matmul(rotmat(-n1:n1,nn1:nn2),cxit(-n:n,n,1:2,i))
                  endif
               enddo
               if(tranop) then
                  do i=1,nrhs
                     if(rhslist(i)) then
                        cxit_t(-n1:n1,n,1:2,i)&
                        =matmul(rotmat(-n1:n1,nn1:nn2),cxit_t(-n:n,n,1:2,i))
                     endif
                  enddo
               endif
            enddo
         endif
!
!  axial translation to target
!
         do m=-nmin,nmin
            m1=max(1,abs(m))
            offset=moffset(m,nodri,nodrj)
            blocksize=(nodri-m1+1)*(nodrj-m1+1)*2
            atc(m1:nodri,m1:nodrj,1:2)= &
                reshape(tranmat(offset+1:offset+blocksize),(/nodri-m1+1,nodrj-m1+1,2/))
            do i=1,nrhs
               if(rhslist(i)) then
                  do p=1,2
                     cxjt(m,m1:nodri,p,i) &
                        =matmul(atc(m1:nodri,m1:nodrj,p),cxjt(m,m1:nodrj,p,i))
                  enddo
               endif
            enddo
            if(tranop) then
               do i=1,nrhs
                  if(rhslist(i)) then
                     do p=1,2
                        cxjt_t(m,m1:nodri,p,i) &
                           =matmul(atc(m1:nodri,m1:nodrj,p),cxjt_t(m,m1:nodrj,p,i))
                     enddo
                  endif
               enddo
            endif
            if(twoop) then
               do i=1,nrhs
                  if(rhslist(i)) then
                     do p=1,2
                        cxit(m,m1:nodrj,p,i) &
                           =matmul(cxit(m,m1:nodri,p,i),atc(m1:nodri,m1:nodrj,p))
                     enddo
                  endif
               enddo
            endif
            if(twoop.and.tranop) then
               do i=1,nrhs
                  if(rhslist(i)) then
                     do p=1,2
                        cxit_t(m,m1:nodrj,p,i) &
                           =matmul(cxit_t(m,m1:nodri,p,i),atc(m1:nodri,m1:nodrj,p))
                     enddo
                  endif
               enddo
            endif
         enddo
!
!  rotation back to original frame
!
         do n=1,nodri
            nn1=n*(n+1)-n
            nn2=nn1+(2*n+1)-1
            n1=min(n,nodrj)
            do i=1,nrhs
               if(rhslist(i)) then
                  cxjt(-n:n,n,1,i)=matmul(cxjt(-n1:n1,n,1,i),rotmat(-n1:n1,nn1:nn2))
                  cxjt(-n:n,n,2,i)=matmul(cxjt(-n1:n1,n,2,i),rotmat(-n1:n1,nn1:nn2))
               endif
            enddo
            if(tranop) then
               do i=1,nrhs
                  if(rhslist(i)) then
                     cxjt_t(-n:n,n,1,i) &
                        =matmul(cxjt_t(-n1:n1,n,1,i),rotmat(-n1:n1,nn1:nn2))
                     cxjt_t(-n:n,n,2,i) &
                        =matmul(cxjt_t(-n1:n1,n,2,i),rotmat(-n1:n1,nn1:nn2))
                  endif
               enddo
            endif
         enddo
         if(twoop) then
            do n=1,nodrj
               nn1=n*(n+1)-n
               nn2=nn1+(2*n+1)-1
               n1=min(n,nodri)
               do i=1,nrhs
                  if(rhslist(i)) then
                     cxit(-n:n,n,1,i)=matmul(cxit(-n1:n1,n,1,i),rotmat(-n1:n1,nn1:nn2))
                     cxit(-n:n,n,2,i)=matmul(cxit(-n1:n1,n,2,i),rotmat(-n1:n1,nn1:nn2))
                  endif
               enddo
               if(tranop) then
                  do i=1,nrhs
                     if(rhslist(i)) then
                        cxit_t(-n:n,n,1,i)=matmul(cxit_t(-n1:n1,n,1,i),rotmat(-n1:n1,nn1:nn2))
                        cxit_t(-n:n,n,2,i)=matmul(cxit_t(-n1:n1,n,2,i),rotmat(-n1:n1,nn1:nn2))
                     endif
                  enddo
               endif
            enddo
         endif
!
!  transfer into output arrays
!
! regular
         do i=1,nrhs
            if(rhslist(i)) then
               cyi(0,1:nodri,1:2,i)=cyi(0,1:nodri,1:2,i)+cxjt(0,1:nodri,1:2,i)
               do m=1,nodri
                  cyi(m,m:nodri,1:2,i)=cyi(m,m:nodri,1:2,i) &
                     +cxjt(m,m:nodri,1:2,i)*ephimat(-m)
                  cyi(m+1:nodri+1,m,1:2,i)=cyi(m+1:nodri+1,m,1:2,i) &
                     +cxjt(-m,m:nodri,1:2,i)*ephimat(m)
               enddo
            endif
         enddo
         deallocate(cxjt)
! double operation
         if(twoop) then
            do i=1,nrhs
               if(rhslist(i)) then
                  cyj(0,1:nodrj,1:2,i)=cyj(0,1:nodrj,1:2,i)+cxit(0,1:nodrj,1:2,i)
                  do m=1,nodrj
                     cyj(m,m:nodrj,1:2,i)=cyj(m,m:nodrj,1:2,i) &
                        +cxit(-m,m:nodrj,1:2,i)*monen(m)*ephimat(-m)
                     cyj(m+1:nodrj+1,m,1:2,i)=cyj(m+1:nodrj+1,m,1:2,i) &
                        +cxit(m,m:nodrj,1:2,i)*monen(m)*ephimat(m)
                  enddo
               endif
            enddo
            deallocate(cxit)
         endif
! regular transpose
         if(tranop) then
            do i=1,nrhs
               if(rhslist(i)) then
                  cyi_t(0,1:nodri,1:2,i)=cyi_t(0,1:nodri,1:2,i)+cxjt_t(0,1:nodri,1:2,i)
                  do m=1,nodri
                     cyi_t(m,m:nodri,1:2,i)=cyi_t(m,m:nodri,1:2,i)  &
                         +cxjt_t(-m,m:nodri,1:2,i)*monen(m)*ephimat(m)
                     cyi_t(m+1:nodri+1,m,1:2,i)=cyi_t(m+1:nodri+1,m,1:2,i) &
                         +cxjt_t(m,m:nodri,1:2,i)*monen(m)*ephimat(-m)
                  enddo
               endif
            enddo
            deallocate(cxjt_t)
         endif
! double transpose
         if(twoop.and.tranop) then
            do i=1,nrhs
               if(rhslist(i)) then
                  cyj_t(0,1:nodrj,1:2,i)=cyj_t(0,1:nodrj,1:2,i)+cxit_t(0,1:nodrj,1:2,i)
                  do m=1,nodrj
                     cyj_t(m,m:nodrj,1:2,i)=cyj_t(m,m:nodrj,1:2,i)  &
                         +cxit_t(m,m:nodrj,1:2,i)*ephimat(m)
                     cyj_t(m+1:nodrj+1,m,1:2,i)=cyj_t(m+1:nodrj+1,m,1:2,i) &
                         +cxit_t(-m,m:nodrj,1:2,i)*ephimat(-m)
                  enddo
               endif
            enddo
            deallocate(cxit_t)
         endif

         deallocate(rotmat,tranmat,ephimat)

         end subroutine rottran
!
!  GB coefficients for sphere-centered expansions, obtained via translation
!
!  last revised: 15 January 2011
!  april 2012: lr formulation
!  february 2013:  mpi comm option
!
         subroutine spheregaussianbeamcoef(nsphere,neqns,nodr,alpha,beta,cbeam, &
                    rpos,hostsphere,numberfieldexp,rimed,rbeam,epstran,pmnp, &
                    mpi_comm)
         use specialfuncs
         use mpidefs
         implicit none
         integer :: n,nsphere,i,nodr(nsphere),nblk,noff,nodrgb,neqns, &
                    hostsphere(nsphere),numberfieldexp(nsphere),j,mpicomm,&
                    rank,numprocs,task,nsend
         integer, optional :: mpi_comm
         real(8) :: alpha,beta,rpos(3,nsphere),rmax,rbeam(3),xib(3),rib, &
                    cbeam,epstran
         complex(8) :: pmnp(neqns*2),rimed(2),rimed0
         complex(8), allocatable :: pmnp0(:,:,:,:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)

         nodrgb=0
         rmax=0.d0
         rimed0=2.d0/(1.d0/rimed(1)+1.d0/rimed(2))
         do i=1,nsphere
            if(hostsphere(i).eq.0) then
               xib(:)=rpos(:,i)-rbeam(:)
               rib=sqrt(dot_product(xib,xib))
               rmax=max(rmax,rib)
               call tranordertest(rib,rimed0,nodr(i),epstran,n)
               nodrgb=max(n,nodrgb)
            endif
         enddo
         allocate(pmnp0(0:nodrgb+1,nodrgb,2,2))
         call gaussianbeamcoef(alpha,beta,cbeam,nodrgb,pmnp0)
         pmnp=0.d0
         noff=0
         task=0
         do i=1,nsphere
            do j=1,numberfieldexp(i)
               nblk=2*nodr(i)*(nodr(i)+2)*2
               if(hostsphere(i).eq.0.and.j.eq.1) then
                  task=task+1
                  if(mod(task,numprocs).eq.rank) then
                     xib(:)=rpos(:,i)-rbeam(:)
                     call rottran(nodrj=nodrgb,nodri=nodr(i),translation_vector=xib, &
                          refractive_index=rimed,vswf_type=1, &
                          cxj=pmnp0(0:nodrgb+1,1:nodrgb,1:2,1:2), &
                          cyi=pmnp(noff+1:noff+nblk), &
                          number_rhs=2)
                  endif
               endif
               noff=noff+nblk
            enddo
         enddo
         if(numprocs.gt.1) then
            nsend=neqns*2
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=pmnp, &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum, &
                    mpi_comm=mpicomm)
         endif
         deallocate(pmnp0)
         end subroutine spheregaussianbeamcoef
!
!  rotation of expansion coefficients amn by euler angles alpha,beta,gamma
!  idir=1: forward rotation, idir=-1, reverse rotation.
!
!
!  last revised: 15 January 2011
!
         subroutine rotvec(alpha,beta,gamma,nmax,mmax,amn,idir)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nmax,mmax,idir,k,n,m,in,kmax,ka,na,im,m1
         real(8) :: dc(-nmax-1:nmax+1,-nmax-1:nmax+1),dk0(-nmax-1:nmax+1), &
                    dk01(-nmax-1:nmax+1),sbe,cbe,sbe2,cbe2,sben,dkt, &
                    fmn,dkm0,dkm1,alpha,beta,gamma
         complex(8) :: ealpha,amn(0:nmax+1,nmax,2),ealpham(-nmax:nmax), &
                       amnt(2,-nmax:nmax),a,b,ci,egamma,egammam(-nmax:nmax)
         data ci/(0.d0,1.d0)/
         call init(nmax)
         dc=0.d0
         dk01=0.d0
         dk0=0.d0
         ealpha=cdexp(ci*alpha)
         egamma=cdexp(ci*gamma)
         cbe=cos(beta)
         sbe=sqrt((1.d0+cbe)*(1.d0-cbe))
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         call ephicoef(ealpha,nmax,ealpham)
         call ephicoef(egamma,nmax,egammam)
         in=1
         dk0(0)=1.d0
         sben=1.d0
         dk01(0)=0.d0
         do n=1,nmax
            kmax=min(n,mmax)
            do k=-kmax,kmax
               if(k.le.-1) then
                  ka=n+1
                  na=-k
               else
                  ka=k
                  na=n
               endif
               if(idir.eq.1) then
                  amnt(1,k)=amn(ka,na,1)*ealpham(k)
                  amnt(2,k)=amn(ka,na,2)*ealpham(k)
               else
                  amnt(1,-k)=amn(ka,na,1)*egammam(k)
                  amnt(2,-k)=amn(ka,na,2)*egammam(k)
               endif
            enddo
            in=-in
            sben=sben*sbe/2.d0
            dk0(n)=in*sben*bcof(n,n)
            dk0(-n)=in*dk0(n)
            dk01(n)=0.d0
            dk01(-n)=0.d0
            dc(0,n)=dk0(n)
            dc(0,-n)=dk0(-n)
            do k=-n+1,n-1
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt) &
                     /(fnr(n+k)*fnr(n-k))
               dc(0,k)=dk0(k)
            enddo
            im=1
            do m=1,kmax
               im=-im
               fmn=1./fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.
               do k=-n,n
                  dkm1=dkm0
                  dkm0=dc(m1,k)
                  dc(m,k)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                         -fnr(n-k)*fnr(n+k+1)*sbe2*dc(m1,k+1) &
                         -k*sbe*dc(m1,k))*fmn
                  dc(-m,-k)=dc(m,k)*(-1)**(k)*im
               enddo
            enddo
            do m=-n,n
               if(m.le.-1) then
                  ka=n+1
                  na=-m
               else
                  ka=m
                  na=n
               endif
               a=0.
               b=0.
               do k=-kmax,kmax
                  a=a+dc(-k,-m)*amnt(1,k)
                  b=b+dc(-k,-m)*amnt(2,k)
               enddo
               if(idir.eq.1) then
                  amn(ka,na,1)=a*egammam(m)
                  amn(ka,na,2)=b*egammam(m)
               else
                  amn(ka,na,1)=a*ealpham(m)
                  amn(ka,na,2)=b*ealpham(m)
               endif
            enddo
         enddo
         end subroutine rotvec

      end module translation
!
! scatprops module: various subroutines for calculation of observables from the solution
!
!
!  last revised: 15 January 2011
!
      module scatprops
      implicit none
      contains
!
!  determination of maximum orders for target--based expansions
!
!
!  last revised: 15 January 2011
!
         subroutine tranorders(nsphere,nodr,rpos,ri,hostsphere,eps,ntran,nodrt)
         use numconstants
         use specialfuncs
         use translation
         implicit none
         integer :: nsphere,nodr(nsphere),nodrt,ntran(nsphere),i, &
                    hostsphere(nsphere),host,nodrmax
         real(8) :: rpos(3,nsphere),r,eps,rpos0(3),xi0(3)
         complex(8) :: ri(2,0:nsphere),ri0
         nodrt=0
         nodrmax=maxval(nodr)
         do i=1,nsphere
            host=hostsphere(i)
            ri0=2.d0/(1.d0/ri(1,host)+1.d0/ri(2,host))
            if(host.eq.0) then
               rpos0=0.
            else
               rpos0(:)=rpos(:,host)
            endif
            xi0(:)=rpos(:,i)-rpos0(:)
            r=sqrt(dot_product(xi0,xi0))
            call tranordertest(r,ri0,nodr(i),eps,ntran(i))
            if(host.eq.0) nodrt=max(nodrt,ntran(i),nodrmax)
         enddo
         end subroutine tranorders
!
!  translation of sphere-based expansions to common target origin
!
!
!  last revised: 15 January 2011
!  april 2012: l/r formulation: amnp is in l/r basis, amnp0 is in e/m basis
!  february 2013: number rhs, mpi comm options.   This is for general sphere configurations.
!
         subroutine amncommonorigin(nsphere,nodr,ntran,nodrt,rpos,hostsphere, &
                    numberfieldexp,rimed,amnp,amnp0,number_rhs,mpi_comm)
         use specialfuncs
         use mpidefs
         use translation
         implicit none
         integer :: nsphere,nodr(nsphere),nodrt,i,m,n,nblk,ntran(nsphere),noff, &
                    hostsphere(nsphere),numberfieldexp(nsphere),ntrani,nrhs,task, &
                    rank,numprocs,nsend,proc,mpicomm
         integer, optional :: number_rhs,mpi_comm
         real(8) :: rpos(3,nsphere),xij(3)
         complex(8) :: amnp(*),amnp0(0:nodrt+1,nodrt,2,*),rimed(2)
         complex(8), allocatable :: amnpt(:,:,:,:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         amnp0(:,:,:,1:nrhs)=(0.d0,0.d0)
         noff=0
         task=0
         do i=1,nsphere
            nblk=nodr(i)*(nodr(i)+2)*2*nrhs
            if(hostsphere(i).eq.0) then
               task=task+1
               proc=mod(task,numprocs)
               if(proc.eq.rank) then
                  ntrani=min(nodrt,ntran(i))
                  allocate(amnpt(0:ntrani+1,ntrani,2,nrhs))
                  amnpt=(0.d0,0.d0)
                  xij=-rpos(:,i)
                  call rottran(nodrj=nodr(i),nodri=ntrani,translation_vector=xij, &
                       refractive_index=rimed,cxj=amnp(noff+1:noff+nblk), &
                       cyi=amnpt,vswf_type=1,number_rhs=nrhs)
                  do n=1,ntrani
                     do m=0,ntrani+1
                        amnp0(m,n,1,1:nrhs)=amnp0(m,n,1,1:nrhs) &
                             +amnpt(m,n,1,1:nrhs)+amnpt(m,n,2,1:nrhs)
                        amnp0(m,n,2,1:nrhs)=amnp0(m,n,2,1:nrhs) &
                             +amnpt(m,n,1,1:nrhs)-amnpt(m,n,2,1:nrhs)
                     enddo
                  enddo
                  deallocate(amnpt)
               endif
            endif
            noff=noff+nblk*numberfieldexp(i)
         enddo
         if(numprocs.gt.1) then
            nsend=2*nodrt*(nodrt+2)*nrhs
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=amnp0, &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
         end subroutine amncommonorigin
!
!  general efficiency factor calcuation
!  L/R formulation
!  April 2012
!  march 2013: something is always changing in this one
!
         subroutine lrsphereqeff(ri0,nodr,npol,xsp,anp,gnpinc,gnp,qe,qs,qa)
         use miecoefdata
         use spheredata
         use translation
         implicit none
         integer :: nodr,npol,m,n,p,p1,p2,k,ma,na,s,t,ss,st
         real(8) :: qe(2*npol-1),qa(2*npol-1),qs(2*npol-1),const,xsp,qi(2*npol-1)
         complex(8) :: anp(0:nodr+1,nodr,2,npol),gnp(0:nodr+1,nodr,2,npol), &
                       gnpinc(0:nodr+1,nodr,2,npol), &
                       psi(0:nodr,2),xi(0:nodr,2),psip(2),xip(2),xri(2),&
                       rib,aamat(2,2),ggmat(2,2),agmat(2,2), &
                       gamat(2,2),cterm,ri0(2)
         qe=0.d0
         qa=0.d0
         qs=0.d0
         qi=0.d0
         xri=xsp*ri0
         rib=2.d0/(1.d0/ri0(1)+1.d0/ri0(2))
         do p=1,2
            call cricbessel(nodr,xri(p),psi(0,p))
            call crichankel(nodr,xri(p),xi(0,p))
         enddo
         do n=1,nodr
            do s=1,2
               psip(s)=psi(n-1,s)-n*psi(n,s)/xri(s)
               xip(s)=xi(n-1,s)-n*xi(n,s)/xri(s)
            enddo
            do s=1,2
               ss=(-1)**s
               do t=1,2
                  st=(-1)**t
                  cterm=dcmplx(0.d0,1.d0)*conjg(1./ri0(t))/ri0(s)
                  aamat(s,t)=cterm*(xip(s)*conjg(xi(n,t)) &
                        -ss*st*xi(n,s)*conjg(xip(t)))
                  ggmat(s,t)=cterm*(psip(s)*conjg(psi(n,t)) &
                        -ss*st*psi(n,s)*conjg(psip(t)))
                  agmat(s,t)=cterm*(xip(s)*conjg(psi(n,t)) &
                        -ss*st*xi(n,s)*conjg(psip(t)))
                  gamat(s,t)=cterm*(psip(s)*conjg(xi(n,t)) &
                        -ss*st*psi(n,s)*conjg(xip(t)))
               enddo
            enddo
            do m=-n,n
               if(m.le.-1) then
                  ma=n+1
                  na=-m
               else
                  ma=m
                  na=n
               endif
               do p1=1,npol
                  do p2=1,npol
                     if(p1.eq.1.and.p2.eq.1) then
                        k=1
                        const=1.d0
                     elseif(p1.eq.2.and.p2.eq.2) then
                        k=2
                        const=1.d0
                     else
                        k=3
                        const=.5d0
                     endif
                     do s=1,2
                        do t=1,2
                           qi(k)=qi(k)+const*ggmat(s,t)*gnp(ma,na,s,p1) &
                             *conjg(rib*gnp(ma,na,t,p2))
                           qa(k)=qa(k)+const &
                              *(aamat(s,t)*anp(ma,na,s,p1)*conjg(rib*anp(ma,na,t,p2)) &
                              +ggmat(s,t)*gnp(ma,na,s,p1)*conjg(rib*gnp(ma,na,t,p2)) &
                              +agmat(s,t)*anp(ma,na,s,p1)*conjg(gnp(ma,na,t,p2)) &
                               *(conjg(rib)+(-1)**(s+t)*rib))
                           qs(k)=qs(k)-const &
                              *aamat(s,t)*anp(ma,na,s,p1)*conjg(rib*anp(ma,na,t,p2))
                           qe(k)=qe(k)+const &
                              *(agmat(s,t)*anp(ma,na,s,p1)*conjg(gnpinc(ma,na,t,p2)) &
                              *(conjg(rib)+(-1)**(s+t)*rib))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         do k=1,1+2*(npol-1)
            qi(k)=qi(k)*2./xsp/xsp
            qe(k)=qe(k)*2./xsp/xsp
            qa(k)=qa(k)*2./xsp/xsp
            qs(k)=qs(k)*2./xsp/xsp
         enddo

         end subroutine lrsphereqeff
!
! calling routine for efficiency calculation
! april 2012: lr formulation
! february 2013:  number of rhs, mpi comm options added.
!
         subroutine qefficiencyfactors(nsphere,neqns,nodr,npol,xsp,hostsphere, &
                    numberfieldexp,amnp,gmnp0,qext,qabs,qsca, &
                    number_rhs,mpi_comm)
         use mpidefs
         use miecoefdata
         use spheredata
         use translation
         implicit none
         integer :: nsphere,i,nodr(nsphere),nblk,noff,neqns, &
                    hostsphere(nsphere),npol,ntot,j,mpicomm, &
                    numberfieldexp(nsphere),b11,b12,nrhs,rank,numprocs,nsend
         integer, optional :: number_rhs,mpi_comm
         real(8) :: xsp(nsphere),qabs(nsphere,2*npol-1,*),qext(nsphere,2*npol-1,*), &
                    qsca(nsphere,2*npol-1,*)
         real(8), allocatable :: qa(:,:),qs(:,:),qe(:,:)
         complex(8) :: amnp(*),gmnp0(*),ri0(2)
         complex(8), allocatable :: amnpi(:,:,:,:),gmnpi(:,:,:,:),fmnpi(:,:,:,:), &
                     gmnpt(:),gmnp(:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         ntot=neqns*nrhs*npol
         allocate(gmnp(ntot),gmnpt(ntot),qa(2*npol-1,nrhs), &
            qe(2*npol-1,nrhs),qs(2*npol-1,nrhs))
         gmnp=0.d0
         gmnpt=0.d0
         call external_to_internal_expansion(neqns,nrhs*npol,amnp,gmnpt, &
              mpi_comm=mpicomm)
         call external_to_external_expansion(neqns,nrhs*npol,amnp,gmnp, &
                 far_field_option=.false.,store_matrix_option=.false., &
                 mpi_comm=mpicomm)
         if(numprocs.gt.1) then
            nsend=ntot
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=gmnpt(1:nsend), &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=gmnp(1:nsend), &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif

         qabs(1:nsphere,1:2*npol-1,1:nrhs)=0.d0
         qsca(1:nsphere,1:2*npol-1,1:nrhs)=0.d0
         qext(1:nsphere,1:2*npol-1,1:nrhs)=0.d0
         gmnp=gmnp+gmnp0(1:ntot)+gmnpt
         noff=0
         do i=1,nsphere
            nblk=2*nodr(i)*(nodr(i)+2)
            allocate(amnpi(0:nodr(i)+1,nodr(i),2,npol), &
                     gmnpi(0:nodr(i)+1,nodr(i),2,npol), &
                     fmnpi(0:nodr(i)+1,nodr(i),2,npol))
            call getspheredataone(hostsphere(i),sphere_refractive_index=ri0)
            b11=noff+1
            do j=1,nrhs
               b12=b11+nblk*npol-1
               amnpi(0:nodr(i)+1,1:nodr(i),1:2,1:npol)=reshape(amnp(b11:b12), &
                     (/nodr(i)+2,nodr(i),2,npol/))
               fmnpi(0:nodr(i)+1,1:nodr(i),1:2,1:npol)=reshape(gmnp0(b11:b12), &
                     (/nodr(i)+2,nodr(i),2,npol/))
               gmnpi(0:nodr(i)+1,1:nodr(i),1:2,1:npol)=reshape(gmnp(b11:b12), &
                     (/nodr(i)+2,nodr(i),2,npol/))
               call lrsphereqeff(ri0,nodr(i),npol,xsp(i),amnpi,fmnpi,gmnpi, &
                     qe(:,j),qs(:,j),qa(:,j))
               b11=b12+1
            enddo
            noff=noff+nblk*npol*nrhs*numberfieldexp(i)
            deallocate(gmnpi,amnpi,fmnpi)
            qext(i,:,1:nrhs)=qe(:,1:nrhs)
            qabs(i,:,1:nrhs)=qa(:,1:nrhs)
            qsca(i,:,1:nrhs)=qs(:,1:nrhs)
         enddo
         end subroutine qefficiencyfactors
!
!  scattering amplitude sa and matrix sm calculation
!
!  original: 15 January 2011
!  revised: 21 February 2011: S11 normalization changed
!  april 2013: moved things around to try to get it to work.
!
         subroutine scatteringmatrix(amn0,nodrt,ct,phi,sa,sm)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nodrt,m,n,p,m1,n1,i,j
         integer, save :: nodrtold
         real(8) :: ct,phi,sm(4,4),cphi,sphi,qsca
         real(8), save :: ctold
         real(8), save, allocatable :: tau(:,:,:)
         complex(8) :: amn0(0:nodrt+1,nodrt,2,2),sa(4),ephi,ephim(-nodrt:nodrt), &
                       ci,cin,a,b,sp(4,4)
         data ci,nodrtold,ctold/(0.d0,1.d0),0,0.d0/
         if(nodrt.ne.nodrtold.or.ct.ne.ctold) then
            if(allocated(tau)) deallocate(tau)
            allocate(tau(0:nodrt+1,nodrt,2))
            call taufunc(ct,nodrt,tau)
         endif
         nodrtold=nodrt
         ctold=ct
         cphi=cos(phi)
         sphi=sin(phi)
         ephi=dcmplx(cphi,sphi)
         call ephicoef(ephi,nodrt,ephim)
         sa=(0.d0,0.d0)
         qsca=0.d0
         do n=1,nodrt
            cin=(-ci)**n
            do m=-n,n
               if(m.le.-1) then
                  m1=n+1
                  n1=-m
               else
                  m1=m
                  n1=n
               endif
               do p=1,2
                  qsca=qsca+amn0(m1,n1,p,1)*dconjg(amn0(m1,n1,p,1)) &
                           + amn0(m1,n1,p,2)*dconjg(amn0(m1,n1,p,2))
                  a=amn0(m1,n1,p,1)*cphi+amn0(m1,n1,p,2)*sphi
                  b=amn0(m1,n1,p,1)*sphi-amn0(m1,n1,p,2)*cphi
                  sa(1)=sa(1)+cin*tau(m1,n1,3-p)*b*ephim(m)
                  sa(2)=sa(2)+ci*cin*tau(m1,n1,p)*a*ephim(m)
                  sa(3)=sa(3)+ci*cin*tau(m1,n1,p)*b*ephim(m)
                  sa(4)=sa(4)+cin*tau(m1,n1,3-p)*a*ephim(m)
               enddo
            enddo
         enddo
         qsca=qsca*2.d0
         do i=1,4
            do j=1,4
               sp(i,j)=sa(i)*dconjg(sa(j))*16.d0/qsca
            enddo
         enddo
         sm(1,1)=sp(1,1)+sp(2,2)+sp(3,3)+sp(4,4)
         sm(1,2)=-sp(1,1)+sp(2,2)-sp(3,3)+sp(4,4)
         sm(2,1)=-sp(1,1)+sp(2,2)+sp(3,3)-sp(4,4)
         sm(2,2)=sp(1,1)+sp(2,2)-sp(3,3)-sp(4,4)
         sm(3,3)=2.*(sp(1,2)+sp(3,4))
         sm(3,4)=2.*dimag(sp(2,1)+sp(4,3))
         sm(4,3)=2.*dimag(sp(1,2)-sp(3,4))
         sm(4,4)=2.*(sp(1,2)-sp(3,4))
         sm(1,3)=2.*(sp(2,3)+sp(1,4))
         sm(3,1)=2.*(sp(2,4)+sp(1,3))
         sm(1,4)=-2.*dimag(sp(2,3)-sp(1,4))
         sm(4,1)=-2.*dimag(sp(4,2)+sp(1,3))
         sm(2,3)=-2.*(sp(2,3)-sp(1,4))
         sm(3,2)=-2.*(sp(2,4)-sp(1,3))
         sm(2,4)=-2.*dimag(sp(2,3)+sp(1,4))
         sm(4,2)=-2.*dimag(sp(4,2)-sp(1,3))
!         do i=1,4
!            do j=1,4
!               if(i.ne.1.or.j.ne.1) then
!                  sm(i,j)=sm(i,j)/sm(1,1)
!               endif
!            enddo
!         enddo
         end subroutine scatteringmatrix
!   c                                                                               c
!   c  subroutine scatexp(amn0,nodrt,nodrg,gmn) computes the expansion coefficients c
!   c  for the spherical harmonic expansion of the scattering phase function from   c
!   c  the scattering coefficients amn0.  For a complete expansion, the max. order  c
!   c  of the phase function expansion (nodrg) will be 2*nodrt, where nodrt is      c
!   c  the max. order of the scattered field expansion.   In this code nodrg is     c
!   c  typically set to 1, so that the subroutine returns the first moments         c
!   c  of the phase function; gmn(1) and gmn(2).                                    c
!   c                                                                               c
!   c  The expansion coefficients are normalized so that gmn(0)=1                   c
!   c                                                                               c
!   c  gmn(1)/3 is the asymmetry parameter.                                         c
!   c                                                                               c
         subroutine s11expansion(amn0,nodrt,mmax,nodrg,gmn)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nodrt,m,n,p,ma,na,mmax,nodrg,w,w1,w2,u,uw,ww1, &
                    l1,l2,ka,la,k,l,q,ik
         real(8) :: vc1(0:nodrt*2+1),vc2(0:nodrt*2+1),g0
         complex(8) :: amn0(0:nodrt+1,nodrt,2,2),gmn(0:nodrg*(nodrg+3)/2), &
                       a(2,2),c,c2
         gmn=(0.d0,0.d0)
         do n=1,nodrt
            l1=max(1,n-nodrg)
            l2=min(nodrt,n+nodrg)
            do l=l1,l2
               c=sqrt(dble((n+n+1)*(l+l+1)))*dcmplx(0.d0,1.d0)**(l-n)
               w2=min(n+l,nodrg)
               call vcfunc(-1,l,1,n,vc2)
               do m=-n,n
                  if(m.le.-1) then
                     ma=n+1
                     na=-m
                  else
                     ma=m
                     na=n
                  endif
                  do k=-l,min(l,m)
                     if(k.le.-1) then
                        ka=l+1
                        la=-k
                     else
                        ka=k
                        la=l
                     endif
                     u=m-k
                     if(u.le.mmax) then
                        ik=(-1)**k
                        c2=ik*c
                        do p=1,2
                           do q=1,2
                              a(p,q)=c2*(amn0(ma,na,p,1)*conjg(amn0(ka,la,q,1)) &
                                    +amn0(ma,na,p,2)*conjg(amn0(ka,la,q,2)))
                           enddo
                        enddo
                        w1=max(abs(n-l),abs(u))
                        w2=min(n+l,nodrg)
                        call vcfunc(-k,l,m,n,vc1)
                        do w=w1,w2
                           uw=(w*(w+1))/2+u
                           do p=1,2
                              if(mod(n+l+w,2).eq.0) then
                                 q=p
                              else
                                 q=3-p
                              endif
                              gmn(uw)=gmn(uw)-vc1(w)*vc2(w)*a(p,q)
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
         g0=dble(gmn(0))
         gmn(0)=1.d0
         do w=1,nodrg
            ww1=(w*(w+1))/2
            gmn(ww1)=dcmplx(dble(gmn(ww1)),0.d0)/g0
            do u=1,min(mmax,w)
               uw=ww1+u
               gmn(uw)=(-1)**u*2.d0*gmn(uw)/g0
            enddo
         enddo
         end subroutine s11expansion
!
!  calculate azimuth--averaged scattering matrix from expansion, for cos(theta) = ct
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  this is currently not used in v. 3.0
!
         subroutine fosmcalc(ntot,s00,s02,sp22,sm22,ct,sm)
         use numconstants
         use specialfuncs
         integer :: ntot,w,ww1
         real(8) :: s00(4,4,0:ntot*2),s02(4,4,0:ntot*2),sp22(4,4,0:ntot*2),sm22(4,4,0:ntot*2), &
                    sm(4,4),dc(-2:2,0:2*ntot*(2*ntot+2)),ct,temp
         call rotcoef(ct,2,2*ntot,dc)
         sm=0.d0
         do w=0,2*ntot
            ww1=w*(w+1)
            sm(:,:)=sm(:,:)+s00(:,:,w)*dc(0,ww1)+s02(:,:,w)*dc(0,ww1+2) &
                   +sp22(:,:,w)*dc(2,ww1+2)+sm22(:,:,w)*dc(-2,ww1+2)
         enddo
         sm=sm/s00(1,1,0)
!
!  a patch
!
         sm(3,1)=-sm(3,1)
         sm(1,3)=-sm(1,3)
         sm(4,3)=-sm(4,3)
         temp=sm(4,1)
         sm(4,1)=sm(4,2)
         sm(4,2)=temp
!         do i=1,4
!            do j=1,4
!               if(i.ne.1.or.j.ne.1) then
!                  sm(i,j)=sm(i,j)/sm(1,1)
!               endif
!            enddo
!         enddo
         end subroutine fosmcalc
!
!  determine the generalized spherical function expansion for the azimuth-averaged scattering matrix
!  corresponding to the target-based scattering field expansion of amnp.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fixed flush call.
!  this is currently not used in v. 3.0
!
         subroutine fosmexpansion(ntot,amnp,s00,s02,sp22,sm22)
         use mpidefs
         use mpidata
         use specialfuncs
         use numconstants
         use spheredata
         integer :: ntot,n,p,m,l,wmin,wmax,m1m,q,m1mq,m1mnpl,w,m1w,fe,fo,i,wtot,j
         integer :: rank,numprocs,nl,nsend,runprintunit
         integer, allocatable :: nlindex(:),nlnum(:)
         real(8) :: s00(4,4,0:ntot*2),s02(4,4,0:ntot*2),sp22(4,4,0:ntot*2),sm22(4,4,0:ntot*2), &
                    cm1p1(0:ntot*2),cm1m1(0:ntot*2),cmmpm(0:ntot*2),cmmm2pm(0:ntot*2), &
                    cmmp2pm(0:ntot*2),sum,nlperproc
         complex(8) :: amnp(0:ntot+1,ntot,2,2),a1(-ntot-2:ntot+2,ntot,2),a2(-ntot-2:ntot+2,ntot,2), &
                       ci,fnl,a1122,a2112,a1p2,a1m2
         data ci/(0.d0,1.d0)/
         call init(2*ntot)
         call getrunparameters(run_print_unit=runprintunit)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
         allocate(nlindex(0:numprocs-1),nlnum(0:numprocs-1))
         nlperproc=dble(ntot*ntot)/dble(numprocs)
         sum=0.
         do i=0,numprocs-1
            nlindex(i)=floor(sum)
            sum=sum+nlperproc
         enddo
         do i=0,numprocs-2
            nlnum(i)=nlindex(i+1)-nlindex(i)
         enddo
         nlnum(numprocs-1)=ntot*ntot-nlindex(numprocs-1)
         if(rank.eq.0) then
            write(runprintunit,'('' SM calc, orders per processor:'',f10.4)') nlperproc
            call flush(runprintunit)
         endif
         a1=(0.d0,0.d0)
         a2=(0.d0,0.d0)
         s00=0.d0
         s02=0.d0
         sp22=0.d0
         sm22=0.d0
         wtot=ntot+ntot
         do n=1,ntot
            do p=1,2
               do m=-n,-1
                  a1(m,n,p)=amnp(n+1,-m,p,1)
                  a2(m,n,p)=amnp(n+1,-m,p,2)
               enddo
               do m=0,n
                  a1(m,n,p)=amnp(m,n,p,1)
                  a2(m,n,p)=amnp(m,n,p,2)
               enddo
            enddo
         enddo
         do nl=nlindex(rank)+1,nlindex(rank)+nlnum(rank)
            n=floor((nl-1)/dble(ntot))+1
            l=mod(nl-1,ntot)+1
            wmin=abs(n-l)
            wmax=n+l
            fnl=sqrt(dble((n+n+1)*(l+l+1)))*ci**(l-n)
            call vcfunc(-1,n,1,l,cm1p1)
            call vcfunc(-1,n,-1,l,cm1m1)
            do m=-min(n,l+2),min(n,l+2)
               m1m=(-1)**m
               if(abs(m).le.l) then
                  call vcfunc(-m,n,m,l,cmmpm)
               else
                  cmmpm=0.d0
               endif
               if(abs(-2+m).le.l) then
                  call vcfunc(-m,n,-2+m,l,cmmm2pm)
               else
                  cmmm2pm=0.d0
               endif
               if(abs(2+m).le.l) then
                  call vcfunc(-m,n,2+m,l,cmmp2pm)
               else
                  cmmp2pm=0.d0
               endif
               do p=1,2
                  do q=1,2
                     m1mq=(-1)**(m+q)
                     m1mnpl=(-1)**(m+n+p+l)
                     a1122=(a1(m,n,p)*conjg(a1(m,l,q)) + a2(m,n,p)*conjg(a2(m,l,q)))
                     a2112=(a2(m,n,p)*conjg(a1(m,l,q)) - a1(m,n,p)*conjg(a2(m,l,q)))
                     a1p2=(a1(m,n,p)+ci*a2(m,n,p))*conjg(a1(m-2,l,q)-ci*a2(m-2,l,q))
                     a1m2=(a1(m,n,p)-ci*a2(m,n,p))*conjg(a1(m+2,l,q)+ci*a2(m+2,l,q))
                     do w=wmin,wmax
                        m1w=(-1)**w
                        if(mod(n+l+w+p+q,2).eq.0) then
                           fe=1
                           fo=0
                        else
                           fe=0
                           fo=1
                        endif
                        s00(1,1,w) = s00(1,1,w)-(m1m*fe*fnl*a1122*cm1p1(w)*cmmpm(w))/2.
                        s00(3,2,w) = s00(3,2,w)+ (ci/2.*m1m*fnl*fo*a1122*cm1p1(w)*cmmpm(w))
                        s00(4,2,w) = s00(4,2,w)+ dimag(-ci/2.*m1m*fnl*fo*a1122*cm1p1(w)*cmmpm(w))
                        s00(1,4,w) = s00(1,4,w)+ dimag(m1m*fe*fnl*(-a2112)*cm1p1(w)*cmmpm(w))/2.
                        s00(2,3,w) = s00(2,3,w)+ (m1m*fe*fnl*(-a2112)*cm1p1(w)*cmmpm(w))/2.
                        s00(4,3,w) = s00(4,3,w)+ dimag(ci/2.*m1m*fnl*fo*a2112*cm1p1(w)*cmmpm(w))
                        s00(4,4,w) = s00(4,4,w)+ (ci/2.*m1m*fnl*fo*a2112*cm1p1(w)*cmmpm(w))

                        if(w.lt.2) cycle

                        s02(2,1,w) = s02(2,1,w)-(m1mq*a1122*fe*fnl*cm1m1(w)*cmmpm(w))/2.
                        s02(3,1,w) = s02(3,1,w)+ (-ci/2.*m1mq*a1122*fnl*fo*cm1m1(w)*cmmpm(w))
                        s02(4,1,w) = s02(4,1,w)+ dimag(ci/2.*m1mq*a1122*fnl*fo*cm1m1(w)*cmmpm(w))
                        s02(1,3,w) = s02(1,3,w)-(m1mq*a2112*fe*fnl*cm1m1(w)*cmmpm(w))/2.
                        s02(2,4,w) = s02(2,4,w)-dimag(m1mq*a2112*fe*fnl*cm1m1(w)*cmmpm(w))/2.
                        s02(3,3,w) = s02(3,3,w)+ (ci/2.*m1mq*a2112*fnl*fo*cm1m1(w)*cmmpm(w))
                        s02(3,4,w) = s02(3,4,w)+ dimag(-ci/2.*m1mq*a2112*fnl*fo*cm1m1(w)*cmmpm(w))

                        s02(1,2,w) = s02(1,2,w)-(m1m*a1p2*fe*fnl*cm1p1(w)*cmmm2pm(w))/4.
                        s02(1,3,w) = s02(1,3,w)+ (-ci/4.*m1m*a1p2*fe*fnl*cm1p1(w)*cmmm2pm(w))
                        s02(2,4,w) = s02(2,4,w)+ dimag(-ci/4.*m1m*a1p2*fe*fnl*cm1p1(w)*cmmm2pm(w))
                        s02(3,1,w) = s02(3,1,w)+ (ci/4.*m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))
                        s02(4,1,w) = s02(4,1,w)+ dimag(-ci/4.*m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))
                        s02(4,3,w) = s02(4,3,w)+ dimag(m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))/4.
                        s02(4,4,w) = s02(4,4,w)+ (m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))/4.

                        sm22(1,4,w) = sm22(1,4,w)+ dimag(-ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sm22(2,2,w) = sm22(2,2,w)-(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sm22(2,3,w) = sm22(2,3,w)+ (-ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sm22(3,2,w) = sm22(3,2,w)+ (ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sm22(3,3,w) = sm22(3,3,w)-(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sm22(3,4,w) = sm22(3,4,w)+ dimag(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sm22(4,2,w) = sm22(4,2,w)+ dimag(-ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))

                        sp22(1,4,w) = sp22(1,4,w)+ dimag(-ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sp22(2,2,w) = sp22(2,2,w)-(m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sp22(2,3,w) = sp22(2,3,w)+ (-ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sp22(3,2,w) = sp22(3,2,w)+ (-ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sp22(3,3,w) = sp22(3,3,w)+ (m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sp22(3,4,w) = sp22(3,4,w)-dimag(m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sp22(4,2,w) = sp22(4,2,w)+ dimag(ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))

                        s02(1,2,w) = s02(1,2,w)-(m1m*a1m2*fe*fnl*cm1p1(w)*cmmp2pm(w))/4.
                        s02(1,3,w) = s02(1,3,w)+ (ci/4.*m1m*a1m2*fe*fnl*cm1p1(w)*cmmp2pm(w))
                        s02(2,4,w) = s02(2,4,w)+ dimag(ci/4.*m1m*a1m2*fe*fnl*cm1p1(w)*cmmp2pm(w))
                        s02(3,1,w) = s02(3,1,w)+ (ci/4.*m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))
                        s02(4,1,w) = s02(4,1,w)+ dimag(-ci/4.*m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))
                        s02(4,3,w) = s02(4,3,w)-dimag(m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))/4.
                        s02(4,4,w) = s02(4,4,w)-(m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))/4.

                        sm22(1,4,w) = sm22(1,4,w)+ dimag(ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sm22(2,2,w) = sm22(2,2,w)-(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sm22(2,3,w) = sm22(2,3,w)+ (ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sm22(3,2,w) = sm22(3,2,w)+ (-ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sm22(3,3,w) = sm22(3,3,w)-(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sm22(3,4,w) = sm22(3,4,w)+ dimag(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sm22(4,2,w) = sm22(4,2,w)+ dimag(ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))

                        sp22(1,4,w) = sp22(1,4,w)+ dimag(ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sp22(2,2,w) = sp22(2,2,w)-(m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sp22(2,3,w) = sp22(2,3,w)+ (ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sp22(3,2,w) = sp22(3,2,w)+ (ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sp22(3,3,w) = sp22(3,3,w)+ (m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sp22(3,4,w) = sp22(3,4,w)-dimag(m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sp22(4,2,w) = sp22(4,2,w)+ dimag(-ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         call mstm_mpi(mpi_command='barrier')
         nsend=4*4*(2*ntot+1)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=s00,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=s02,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=sp22,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=sm22,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum)
!
!  a patch
!
         do i=3,4
            do j=1,i
               s00(j,i,0:wtot)=-s00(j,i,0:wtot)
               s02(j,i,0:wtot)=-s02(j,i,0:wtot)
               sm22(j,i,0:wtot)=-sm22(j,i,0:wtot)
               sp22(j,i,0:wtot)=-sp22(j,i,0:wtot)
            enddo
         enddo
         deallocate(nlindex,nlnum)
         end subroutine fosmexpansion
!
!  compute the coefficients for the GSF expansion of the random orientation
!  scattering matrix.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  January 2012: added computation of coherent field average
!  February 2013:  added number processors option.
!
         subroutine ranorientscatmatrix(xv,nsphere,nodr,nodrw,cbeam,tmatrixfile,&
                    sm,smcf,qext,qabs,qsca,number_processors)
         use mpidefs
         use mpidata
         use intrinsics
         use specialfuncs
         use spheredata
         use numconstants
         implicit none
         integer :: nodr,nodrw,nodr2,m,n,p,k,l,q,t,v,u,w,nblk,kl,mn,nn1,tn, &
                    lmax,ll1,tvl,ku,k1,ns,ik,ik1,m1,nu,n1s,n1e,nu1,p1,n1max, &
                    in,n1,i,lt,kt,qt,nt,mt,ikm,klm,mnm,nodrt,nsphere, &
                    rank,iunit,numprocs,j,rhscondition, &
                    nodrrhs,nblkrhs,nstop,numprocscalc,orig_group,new_group, &
                    new_comm,new_rank
         integer, optional :: number_processors
         integer, allocatable :: group_list(:)
         real(8) :: sm(4,4,0:nodrw),fl,vc(0:4*nodr+2),xv,fl2,cbeam,gbn,&
                    qext(nsphere),qabs(nsphere),qsca(nsphere),qel,qal,qsl, &
                    fc(4),time1,time2,qsca0,smcf(4,4,0:nodrw)
         complex(8) :: ci,cin,a
         complex(8) :: aw(0:2,-1:1,0:nodrw),bw(0:2,-1:1,0:nodrw),cw(0:nodrw), &
                       dw(0:nodrw),pp(nodr,2,2),bm(2,nodr*(nodr+2),2), &
                       am(2,nodr+1,2),fm(3,nodr,2,nodr,2),bmcf(2,nodr*(nodr+2),2), &
                       amcf(2,nodr+1,2),fmcf(3,nodr,2,nodr,2),awcf(0:2,-1:1,0:nodrw), &
                       bwcf(0:2,-1:1,0:nodrw),cwcf(0:nodrw),dwcf(0:nodrw)
         complex(8), allocatable :: dm(:,:,:,:,:,:),dmcf(:,:,:,:,:,:)
         complex(4), allocatable :: tc(:,:,:,:)
         integer :: nblkw,wv,sizedm,sizetm
         integer, allocatable :: windex(:),vindex(:),wvindex(:),wvnum(:)
         real(8) :: wvperproc,sum
         character*30 :: tmatrixfile
         data ci/(0.d0,1.d0)/
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
         if(present(number_processors)) then
            numprocscalc=min(number_processors,numprocs)
         else
            numprocscalc=numprocs
         endif
         allocate(group_list(numprocscalc))
         group_list=(/(i,i=0,numprocscalc-1)/)
         call mstm_mpi(mpi_command='group',mpi_group=orig_group)
         call mstm_mpi(mpi_command='incl',mpi_group=orig_group, &
            mpi_size=numprocscalc,mpi_new_group_list=group_list, &
            mpi_new_group=new_group)
         call mstm_mpi(mpi_command='create',mpi_group=new_group, &
            mpi_new_comm=new_comm)
         if(rank.le.numprocscalc-1) then

            call mstm_mpi(mpi_command='rank',mpi_comm=new_comm, &
               mpi_rank=new_rank)
            if(rank.eq.0) time1=mytime()
!
!  read the T matrix from the file
!
            if(new_rank.eq.0) then
               call getrunparameters(run_print_unit=iunit)
               open(3,file=tmatrixfile)
               read(3,*) rhscondition,nodrt,nodrrhs
               read(3,*) nt,xv
            endif
            nodrt=nodr
            nblk=nodr*(nodr+2)
            nodrrhs=nodr
            nblkrhs=nodrrhs*(nodrrhs+2)
            sizetm=4*nblk*nblkrhs
            allocate(tc(2,nblk,2,nblkrhs))
            tc=(0.,0.)
            if(new_rank.eq.0) then
               qext=0.d0
               qabs=0.d0
               qsca=0.d0
               do l=1,nodrrhs
                  if(rhscondition.eq.0) then
                     nstop=l
                  else
                     nstop=nodr
                  endif
                  gbn=dexp(-((dble(l)+.5d0)*cbeam)**2.)
                  do k=-l,l
                     kl=l*(l+1)+k
                     klm=l*(l+1)-k
                     do q=1,2
                        read(3,*) lt,kt,qt
                        do n=1,nstop
                           do m=-n,n
                              mn=n*(n+1)+m
                              mnm=n*(n+1)-m
                              read(3,*) nt,mt,fc
                              tc(1,mn,q,kl)=cmplx(fc(1),fc(2))
                              tc(2,mn,q,kl)=cmplx(fc(3),fc(4))
                              if(n.lt.l.and.rhscondition.eq.0) then
                                 ikm=(-1)**(m+k)
                                 do p=1,2
                                    tc(q,klm,p,mnm)=tc(p,mn,q,kl)*ikm
                                 enddo
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
                  do i=1,nsphere
                     read(3,*) n,qel,qal,qsl
                     qext(i)=qext(i)+qel*gbn*gbn
                     qabs(i)=qabs(i)+qal*gbn*gbn
                     qsca(i)=qsca(i)+qsl*gbn*gbn
                  enddo
               enddo
               close(3)
            endif
!
!  send to the other processors
!
            if(numprocscalc.gt.1) then
               call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qext, &
                    mpi_number=nsphere,mpi_rank=0,mpi_comm=new_comm)
               call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qabs, &
                    mpi_number=nsphere,mpi_rank=0,mpi_comm=new_comm)
               call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qsca, &
                    mpi_number=nsphere,mpi_rank=0,mpi_comm=new_comm)
               call mstm_mpi(mpi_command='bcast',mpi_send_buf_c=tc, &
                    mpi_number=sizetm,mpi_rank=0,mpi_comm=new_comm)
            endif

            allocate(dm(-nodr-1:nodr+1,3,nodr,2,nodr,2),dmcf(-nodr-1:nodr+1,3,nodr,2,nodr,2))
            if(new_rank.eq.0) then
               time2=mytime()-time1
               call timewrite(iunit,' t matrix read time:',time2)
               time1=mytime()
            endif
            nodr2=nodr+nodr
            nblk=nodr*(nodr+2)
            dm=(0.d0,0.d0)
            dmcf=(0.d0,0.d0)
            sizedm=size(dm)
            call init(nodr2)
!
!  compute the GB modified T matrix
!
            do n=1,nodrrhs
               gbn=dexp(-((dble(n)+.5d0)*cbeam)**2.)
               cin=ci**(n+1)
               pp(n,1,1) =-.5d0*cin*fnr(n+n+1)*gbn
               pp(n,2,1) =-pp(n,1,1)
               pp(n,1,2)=-pp(n,1,1)
               pp(n,2,2)=pp(n,2,1)
            enddo
            do n=1,nodr
               nn1=n*(n+1)
               do m=-n,n
                  mn=nn1+m
                  do p=1,2
                     do l=1,nodrrhs
                        do k=-l,l
                           kl=l*(l+1)+k
                           a=tc(p,mn,1,kl)
                           tc(p,mn,1,kl)=tc(p,mn,1,kl)*pp(l,1,1)&
                              +tc(p,mn,2,kl)*pp(l,1,2)
                           tc(p,mn,2,kl)=a*pp(l,2,1)+tc(p,mn,2,kl)*pp(l,2,2)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
!
!  determine the distribution of work load among the processors
!
            nblkw=nodr2*(nodr2+2)+1
            allocate(windex(nblkw),vindex(nblkw),wvindex(0:numprocscalc-1), &
               wvnum(0:numprocscalc-1))
            w=0
            do n=0,nodr2
               do m=-n,n
                  w=w+1
                  windex(w)=n
                  vindex(w)=m
               enddo
            enddo
            wvperproc=dble(nblkw)/dble(numprocscalc)
            sum=0.
            do i=0,numprocscalc-1
               wvindex(i)=floor(sum)
               sum=sum+wvperproc
            enddo
            do i=0,numprocscalc-2
               wvnum(i)=wvindex(i+1)-wvindex(i)
            enddo
            wvnum(numprocscalc-1)=nblkw-wvindex(numprocscalc-1)
            if(new_rank.eq.0) then
               write(iunit,'('' d matrix calculation, order+degree per proc.:'',f9.2)') &
                   wvperproc
               call flush(iunit)
            endif
!
!  the big loop
!
            do wv=wvindex(rank)+1,wvindex(rank)+wvnum(rank)
               w=windex(wv)
               v=vindex(wv)
               bm=(0.d0,0.d0)
               bmcf=(0.d0,0.d0)
               do n=1,nodr
                  nn1=n*(n+1)
                  do l=max(1,abs(w-n)),min(nodr,w+n)
                     am(1,l,1)=0.d0
                     am(1,l,2)=0.d0
                     am(2,l,1)=0.d0
                     am(2,l,2)=0.d0
                     amcf(1,l,1)=0.d0
                     amcf(1,l,2)=0.d0
                     amcf(2,l,1)=0.d0
                     amcf(2,l,2)=0.d0
                  enddo
                  do t=-n,n
                     tn=nn1+t
                     lmax=min(nodrrhs,w+n)
                     call vcfunc(v,w,-t,n,vc)
                     do l=max(1,abs(v-t),abs(n-w)),lmax
                        ll1=l*(l+1)
                        tvl=ll1+t-v
                        do k=1,2
                           do p=1,2
                              am(k,l,p)=am(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                              if(l.eq.n.and.v.eq.0) then
                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
                  do m=-n,n
                     mn=nn1+m
                     do k=1,2
                        u=m-(-3+2*k)
                        if(abs(u).le.w) then
                           lmax=min(nodrrhs,w+n)
                           call vcfunc(-u,w,m,n,vc)
                           do l=max(1,abs(w-n)),lmax
                              fl=-(-1)**l*vc(l)/dble(l+l+1)
                              do p=1,2
                                 bm(k,mn,p)=bm(k,mn,p)+am(k,l,p)*fl
                                 if(v.eq.0) then
                                    bmcf(k,mn,p)=bmcf(k,mn,p)+amcf(k,l,p)*fl
                                 endif
                              enddo
                           enddo
                        endif
                     enddo
                  enddo
               enddo
               do u=-min(w,nodr+1),min(w,nodr+1)
                  do ku=1,3
                     if(ku.eq.1) then
                        k=-1
                        k1=-1
                     elseif(ku.eq.2) then
                        k=1
                        k1=1
                     else
                        k=1
                        k1=-1
                     endif
                     m=u+k
                     ns=max(1,abs(m))
                     ik=(k+1)/2+1
                     ik1=(k1+1)/2+1
                     m1=u+k1
                     do n=ns,nodr
                        nu=n*(n+1)+m
                        n1s=max(1,abs(m1),n-nodrw)
                        n1e=min(nodr,n+nodrw)
                        do n1=n1s,n1e
                           cin=ci**(n-n1)
                           nu1=n1*(n1+1)+m1
                           fl=-fnr(n+n+1)*fnr(n1+n1+1)*dble(w+w+1)
                           do p=1,2
                              do p1=1,2
                                 a=bm(ik,nu,p)*cin*fl*conjg(bm(ik1,nu1,p1))
                                 dm(u,ku,n,p,n1,p1)=dm(u,ku,n,p,n1,p1)+a
                                 if(v.eq.0) then
                                    a=bmcf(ik,nu,p)*cin*fl*conjg(bmcf(ik1,nu1,p1))
                                    dmcf(u,ku,n,p,n1,p1)=dmcf(u,ku,n,p,n1,p1)+a
                                 endif
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            deallocate(tc)

            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=dm,&
                 mpi_number=sizedm,mpi_operation=mstm_mpi_sum,mpi_comm=new_comm)
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=dmcf,&
                 mpi_number=sizedm,mpi_operation=mstm_mpi_sum,mpi_comm=new_comm)
            if(new_rank.eq.0) then
               time2=mytime()-time1
               call timewrite(iunit,' d matrix time:',time2)
               time1=mytime()
            endif
!
!  compute the expansion coefficients
!
            aw=0.d0
            bw=0.d0
            cw=0.d0
            dw=0.d0
            awcf=0.d0
            bwcf=0.d0
            cwcf=0.d0
            dwcf=0.d0
            do w=0,nodrw
               do n=1,nodr
                  n1s=max(1,abs(n-w))
                  n1e=min(nodr,n+w)
                  do n1=n1s,n1e
                     do k=1,3
                        do p=1,2
                           do p1=1,2
                              fm(k,n,p,n1,p1)=0.
                              fmcf(k,n,p,n1,p1)=0.
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               do u=-nodr-1,nodr+1
                  do k=-1,1,2
                     m=u+k
                     ik=(k+1)/2+1
                     ns=max(1,abs(m))
                     do n=ns,nodr
                        n1max=min(w+n,nodr)
                        call vcfunc(m,n,0,w,vc)
                        do n1=ns,nodr
                           if((n+n1.lt.w).or.(abs(n-n1).gt.w)) cycle
                           fl=-(-1)**n*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                           do p=1,2
                              do p1=1,2
                                 fm(ik,n,p,n1,p1)=fm(ik,n,p,n1,p1) &
                                    +dm(u,ik,n,p,n1,p1)*fl
                                 fmcf(ik,n,p,n1,p1)=fmcf(ik,n,p,n1,p1) &
                                    +dmcf(u,ik,n,p,n1,p1)*fl
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  if(w.lt.2) cycle
                  m=u+1
                  m1=u-1
                  ns=max(1,abs(m))
                  n1s=max(1,abs(m1))
                  do n=ns,nodr
                     n1max=min(w+n,nodr)
                     call vcfunc(m,n,-2,w,vc)
                     do n1=n1s,nodr
                        if((n+n1.lt.w).or.(abs(n-n1).gt.w)) cycle
                        fl=-(-1)**n*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                        do p=1,2
                           do p1=1,2
                              fm(3,n,p,n1,p1)=fm(3,n,p,n1,p1) &
                                 +dm(u,3,n,p,n1,p1)*fl
                              fmcf(3,n,p,n1,p1)=fmcf(3,n,p,n1,p1) &
                                 +dmcf(u,3,n,p,n1,p1)*fl
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               do n=1,nodr
                  n1s=max(1,abs(n-w))
                  n1e=min(nodr,n+w)
                  in=(-1)**n
                  n1max=min(w+n,nodr)
                  call vcfunc(1,n,0,w,vc)
                  do n1=n1s,n1e
                     fl=2.d0*in*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     i=mod(n+n1-w,2)+1
                     do p=1,2
                        p1=(2-i)*p+(i-1)*(3-p)
                        do k=-1,1,2
                           ik=(k+1)/2+1
                           aw(0,k,w)=aw(0,k,w)+fm(ik,n,p,n1,p1)*fl
                           bw(0,k,w)=bw(0,k,w)+fm(ik,n,p,n1,3-p1)*fl
                           awcf(0,k,w)=awcf(0,k,w)+fmcf(ik,n,p,n1,p1)*fl
                           bwcf(0,k,w)=bwcf(0,k,w)+fmcf(ik,n,p,n1,3-p1)*fl
                        enddo
                        bw(2,0,w)=bw(2,0,w)+fm(3,n,p,n1,3-p1)*fl
                        aw(2,0,w)=aw(2,0,w)+fm(3,n,p,n1,p1)*fl
                        bwcf(2,0,w)=bwcf(2,0,w)+fmcf(3,n,p,n1,3-p1)*fl
                        awcf(2,0,w)=awcf(2,0,w)+fmcf(3,n,p,n1,p1)*fl
                     enddo
                  enddo
                  if(w.lt.2) cycle
                  call vcfunc(1,n,-2,w,vc)
                  do n1=n1s,n1e
                     fl=2.d0*in*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     i=mod(n+n1-w,2)+1
                     do p=1,2
                        p1=(2-i)*p+(i-1)*(3-p)
                        do k=-1,1,2
                           ik=(k+1)/2+1
                           aw(2,k,w)=aw(2,k,w)+fm(ik,n,p,n1,p1)*fl*(-1)**p1
                           bw(2,k,w)=bw(2,k,w)+fm(ik,n,p,n1,3-p1)*fl*(-1)**(3-p1)
                           awcf(2,k,w)=awcf(2,k,w)+fmcf(ik,n,p,n1,p1)*fl*(-1)**p1
                           bwcf(2,k,w)=bwcf(2,k,w) &
                              +fmcf(ik,n,p,n1,3-p1)*fl*(-1)**(3-p1)
                        enddo
                     enddo
                     fl2=2.*(-1)**(n1+w)*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     do p=1,2
                        do p1=1,2
                           cw(w)=cw(w)+fm(3,n,p,n1,p1)*fl*(-1)**p1
                           dw(w)=dw(w)+fm(3,n,p,n1,p1)*fl2*(-1)**p
                           cwcf(w)=cwcf(w)+fmcf(3,n,p,n1,p1)*fl*(-1)**p1
                           dwcf(w)=dwcf(w)+fmcf(3,n,p,n1,p1)*fl2*(-1)**p
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do w=0,nodrw
               do k=-1,1
                  do i=0,2
                     aw(i,k,w)=aw(i,k,w)*2./xv/xv
                     bw(i,k,w)=bw(i,k,w)*2./xv/xv
                     awcf(i,k,w)=awcf(i,k,w)*2./xv/xv
                     bwcf(i,k,w)=bwcf(i,k,w)*2./xv/xv
                  enddo
               enddo
               cw(w)=cw(w)*2./xv/xv
               dw(w)=dw(w)*2./xv/xv
               cwcf(w)=cwcf(w)*2./xv/xv
               dwcf(w)=dwcf(w)*2./xv/xv
            enddo
            do n=0,nodrw
               sm(1,1,n)=aw(0,-1,n)+aw(0,1,n)
               sm(1,2,n)=aw(2,-1,n)+aw(2,1,n)
               sm(1,3,n)=2.d0*dimag(aw(2,0,n))
               sm(1,4,n)=aw(0,1,n)-aw(0,-1,n)
               sm(2,2,n)=dw(n)
               sm(2,3,n)=dimag(dw(n))
               sm(2,4,n)=aw(2,1,n)-aw(2,-1,n)
               sm(3,1,n)=-dimag(bw(2,-1,n)+bw(2,1,n))
               sm(3,2,n)=dimag(cw(n))
               sm(3,3,n)=cw(n)
               sm(3,4,n)=dimag(bw(2,-1,n)-bw(2,1,n))
               sm(4,1,n)=bw(0,-1,n)+bw(0,1,n)
               sm(4,2,n)=2.*bw(2,0,n)
               sm(4,4,n)=bw(0,1,n)-bw(0,-1,n)

               smcf(1,1,n)=awcf(0,-1,n)+awcf(0,1,n)
               smcf(1,2,n)=awcf(2,-1,n)+awcf(2,1,n)
               smcf(1,3,n)=2.d0*dimag(awcf(2,0,n))
               smcf(1,4,n)=awcf(0,1,n)-awcf(0,-1,n)
               smcf(2,2,n)=dwcf(n)
               smcf(2,3,n)=dimag(dwcf(n))
               smcf(2,4,n)=awcf(2,1,n)-awcf(2,-1,n)
               smcf(3,1,n)=-dimag(bwcf(2,-1,n)+bwcf(2,1,n))
               smcf(3,2,n)=dimag(cwcf(n))
               smcf(3,3,n)=cwcf(n)
               smcf(3,4,n)=dimag(bwcf(2,-1,n)-bwcf(2,1,n))
               smcf(4,1,n)=bwcf(0,-1,n)+bwcf(0,1,n)
               smcf(4,2,n)=2.*bwcf(2,0,n)
               smcf(4,4,n)=bwcf(0,1,n)-bwcf(0,-1,n)

            enddo
!
!  normalization
!
            qsca0=sm(1,1,0)
            do n=0,nodrw
               do i=1,4
                  do j=1,4
                     sm(i,j,n)=sm(i,j,n)/qsca0
                     smcf(i,j,n)=smcf(i,j,n)/qsca0
                  enddo
               enddo
            enddo
            if(new_rank.eq.0) then
               time2=mytime()-time1
               call timewrite(iunit,' scat matrix coef time:',time2)
            endif
            deallocate(windex,vindex,wvindex,wvnum,dm)
         endif
         call mstm_mpi(mpi_command='barrier')
         end subroutine ranorientscatmatrix
!
!  calculation of the RO scattering matrix from the GSF expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
         subroutine ranorienscatmatrixcalc(ct,smc,nodrexp,sm)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nodrexp,n,nn0,nnp2,nnm2
         real(8) :: smc(4,4,0:nodrexp),sm(4,4),dc(-2:2,0:nodrexp*(nodrexp+2)), &
                    ct
!
!     dc is the normalized generalized spherical function
!     dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
!     where D^k_{mn} is defined in M&M JOSA 96
!
         call rotcoef(ct,2,nodrexp,dc)
         sm=0.d0
         do n=0,nodrexp
            nn0=n*(n+1)
            nnp2=nn0+2
            nnm2=nn0-2
            sm(1,1)=sm(1,1)+dc(0,nn0)*smc(1,1,n)
            sm(1,4)=sm(1,4)+dc(0,nn0)*smc(1,4,n)
            sm(4,1)=sm(4,1)+dc(0,nn0)*smc(4,1,n)
            sm(4,4)=sm(4,4)+dc(0,nn0)*smc(4,4,n)
            if(n.ge.2) then
               sm(1,2)=sm(1,2)+dc(2,nn0)*smc(1,2,n)
               sm(2,4)=sm(2,4)+dc(2,nn0)*smc(2,4,n)
               sm(3,4)=sm(3,4)+dc(2,nn0)*smc(3,4,n)
               sm(1,3)=sm(1,3)+dc(2,nn0)*smc(1,3,n)
               sm(3,1)=sm(1,3)+dc(2,nn0)*smc(3,1,n)
               sm(4,2)=sm(4,2)+dc(2,nn0)*smc(4,2,n)
               sm(2,2)=sm(2,2)+dc(2,nnm2)*smc(2,2,n)+dc(2,nnp2)*smc(3,3,n)
               sm(2,3)=sm(2,3)+dc(2,nnp2)*smc(2,3,n)+dc(2,nnp2)*smc(3,2,n)
               sm(3,3)=sm(3,3)-dc(2,nnm2)*smc(2,2,n)+dc(2,nnp2)*smc(3,3,n)
               sm(3,2)=sm(3,2)+dc(2,nnp2)*smc(2,3,n)-dc(2,nnp2)*smc(3,2,n)
            endif
         enddo
         sm(2,1)=sm(1,2)
         sm(4,3)=-sm(3,4)
!
!  discontiued scaling option: now done in main program
!
!            if(iscale.eq.1) then
!               do j=1,4
!                  do k=j,4
!                     if(j.ne.1.or.k.ne.1) then
!                        sm(j,k,i)=sm(j,k,i)/sm(1,1,i)
!                     endif
!                  enddo
!               enddo
!            endif
!
!    here are the VV and HH differential cross sections
!
!            gvv=.25*(sm(1,1)+sm(2,2)-2.*sm(1,2))
!            ghh=.25*(sm(1,1)+sm(2,2)+2.*sm(1,2))
!
         return
         end subroutine ranorienscatmatrixcalc
      end module scatprops
!
! module nearfield contains local data and subroutines for near field calculation
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  april 2012: lr formulation throughout
!
      module nearfield
      implicit none
      integer, private :: axialinc,ndimpw,nodrpwmax,naeqns
      integer, private, allocatable :: nblk_nf(:),noff_nf(:)
      real(8), private :: rplotmax
      complex(8), allocatable, private :: amnp_nf(:),fmnp_nf(:),pmnp_nf(:), &
               amn3mp_nf(:),fmn3mp_nf(:),pmn3mp_nf(:)
      contains

         subroutine packcoefficient(nodr,ain,aout)
         implicit none
         integer :: nodr,n,p,nn1
         complex(8) :: ain(0:nodr+1,nodr,2),aout(2,nodr*(nodr+2))
         do n=1,nodr
            nn1=n*(n+1)
            do p=1,2
               aout(p,nn1-n:nn1-1)=ain(n+1,n:1:-1,p)
               aout(p,nn1:nn1+n)=ain(0:n,n,p)
            enddo
         enddo
         end subroutine packcoefficient

         subroutine unpackcoefficient(nodr,ain,aout)
         implicit none
         integer :: nodr,n,p,nn1
         complex(8) :: aout(0:nodr+1,nodr,2),ain(2,nodr*(nodr+2))
         do n=1,nodr
            nn1=n*(n+1)
            do p=1,2
               aout(n+1,n:1:-1,p)=ain(p,nn1-n:nn1-1)
               aout(0:n,n,p)=ain(p,nn1:nn1+n)
            enddo
         enddo
         end subroutine unpackcoefficient
!
!  nearfieldspherepart calculates the field at point xg generated by the spheres
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  april 2012: lr formulation: medium is oa by default
!
         subroutine nearfieldspherepart(xg,nsphere,xsp,rpos,ri,&
                    hostsphere,nodr,insphere,efield,hfield)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nsphere,nodr(nsphere),i,insphere, &
                    hostsphere(nsphere),b11,b12,j
         real(8) :: xg(3),xsp(nsphere),rpos(3,nsphere),x(3),r,xspmin
         complex(8) :: ri(2,0:nsphere),efield(3),hfield(3),rib
         complex(8), allocatable :: vwh(:,:)

!
!  find if the point is inside a sphere
!
         insphere=0
         xspmin=1.d10
         do i=1,nsphere
            x=xg(:)-rpos(:,i)
            r=sqrt(dot_product(x,x))
            if(r.le.xsp(i).and.xsp(i).lt.xspmin) then
               xspmin=xsp(i)
               insphere=i
            endif
         enddo
!
!  do the calculations
!
         if(insphere.eq.0) then
!
!  outside all spheres: field = external scattered
!
            efield=0.
            hfield=0.
            do i=1,nsphere
               if(hostsphere(i).eq.0) then
                  allocate(vwh(3,nblk_nf(i)))
                  x=xg(:)-rpos(:,i)
                  call vwhcalc(x,ri(:,0),nodr(i),3,vwh)
                  b11=noff_nf(i)+1
                  b12=noff_nf(i)+nblk_nf(i)
                  efield(:)=efield(:)+matmul(vwh(:,1:nblk_nf(i)),amnp_nf(b11:b12))
                  hfield(:)=hfield(:)+matmul(vwh(:,1:nblk_nf(i)),amn3mp_nf(b11:b12)) &
                     *ri(1,0)/dcmplx(0.d0,1.d0)
                  deallocate(vwh)
               endif
            enddo
         else
!
!  inside a sphere: field = regular + outgoing parts
!
            i=insphere
            rib=2.d0/(1.d0/ri(1,i)+1.d0/ri(2,i))
            allocate(vwh(3,nblk_nf(i)))
            x=xg(:)-rpos(:,i)
            call vwhcalc(x,ri(:,i),nodr(i),1,vwh)
            b11=noff_nf(i)+1
            b12=noff_nf(i)+nblk_nf(i)
            efield(:)=matmul(vwh(:,1:nblk_nf(i)), &
                fmnp_nf(b11:b12))
            hfield(:)=matmul(vwh(:,1:nblk_nf(i)), &
                fmn3mp_nf(b11:b12))*rib/dcmplx(0.d0,1.d0)
            deallocate(vwh)
            do j=1,nsphere
               if(hostsphere(j).eq.i) then
                  allocate(vwh(3,nblk_nf(j)))
                  x=xg(:)-rpos(:,j)
                  call vwhcalc(x,ri(:,i),nodr(j),3,vwh)
                  b11=noff_nf(j)+1
                  b12=noff_nf(j)+nblk_nf(j)
                  efield(:)=efield(:)+matmul(vwh(:,1:nblk_nf(j)), &
                      amnp_nf(b11:b12))
                  hfield(:)=hfield(:)+matmul(vwh(:,1:nblk_nf(j)), &
                      amn3mp_nf(b11:b12))*rib/dcmplx(0.d0,1.d0)
                  deallocate(vwh)
               endif
            enddo
         endif
         end subroutine nearfieldspherepart
!
!  nearfieldincidentpart calculates the incident field at point xg using a regular
!  vswh expansion
!
!
!  last revised: 15 January 2011
!  april 2012: lr formulation.
!
         subroutine nearfieldincidentpart(xg,nodrpw,rimedium,efield,hfield)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nblkpw,nodrpw,nodrc
         real(8) :: xg(3),rdist
         complex(8) :: vwhpw(3,nodrpw*(nodrpw+2)*2),vwhpwaxial(3,4*nodrpw), &
                       efield(3),hfield(3),rimedium(2),rib
         rib=2.d0/(1.d0/rimedium(1)+1.d0/rimedium(2))
         rdist=sqrt(dot_product(xg,xg))
         nodrc=ceiling(rdist+4.*rdist**0.3333+2.)
         nodrc=min(nodrpw,nodrpwmax,nodrc)
!
!  oblique incidence: use the full expansion
!
         if(axialinc.eq.0) then
            call vwhcalc(xg,rimedium,nodrc,1,vwhpw)
            nblkpw=nodrc*(nodrc+2)*2
            efield(:)=matmul(vwhpw(:,1:nblkpw),pmnp_nf(1:nblkpw))
            hfield(:)=matmul(vwhpw(:,1:nblkpw),pmn3mp_nf(1:nblkpw)) &
                    *rib/dcmplx(0.d0,1.d0)
         else
!
!  axial incidence: use the shortcut
!
            call vwhaxialcalc(xg,rimedium,nodrc,1,vwhpwaxial)
            nblkpw=4*nodrc
            efield(:)=matmul(vwhpwaxial(:,1:nblkpw),pmnp_nf(1:nblkpw))
            hfield(:)=matmul(vwhpwaxial(:,1:nblkpw),pmn3mp_nf(1:nblkpw)) &
               *rib/dcmplx(0.d0,1.d0)
         endif
         end subroutine nearfieldincidentpart
!
!  nearfieldincidentcoef generates the reshaped array of incident field coefficients
!
!
!  last revised: 15 January 2011
!  april 2012: lr formulation
!
         subroutine nearfieldincidentcoef(nodrpw,alpha,beta,gamma,cbeam)
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: m,n,p,nn1,mn,mnp,nodrpw,sp
         real (8) :: alpha,beta,cbeam,gamma,cgamma,sgamma
         complex(8), allocatable :: pmnp0(:,:,:,:)
         allocate(pmnp0(0:nodrpw+1,nodrpw,2,2))
         if(allocated(pmnp_nf)) deallocate(pmnp_nf,pmn3mp_nf)
         if(beta.ne.0.d0) then
            axialinc=0
            ndimpw=2*nodrpw*(nodrpw+2)
         else
            axialinc=1
            ndimpw=4*nodrpw
         endif
         allocate(pmnp_nf(ndimpw),pmn3mp_nf(ndimpw))
         if(cbeam.eq.0.d0) then
            call planewavecoef(alpha,beta,nodrpw,pmnp0)
         else
            call gaussianbeamcoef(alpha,beta,cbeam,nodrpw,pmnp0)
         endif
         cgamma=cos(gamma)
         sgamma=sin(gamma)
         if(axialinc.eq.0) then
            do n=1,nodrpw
               nn1=n*(n+1)
               do p=1,2
                  sp=-(-1)**p
                  do m=-n,-1
                     mn=nn1+m
                     mnp=2*(mn-1)+p
                     pmnp_nf(mnp)=pmnp0(n+1,-m,p,1)*cgamma+pmnp0(n+1,-m,p,2)*sgamma
                     pmn3mp_nf(mnp)=sp*(pmnp0(n+1,-m,p,1)*cgamma+pmnp0(n+1,-m,p,2)*sgamma)
                  enddo
                  do m=0,n
                     mn=nn1+m
                     mnp=2*(mn-1)+p
                     pmnp_nf(mnp)=pmnp0(m,n,p,1)*cgamma+pmnp0(m,n,p,2)*sgamma
                     pmn3mp_nf(mnp)=sp*(pmnp0(m,n,p,1)*cgamma+pmnp0(m,n,p,2)*sgamma)
                  enddo
               enddo
            enddo
         else
            do n=1,nodrpw
               do p=1,2
                  sp=-(-1)**p
                  mnp=4*(n-1)+p
                  pmnp_nf(mnp)=pmnp0(n+1,1,p,1)*cgamma+pmnp0(n+1,1,p,2)*sgamma
                  pmn3mp_nf(mnp)=sp*(pmnp0(n+1,1,p,1)*cgamma+pmnp0(n+1,1,p,2)*sgamma)
                  pmnp_nf(mnp+2)=pmnp0(1,n,p,1)*cgamma+pmnp0(1,n,p,2)*sgamma
                  pmn3mp_nf(mnp+2)=sp*(pmnp0(1,n,p,1)*cgamma+pmnp0(1,n,p,2)*sgamma)
               enddo
            enddo
         endif
         deallocate(pmnp0)
         end subroutine nearfieldincidentcoef
!
!  nearfieldpointcalc: if newcalc = 1, generates the reshaped incident, scattered, and
!                      internal field coefficients, and returns with newcalc=0
!                      if newcalc = 0, generates the field at point xg
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  april 2012: lr formulation: oa is implied everywhere
!
         subroutine nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri, &
                    hostsphere,numberfieldexp,amnp,gamma,epspw,xg,newcalc,efield,hfield)
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: nsphere,neqns,nodr(nsphere),i,nodrpw,newcalc, &
                    insphere,hostsphere(nsphere),numberfieldexp(nsphere), &
                    b11,b12,b21,b22,noffi,noi,noi2,j,n,jp1,nbi
         real (8) :: alpha,beta,cbeam,xsp(nsphere),rpos(3,nsphere),xg(3), &
                     gamma,epspw,rplot,cgamma,sgamma
         complex(8) :: amnp(neqns*2),ri(2,0:nsphere),efield(3),hfield(3),einc(3),hinc(3)
         complex(8), allocatable :: amnp1(:,:,:),amnpt(:,:),fmnpt(:,:), &
                     fmnp1(:,:,:)
!
!  initialization operations: newcalc=1
!
         if(newcalc.eq.1) then
            naeqns=0.
            do i=1,nsphere
               naeqns=naeqns+nodr(i)*(nodr(i)+2)*2
            enddo
            if(allocated(amnp_nf)) deallocate(amnp_nf,fmnp_nf,amn3mp_nf,fmn3mp_nf, &
                     noff_nf,nblk_nf)
            allocate(amnp_nf(naeqns),fmnp_nf(naeqns),amn3mp_nf(naeqns),fmn3mp_nf(naeqns), &
                     noff_nf(nsphere),nblk_nf(nsphere))
            noff_nf(1)=0
            do i=1,nsphere
               nblk_nf(i)=2*nodr(i)*(nodr(i)+2)
               if(i.lt.nsphere) noff_nf(i+1)=noff_nf(i)+nblk_nf(i)
            enddo
            cgamma=cos(gamma)
            sgamma=sin(gamma)
            noffi=0
            do i=1,nsphere
               noi=nodr(i)
               noi2=noi+2
               nbi=noi*noi2
               b11=noffi+1
               b12=noffi+nblk_nf(i)
               allocate(amnp1(0:noi+1,noi,2),amnpt(2,nbi), &
                        fmnpt(2,nbi),fmnp1(0:noi+1,noi,2))
               amnp1(0:noi+1,1:noi,1:2) &
                      =cgamma*reshape(amnp(b11:b12),(/noi2,noi,2/)) &
                      +sgamma*reshape(amnp(b11+nblk_nf(i):b12+nblk_nf(i)), &
                          (/noi2,noi,2/))
               if(numberfieldexp(i).eq.1) then
                  call onemiecoeffmult(i,noi,amnp1,fmnp1,'c')
               else
                  b21=b11+nblk_nf(i)*2
                  b22=b12+nblk_nf(i)*2
                  fmnp1(0:noi+1,1:noi,1:2) &
                      =cgamma*reshape(amnp(b21:b22),(/noi2,noi,2/)) &
                      +sgamma*reshape(amnp(b21+nblk_nf(i):b22+nblk_nf(i)), &
                          (/noi2,noi,2/))
               endif
               noffi=noffi+nblk_nf(i)*numberfieldexp(i)*2
               call packcoefficient(noi,amnp1,amnpt)
               call packcoefficient(noi,fmnp1,fmnpt)
               do n=1,nbi
                  j=noff_nf(i)+n+n-1
                  jp1=j+1
                  amnp_nf(j)=amnpt(1,n)
                  fmnp_nf(j)=fmnpt(1,n)
                  amnp_nf(jp1)=amnpt(2,n)
                  fmnp_nf(jp1)=fmnpt(2,n)
                  amn3mp_nf(j)=amnpt(1,n)
                  fmn3mp_nf(j)=fmnpt(1,n)
                  amn3mp_nf(jp1)=-amnpt(2,n)
                  fmn3mp_nf(jp1)=-fmnpt(2,n)
               enddo
               deallocate(amnp1,fmnp1,amnpt,fmnpt)
            enddo
            rplot=sqrt(dot_product(xg,xg))
            rplotmax=rplot
            call planewavetruncationorder(rplot,ri(:,0),epspw,nodrpw)
            nodrpwmax=nodrpw
            call nearfieldincidentcoef(nodrpw,alpha,beta,gamma,cbeam)
            newcalc=0
            return
         endif
!
!  point calculation operations: newcalc=0
!  first determine the required order of the incident field expansion, and recalculate
!  field coefficients, if necessary
!
         rplot=sqrt(dot_product(xg,xg))
         if(rplot.gt.rplotmax) then
            rplotmax=rplot
            call planewavetruncationorder(rplot,ri(:,0),epspw,n)
            if(n.gt.nodrpwmax) then
               nodrpwmax=n
               nodrpw=n
               call nearfieldincidentcoef(nodrpw,alpha,beta,gamma,cbeam)
            endif
         endif
         efield=0.d0
         hfield=0.d0
!
!  calculate the sphere contribution to the field
!
         call nearfieldspherepart(xg,nsphere,xsp,rpos,ri,hostsphere,&
                    nodr,insphere,efield,hfield)
!
!  if the point is external to the spheres, calculate the incident field
!
         if(insphere.eq.0) then
            call nearfieldincidentpart(xg,nodrpwmax,ri(1,0),einc,hinc)
!
! a temporary patch for alpha=beta=0
!
!            einc=(/1.d0,0.d0,0.d0/)*cdexp(dcmplx(0.d0,1.d0)*xg(3))
!            hinc=(/0.d0,1.d0,0.d0/)*cdexp(dcmplx(0.d0,1.d0)*xg(3))
!            einc=(/cos(beta)*cos(alpha + gamma),sin(alpha + gamma), &
!               -(cos(alpha + gamma)*sin(beta))/) &
!               *exp(dcmplx(0.d0,1.d0)*ri(1,0) &
!               *((xg(1)*cos(alpha)+xg(2)*sin(alpha))*sin(beta) &
!               +xg(3)*cos(beta)))
!            hinc=ri(1,0)*(/-((cos(alpha + gamma)*sin(alpha) &
!               *sin(beta)**2 + cos(beta)*sin(alpha + gamma))), &
!               cos(alpha + gamma)*(cos(beta)**2 + cos(alpha)*sin(beta)**2), &
!               sin(beta)*(-(cos(beta)*cos(alpha + gamma)  &
!               *sin(alpha)) + cos(alpha)*sin(alpha + gamma))/) &
!               *exp(dcmplx(0.d0,1.d0)*ri(1,0) &
!               *((xg(1)*cos(alpha)+xg(2)*sin(alpha))*sin(beta) &
!               +xg(3)*cos(beta)))
            efield=efield+einc
            hfield=hfield+hinc
         endif
         end subroutine nearfieldpointcalc
!
!  nearfieldgridcalc is an MPI--enabled subroutine for calculating field points on a
!  rectangular grid.   Writes the data to nfoutunit.
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  changed so that input positions are defined relative to sphere position file origin, and
!  not the gb focal point.
!  february 2013:  reorganized mpi task sharing.
!  optional output arrays
!
         subroutine nearfieldgridcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,  &
                    hostsphere,numberfieldexp,amnp,nfplane,nfplanepos0,nfplanevert0, &
                    gbfocus,deltax,gamma,nfoutunit,epspw,nfoutdata,runprintunit, &
                    ef_array,hf_array)
         use mpidefs
         use mpidata
         use intrinsics
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: nsphere,neqns,nodr(nsphere),nfplane,runprintunit,npoints1,npoints2, &
                    npoints,i,j,k,np23,gcoord(3),rank,numprocs,npoints1by5,nsp, &
                    nfoutunit,nfoutdata,newcalc,hostsphere(nsphere), &
                    numberfieldexp(nsphere),task,nsend,noff,nppp,k1
         real (8) :: alpha,beta,cbeam,xsp(nsphere),rpos(3,nsphere),nfplanepos,&
                     nfplanevert(2,2),xg(3),xgp(3),deltax,gamma,epspw, &
                     time1,xplot(3,nsphere),xi0,ri0,esquare,xgpmax(3),rplot, &
                     gbfocus(3),nfplanepos0,nfplanevert0(2,2)
         complex(8) :: amnp(neqns*2),ri(2,0:nsphere),efield(3),hfield(3)
         complex(8), optional :: ef_array(3,*),hf_array(3,*)
         complex(8), allocatable :: efvec(:,:),hfvec(:,:),efvec0(:,:),hfvec0(:,:)
         logical :: ef_array_present,hf_array_present
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         ef_array_present=present(ef_array)
         hf_array_present=present(hf_array)
!
!  determine the plane
!
         if(nfplane.eq.1) then
            gcoord=(/2,3,1/)
         elseif(nfplane.eq.2) then
            gcoord=(/3,1,2/)
         else
            gcoord=(/1,2,3/)
         endif
!
!  shift the coordinates to gb focal origin
!
         nfplanevert(1,1)=nfplanevert0(1,1)-gbfocus(gcoord(1))
         nfplanevert(1,2)=nfplanevert0(1,2)-gbfocus(gcoord(1))
         nfplanevert(2,1)=nfplanevert0(2,1)-gbfocus(gcoord(2))
         nfplanevert(2,2)=nfplanevert0(2,2)-gbfocus(gcoord(2))
         nfplanepos=nfplanepos0-gbfocus(gcoord(3))
         xg(gcoord(3))=nfplanepos
!
!  determine the number of points
!
         npoints1=nint((nfplanevert(1,2)-nfplanevert(1,1))/deltax)+1
         npoints2=nint((nfplanevert(2,2)-nfplanevert(2,1))/deltax)+1
         npoints=npoints1*npoints2
         nppp=ceiling(dble(npoints)/dble(numprocs))
         allocate(efvec(3,nppp),hfvec(3,nppp))
         if(rank.eq.0) then
            allocate(efvec0(3,nppp*numprocs),hfvec0(3,nppp*numprocs))
         endif
!
!  find the maximum point-to-target origin distance and initialize the field calculation
!
         xgp(3)=nfplanepos
         rplotmax=0.d0
         xgpmax=0.d0
         do i=1,npoints1
            xgp(1)=nfplanevert(1,1)+deltax*dble(i-1)
            do j=1,npoints2
               xgp(2)=nfplanevert(2,1)+deltax*dble(j-1)
               rplot=sqrt(dot_product(xgp,xgp))
               if(rplot.gt.rplotmax) then
                  rplotmax=rplot
                  xgpmax=xgp
               endif
            enddo
         enddo
!
! initialize the calculation: xgpmax is the maximum displacement from origin
!
         newcalc=1
         call nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri, &
                    hostsphere,numberfieldexp,amnp,gamma,epspw,xgpmax,newcalc, &
                    efield,hfield)
!
!  determine the intersecting spheres
!
         nsp=0
         do i=1,nsphere
            xi0=abs(rpos(gcoord(3),i)-xg(gcoord(3)))
            if(xi0.le.xsp(i)) then
               nsp=nsp+1
               xplot(1,nsp)=rpos(gcoord(1),i)+gbfocus(gcoord(1))
               xplot(2,nsp)=rpos(gcoord(2),i)+gbfocus(gcoord(2))
               ri0=xsp(i)*xsp(i)-xi0*xi0
               if(ri0.ne.0.) ri0=sqrt(ri0)
               xplot(3,nsp)=ri0
            endif
         enddo
!
!  report to runprintunit
!
         if(rank.eq.0.and.runprintunit.ne.0) then
            write(runprintunit,'('' near field calculations'')')
            write(runprintunit,'('' plane, position:'',i5,f9.3)') nfplane,nfplanepos0
            write(runprintunit,'('' rectangular plot vertices:'')')
            write(runprintunit,'('' min:'',3f9.3)') nfplanevert0(1:2,1)
            write(runprintunit,'('' max:'',3f9.3)') nfplanevert0(1:2,2)
            write(runprintunit,'('' number of plotting points, step size:'',i8,f8.3)') npoints, deltax
            write(runprintunit,'('' max plane wave order:'',i5)') nodrpwmax
         endif
!
!  determine the distribution of work among the processors
!
         np23=3*npoints2
         npoints1by5=int(npoints1/5.+.5)
!
!  do the calculations and write the results to the file
!
         if(rank.eq.0) then
            if(nfoutunit.ne.0) then
               write(nfoutunit,*) npoints1,npoints2
               write(nfoutunit,*) nsp
               do i=1,nsp
                  write(nfoutunit,'(3e13.5)') xplot(1,i),xplot(2,i),xplot(3,i)
               enddo
            endif
            time1=mytime()
         endif
         xg(gcoord(3))=nfplanepos
         newcalc=0
         efvec=0.d0
         hfvec=0.d0
         task=-1
         k=0
         do i=1,npoints1
            xg(gcoord(1))=nfplanevert(1,1)+deltax*dble(i-1)
            xgp(gcoord(1))=nfplanevert0(1,1)+deltax*dble(i-1)
            do j=1,npoints2
               task=task+1
               if(mod(task,numprocs).eq.rank) then
                  k=k+1
                  xg(gcoord(2))=nfplanevert(2,1)+deltax*dble(j-1)
                  call nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,&
                       xsp,rpos,ri,hostsphere,numberfieldexp,amnp,gamma,epspw, &
                       xg,newcalc,efvec(:,k),hfvec(:,k))
               endif
            enddo
         enddo
         nsend=3*nppp
         call mstm_mpi(mpi_command='gather',mpi_recv_buf_dc=efvec0, &
              mpi_send_buf_dc=efvec(1:3,1:nppp),mpi_number=nsend,mpi_rank=0)
         call mstm_mpi(mpi_command='gather',mpi_recv_buf_dc=hfvec0, &
              mpi_send_buf_dc=hfvec(1:3,1:nppp),mpi_number=nsend,mpi_rank=0)
         if(rank.eq.0) then
            k=0
            k1=0
            task=-1
            do i=1,npoints1
               xgp(gcoord(1))=nfplanevert0(1,1)+deltax*dble(i-1)
               do j=1,npoints2
                  xgp(gcoord(2))=nfplanevert0(2,1)+deltax*dble(j-1)
                  task=task+1
                  noff=mod(task,numprocs)*nppp
                  if(noff.eq.0) then
                     k=k+1
                  endif
                  efield(:)=efvec0(:,noff+k)
                  hfield(:)=hfvec0(:,noff+k)
                  if(nfoutunit.ne.0) then
                     if(nfoutdata.eq.0) then
                        esquare=dot_product(efield(:),efield(:))
                        write(nfoutunit,'(2f9.4,e13.5)')  &
                           xgp(gcoord(1)),xgp(gcoord(2)),esquare
                     elseif(nfoutdata.eq.1) then
                        write(nfoutunit,'(2f9.4,6e13.5)')  &
                           xgp(gcoord(1)),xgp(gcoord(2)),efield(:)
                     else
                        write(nfoutunit,'(2f9.4,12e13.5)') &
                           xgp(gcoord(1)),xgp(gcoord(2)),efield(:),hfield(:)
                     endif
                  endif
                  k1=(j-1)*npoints1+i
                  if(ef_array_present) then
                     ef_array(:,k1)=efield(:)
                  endif
                  if(hf_array_present) then
                     hf_array(:,k1)=hfield(:)
                  endif
               enddo
            enddo
         endif
         end subroutine nearfieldgridcalc
      end module nearfield
!
!  module solver: subroutines for solving interaction equations for fixed orientation
!  and T matrix problems
!
!
!  last revised: 15 January 2011
!                february 2013
!
      module solver
      implicit none

      contains

!
!  tmatrixsoln: calculation of T matrix via solution of interaction equations for
!  a generalized plane wave expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: call for sphereqeff changed
!  april 2012: lr formulation.   required changing rhs.
!  february 2013: completely rewritten.    groups of solutions are calculated simultaneously by
!   all processors for each rhs order.    The processer groups are set by the option mpi comm.
! march 2013: max mb per array optional input arg.   Limit number rhs.
!
         subroutine tmatrixsoln(neqns,nsphere,nodr,ntran,nodrt,xsp,rpos,hostsphere,numberfieldexp, &
                    rimedium,epssoln,epscon,niter,calctmatrix,tmatrixfile,fftranpresent,niterstep, &
                    excitedsphere,qext,qabs,qsca,istat,mpi_comm,max_mb_per_array)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         use miecoefdata
         use spheredata
         use translation
         use scatprops
         implicit none
         integer :: iter,niter,neqns,nsphere,nodr(nsphere),ntran(nsphere),nodrt, &
                    nodrmax,i,istat,m,n,k,l,q,noff,nblk,ma,na,ka,la,iunit, &
                    rank,calctmatrix,lt,kt,qt,nt,mt,it,nodrtt,lstart,numsolns, &
                    lold,lstarta(1),numprocs,istart,task, &
                    niterstep,hostsphere(nsphere),numberfieldexp(nsphere), &
                    excitedsphere,numextspheres,nodrrhs,nstop,i0,ntran0(nsphere),kk1, &
                    nrhs,mpicomm,nsend,nrhsmax,narraymax,maxmbperarray,nbatch,kqstop, &
                    batch,kqstart,kq
         integer, optional :: mpi_comm,max_mb_per_array
         logical :: fftranpresent
         real(8) :: err,qext(nsphere),qabs(nsphere),qsca(nsphere),xsp(nsphere),xv, &
                    rpos(3,nsphere),xij(3),ftemp, &
                    qexttot,qabstot,qscatot,qextold(1),qscaold(1),errqe,errqs, &
                    time1,time2,epssoln,epscon,timeorder,&
                    at1,at2,at3,at4,rpos0(3)
         real(8) :: qextl(nsphere),qabsl(nsphere),qscal(nsphere)
         real(8), allocatable :: qextlt(:,:),qabslt(:,:),qscalt(:,:)
         complex(8) :: rimedium(2)
         complex(8), allocatable :: pmnp0(:,:,:,:),amnp0(:,:,:,:), &
                     amnp(:),pmnp(:),pmnpan(:)
         character*30 :: tmatrixfile
         character*4 :: timeunit
         data istart/1/
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(mpicomm.eq.mpi_comm_null) return
         if(present(max_mb_per_array)) then
            maxmbperarray=max_mb_per_array
         else
            maxmbperarray=1000
         endif
         call getrunparameters(run_print_unit=iunit)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         numextspheres=0
         nodrmax=maxval(nodr)
         if(excitedsphere.eq.0) then
            nodrrhs=nodrt
            rpos0=0.
            xv=0.d0
            do i=1,nsphere
               if(hostsphere(i).eq.0) then
                  xv=xv+xsp(i)**3.d0
                  numextspheres=numextspheres+1
               endif
            enddo
            xv=xv**(1.d0/3.d0)
         else
            nodrrhs=nodr(excitedsphere)
            rpos0=rpos(:,excitedsphere)
            xv=xsp(excitedsphere)
         endif
         qext=0.d0
         qabs=0.d0
         qsca=0.d0
         qextold=0.d0
         qscaold=0.d0
         lstart=1
!
!  perform T matrix file operations as needed
!
         if(rank.eq.0) then
            if(calctmatrix.eq.1) then
               open(3,file=tmatrixfile)
               write(3,'(3i4)') excitedsphere,nodrt,nodrrhs
               write(3,'(i6,e13.5)') nsphere,xv
            else
               open(3,file=tmatrixfile)
               if(rank.eq.0) write(iunit,'('' finding end of record to file '',a)') &
                  tmatrixfile
               read(3,*) i0,nodrtt,i0
               read(3,*) i0,ftemp
               do l=1,nodrt
                  if(excitedsphere.eq.0) then
                     nstop=l
                  else
                     nstop=nodrtt
                  endif
                  do k=-l,l
                     do q=1,2
                        read(3,'(3i5)',end=20,err=20) lt,kt,qt
                        do n=1,nstop
                           do m=-n,n
                              read(3,'(2i5,4e17.9)',end=20,err=20) nt,mt,at1,at2,at3,at4
                           enddo
                        enddo
                     enddo
                  enddo
                  do i=1,nsphere
                     read(3,'(i5,3e17.9)',end=20,err=20) it,qextl(i),qabsl(i),qscal(i)
                  enddo
                  qext=qext+qextl
                  qabs=qabs+qabsl
                  qsca=qsca+qscal
               enddo
20             close(3)
               open(3,file=tmatrixfile)
               qextold(1)=0.d0
               qabstot=0.d0
               do i=1,nsphere
                  qextold=qextold+qext(i)*xsp(i)*xsp(i)/xv/xv
                  qabstot=qabstot+qabs(i)*xsp(i)*xsp(i)/xv/xv
               enddo
               qscaold(1)=qextold(1)-qabstot
               lstart=lt
               read(3,*) i0,nodrtt,i0
               read(3,*) i0,ftemp
               if(rank.eq.0) then
                  write(iunit,'('' calculations begin with order '',i5)') lstart
                  call flush(iunit)
               endif
               do l=1,lstart-1
                  if(excitedsphere.eq.0) then
                     nstop=l
                  else
                     nstop=nodrtt
                  endif
                  do k=-l,l
                     do q=1,2
                        read(3,'(3i5)') lt,kt,qt
                        do n=1,nstop
                           do m=-n,n
                              read(3,'(2i5,4e17.9)') nt,mt,at1,at2,at3,at4
                           enddo
                        enddo
                     enddo
                  enddo
                  do i=1,nsphere
                     read(3,'(i5,3e17.9)') it,at1,at2,at3
                  enddo
               enddo
            endif
         endif
         if(calctmatrix.ne.1) then
            lstarta(1)=lstart
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_i=lstarta,mpi_number=1, &
                mpi_rank=0,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qextold,mpi_number=1, &
                mpi_rank=0,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qscaold,mpi_number=1, &
                mpi_rank=0,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qext,mpi_number=nsphere, &
                mpi_rank=0,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qabs,mpi_number=nsphere, &
                mpi_rank=0,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=qsca,mpi_number=nsphere, &
                mpi_rank=0,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
            lstart=lstarta(1)
         endif
         numsolns=2*nodrrhs*(nodrrhs+2)
         qextl=0.d0
         qabsl=0.d0
         qscal=0.d0
         lold=0
!
!  begin the loop over RHS of the interaction equations.
!
         do l=lstart,nodrt
            if(rank.eq.0) timeorder=mytime()
            nrhsmax=2*(2*l+1)
            narraymax=nrhsmax*neqns
            nbatch=ceiling(dble(narraymax)*16.*dble(numprocs)/(dble(maxmbperarray)*1.d6))
            kqstop=0
            qextl=0.
            qabsl=0.
            qscal=0.
            do batch=1,nbatch
               kqstart=kqstop+1
               if(kqstart.gt.nrhsmax) then
                  cycle
               endif
               kqstop=nint(batch*nrhsmax/dble(nbatch))
               kqstop=min(kqstop,nrhsmax)
               kqstop=max(kqstop,kqstart)
               nrhs=kqstop-kqstart+1
               allocate(pmnp0(0:l+1,l,2,nrhs),amnp(neqns*nrhs), &
                  pmnp(neqns*nrhs),pmnpan(neqns*nrhs), &
                  qextlt(nsphere,nrhs),qabslt(nsphere,nrhs), &
                  qscalt(nsphere,nrhs))
!
!  evaluate the nrhs for the current l, and store in the rhs vector
!
               do kq=kqstart,kqstop
                  q=mod(kq-1,2)+1
                  k=-l+floor((kq-1)/2.d0)
                  if(k.le.-1) then
                     ka=l+1
                     la=-k
                  else
                     ka=k
                     la=l
                  endif
                  kk1=kq-kqstart+1
                  pmnp0(:,:,:,kk1)=0.d0
!
!  the t matrix is te-tm based; hence the following two lines
!
                  pmnp0(ka,la,1,kk1)=.5d0
                  pmnp0(ka,la,2,kk1)=-.5d0*(-1)**q
               enddo
!
!  calculate the sphere-based rhs vectors
!
               noff=0
               pmnp=0.d0
               pmnpan=0.d0
               task=0
               do i=1,nsphere
                  nblk=nodr(i)*(nodr(i)+2)*2*nrhs
                  if(hostsphere(i).eq.0.and. &
                     (excitedsphere.eq.0.or.excitedsphere.eq.i)) then
                     task=task+1
                     if(mod(task,numprocs).eq.rank) then
                        xij=rpos(:,i)-rpos0
                        call rottran(nodrj=l,nodri=nodr(i),translation_vector=xij, &
                           refractive_index=rimedium,vswf_type=1, &
                           cxj=pmnp0,cyi=pmnp(noff+1:noff+nblk),number_rhs=nrhs)
                     endif
                  endif
                  noff=noff+nblk*numberfieldexp(i)
               enddo
               deallocate(pmnp0)
               if(numprocs.gt.1) then
                  nsend=neqns*nrhs
                  call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=pmnp, &
                       mpi_number=nsend,mpi_operation=mstm_mpi_sum)
               endif
               call multmiecoeffmult(neqns,nrhs,1,pmnp,pmnpan)
               amnp=pmnpan
!
!  call the solver
!
               if(fftranpresent) then
                  call cbicgff(neqns,nsphere,niter,epssoln,pmnpan,amnp,0, &
                       niterstep,iter,err,number_rhs=nrhs,mpi_comm=mpicomm)
               else
                  call cbicg(neqns,nsphere,niter,epssoln,pmnpan,amnp,0, &
                       iter,err,number_rhs=nrhs,mpi_comm=mpicomm)
               endif
               if(iter.gt.niter.or.err.gt.epssoln) istat=1
!
! calculate the efficiency factors for the l column order
!
               call qefficiencyfactors(nsphere,neqns,nodr,1,xsp,hostsphere, &
                      numberfieldexp,amnp,pmnp,qextlt,qabslt,qscalt, &
                      number_rhs=nrhs,mpi_comm=mpicomm)
               do i=1,nsphere
                  qextl(i)=qextl(i)+sum(qextlt(i,1:nrhs))
                  qabsl(i)=qabsl(i)+sum(qabslt(i,1:nrhs))
                  qscal(i)=qscal(i)+sum(qscalt(i,1:nrhs))
               enddo
               qext=qext+qextl
               qabs=qabs+qabsl
               qsca=qsca+qscal
!
!  compute the target-based expansion
!
               if(excitedsphere.eq.0) then
                  ntran0=l
                  ntran0=min(ntran0,ntran)
                  nstop=l
               else
                  ntran0=ntran
                  nstop=nodrt
               endif
               deallocate(pmnp,pmnpan,qextlt,qabslt,qscalt)
               allocate(amnp0(0:nstop+1,nstop,2,nrhs))
               amnp0=0.d0
               call amncommonorigin(nsphere,nodr,ntran0,nstop,rpos, &
                      hostsphere,numberfieldexp,rimedium,amnp,amnp0,&
                      number_rhs=nrhs,mpi_comm=mpicomm)
!
!  write results, check for convergence
!
               if(rank.eq.0) then
                  do kq=kqstart,kqstop
                     q=mod(kq-1,2)+1
                     k=-l+floor((kq-1)/2.d0)
                     kk1=kq-kqstart+1
                     write(3,'(3i5)') l,k,q
                     do n=1,nstop
                        do m=-n,n
                           if(m.le.-1) then
                              ma=n+1
                              na=-m
                           else
                              ma=m
                              na=n
                           endif
                           write(3,'(2i5,4e17.9)') n,m,amnp0(ma,na,1,kk1), &
                                                   amnp0(ma,na,2,kk1)
                        enddo
                     enddo
                  enddo
                  if(batch.eq.nbatch) then
                     if(istart.eq.1) then
                        time1=(mytime()-timeorder)/dble(nrhsmax)
                        call timewrite(iunit,' time per solution:',time1)
                        time2=time1*dble(numsolns-2*lstart*(lstart+2))
                        call timewrite(iunit,' estimated t matrix calculation time:',time2)
                        write(iunit,'(''  n   # its  qext         qabs'',&
                            &''         qsca      error     est. time rem.'')')
                        call flush(iunit)
                        istart=0
                     endif
                     timeorder=(mytime()-timeorder)/dble(nrhsmax)
                     time2=timeorder*dble(numsolns-2*l*(l+2))
                     if(time2.gt.3600.d0) then
                        time2=time2/3600.d0
                        timeunit=' hrs'
                     elseif(time2.gt.60.d0) then
                        time2=time2/60.d0
                        timeunit=' min'
                     else
                        timeunit=' sec'
                     endif
                  endif
               endif
               deallocate(amnp0,amnp)
!
! end of batch loop
!
            enddo

            qexttot=0.d0
            qabstot=0.d0
            do i=1,nsphere
               qexttot=qexttot+qext(i)*xsp(i)*xsp(i)/xv/xv
               qabstot=qabstot+qabs(i)*xsp(i)*xsp(i)/xv/xv
               if(rank.eq.0) then
                  write(3,'(i5,3e17.9)') i,qextl(i),qabsl(i),qscal(i)
               endif
            enddo
            qscatot=qexttot-qabstot
            if(excitedsphere.ne.0) qscatot=qext(excitedsphere)
            errqe=qexttot-qextold(1)
            errqs=qscatot-qscaold(1)
            err=max(errqe,errqs)
            if(rank.eq.0) then
               write(iunit,'(i4,i5,4e13.5,f8.2,a4)') l,iter,qexttot,qabstot, &
                  qscatot,err,time2,timeunit
               call flush(iunit)
            endif
            qextold(1)=qexttot
            qscaold(1)=qscatot
            if(err.le.epscon.and.excitedsphere.eq.0) then
               exit
            endif
         enddo
!
!  solution has converged
!
         if(rank.eq.0) then
            close(3)
            if(l.lt.nodrt) then
               nodrt=l
               if(rank.eq.0) then
                  write(iunit,'('' T matrix converged, order:'',i5)') nodrt
                  call flush(iunit)
               endif
               open(3,file=tmatrixfile,form='formatted',access='direct',recl=12)
               write(3,'(3i4)',rec=1) excitedsphere,nodrt,nodrt
               close(3)
            elseif(excitedsphere.eq.0) then
               if(rank.eq.0) then
                  write(iunit,'('' T matrix did not converge to set epsilon'')')
                  call flush(iunit)
               endif
            endif
         endif
         lstarta(1)=nodrt
         call mstm_mpi(mpi_command='bcast',mpi_send_buf_i=lstarta,mpi_number=1, &
         mpi_rank=0,mpi_comm=mpicomm)
         nodrt=lstarta(1)
         end subroutine tmatrixsoln
!
!  solution of interaction equations for a fixed orientation
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: modification of efficiency calculation, to calculate
!           polarized components
!  30 March 2011: took out gbfocus argument: this is not needed since positions are defined
!  relative to the gb focus.
!  20 April 2011: used 2-group MPI formulation
!  October 2011: adapted to far field approximation.
!  December 2011: changed efficiency factor calculation, adapted to generalized sphere
!                 configuration.
! february 2013: number rhs and mpi comm options added, completely rewritten.
!
!
         subroutine fixedorsoln(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,&
                    hostsphere,numberfieldexp,rimed,eps,epstran,niter,amnp,qext,qabs,qsca, &
                    maxerr,maxiter,iterwrite,fftranpresent,niterstep,istat, &
                    far_field_coefficient,mpi_comm)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         use miecoefdata
         use translation
         use scatprops
         implicit none
         integer :: iter,niter,neqns,nodrmax,nsphere,istat,rank,maxiter,iterwrite,&
                    nodr(nsphere),numprocs,niterstep,hostsphere(nsphere), &
                    numberfieldexp(nsphere),mpicomm
         integer, optional :: mpi_comm
         logical :: fftranpresent
         real(8) :: alpha,beta,eps,err,qext(nsphere,3),maxerr,&
                    qabs(nsphere,3),qsca(nsphere,3),cbeam,gbfocus(3),epstran
         real(8) :: xsp(nsphere), rpos(3,nsphere)
         complex(8) :: amnp(neqns*2),rimed(2)
         complex(8), optional :: far_field_coefficient(neqns*2)
         complex(8), allocatable :: pmnp(:),pmnpan(:),pmnp0(:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(mpicomm.eq.mpi_comm_null) return
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nodrmax=maxval(nodr)
         allocate(pmnp(neqns*2),pmnp0(neqns*2))
         gbfocus=0.d0
         if(cbeam.eq.0.d0) then
            call sphereplanewavecoef(nsphere,neqns,nodr,nodrmax,alpha,beta,rpos, &
                  hostsphere,numberfieldexp,rimed,pmnp0)
         else
            call spheregaussianbeamcoef(nsphere,neqns,nodr,alpha,beta,cbeam, &
                    rpos,hostsphere,numberfieldexp,rimed,gbfocus,epstran, &
                    pmnp0,mpi_comm=mpicomm)
         endif
         istat=0
         maxiter=0
         maxerr=0.
         if(present(far_field_coefficient)) then
            pmnp=pmnp0+far_field_coefficient
         else
            pmnp=pmnp0
         endif
!
!  calculate the two solutions
!
         allocate(pmnpan(neqns*2))
         call multmiecoeffmult(neqns,2,1,pmnp,pmnpan)
         amnp=pmnpan
         if(fftranpresent) then
            call cbicgff(neqns,nsphere,niter,eps,pmnpan,amnp,iterwrite, &
                         niterstep,iter,err,number_rhs=2,mpi_comm=mpicomm)
         else
            if(niter.ne.0) then
               call cbicg(neqns,nsphere,niter,eps,pmnpan,amnp,iterwrite, &
                         iter,err,number_rhs=2,mpi_comm=mpicomm)
            endif
         endif
         maxiter=max(iter,maxiter)
         maxerr=max(err,maxerr)
         if(iter.gt.niter.or.err.gt.eps) istat=1
         deallocate(pmnpan)
!
!  efficiency factor calculations
!
         call qefficiencyfactors(nsphere,neqns,nodr,2,xsp,hostsphere, &
                    numberfieldexp,amnp,pmnp,qext,qabs,qsca, &
                    mpi_comm=mpicomm)
         deallocate(pmnp)
         end subroutine fixedorsoln
!
! hybrid bcgm, using far field translation
! november 2011
! february 2013: number rhs option added, completely rewritten.
!
!
         subroutine cbicgff(neqns,nsphere,niter,eps,pnp,anp,iterwrite,&
                    niterstep,iter,errmax,number_rhs,mpi_comm)
         use mpidefs
         use mpidata
         use intrinsics
         use spheredata
         use miecoefdata
         use numconstants
         use specialfuncs
         use translation
         implicit none
         integer :: neqns,niter,iter,nsphere,writetime,&
                    rank,iunit,iterwrite,ntot,numprocs,itermax, &
                    istep,niterstep,nrhs,i,neval,mpicomm
         integer, optional :: number_rhs,mpi_comm
         real(8) :: eps,epsstep,errstep,errmin,errmax
         real(8), allocatable :: err(:),enorm(:)
         complex(8)  :: pnp(*),anp(*)
         complex(8), allocatable :: gnp(:),gnpold(:),pgnp(:)
         logical, allocatable :: rhslist(:)
         data writetime/0/
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         ntot=neqns*nrhs
         if(allocated(gnp)) deallocate(gnp,gnpold,pgnp,rhslist,err,enorm)
         allocate(gnp(neqns*nrhs),gnpold(neqns*nrhs),pgnp(neqns*nrhs),&
                  rhslist(nrhs),err(nrhs),enorm(nrhs))
         if(rank.eq.0) then
            call getrunparameters(run_print_unit=iunit)
         endif
         rhslist=.true.
         call dotproduct(neqns,nrhs,pnp,enorm)
         do i=1,nrhs
            if(enorm(i).eq.0.d0) rhslist(i)=.false.
         enddo
         if(.not.any(rhslist)) return

         err=0.d0
         iter=0
         gnpold=0.d0
         gnp=0.d0
         call  farfieldtranslationerror(neqns,anp(1:ntot),gnp,&
               number_rhs=nrhs,rhs_list=rhslist,mpi_comm=mpicomm)
         call multmiecoeffmult(neqns,nrhs,1,gnp,gnp,rhs_list=rhslist)
         iter=0
         epsstep=.01
         istep=0
         do
            istep=istep+1
            gnpold=gnp
            pgnp=pnp(1:ntot)+gnp
            call cbicg(neqns,nsphere,niterstep,epsstep,pgnp,anp(1:ntot),0,itermax,&
                 errstep,number_rhs=nrhs,rhs_list=rhslist,mpi_comm=mpicomm)
            iter=iter+min(itermax,niterstep)
            gnp=0.d0
            call  farfieldtranslationerror(neqns,anp(1:ntot),gnp, &
                  number_rhs=nrhs,rhs_list=rhslist,mpi_comm=mpicomm)
            call multmiecoeffmult(neqns,nrhs,1,gnp,gnp,rhs_list=rhslist)
            pgnp=gnp-gnpold
            call dotproduct(neqns,nrhs,pgnp,err)
            do i=1,nrhs
               if(rhslist(i)) then
                  err(i)=err(i)/enorm(i)
               endif
            enddo
            errmin=minval(err,mask=rhslist)
            errmax=maxval(err,mask=rhslist)
            neval=sum((/(1,i=1,nrhs)/),mask=rhslist)
            if(rank.eq.0.and.iterwrite.eq.1) then
               write(iunit,'(3i5,3e13.5)') &
                  istep,iter,neval,errstep,errmin,errmax
               call flush(iunit)
            endif
            epsstep=max(errmax,eps)
            if((errmax.lt.eps.and.epsstep.le.eps).or.iter.gt.niter) exit
         enddo
         end subroutine cbicgff
!
! iteration solver
! generalized complex biconjugate gradient method
! original code: Piotr Flatau, although not much remains.
! specialized to the multiple sphere problem
!
!
!  last revised: 15 January 2011
!  october 2011: translation calls modified
!  february 2013: number rhs option added, completely rewritten.
!
         subroutine cbicg(neqns,nsphere,niter,eps,pnp,anp,iterwrite,iter,errmax,&
                     number_rhs,rhs_list,mpi_comm)
         use mpidefs
         use mpidata
         use intrinsics
         use spheredata
         use miecoefdata
         use numconstants
         use specialfuncs
         use translation
         implicit none
         integer :: neqns,niter,iter,nsphere,writetime,nodr(nsphere),&
                    rank,iunit,iterwrite,nblk,nodri,j, &
                    numprocs,nrhs,i,neval,ntot,numberfieldexp(nsphere),&
                    n1,n2,k,mpicomm
         integer, optional :: number_rhs,mpi_comm
         real(8) :: eps,time1,time2,errmin,errmax
         real(8), allocatable :: err(:),enorm(:)
         complex(8)  :: pnp(*),anp(*)
         complex(8), allocatable :: cr(:),cp(:),cw(:),cq(:),cap(:),caw(:), &
                       capt(:),cawt(:),cak(:),csk(:),cbk(:),csk2(:)
         logical, allocatable :: rhslist(:)
         logical, optional :: rhs_list(*)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         ntot=neqns*nrhs
         if(allocated(cr)) deallocate(cr,cp,cw,cq,cap,caw,capt,cawt, &
                      err,enorm,cak,csk,cbk,csk2,rhslist)
         allocate(cr(ntot),cp(ntot),cw(ntot),cq(ntot), &
               cap(ntot),caw(ntot),capt(ntot), &
               cawt(ntot),err(nrhs),enorm(nrhs),cak(nrhs),csk(nrhs), &
               cbk(nrhs),csk2(nrhs),rhslist(nrhs))
         data writetime/0/
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call getmiedataall(sphere_order=nodr,number_field_expansions=numberfieldexp)
         if(rank.eq.0) then
            call getrunparameters(run_print_unit=iunit)
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list(1:nrhs)
         else
            rhslist=.true.
         endif
         iter=0
         errmax=0.
         call dotproduct(neqns,nrhs,pnp,enorm)
         do i=1,nrhs
            if(enorm(i).eq.0.d0) rhslist(i)=.false.
         enddo
         if(.not.any(rhslist)) return
!
!  setting niter < 0 runs the following simple order--of--scattering solution
!
         if(niter.lt.0) then
            cp=anp(1:ntot)
            do iter=1,-niter
               cr=0.
               if(iter.eq.1) then
                  call sphereinteraction(neqns,nrhs,cp,cr, &
                       initial_run=.true.,rhs_list=rhslist, &
                       mpi_comm=mpicomm)
               else
                  call sphereinteraction(neqns,nrhs,cp,cr, &
                       rhs_list=rhslist,mpi_comm=mpicomm)
               endif
               anp(1:ntot)=anp(1:ntot)+cr(1:ntot)
               cp=cr
               call dotproduct(neqns,nrhs,cr,err)
               do i=1,nrhs
                  if(rhslist(i)) then
                    err(i)=err(i)/enorm(i)
                    if(err(i).lt. eps) then
                       rhslist(i)=.false.
                    endif
                  endif
               enddo
               if(.not.any(rhslist)) then
                  exit
               endif
               if(rank.eq.0.and.iter.eq.1.and.writetime.eq.0) then
                  time2=mytime()-time1
                  call timewrite(iunit,' time per iteration:',time2)
                  writetime=1
               endif
               if(rank.eq.0.and.iterwrite.eq.1) then
                  errmin=minval(err,mask=rhslist)
                  errmax=maxval(err,mask=rhslist)
                  neval=sum((/(1,i=1,nrhs)/),mask=rhslist)
                  write(iunit,'('' iter,nrhs,errmin,errmax:'',2i5,2e13.5)') &
                     iter,neval,errmin,errmax
                  call flush(iunit)
               endif
            enddo
            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
            return
         endif
!
! the following is the implementation of the complex biconjugate gradient
! iteration scheme
!
         cr=0.d0
         call sphereinteraction(neqns,nrhs,anp(1:ntot),cr,initial_run=.true., &
              rhs_list=rhslist,mpi_comm=mpicomm)
         n1=1
         csk=0.
         do i=1,nsphere
            nodri=nodr(i)
            nblk=2*nodri*(nodri+2)
            do j=1,numberfieldexp(i)
               do k=1,nrhs
                  if(rhslist(k)) then
                     n2=n1+nblk-1
                     cr(n1:n2)=pnp(n1:n2)-anp(n1:n2)+cr(n1:n2)
                     cq(n1:n2)=conjg(cr(n1:n2))
                     cw(n1:n2)=cq(n1:n2)
                     cp(n1:n2)=cr(n1:n2)
                     csk(k)=csk(k)+dot_product(conjg(cr(n1:n2)),cr(n1:n2))
                  endif
                  n1=n1+nblk
               enddo
            enddo
         enddo
         if(maxval(cdabs(csk)).eq.0.d0) return
!
!  here starts the main iteration loop
!
         do iter=1,niter
            n1=1
            do i=1,nsphere
               nodri=nodr(i)
               nblk=2*nodri*(nodri+2)
               do j=1,numberfieldexp(i)
                  do k=1,nrhs
                     if(rhslist(k)) then
                        n2=n1+nblk-1
                        cak(k)=(0.d0,0.d0)
                        cawt(n1:n2)=(0.d0,0.d0)
                        capt(n1:n2)=(0.d0,0.d0)
                        cap(n1:n2)=0.d0
                        caw(n1:n2)=0.d0
                     endif
                     n1=n1+nblk
                  enddo
               enddo
            enddo
            if(rank.eq.0) then
               if(writetime.eq.0) time1=mytime()
            endif
            call sphereinteraction(neqns,nrhs,cp,cap,cw,caw, &
                 rhs_list=rhslist,mpi_comm=mpicomm)
            n1=1
            do i=1,nsphere
               nodri=nodr(i)
               nblk=2*nodri*(nodri+2)
               do j=1,numberfieldexp(i)
                  do k=1,nrhs
                     if(rhslist(k)) then
                        n2=n1+nblk-1
                        cap(n1:n2)=cp(n1:n2)-cap(n1:n2)
                        caw(n1:n2)=cw(n1:n2)-caw(n1:n2)
                        cak(k)=cak(k)+dot_product(cw(n1:n2),cap(n1:n2))
                     endif
                     n1=n1+nblk
                  enddo
               enddo
            enddo
            do k=1,nrhs
               if(rhslist(k)) then
                  if(cdabs(cak(k)).eq.0.d0) then
                     rhslist(k)=.false.
                  else
                     cak(k)=csk(k)/cak(k)
                  endif
               endif
            enddo
            n1=1
            err=0.
            csk2=0.
            do i=1,nsphere
               nodri=nodr(i)
               nblk=2*nodri*(nodri+2)
               do j=1,numberfieldexp(i)
                  do k=1,nrhs
                     if(rhslist(k)) then
                        n2=n1+nblk-1
                        anp(n1:n2)=anp(n1:n2)+cak(k)*cp(n1:n2)
                        cr(n1:n2)=cr(n1:n2)-cak(k)*cap(n1:n2)
                        cq(n1:n2)=cq(n1:n2)-conjg(cak(k))*caw(n1:n2)
                        csk2(k)=csk2(k)+dot_product(cq(n1:n2),cr(n1:n2))
                        err(k)=err(k)+dot_product(cr(n1:n2),cr(n1:n2))
                     endif
                     n1=n1+nblk
                  enddo
               enddo
            enddo
            do k=1,nrhs
               if(rhslist(k)) then
                  err(k)=err(k)/enorm(k)
                  if(err(k).lt. eps.or.cdabs(csk(k)).eq.0.d0) then
                     rhslist(k)=.false.
                  else
                     cbk(k)=csk2(k)/csk(k)
                     csk(k)=csk2(k)
                  endif
               endif
            enddo
            if(.not.any(rhslist)) then
               return
            endif
            n1=1
            do i=1,nsphere
               nodri=nodr(i)
               nblk=2*nodri*(nodri+2)
               do j=1,numberfieldexp(i)
                  do k=1,nrhs
                     if(rhslist(k)) then
                        n2=n1+nblk-1
                        cp(n1:n2)=cr(n1:n2)+cbk(k)*cp(n1:n2)
                        cw(n1:n2)=cq(n1:n2)+conjg(cbk(k))*cw(n1:n2)
                     endif
                     n1=n1+nblk
                  enddo
               enddo
            enddo
            if(rank.eq.0.and.iter.eq.1.and.writetime.eq.0) then
               time2=mytime()-time1
               call timewrite(iunit,' time per iteration:',time2)
               writetime=1
            endif
            errmin=minval(err,mask=rhslist)
            errmax=maxval(err,mask=rhslist)
            neval=sum((/(1,i=1,nrhs)/),mask=rhslist)
            if(rank.eq.0.and.iterwrite.eq.1) then
               write(iunit,'('' iter,nrhs,errmin,errmax:'',2i5,2e13.5)') &
                  iter,neval,errmin,errmax
               call flush(iunit)
            endif
         enddo
         end subroutine cbicg

      end module solver
