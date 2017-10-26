      subroutine airesint (int1, int2, int3, t0, sig, omega,
     &                 phase, intens, freq,gam, t, delt,tw)
c
      double precision int1, int2, int3, t0, sig, omega
      double precision phase, intens, freq, t, delt,gam
      double precision env, tmp, tw
c
c      write(*,*)t,env(t0, sig,omega,intens,freq,gam,phase,t,tw)
      tmp = int1
      int1 = int1 -
     &(env(t0, sig,omega,intens,freq,gam,phase,t-delt,tw)+
     &env(t0, sig, omega, intens, freq,gam, phase, t,tw))
     &        * 0.25d0 * delt
      int2 = int2 + (tmp + int1) * 0.5d0 * delt
      int3 = int3 + (tmp**2 + int1**2) * 0.5d0 * delt
      return
      end      
      double precision function env (t0, sig, omega,
     &                 intens, freq,gam, phase, tim,tp)
c
      double precision t0, sig, omega, intens, freq, phase, tim
      double precision field1, pi,gam,tp,conv
csig * 0.600561204d0

      conv= sig * 0.600561204d0

c    tw= (t0+336(impulsion complete))
      pi = 3.141592654d0
      

c      if ((tim-t0).le.tp) then
c          env=dsin(pi*(tim-t0)/tp)**2*(-field1(t0, sig, freq, intens,
c     &    tim,phase,betatp))
c      else
c          env=0.d0
c      end if
c      write(*,*)t0, sig, freq, intens,
c     &    tim,phase,beta
c      write(*,*)tim,"TEST"
c      if (tim-t0.gt.12000)then
c        write(*,*)"Here is the problem"
c        write(*,*) dexp(-((tim-tp)/conv)*
c     &                ((tim-tp)/conv))*intens
c      endif
         write(*,*)phase
      if ((tim-t0).le.tp) then
c          write(*,*)"In the max ",tim
          env=(-intens * 
     &(dcos(freq * (tim-t0))+gam*dcos(2.d0*freq*(tim-t0)+pi*phase)))


c	  field1(t0, sig, freq, intens,
c     &    tim,phase,tp)
      else  
c          write(*,*)"In the expo ",tim
          env= dexp(-((tim-tp)/conv)*
     &                ((tim-tp)/conv))*
c     &                intens 
c     &                (field1(t0, sig, freq, intens,
c     &    tim,phase,tp))
     &    ( -intens * 
     & (dcos(freq * (tim-t0))+gam*dcos(2.d0*freq*(tim-t0)+pi*phase)))



c      else if ((tim-t0).gt.(12000.d0)) then 
c         write(*,*)tim,"In the 0 "
c         env=0.d0
      endif
c          write(*,*)freq,tim,env
        
c      write(*,*)tim

      return
      end
      
      double precision function field1(t0,sig,freq,intens,temps,
     &                          phase,tp)
c
      double precision t0, sig, freq, intens, temps
      double precision conv,phase,pi,tp
c
c      write(*,*)"test",freq,temps,t0,phase,tp
      pi = 3.141592654d0
      conv = sig * 0.600561204d0
      field1 = -intens * dcos( freq*((temps-t0)-pi*phase) )
c         field1 = intens
c      write(*,*)"field1 ",temps,field1
      
      return
      end

c      double precision function env (t0, sig, omega,
c     &                 intens, freq,gam, phase, tim,tw)
c
c      double precision t0, sig, omega, intens, freq, phase, tim
c      double precision field1, pi,gam,tw
c
c    tw= (t0+336(impulsion complete))
c      pi = 3.141592654d0
      
c       if(((tim-t0).ge.(t0-tw)).and.(tim.le.tw))then
c      env = -field1(t0, sig, omega, intens, tim)*
c     &(dcos(freq * (tim-t0))+gam*dcos(2.d0*freq*(tim-t0)+pi/phase))
c       else
c          env=0.d0
c      end if
c      return
c      end

c      double precision function field1 (t0, sig, omega, intens, temps)
c
c      double precision t0, sig, omega, intens, temps
c      double precision conv
cc
c      conv = sig * 0.600561204d0
c      
c         field1 = intens * dcos(omega*temps)
c      
c      return
c      end




      subroutine asympt(t, tf, zpsik1, zpsik2, zfi1, zfi2,
     &                  zcutA, npbig, npun, xmu,
     &                  int1tf, int2tf, int3tf,
     &                  int1t0, int2t0, int3t0, dr,
     &                 zwork1, zwork2, zwork3, 
     &                  worka, workb,tablea, tableb)
c

      integer npbig, npun
      double complex zpsik1(npbig), zpsik2(npbig), zcutA(npbig)
      double complex zfi1(npun), zfi2(npun), zwork3(npbig)

      double precision int1tf, int2tf, int3tf
      double precision int1t0, int2t0, int3t0, dr
      double complex zwork1(npbig), zwork2(npbig), zwork(npbig)
c
      double precision work2(4*4096), xnorm1, xnorm2
c
      double precision xmu
c      double precision work(4*npbig+15)
      double precision cte, phase, phase1, phase2, phase3
      double precision dk, A, pi, rk, t, ctet, tf, scale
      integer iA, ll, l, i, isign
      double complex zexpph, zexph1, zexph2, zexph3, ztmp1, ztmp2
      double complex zsqr2
      double precision worka(2*npun), workb(2*npbig)
      double precision tablea(2*npun+30), tableb(2*npbig+30)   
c
      pi = 3.141592654d0
      cte = 0.5d0/xmu
c      write(*,*)"int",int1tf,int1t0
      phase1 = int1tf - int1t0
      phase2 = int1tf*(tf-t) - int2tf + int2t0
      phase3 = (tf-t)*int1tf**2 - 2.d0*int1tf*(int2tf-int2t0)
     &         + int3tf -int3t0
c
      zexph3 = cdexp(dcmplx(0.d0,-cte*phase3))
c
      dk = 2.d0*pi/(npbig*dr)
      A = phase1/dk
      iA = idnint(A)

c
      call ZVEA(npun,zfi1(1),1,zfi2(1),1,zwork1(1),1)
      call ZVES(npun,zfi1(1),1,zfi2(1),1,zwork2(1),1)
c
      do i = npun+1, npbig
         zwork1(i) = dcmplx(0.d0,0.d0)
         zwork2(i) = dcmplx(0.d0,0.d0)
      end do
c
      call ZVEM(npbig,zcutA(1),1,zwork1(1),1,zwork1(1),1)
      call ZVEM(npbig,zcutA(1),1,zwork2(1),1,zwork2(1),1)
c
      scale = dr/dsqrt(2.d0*pi)
      isign = -1
c
      
      
c      call zzfft(isign, npbig, scale, zwork1(1),zwork1(1),tableb(1),
c     & workb(1),0)

c       call zzfft(isign, npbig, scale, zwork2(1),zwork2(1),tableb(1),
c     & workb(1),0)
      
c      call dcfft(zwork1(1), npbig, isign, scale, work(1))
c      call dcfft(zwork2(1), npbig, isign, scale, work(1))
c
      if(iA.ne.0) then
         call decalk(zwork1(1),zwork(1),iA,npbig)
         call decalk(zwork2(1),zwork1(1),-iA,npbig)
      else
         call ZCOPY(npbig,zwork1(1),1,zwork(1),1)
         call ZCOPY(npbig,zwork2(1),1,zwork1(1),1)
      end if
c
      do i = 0, npbig/2-1
         rk = i*dk
         phase = -2.d0*cte*rk*phase2
         zexph1 =  cdexp(dcmplx(0.d0,phase))
         zexph2 =  dconjg(zexph1)
         ztmp1 = zexph1*zwork(i+1)
         ztmp2 = zexph2*zwork1(i+1)
         zwork1(i+1) = ztmp2
         zwork(i+1) = ztmp1
      end do
c
      do i = npbig/2, npbig-1
         rk = (i-npbig)*dk
         phase = -2.d0*cte*rk*phase2
         zexph1 =  cdexp(dcmplx(0.d0,phase))
         zexph2 =  dconjg(zexph1)
         ztmp1 = zexph1*zwork(i+1)
         ztmp2 = zexph2*zwork1(i+1)
         zwork1(i+1) = ztmp2
         zwork(i+1) = ztmp1
      end do
c
      call ZVES(npbig,zwork(1),1,zwork1(1),1,zwork2(1),1)
      call ZVEA(npbig,zwork(1),1,zwork1(1),1,zwork1(1),1)
c
      ctet = cte*(t-tf)
      do i = 0, npbig/2-1
         rk = i*dk
         phase = ctet*(rk**2)
         zexpph =  cdexp(dcmplx(0.d0,phase))
         ztmp1 = zexpph*zwork1(i+1)
         ztmp2 = zexpph*zwork2(i+1)
         zwork1(i+1) = ztmp1
         zwork2(i+1) = ztmp2
      end do
      do i = npbig/2, npbig-1
         rk = (i-npbig)*dk
         phase = ctet*(rk**2)
         zexpph =  cdexp(dcmplx(0.d0,phase))
         ztmp1 = zexpph*zwork1(i+1)
         ztmp2 = zexpph*zwork2(i+1)
         zwork1(i+1) = ztmp1
         zwork2(i+1) = ztmp2
      end do
c
      zsqr2 = 0.5d0 * zexph3
c
      call ZAXPY(npbig,zsqr2,zwork1(1),1,zpsik1(1),1)
      call ZAXPY(npbig,zsqr2,zwork2(1),1,zpsik2(1),1)
c
1004  format(4(e16.8e3,2x))
      return
      end

      subroutine decalk(ztab1,ztab2,id,n)
c
      integer id, n, idd
      double complex ztab1(0:n-1),ztab2(0:n-1)
c
      if(id.lt.0) then
            idd = -id
            call ZCOPY(n-idd,ztab1(0),1,ztab2(idd),1)
            call ZCOPY(idd,ztab1(n-idd),1,ztab2(0),1)
      else
            call ZCOPY(n-id,ztab1(id),1,ztab2(0),1)
            call ZCOPY(id,ztab1(0),1,ztab2(n-id),1)
      end if
c
c      write(*,*) id, n
      return
      end

      subroutine calczcut(zcutA, zcutI, rdeb, delr, npbig,
     &                    r0cut, scut)
c
      integer npbig
      double complex zcutA(npbig), zcutI(npbig)
      double precision rdeb, delr, r0cut, scut, pi
c
      double precision r, r1, tmp
      integer n
c
      r1 = r0cut + scut
      pi = dacos(-1.d0)
      r = rdeb - delr
      do n = 1, npbig
         r = r + delr
         if(r.le.r0cut) then
            tmp =0.d0
         else if((r.gt.r0cut).and.(r.le.r1)) then
            tmp = (dsin(0.5d0*(r-r0cut)*pi/scut))**2
         else
            tmp = 1.d0
         end if
c     **************************************************************
c         tmp = 1.d0 + dexp((r-r0cut)/scut)
c         tmp = 1.d0/tmp
c     **************************************************************
         zcutI(n) = dcmplx(1.d0 - tmp, 0.d0)
         zcutA(n) = dcmplx(tmp, 0.d0)
c     **************************************************************
c         zcutA(n) = dcmplx(1.d0 - tmp, 0.d0)
c         zcutI(n) = dcmplx(tmp, 0.d0)
c     **************************************************************
c      write(15,1003) r, dreal(zcutI(n)), dreal(zcutA(n))
      end do
1003  format(3(e16.8e3,2x))
      return
      end

