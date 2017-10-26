      program wpack
      
c     Program developped by H.Abou-Rachid and Catherine Levebre May 2001
c     propagation de paquets d'ondes par split-operator
c     **************************************************************************
      character*24 fichier1, fichier2, fichier3, fichier4, fichier5
      logical writepot /.false./
      logical Etatpro
      integer idim, ndim
      parameter (idim=4096,ndim=4*idim)
      parameter(limit=1024*8)
      double complex chi1(idim), chi2(idim),psik1(ndim),psik2(ndim)
      double complex chi11(idim),chi22(idim)
      double complex zetdt(idim)
      double precision v1(idim),v2(idim),xmu12(idim)
      double precision morse, r0cut, scut, vibfunc(0:18,idim)
      double precision delr, rdeb, rmax, delt, tf, seuil,tw
      double precision t0, sig, omega, intens, freq,gama,phase
      integer npos, ntemps, npbig
c
      double complex zwork1(ndim), zwork2(ndim), 
     &zwork3(ndim), zwork4(ndim)
      double complex cun, cim, cnul, zcutA(ndim), zcutI(ndim)
      double precision r, xmu, pi, rc0, p0, alpha, t
      double precision int1tf, int2tf, int3tf
      double precision int1t0, int2t0, int3t0
      double precision projreal, projimag, dissprob, champ, env
      double precision w1(idim), w2(idim), proj(idim)
      double precision xnorm1, xnorm2, work1(ndim), work2(ndim)
      double precision worka(2*idim), workb(2*ndim)
      double precision tablea(2*idim+30),tableb(2*ndim+30)
      double precision xnormk1, xnormk2, work3(ndim), work4(ndim)
c      double precision wspos(4*idim+15),wsbig(4*ndim+15)
c      double complex wspos(2*idim+30),wsbig(2*ndim+30)
      double precision dk, cte, xk, spk
      double precision dtper, timeper, evbyau, normedeb, period
      double precision delta, xmue, thet, vp1, vp2
      double complex cwtemp1, cwtemp2
      integer k, m
      integer i, j, l, nu, nudep
      external time
      integer time, cput1, cput2, metint 
      double precision epsabs,epsrel
       integer key, limit, neval, ier, iord(limit),last
       double precision result,abserr, alist(limit), blist(limit)
       double precision rlist(limit),elist(limit)
       external efield2
       double precision field2
       double precision calcfreq,freqt,tc
c        call dqage(efield2,t0,tf,epsabs,epsrel,key,limit,result,abserr,
c     &  neval,ier,alist,blist,rlist,elist,iord,last)
      namelist /iofile/
     & fichier1, fichier2, fichier3, fichier4, fichier5,
     & delt, tf, t0, sig, omega, intens, freq,phase,
     & npos, npbig, rdeb, rmax, seuil, nudep,gama,tw,
     & r0cut, scut, xmu, Etatpro,metint,beta,tc
c
c
      common // t0,tw,omega,intens,phase,gama,sig,freq,beta

      cput1 = time()
c
      fichier1='fch1.dat'
      fichier2='fch2.dat'
      fichier3='fch3.dat'
      fichier4='fch4.dat'
      fichier5='fch5.dat'
 
      write(*,*) 'start'
      t0 = 0.d0
      sig = 1.d0
      omega = 1.d0
      intens = 5.33d-2
      freq = .1381962d0
      phase=1.d0
      gama=0.5d0
      npos = 1024
      npbig = 4096
      rdeb = 2.5d0
      rmax = 60.d0
      delt = 1.778528d-01
      tf = 400.d0
      tw=500.d0
        epsabs=0.00E+00
        epsrel=1.0E-03
        key=6
c     ******************************************************************
      nu = 0
      xmu = 918.074d0
c     ******************************************************************
      read(5,iofile)
c
      cun = dcmplx(1.0d0,0.d0)
      cim = dcmplx(0.d0,1.0d0)
      cnul = dcmplx(0.d0,0.d0)
      pi = 3.141592654d0
      evbyau = 27.212d0
      rc0 = 1.3989d0
      p0 = 0.0d0
      alpha = 9.2085226d0
      freq=freq/(27.212*8065)
      final_freq=final_freq/(27.212*8065)
c      alpha = 13.019d0
      delr = (rmax-rdeb)/(npos-1)
      dtper = 2.d0*pi/freq/(1024.d0)
      period = 2.d0*pi/freq
c
      ntemps = idnint(tf/delt+0.5d0)
c       write(*,*) 'OK'
      open(99,file="field.dat",status='unknown')
      open(98,file="field2.dat",status='unknown')
      open(97,file="freq.dat",status='unknown')
      
c      call zzfft(0,npos,1.0d0,chi11(1),chi11(1),tablea(1),worka(1),0)
c      call zzfft(0,npbig, 1.0d0,chi22(1),chi22(1),tableb(1),workb(1),0)
c       write(*,*) 'OK'
c
      call eval(chi1(1), chi2(1), psik1(1), psik2(1), delr, rdeb,
     &          p0, rc0, alpha, npos, npbig)
c
      call pot_spec(v1(1), v2(1), xmu12(1), npos, delr, rdeb)
c
      call calczcut(zcutA(1), zcutI(1), rdeb, delr, npbig,
     &              r0cut, scut)
c
      call zexptdt(zetdt(1), npos, delr, xmu, delt)
c
      do nu = 0, 18
         r = rdeb - delr
         do l = 1, npos
            r = r + delr
            vibfunc(nu,l) = morse(1.026202d-1,0.72d0,xmu,2.d0,r,nu)
         end do
      end do
c
      if (Etatpro) then
         do l = 1, npos
            chi1(l) = vibfunc(nudep,l)
         end do
      end if
c
      do l = 1, npos
         work1(l) = cdabs(chi1(l))**2
      end do
      call simps(work1(1), normedeb, delr, npos)
c
      do l = 1, npos
         chi1(l) = chi1(l)/dsqrt(normedeb)
         work1(l) = cdabs(chi1(l))**2
      end do
      call simps(work1(1), normedeb, delr, npos)
c
      t = 0.d0
      int1tf = 0.d0
      int2tf = 0.d0
      int3tf = 0.d0
      timeper = 0.d0
      do while(timeper.lt.tf) 
        timeper = timeper + dtper
	freqt=calcfreq(freq,beta,timeper,t0,tc)

c  INTEGRATION BY TRAPEZOIDE RULES
        if (metint.eq.0) then
c
c	 write(*,*)"freqt",timeper,int1tf
         call airesint (int1tf, int2tf, int3tf, t0, sig, omega,
     &                  phase, intens, freqt,gama, timeper, dtper,tw)
         else
c  INTEGRATION BY QuADRATURE RULES
       call dqage(efield2,timeper-dtper,timeper,epsabs,epsrel,key,
     & limit,result,
     &  abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if (ier.gt.0) then
	    write(*,*)"WARNING!!! ERROR IN INTEGRATION 
     &      (final time)",timeper,field2(timeper),ier
          endif 
	   
            temp=int1tf
            int1tf=int1tf-result/2
            int2tf=int2tf + (temp+int1tf)*0.5d0*dtper
            int3tf=int3tf + (temp**2 +int1tf**2)*0.5d0*dtper
          end if
      end do
      freqt=calcfreq(freq,beta,tf,t0,tc)
      if (metint.eq.0) then
      call airesint(int1tf, int2tf, int3tf, t0, sig, omega,
     &              phase, intens, freqt,gama, tf, tf-timeper,tw)
      else
        call dqage(efield2,t-timeper,tf,epsabs,epsrel,key,limit,
     &  result,
     &  abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if (ier.gt.0) then
	    write(*,*)"WARNING!!! ERROR IN INTEGRATION 
     &      (final time last...)",timeper,field2(timeper),ier
          end if
        temp=int1tf
        int1tf=int1tf-result/2
	int2tf=int2tf + (temp+int1tf)*0.5d0*dtper
	int3tf=int3tf + (temp**2 +int1tf**2)*0.5d0*dtper
       endif
      
      int1t0 = 0.d0
      int2t0 = 0.d0
      int3t0 = 0.d0
      t = 0.d0
      timeper = 0.d0
      m = 0
      k = 0

      do i = 1, ntemps
         t = t + delt
c         write(*,*)field3(t)
         if (metint.eq.0) then
         do while(timeper.lt.t)
            timeper = timeper + dtper
            freqt=calcfreq(freq,beta,timeper,t0,tc)
c	    write(*,*)timeper,freqt,beta,t0
            call airesint(int1t0, int2t0, int3t0, t0, sig, omega,
     &                 phase, intens, freqt,gama, timeper, dtper,tw)
        end do
            freqt=calcfreq(freq,beta,t,t0,tc)
         call airesint(int1t0, int2t0, int3t0, t0, sig, omega,
     &                 phase, intens, freqt,gama, t, t-timeper,tw)
        else if (metint.eq.1) then 
	do while(timeper.lt.t)
            timeper = timeper + dtper
       call dqage(efield2,timeper-dtper,timeper,epsabs,epsrel,key,limit,
     & result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if (ier.gt.0) then
	     write(*,*)"WARNING!!! ERROR IN INTEGRATION 
     & (current time)",ier
      end if
        temp=int1t0
        int1t0=int1t0-result/2
	int2t0=int2t0 + (temp+int1t0)*0.5d0*dtper
	int3t0=int3t0 + (temp**2 +int1t0**2)*0.5d0*dtper
        end do
	
         call dqage(efield2,timeper,t,epsabs,epsrel,key,limit,
     & result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     
      if (ier.gt.0) then
	     write(*,*)"WARNING!!! ERROR IN INTEGRATION 
     & (current time)",ier
      end if
     
        temp=int1t0
        int1t0=int1t0-result/2
	int2t0=int2t0 + (temp+int1t0)*0.5d0*dtper
	int3t0=int3t0 + (temp**2 +int1t0**2)*0.5d0*dtper

      endif

         timeper = t

c    INTEGRATION BY QUADPACK
c    funtion efield2 is field2 but called external
c    t0 initial time
c    t is current time (upper limit of integration)
c    the key=6 give good result
c    if ier < 0 ; error ...
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value
c                             of limit.
c                             however, if this yields no improvement it
c                             is rather advised to analyze the integrand
c                             in order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined(e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1) ,
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c       write(*,*)t,result/2.d0,int1t0
c       qint1tf=result/2.d0
c      qint2tf = int2 + (tmp + int1) * 0.5d0 * delt
c      qint3tf = int3 + (tmp**2 + int1**2) * 0.5d0 * delt

            freqt=calcfreq(freq,beta,t,t0,tc)
         write(97,999)t,freqt*(27.212*8065)

c            write(*,*)t,freqt
         champ = env (t0, sig, omega,
     &            intens, freqt,gama, phase,t - 0.5d0*delt,tw)
c
         write(99,999)t,field2(t-0.5d0*delt)
         write(98,999)t,champ
999   format(2(e16.8e3,e20.8E5,2x))
         call splitop(chi1(1), chi2(1), zetdt(1), v1(1), v2(1),
     &                xmu12(1), npos, champ, delr, xmu, delt,
     &                 tablea(1), worka(1))
c
c       write(*,*)"testtest2"
         if( (cdabs(chi1(npos)).gt.seuil).or.
     &                    (cdabs(chi2(npos)).gt.seuil) )then
            call asympt(t, tf, psik1(1), psik2(1), chi1(1), chi2(1),
     &                  zcutA(1), npbig, npos, xmu,
     &                  int1tf, int2tf, int3tf,
     &                  int1t0, int2t0, int3t0, delr,
     &                  zwork1(1), zwork2(1), zwork3(1),
     &                  worka(1), workb(1), tablea(1), tableb(1))
            call ZVEM(npos,chi1(1),1,zcutI(1),1,chi1(1),1)
            call ZVEM(npos,chi2(1),1,zcutI(1),1,chi2(1),1)
         endif
c       write(*,*)"testtest1"
c
         if(mod(i,10).eq.0)then
            dissprob = 0.d0
            do nu = 0, 18
               do j = 1, npos
                  proj(j) = vibfunc(nu,j)*dreal(chi1(j))
               end do
               call simps(proj(1),projreal,delr,npos)
               do j = 1, npos
                  proj(j) = vibfunc(nu,j)*dimag(chi1(j))
               end do
               call simps(proj(1),projimag,delr,npos)
c               if(nu.lt.11)then
c                  write(10+nu,1003) i*delt, projreal, projimag
c               end if
               dissprob = dissprob + projreal**2 + projimag**2
            end do
            write(29,1002) t, dissprob
         endif
c
         if (((t.gt.0.d0).and.(t.le.0.d0+delt)).or.
     &((t.gt.1.d0*period/4.d0).and.(t.le.1.d0*period/4.d0+delt)).or.
     &((t.gt.2.d0*period/4.d0).and.(t.le.2.d0*period/4.d0+delt)).or.
     &((t.gt.3.d0*period/4.d0).and.(t.le.3.d0*period/4.d0+delt)).or.
     &((t.gt.4.d0*period/4.d0).and.(t.le.4.d0*period/4.d0+delt)).or.
c
     &((t.gt.76.d0*period/4.d0).and.(t.le.76.d0*period/4.d0+delt)).or.
     &((t.gt.77.d0*period/4.d0).and.(t.le.77.d0*period/4.d0+delt)).or.
     &((t.gt.78.d0*period/4.d0).and.(t.le.78.d0*period/4.d0+delt)).or.
     &((t.gt.79.d0*period/4.d0).and.(t.le.79.d0*period/4.d0+delt)).or.
     &((t.gt.80.d0*period/4.d0).and.(t.le.80.d0*period/4.d0+delt)))then
            m = m + 1
c            write(30+m,1002) t, t/period
            r = rdeb - delr
            do j = 1, npos
               r = r + delr
               xmue = xmu12(j) * champ
	       delta = (v2(j) - v1(j))**2 + (2.d0*xmue)**2
	       delta = dsqrt(delta)
	       vp1 = (v2(j) + v1(j) - delta)*0.5d0
	       vp2 = (v2(j) + v1(j) + delta)*0.5d0
c               write(30+m,1007) r, dsqrt(dreal(chi1(j))**2+
c     &  dimag(chi1(j))**2),
c     & dsqrt(dreal(chi2(j))**2+ dimag(chi2(j))**2), vp1, vp2
            end do
         end if
c
         if (((t.gt.period).and.(t.le.period+delt)).or.
     &((t.gt.2.d0*period).and.(t.le.2.d0*period+delt)).or.
     &((t.gt.3.d0*period).and.(t.le.3.d0*period+delt)).or.
     &((t.gt.4.d0*period).and.(t.le.4.d0*period+delt)).or.
     &((t.gt.5.d0*period).and.(t.le.5.d0*period+delt)).or.
     &((t.gt.6.d0*period).and.(t.le.6.d0*period+delt)).or.
     &((t.gt.7.d0*period).and.(t.le.7.d0*period+delt)).or.
     &((t.gt.8.d0*period).and.(t.le.8.d0*period+delt)).or.
     &((t.gt.9.d0*period).and.(t.le.9.d0*period+delt)).or.
     &((t.gt.10.d0*period).and.(t.le.10.d0*period+delt)).or.
     &((t.gt.11.d0*period).and.(t.le.11.d0*period+delt)).or.
     &((t.gt.12.d0*period).and.(t.le.12.d0*period+delt)).or.
     &((t.gt.13.d0*period).and.(t.le.13.d0*period+delt)).or.
     &((t.gt.14.d0*period).and.(t.le.14.d0*period+delt)).or.
     &((t.gt.15.d0*period).and.(t.le.15.d0*period+delt)).or.
     &((t.gt.16.d0*period).and.(t.le.16.d0*period+delt)).or.
     &((t.gt.17.d0*period).and.(t.le.17.d0*period+delt)).or.
     &((t.gt.18.d0*period).and.(t.le.18.d0*period+delt)).or.
     &((t.gt.19.d0*period).and.(t.le.19.d0*period+delt)).or.
     &((t.gt.20.d0*period).and.(t.le.20.d0*period+delt)))then
            k = k + 1
c            write(50+k,1002) t, t/period
            dk = 2.d0*pi/(npbig*delr)
            cte = 0.5d0/xmu
c
            do j = npbig/2 +2, npbig
               xk = dk*(j-1 - npbig)
               work1(j) = cdabs(psik1(j))**2
               work2(j) = cdabs(psik2(j))**2
               work1(j) = work1(j)*(-xmu/(xk*evbyau))
               work2(j) = work2(j)*(-xmu/(xk*evbyau))
               spk = work1(j) + work2(j)
               xk = -xk*xk*0.5/xmu
c               write(50+k,1004) xk*evbyau, work1(j), work2(j), spk
            end do
c
            do j = 1, npbig/2 +1
               xk = dk*(j-1)
               work1(j) = cdabs(psik1(j))**2
               work2(j) = cdabs(psik2(j))**2
               if (xk.eq.0.d0) then
                  work1(j) = 0.d0
                  work2(j) = 0.d0
               else
                  work1(j) = work1(j)*(xmu/(xk*evbyau))
                  work2(j) = work2(j)*(xmu/(xk*evbyau))
               end if
               spk = work1(j) + work2(j)
               xk = xk*xk*0.5/xmu
c               write(50+k,1004) xk*evbyau, work1(j), work2(j), spk
            end do
         end if
c
c         write(9,*) 100.d0*float(i)/ntemps,'% done'
c         backspace(9)
      end do
c
      call asympt(tf, tf, psik1(1), psik2(1), chi1(1), chi2(1),
     &            zcutA(1), npbig, npos, xmu,
     &            int1tf, int2tf, int3tf,
     &            int1t0, int2t0, int3t0, delr,
     &            zwork1(1), zwork2(1), zwork3(1),
     &            worka(1), workb(1), tablea(1), tableb(1))
      call ZVEM(npos,chi1(1),1,zcutI(1),1,chi1(1),1)
      call ZVEM(npos,chi2(1),1,zcutI(1),1,chi2(1),1)
c
      open(1,file=fichier1,status='unknown')
      open(2,file=fichier2,status='unknown')
      r = rdeb - delr
      do l = 1, npos
	 r = r + delr
	 w1(l) = cdabs(chi1(l))**2
	 w2(l) = cdabs(chi2(l))**2
	 write(1,1004) dreal(chi1(l)),dimag(chi1(l)),
     $              dreal(chi2(l)),dimag(chi2(l))
	 write(2,1003) r, w1(l), w2(l)
      end do
      close(1,status='keep')
      close(2,status='keep')
c
      open(1,file=fichier3,status='unknown')
      open(2,file=fichier4,status='unknown')
      open(3,file=fichier5,status='unknown')
      dk = 2.d0*pi/(npbig*delr)
      cte = 0.5d0/xmu
c
      do i = npbig/2 +2, npbig
         xk = dk*(i-1 - npbig)
         work1(i) = cdabs(psik1(i))**2
         work2(i) = cdabs(psik2(i))**2
         work3(i) = work1(i)
         work4(i) = work2(i)
         work1(i) = work1(i)*(-xmu/(xk*evbyau))
         work2(i) = work2(i)*(-xmu/(xk*evbyau))
         spk = work1(i) + work2(i)
         xk = -xk*xk*0.5/xmu
         write(1,1002) xk*evbyau, work1(i)
         write(2,1002) xk*evbyau, work2(i)
         write(3,1002) xk*evbyau, spk
      end do
c
      do i = 1,npbig/2 +1
         xk = dk*(i-1)
         work1(i) = cdabs(psik1(i))**2
         work2(i) = cdabs(psik2(i))**2
         work3(i) = work1(i)
         work4(i) = work2(i)
         if (xk.eq.0.d0) then
            work1(i) = 0.d0
            work2(i) = 0.d0
         else
            work1(i) = work1(i)*(xmu/(xk*evbyau))
            work2(i) = work2(i)*(xmu/(xk*evbyau))
         end if
         spk = work1(i) + work2(i)
         xk = xk*xk*0.5/xmu
         write(1,1002) xk*evbyau, work1(i)
         write(2,1002) xk*evbyau, work2(i)
         write(3,1002) xk*evbyau, spk
      end do
c
      close(1)
      close(2)
      close(3)
      close(99)
      close(98)
      close(97)
c
      call simps(w1(1), xnorm1, delr, npos)
      call simps(w2(1), xnorm2, delr, npos)
      write(*,*)
      write(*,*) 'Niveau vibrationnel de depart'
      write(*,*) nu
      write(*,*)
      write(*,*) 'NORME de depart'
      write(*,*) normedeb
      write(*,*)
      write(*,*) 'NORMES population restante'
      write(*,*) '    SIGMA g           SIGMA u           TOTAL'
      write(*,1003) xnorm1, xnorm2, xnorm1 + xnorm2
      call simps(work3(1), xnormk1, dk, npbig)
      call simps(work4(1), xnormk2, dk, npbig)
      write(*,*)
      write(*,*) 'NORMES population dissociee'
      write(*,*) '    SIGMA g           SIGMA u           TOTAL'
      write(*,1003) xnormk1, xnormk2, xnormk1 + xnormk2
      write(*,*)
      write(*,*)
      write(*,*) '1 - NORMES population restante:', 1 - xnorm1 - xnorm2
      write(*,*)
      write(*,*) 'Erreur sur normes:',
     &           1 - xnormk1 - xnormk2 - xnorm1 - xnorm2
c

1002  format(2(e16.8e3,2x))
1003  format(3(e16.8e3,2x))
1004  format(4(e16.8e3,2x))
1005  format(5(e16.8e3,2x))
1006  format(6(e16.8e3,2x))
1007  format(7(e16.8e3,2x))
      write(*,*) 'end'
c
      cput2 = time() - cput1
c
      write(*,*) idint(cput2/3600.d0),' hr',
     &           idint((cput2 - mod(cput2,60) -
     &                  3600*idint(cput2/3600.d0))/60.d0),' min',
     &           mod(cput2,60),' sec'
c
      stop
      end

      subroutine simps(func, vint, delti, npl)
c
      integer j, npl
      double precision func(npl), vint, delti
      double precision s
c
      s = -.5d0 * (func(1) + func(npl))
      do j = 1, npl
         s = s + func(j)
      end do
      vint = delti * s
      return
      end

      function dargz(c)
      implicit double complex (a-c)
      implicit double precision (d-h,o-z)
      dargz=datan2(dimag(c),dreal(c))
      return
      end

c      double precision function field3(t)
c      double precision t
c      call field2(field,t,t0,tw,freq,omega,intens,phase,
c     &gama,siq)
c      field3=field
c      write(*,*)omega
c      return
c      end
      double precision function field2(t)
      double precision conv,pi,t,t0,tw,omega,intens,phase,gama,sig,freq
     &freqtime,beta,tc
     
      
      common // t0,tw,omega,intens,phase,gama,sig,freq,beta,tc
      conv = sig * 0.600561204d0
      pi=3.141592654d0
      freqtime=calcfreq(freq,beta,t,t0,tc)

            if ((t-t0).le.tw) then
c          write(*,*)"In the max ",tim
          field2=(-intens * dcos( freqtime*((t-t0)-pi/phase) ))
c	  field1(t0, sig, freq, intens,
c     &    tim,phase,tp)
      else  
c          write(*,*)"In the expo ",tim
          field2= dexp(-((t-tw)/conv)*
     &                ((t-tw)/conv))*
c     &                intens 
c     &                (field1(t0, sig, freq, intens,
c     &    tim,phase,tp))
     &    ( -intens * dcos( freqtime*((t-t0)-pi/phase) ))
c      else if ((tim-t0).gt.(12000.d0)) then 
c         write(*,*)tim,"In the 0 "
c         env=0.d0
      endif

      
c      if(((t-t0).ge.(t0-tw)).and.(t.le.tw))then
c      field2 = -intens * dcos(omega*(t)) *
c     &(dcos( freqtime* (t-t0))+gama*dcos(2.d0*freqtime*(t-t0)+
c     &pi/phase))
c       else
c          field2=0.d0
c      end if
c      write(*,*)t,t0,tw,omega,intens,phase,gama,sig,freq,field2
      return
      end

c   efield2 is EXTERNAL, so WARNING!!!
c    
c
c
c
c
c
      double precision function efield2(t)
      double precision conv,pi,t,t0,tw,omega,intens,phase,gama,sig,freq
     &freqtime,beta,tc
      common // t0,tw,omega,intens,phase,gama,sig,freq,beta,tc
      conv = sig * 0.600561204d0
      pi=3.141592654d0
      
      freqtime=calcfreq(freq,beta,t,t0,tc)
      
      if ((t-t0).le.tw) then
c          write(*,*)"In the max ",tim
          efield2=(-intens * dcos( freqtime*((t-t0)-pi/phase) ))
c	  field1(t0, sig, freq, intens,
c     &    tim,phase,tp)
      else  
c          write(*,*)"In the expo ",tim
          efield2= dexp(-((t-tw)/conv)*
     &                ((t-tw)/conv))*
c     &                intens 
c     &                (field1(t0, sig, freq, intens,
c     &    tim,phase,tp))
     &    ( -intens * dcos( freqtime*((t-t0)-pi/phase) ))
c      else if ((tim-t0).gt.(12000.d0)) then 
c         write(*,*)tim,"In the 0 "
c         env=0.d0
      endif

c      if(((t-t0).ge.(t0-tw)).and.(t.le.tw))then
c      efield2 = -intens * dcos(omega*(t)) *
c     &(dcos( freqtime* (t-t0))+gama*dcos(2.d0*freqtime*(t-t0)+
c     &pi/phase))
c       else
c          efield2=0.d0
c      end if
      return
      end

       
      double precision function calcfreq (freq0,beta,t,t0,tc)
      double precision freq0,beta,t,t0,tc

      if ((t-t0).le.tc) then
        calcfreq=freq0+beta*(t-t0)**2
      else
        calcfreq=freq0+beta*(tc-t0)**2
      endif
      return

      end
