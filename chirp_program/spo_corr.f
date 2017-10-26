      subroutine splitop(cw1, cw2, zetdt, v1, v2, xmu12,
     &                   npos, champ, delr, xmu, delt,
     &                   tablea, worka)
c
      integer npos

      double precision delr, champ, xmu, delt
      double complex cw1(npos),cw2(npos), zetdt(npos)
      double precision worka(2*npos), tablea(2*npos+30)
      double precision v1(npos),v2(npos),xmu12(npos)
c      double complex work(4*npos+15)
c
      double complex cun, cim, cnul, cwtemp, cphase1, cphase2
      double complex ctemp1, ctemp2
      double precision xmue, vp1, vp2, delta
      double precision xk1, xkn, arg, thet
      integer ll

      cun=dcmplx(1.0d0,0.d0)
      cim=dcmplx(0.d0,1.0d0)
      cnul=dcmplx(0.d0,0.d0)
      do ll = 1, npos
         xmue = xmu12(ll) * champ
	 delta = (v2(ll) - v1(ll))**2 + (2.d0*xmue)**2
	 delta = dsqrt(delta)
c	 if(dabs(xmue).gt.1.d-16)then
            thet = .5d0*datan((2.d0*xmue)/(v2(ll)-v1(ll)))
c         else
c            thet = 0.d0
c         end if
	 vp1 = (v2(ll) + v1(ll) - delta)*0.5d0
	 vp2 = (v2(ll) + v1(ll) + delta)*0.5d0
         cwtemp  = dcos(thet)*cw1(ll)-dsin(thet)*cw2(ll)
         cw2(ll) = dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
         cphase1 = cdexp(-cim*vp1*delt/2.d0)
         cphase2 = cdexp(-cim*vp2*delt/2.d0)
         cw1(ll) = cw1(ll)*cphase1
         cw2(ll) = cw2(ll)*cphase2
         cwtemp  = dcos(thet)*cw1(ll)+dsin(thet)*cw2(ll)
         cw2(ll) =-dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
      end do
c      c un commentaire
c     call zzfft(0, npos, 0.0, DUMMY, DUMMY, tablea(1), DUMMY, 0)
c      call zzfft(-1,npos,1.d0, cw1(1), cw1(1), tablea(1), worka(1),0)
c      call zzfft(-1,npos,1.d0,cw2(1), cw2(1), tablea(1), worka(1),0)
c
      call ZVEM(npos, zetdt(1), 1, cw1(1), 1, cw1(1), 1)
      call ZVEM(npos, zetdt(1), 1, cw2(1), 1, cw2(1), 1)
c
c      call zzfft(1,npos,1.d0/npos,cw1(1),cw1(1),tablea(1),worka(1),0)
c      call zzfft(1,npos,1.d0/npos,cw2(1),cw2(1),tablea(1),worka(1),0)
c
      do ll = 1, npos
         xmue = xmu12(ll) * champ
	 delta = (v2(ll) - v1(ll))**2 + (2.d0*xmue)**2
	 delta = dsqrt(delta)
c	 if(dabs(xmue).gt.1.d-16)then
	    thet = .5d0*datan((2.d0*xmue)/(v2(ll)-v1(ll)))
c         else
c            thet = 0.d0
c         end if
	 vp1 = (v2(ll) + v1(ll) - delta)*0.5d0
	 vp2 = (v2(ll) + v1(ll) + delta)*0.5d0
         cwtemp  = dcos(thet)*cw1(ll)-dsin(thet)*cw2(ll)
         cw2(ll) = dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
         cphase1 = cdexp(-cim*vp1*delt/2.d0)
         cphase2 = cdexp(-cim*vp2*delt/2.d0)
         cw1(ll) = cw1(ll)*cphase1
         cw2(ll) = cw2(ll)*cphase2
         cwtemp  = dcos(thet)*cw1(ll)+dsin(thet)*cw2(ll)
         cw2(ll) =-dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
      end do
      return
      end

      subroutine zexptdt(etdt, npos, delr, xmu, delt)
c
      integer npos
      double complex etdt(npos)
      double precision xmu, delr, delt
      double complex carg, cim
      double precision pi, xk1, xkn, arg
      integer ll
c
      pi = 3.141592654d0
      cim = dcmplx(0.d0, 1.d0)
      xk1 = 2.d0*pi/(delr*npos)
      do ll = 1, npos
         if(ll.le.npos/2)then
            xkn = (ll-1) * xk1
         else
            xkn = -(npos-ll+1) * xk1
         endif
	 arg = ((xkn*xkn)/(2.d0*xmu)) * delt
	 etdt(ll) = cdexp(-cim*arg)
      end do
      return
      end
