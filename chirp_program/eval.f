      subroutine eval(cw1, cw2, psik1, psik2, delr, rdeb,
     &                p0, rc0, alpha, npos, npbig)
c
      double complex cw1(npos), cw2(npos)
      double complex psik1(npbig), psik2(npbig)
      double precision delr, rdeb, p0, rc0, alpha
      integer npos, npbig
c
      double precision pi, r
      double complex cnul, cim, cpoi, cval, arg
      integer l
c
      cnul = dcmplx(0.d0,0.d0)
      cim = dcmplx(0.d0,1.d0)
      pi = 3.141592654d0
      cpoi = cdsqrt(cdsqrt(dcmplx(2.d0*alpha/pi,0.d0)))
      r = rdeb-delr
      do l = 1, npos
         r = r + delr
         arg = dcmplx(-alpha*(r-rc0)**2, p0*(r-rc0))
         cval = cpoi*cdexp(arg)
         cw1(l) = cval
         cw2(l) = cnul
      end do
c
      do l = 1, npbig
         psik1(l) = cnul
	 psik2(l) = cnul
      end do
c
      return
      end
