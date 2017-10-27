subroutine travail(E0in,title,pulsetype,wir,phase,le0wattcm2,tc,te,tf,iE0,logfile,t,ntps,ep,npos,v,x,id,dt,nt,xmu12,pot,delr,chi1in,chi2in,massreduite,pbfin)
use spo
use morse1
use gen
use pulse
Use MKL_DFTI
include 'mkl_dfti_examples.fi'

!implicit none
integer , INTENT(IN) :: npos,pulsetype,logfile,le0wattcm2,iE0,ntps,nt,v,id
character(LEN=50) , INTENT(IN) :: title
real(8) , INTENT(IN) ::  wir,phase,tc,te,tf,x(npos),ep(v,npos),xmu12(npos),delr,massreduite,dt
real(8), INTENT(IN) :: E0in
complex(8), INTENT(IN) :: chi1in(npos), chi2in(npos)
real(8), INTENT(OUT) :: pbfin
real(8), INTENT(IN) :: pot(2,npos)


!npos pour la grille de position et d'impulsion, ntps pour la grille de temps

integer :: n,i,j,l,m,maxpot1,minpot2,ideb,ninput,nE0
parameter (ideb=85)
real(8) :: xmin, xmax, requ, diss, lieprobv,E0wattcm2
real(8) :: work1(npos),table1(npos),normedeb, champ(nt),E0,rc0
real(8) :: tablea(npos),worka(npos),projreal, projimag, lieprob 
complex(8) :: chi1(npos),chi2(npos),cun,cim,cnul,zetdt(npos),ctemp(npos),chilie(npos),chi1init(npos)
real(8) :: alpha,p0,rdeb,proj(npos) ,proji(npos),xmue,auto_correl(nt) 
real(8) :: t(nt),delt,t0,pi,omega, dtper, norme,periode,delta,vp1(npos),vp2(npos),sigma,tmax,kmoyen(nt),rmoyen(nt),rmoyenlie(nt),rclapet1(nt),rclapet2(nt),f0
real(8) :: pulset(nt),champir1(nt),periodir,dw
real(8),dimension(npos-ideb) :: vp1reel,vp2reel
real(8) :: time1,time2,time3,btime,ftime
character(LEN=50) :: pbname
character(LEN=50) :: charnum(10)
delt=dt
E0=E0in
!,t2,t3,clock_rate, clock_max
!***********************************************************************
!         Valeurs des paramètres, allocation des variables
!***********************************************************************
! call system_clock ( t1, clock_rate, clock_max )
!call cpu_time ( btime )
!write(*,*) "* Begining of the program ", btime
!write(logfile,*) "Begining of the program ", btime
!call readinputsub(tc,dw,tf)

 cun = dcmplx(1.0d0,0.d0)
 cim = dcmplx(0.d0,1.0d0)
 cnul = dcmplx(0.d0,0.d0)
 pi = MATH_PI
t0=0d0



WRITE(*,*) 'BEGINING OF THE CODE : Thread is ', id



!***********************************************************************
! Mise en place de la grille des positions, potentiels de Morse et
! etats propres associés
!***********************************************************************
 
! call pot_spec(pot(1,:),pot(2,:),xmu12, npos,delr,xmin) !construction des 2 potentiels de H2+ et du moment dipolaire
 ! AJOUTER DANS LE SCAN_E.f90
 ! pour ne pas ouvrir plusieur potentiel.dat
 !open(unit=1,name="potentiel.dat",status='replace') 
!	do i=1,npos
!	  write(1,*)x(i),pot(1,i),pot(2,i) 
!	enddo
! close(1)

write(*,*) 'THIS IS A TEST 11111111111111 !!!!!!!!!!!!! ', iE0


! *******************************************************************  	
!    		Construction de la grille temporelle	
! *******************************************************************
     

!***********************************************************************
!            Calcul de la fonction d'onde initiale
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$OMP PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$OMP DO
!call cpu_time ( btime )
!write(*,*) "Begining of time loop : ", time1
!do ph=0,1
! call eval(chi1, chi2, delr, rdeb, p0, rc0, alpha, npos)
 !do l=1,npos
!	chi1(l)=dcmplx(ep(0,l),0d0)
!	chi2(l)=dcmplx(0d0,0d0)
! enddo
 do j=1,npos
	chi1init(j)=chi1in(j)
 enddo



write(*,*) "LE TEST EST  FONCTIONNEL thread is " , id
!***********************************************************************
! 	      	Normalisation de la fonction d'onde initiale
!***********************************************************************
	do l=1,npos
	worka(l)=(cdabs(chi1in(l)))**2
	enddo

   call simpson( npos,delr,worka, normedeb)
 open(unit=1234567+iE0,name='chi1init.dat',status='replace')
 open(unit=234579+iE0,name='chi2init.dat',status='replace')
   do l = 1, npos
         chi1(l) = chi1in(l)/dsqrt(normedeb)
         chi2(l) = chi2in(l)
	write(1234567+iE0,*)x(l),dreal(chi1(l)),dimag(chi1(l))
	write(234579+iE0,*)x(l),dreal(chi2(l)),dimag(chi2(l))
   enddo
 close(1234567+iE0)
 close(234579+iE0)




!********************************************************************
! 	   		Construction du champ
!********************************************************************
! sigma=periode/(2d0*2.3548d0)
! tmax=2d0*periode
!tc=200.d0*10.d0
!open(unit=2000+iE0,status='replace')
!call champir(champir1,wir,phase,t,ntps,1)
!call pulsechoice(pulsetype,champ,t,ntps,delt,logfile,E0,wir,phase,le0wattcm2,tc,te,tf)
!pseudocarre(champ,t,ntps,delt,ph)
write(*,*) "Temps final : ", tf
    if (le0wattcm2.eq.1) then
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "This is not supossed to be printed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        E0=wattcm22au(E0in)
    end if
    E0wattcm2=au2wattcm2(E0in)
    write(logfile,'(" E0 = ", E16.8 , " u.a. , ",  E16.8 , "  W/cm2" )') E0 , E0wattcm2 
    write(logfile,*) "phase = ", phase/MATH_PI , " Pi"
    write(logfile,'(" tc = ", E16.8 , " a.u. ," , E16.8, " fs" )') tc , tc*tau2fs 
    write(logfile,'(" te = " , E16.8 , " a.u. ," , E16.8, " fs" )') te , te*tau2fs 
    write(logfile,*) "tf = ", tf, "a.u." , tf*tau2fs ," fs"
    write(logfile,*) "nt = ", nt
    write(logfile,*) "dt = ", dt, "a.u." , dt*tau2fs ," fs"

!********************************************************************
!	Application du split operator-Ouverture de la boucle temps 
!********************************************************************
write(*,*) "Juste avant la loop sur le temps , thread : ", id
 do i=1,nt
write(*,*) "Debut de la loop sur le tempsa ", i , " , thread : ", id
!call cpu_time ( time1 )
!write(*,*) "Time 1 (beg of time loop) : ", time1

!    if (mod(i,20).eq.0) then
!        write(*,'(I5," , t = ",F16.8," u.a, PB : ",F6.5,I5,F8.4)')        i,t(i),lieprob,iE0,t(i)-t(i-1)
!        do j=1,npos
            !write((10000*+i,'(I4,7E16.8)')i,t(i),t(i)*2.418884326505E-2,x(j),dabs(chi1(j))**2,dabs(chi2(j))**2,vp1(j),vp2(j)
!            write(90000+i,'(I4,7E16.8)')i,t(i),t(i)*2.418884326505E-2,x(j),abs(chi1(j))**2,abs(chi2(j))**2,vp1(j),vp2(j)

!       end do
!    end if
    if (t(i).lt.tc) then
            pulset(i)=0.d0
            !E0*dexp(-alpha*(t(j)-tc)**2)
    else if (t(i).lt.te) then
            pulset(i)=E0
    else 
            pulset(i)=0.d0
            !E0*dexp(-alpha*(t(j)-te)**2)
    end if
        champ(i)=-dcos(wir*(t(i)-tc)+phase)*pulset(i)
        if (dabs(pulset(i)).lt.1.d-100) then
            pulset(i)=0.d0
        end if
        if (dabs(champ(i)).lt.1.d-100) then
            champ(i)=0.d0
        end if
write(*,*) "fin calcul du champs temps ", i , " , thread : ", id
	do j=1,npos
		worka(j)=x(j)*(cdabs(chi1(j)))**2 !position moyenne du paquet d'onde de l'état fondamental
		tablea(j)=(cdabs(chi1(j)))**2 !densité de probabilité du paquet d'onde de l'état fondamental
	enddo
	call simpson(npos,delr,worka,normedeb)
	call simpson(npos,delr,tablea,norme)
	rmoyen(i)=normedeb/norme
	!kmoyen(i)=normedeb/norme
	write(70000+iE0,*)t(i),rmoyen(i)!r moyen du paquet d'onde total
write(logfile,*) "Rmoyen du iE0 : ", iE0 , " : fort." , 70000+iE0

call cpu_time ( time2 )
!write(*,*) "Time2 : ", time2
! PARTIE TRES LENTE DU PROGRAMME : CE CALCUL PREND PRESQUE 1 SEC PAR APPLICATION
!	do j=1,npos
!		chilie(j)=dcmplx(0d0,0d0)
!		do n=0,v
!			ctemp(j)=ep(n,j)*chi1(j)
!			call simpson(npos,delr,dreal(ctemp),normedeb)
!			call simpson(npos,delr,dimag(ctemp),norme)
!			chilie(j)=chilie(j)+ep(n,j)*(normedeb+cim*norme)
!		enddo
!		worka(j)=x(j)*(cdabs(chilie(j)))**2
!		tablea(j)=(cdabs(chilie(j)))**2
!	enddo
!	call simpson(npos,delr,worka,normedeb)
!	call simpson(npos,delr,tablea,norme)
!	rmoyenlie(i)=normedeb/norme
!	write(30+ph,*)t(i),rmoyenlie(i)!rmoyen de la partie liee du paquet d'onde

	         !if (((t(i).gt.(0.d0)).and.(t(i).le.(0.d0+delt))).or.((t(i).gt.(1.d0*periode/4.d0)).and.(t(i).le.(1.d0*periode/4.d0+delt))).or.((t(i).gt.(2.d0*periode/4.d0)).and.(t(i).le.(2.d0*periode/4.d0+delt))).or.((t(i).gt.(3.d0*periode/4.d0)).and.(t(i).le.(3.d0*periode/4.d0+delt))).or.((t(i).gt.(4.d0*periode/4.d0)).and.(t(i).le.(4.d0*periode/4.d0+delt))).or.((t(i).gt.(12.d0*periode/4.d0)).and.(t(i).le.(12.d0*periode/4.d0+delt))).or.((t(i).gt.(13.d0*periode/4.d0)).and.(t(i).le.(13.d0*periode/4.d0+delt))).or.((t(i).gt.(14.d0*periode/4.d0)).and.(t(i).le.(14.d0*periode/4.d0+delt))).or.((t(i).gt.(15.d0*periode/4.d0)).and.(t(i).le.(15.d0*periode/4.d0+delt))).or.((t(i).gt.(16.d0*periode/4.d0)).and.(t(i).le.(16.d0*periode/4.d0+delt)))) then
	
!	open(unit=2000+i+(ntps+2000)*ph)
!	open(unit=2*ntps+5000+i+(ntps+2000)*ph)
	do j=1,npos
		xmue=xmu12(j)*champ(i)
		delta=(pot(2,j)-pot(1,j))**2+(2d0*xmue)**2
		delta=dsqrt(delta)
		vp1(j)=(pot(2,j)+pot(1,j)-delta)*0.5d0
		vp2(j)=(pot(2,j)+pot(1,j)+delta)*0.5d0
		if (j.gt.ideb) then
			vp1reel(j-ideb)=vp1(j)
			vp2reel(j-ideb)=vp2(j)
		endif
!		write(2000+i+(ntps+2000)*ph,*) x(j), dsqrt((dreal(chi1(j))**2+dimag(chi1(j))**2)),dsqrt((dreal(chi2(j))**2+ dimag(chi2(j))**2))
!construction des potentiels E- et E+
!		write(2*ntps+5000+i+(ntps+2000)*ph,*)x(j), vp1(j), vp2(j)
	enddo
write(*,*) "fin calcul des potentiel dependant du temps, temps : ", i , " , thread : ", id
!	close(2000+i+(ntps+2000)*ph)
!	close(2*ntps+5000+i+(ntps+2000)*ph)
			!endif

! if ((t(i).gt.(11071.75d0)).and.(t(i).le.(11072.5d0)))then
!	do j=1,npos
!		write(51,*)x(j),dreal(chi1(j)),dimag(chi1(j))
!	enddo

! endif
! if ((t(i).gt.(22143.75d0)).and.(t(i).le.(22144.5d0))) then
!	do j=1,npos
!		write(50,*)x(j),dreal(chi1(j)),dimag(chi1(j))
!	enddo
! endif

! do j=1,npos
!	ctemp=chi1init(j)*chi1(j)
!	worka(j)=(cdabs(ctemp(j)))**2 !a revoir (calcul d'autocorrelation)
! enddo
! call simpson (npos,delr,worka,auto_correl(i))
! write(52,*)t(i),auto_correl(i)

 maxpot1=maxloc(vp1reel,1)
 minpot2=minloc(vp2reel,1)
 rclapet1(i)=x(ideb+maxpot1)
 rclapet2(i)=x(ideb+minpot2)
 !write(200000+ph,*)t(i),rclapet1(i),rclapet2(i)
 write(*,*) "avant le split, t : ", i
	call splitop(chi1, chi2, zetdt,pot,xmu12, npos, champ(i), delr, massreduite, delt)
 write(*,*) "apres le split, t : ", i


!********************************************************************
!	    Calcul de probabilité de dissociation
!********************************************************************
	 lieprob = 0.d0
 !   pause 10 
     do n = 1, v
	       do j = 1, npos
                  proj(j) = ep(n,j)*dreal(chi1(j))
                  proji(j) = ep(n,j)*dimag(chi1(j))
                  
	       end do
       call simpson(npos,delr,proj,projreal)
	   call simpson(npos,delr,proji,projimag)
	   lieprobv=(projreal**2 + projimag**2)
	        lieprob = lieprob + lieprobv
!		if (i.eq.1) then		
!			write(49,*)n,lieprobv
!		endif
            end do
            write(10000+iE0,*) t(i), lieprob,champ(i)
write(logfile,*) "Pbount du iE0 : ", iE0 , " : fort." , 10000+iE0
!            write(*,*) t(i), lieprob,champ(i)

!write(*,*) "DIFF 2-1: ", time2-time1
!write(*,*) "DIFF 3-2: ", time3-time2
 enddo
call cpu_time ( time3 )
write(logfile,*) "Time 3 (end of program) : ", time3
!enddo
!call fourier(kmoyen,ntps,1,delt,delt)
    write(*,*) iE0
   write(*,*) "E0 : ", E0, "id : ", id



   !!!! WRITE TITLE IN GNUPLOT
   write(charnum(1),'(E16.8)')E0
   write(charnum(2),'(E16.8)')E0wattcm2

   pbname="set title 'Pbound for E="//trim(charnum(1))//" u.a. , "//trim(charnum(2))//" w/cm2'"
   write(123456,"(A50)") pbname

   !!! SET TERMINAL AND OUTPUT FILE
   write(123456,*) "set terminal png"
   write(charnum(3),"(I5)")10000+iE0
   pbname="set output 'pbound_"//trim(charnum(3))//".png'"
   write(123456,*) pbname
   
   !!! SET PLOTTING DATA
   pbname="plot 'fort."//trim(charnum(3))//"' u 1:2 w l"
   write(123456,*) pbname

   !!! PLOT THE FIELD
   pbname="set title 'Laser field for E="//trim(charnum(1))//" u.a. , "//trim(charnum(2))//" w/cm2'"
   write(123456,"(A50)") pbname
   write(123456,*) "set terminal png"
   pbname="set output 'champ_"//trim(charnum(3))//".png'"
   write(123456,*) pbname
   pbname="plot 'fort."//trim(charnum(3))//"' u 1:3 w l"
   write(123456,*) pbname
   

   !!! PLOT RMOYEN
   write(charnum(4),"(I5)")70000+iE0
   pbname="set title 'rmoyen field for E="//trim(charnum(1))//" u.a. , "//trim(charnum(2))//" w/cm2'"
   write(123456,"(A50)") pbname
   write(123456,*) "set terminal png"
   pbname="set output 'rmoyen_"//trim(charnum(4))//".png'"
   write(123456,*) pbname
   pbname="plot 'fort."//trim(charnum(4))//"' u 1:2 w l"
   write(123456,*) pbname
   pbfin=lieprob
end subroutine travail
!filename='obs_'//job
!open(60,file=filename,status='unknown',form='formatted')
!write (chain(6:8),'(i3.3)') int
!*******************************************************************  
!
!	Calcul du paquet d'onde initial
!
!*******************************************************************


     subroutine eval(cw1, cw2, delr, rdeb,p0, rc0, alpha, npos)

      integer npos
      double complex cw1(npos), cw2(npos)
      double precision delr, rdeb, p0, rc0, alpha
      double precision pi, r
      double complex cnul, cim, cpoi, cval, arg
      integer l

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
 
      return
      end

!***************************************************************
! 		Calcul norme d'une fonction complexe
!***************************************************************
	subroutine simps(func, vint, delti, npl)

      integer j, npl
      double complex func(npl)
      double precision  vint, delti
      
      vint=0d0
      do j = 1, npl-1
         vint=vint+delti*sqrt(cdabs(func(j))**2) 
      end do
      return
      end

!***************************************************************
! 		Integration Simpson
!***************************************************************

SUBROUTINE simpson (N,H,FI,S)
!
! Subroutine for integration over f(x) with the Simpson rule.  FI:
! integrand f(x); H: interval; S: integral.  Copyright (c) Tao Pang 1997.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I
  double precision, INTENT (IN) :: H
  double precision :: S0,S1,S2
  double precision, INTENT (OUT) :: S
  double precision, INTENT (IN), DIMENSION (N) :: FI
!
  S  = 0.0
  S0 = 0.0
  S1 = 0.0
  S2 = 0.0
  DO I = 2, N-1, 2
    S1 = S1+FI(I-1)
    S0 = S0+FI(I)
    S2 = S2+FI(I+1)
  END DO
  S = H*(S1+4.0*S0+S2)/3.0
!
! If N is even, add the last slice separately
!
  IF (MOD(N,2).EQ.0) S = S &
     +H*(5.0*FI(N)+8.0*FI(N-1)-FI(N-2))/12.0
END SUBROUTINE simpson


