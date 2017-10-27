program scanE
    use omp_lib
    Use MKL_DFTI
    use morse1
    use gen
    include 'mpif.h'
    include 'mkl_dfti_examples.fi'
    integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)    
    integer Ethmax,Ethmin
    integer :: nE0,pulsetype,iE0,id,le0wattcm2,logfile,nt,i,v,npos,n,cas
    character(LEN=50) :: title
    real(8), allocatable :: ep(:,:)
    real(8), allocatable :: pot(:,:)
    real(8) :: phase,wir,requ,diss,massreduite,TVIB
    real(8) :: tc,dw,te,tf,t0,period,dt,delr
    real(8) :: E0min,E0max,dE0,E0,norme,alpha,rc0,p0
    real(8) :: xmin,xmax  
    real(8), allocatable :: t(:),x(:),work1(:),xmu12(:)
    real(8), allocatable :: pbfin(:)
    complex(8), allocatable :: chi1(:), chi2(:)


     call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
   print*, 'node', rank, ': Hello world'

    open(5,name="input",status="old")
    namelist /iofile/ title, pulsetype,E0min,E0max,dE0,wir,phase,le0wattcm2,nt,npos
    read(5,iofile)




    allocate(t(nt))
    allocate(work1(npos))
    allocate(chi1(npos))
    allocate(chi2(npos))
    allocate(pot(2,npos))
    allocate(xmu12(npos))
    
    
    xmax=30.d0
    xmin=2.d-3
    
    
rc0 = 1.3989d0 !position du paquet d'onde à t0
p0 = 0.0d0 !impulsion du paquet d'onde à t0
alpha = 13.019d0 !paramètre pour chi1 et chi2

    
    
    
    logfile=987654
open(987654,name="log.dat",status="replace")
v=19
    allocate(ep(v,npos))
    allocate(x(npos))
    requ=2.d0 !valeur de r à l'equilibre = 2*a0 (a0=1 en u.a) 
diss=2.7925d0/27.2d0 !potentiel de dissociation de H2+
massreduite=918.0762887608628d0 !=masse proton/2 en u.a

!write(*,*) "******************************************************************"
!write(*,*) "* Quantum Wave Packet Dynamics of H2+ by split operator technics *"
!write(*,*) "******************************************************************"
!write(*,*) "* Title : " , title 
!write(logfile,*) "Title : " , title
!if (pulsetype.eq.1) then
!    write(*,*) "* Pulsetype = Squarepulse "  
!    write(logfile,*) "Pulsetype = Squarepulse" 
!end if
period=2*MATH_PI/wir
phase=phase*MATH_PI
t0=0.d0
cas=1
TVIB=18.d0/tau2fs
if (cas.eq.1) then 
    tc=TVIB/2.d0-period/2.d0
else if (cas.eq.2) then
    tc=TVIB/2.d0-period/4.d0
else
    write(*,*) "cas must be 1 (trapping) or 2(dissociative)"
    stop 10
end if
!write(*,*) "tc = " , tc , " en fs : " , tc*tau2fs 
!write(*,*) "period = " , period , " en fs : " , period*tau2fs 
!write(*,*) "tvib = " , TVIB , " en fs : " , TVIB*tau2fs 
!pause 10
te=tc+6.d0*period
tf=te+5.d0*period
dt=(tf-t0)/(nt-1)
    write(logfile,*) "dt = ", dt 
!write(*,*) "dt = " , dt, tf, tc, period
do i=1,nt
    t(i)=t0+(i-1)*dt
end do 
delr=(xmax-xmin)/(npos-1)
	do i=1,npos
	 x(i)=xmin+(i-1)*delr !création de l'axe des positions
	enddo

 do n=1,v

	do i=1,npos
 	  ep(n,i)=morse(diss,0.72d0,massreduite,requ,x(i),n-1) !construction des n etats vibrationnels sur la grille
	  work1(i)=(dabs(ep(n,i)))**2
    enddo

	  call simpson(npos,delr,work1,norme)
	  ep(n,:)=ep(n,:)/sqrt(norme) !normalisation des etats vibrationnels
      ! ajouter dans le scan_E.f90
!	open(unit=(10+n),status='replace')
!	do i=1,npos
!		write(10+n,*)x(i),ep(n,i)!ecriture des etats propres d'indice n dans fort.10+n
!	enddo	
!	close(10+n)
 enddo
 call eval(chi1, chi2, delr, xmin, p0, rc0, alpha, npos)
 call pot_spec(pot(1,:),pot(2,:),xmu12, npos,delr,xmin) !construction des 2 potentiels de H2+ et du moment dipolaire

   nE0=((E0max-E0min)/dE0)+1
   Ethmin = (nE0 / size * (rank+1));
   Ethmax = (nE0 / size * (rank + 2) -1 );
   if (Ethmax.gt.nE0+1) then
       Ethmax=nE0+1
   end if
   write(*,*) "nE0 : ", nE0 , "size : " ,size 

    allocate(pbfin(nE0))
   ! write(*,*) "nE0 E0max E0min dE0" , nE0 , E0max , E0min , dE0
    iE0=0
    do while ((iE0+(rank)).le.nE0)  
        !id = omp_get_thread_num ( )
        E0=E0min+(iE0+rank)*dE0
        E0=wattcm22au(E0)
        write(*,*) "E0 : " , E0 , "rank : " , rank, ", iE0 :", iE0 , "iE0+rank : " , iE0+rank
        call travail(E0,title,pulsetype,wir,phase,le0wattcm2,tc,te,tf,iE0,logfile,t,nt,ep,npos,v,x,id,dt,nt,xmu12,pot,delr,chi1,chi2,massreduite,pbfin(iE0))
        write(999654,*)iE0,E0,pbfin(iE0)
        iE0=iE0+size
    end do 
    close(logfile)
   call MPI_FINALIZE(ierror)
end program scanE
