program input
    use readinput
    use gen
    include 'mkl_dfti_examples.fi'
    integer :: choice,pulsetype,n,nt,j
    real(8) :: tc,dw,w,tf,te,E0,alpha,t0,dt,wnm
    real(8), allocatable :: pulse(:) 
    LOGICAL :: file_exists
    INQUIRE(FILE="inputwp.fmt", EXIST=file_exists)
    if (file_exists) then
        call readinputsub(tc,dw,tf)
    else
!        write(*,*) "Which pulse type?"
!        write(*,*) "1 - sinc"
!        write(*,*) "2 - square pulse"
!        read(*,*) pulsetype
        pulsetype=2
        if (pulsetype.eq.2) then
            w=4.55633526d-3*2
            !w=nm2au(w)
!            write (*,*) "enter L (micron)"
!            read (*,*) w
!            w=w*1.d-6
            !w=nm2au(w)
            write(*,*) "w in a.u. : ", w 
            wnm=au2nm(w)
            write(*,*) "w in n.m. : ", wnm 
            write(*,*) "T in a.u. : ", 1/w
            write(*,*) "T in fs : ", tau2fs(1/w)
            stop 
            write (*,*) "enter n (nt=n*1024)"
            read (*,*) n
            nt=n*1024
            write(*,*) "nt = ", nt
            allocate(pulse(nt))
            write (*,*) "enter t0"
            read (*,*) t0
            write (*,*) "enter tc"
            read (*,*) tc
            write(*,*) "Enter te"
            read (*,*) te
            write(*,*) "Enter tf"
            read (*,*) tf
            write(*,*) "Enter alpha"
            read (*,*) alpha
            write(*,*) "Enter E0"
            read (*,*) E0
            dt=(tf-t0)/(nt-1)
            do j=1,nt
                t=t0+(j-1)*dt
                if (t.lt.tc) then
                    pulse(j)=E0*dexp(-alpha*(t-tc)**2)
                else if (t.lt.te) then
                    pulse(j)=E0
                else 
                    pulse(j)=E0*dexp(-alpha*(t-te)**2)
                end if
                write(56789,*)t,pulse(j)
            end do
            go to 300
        end if
        write(*,*) "Enter Tc (u.a.) : "  
        read(*,*) tc
        write(*,*) "Enter dw (cm-1) : "  
        read(*,*) dw
        write(*,*) "Enter tf (u.a.) : "  
        read(*,*) tf
    end if
    100   continue
    write(*,*) "you want to change something? (0 for no)"
    read(*,*) choice
    if (choice.eq.0) then
        write(*,*) "ok, nothing to do"
        go to 200
    else if (choice.eq.1.) then
        write(*,*) "Enter Tc : "  
        read(*,*) tc
    else if (choice.eq.2) then
        write(*,*) "Enter dw (cm-1) : "  
        read(*,*) dw
    else if (choice.eq.3) then
        write(*,*) "Enter tf (u.a.) : "  
        read(*,*) tf    
    end if
    
    go to 100
    200 continue
    if (file_exists) then
        open(10,name="inputwp.fmt",status="replace")
    else 
        open(10,name="inputwp.fmt",status="new")
    end if

    write(10,*) tc 
    write(10,*) dw
    write(10,*) tf
    close(10)
    call readinputsub(tc,dw,tf)
    300 continue
end program input

