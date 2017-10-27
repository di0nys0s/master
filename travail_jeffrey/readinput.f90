module readinput 
contains
subroutine readinputsub(tc,dw,tf)
    real(8) tc,dw,tf

    open(10,name="inputwp.fmt",status="old")
    read(10,*) tc
    write(*,*) "1 - Tc = ", tc
    read(10,*) dw
    write(*,*) "2 - dw = ", dw
    read(10,*) tf
    write(*,*) "3 - tf = ", tf
    close(10)
end subroutine readinputsub
end module readinput
