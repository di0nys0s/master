program test
use gen
include 'mkl_dfti_examples.fi'
real(8) w,E,l,period,wcm1,periodfs
l=2.12d0*5.d0/2.d0
w=micron2au(l)
period=2.d0*MATH_PI/w
wcm1=au2cm1(w)
periodfs=period*tau2fs


write(*,*) "l=2.12 micron , w : " , w , "a.u. , w(cm1) : ", wcm1 , " , period : " , period , " a.u. , " , periodfs , " fs"
write(*,*) "w : " , wcm1 , "cm-1"
write(*,*) "period : " , period , "a.u."
write(*,*) "period/4 : " , period/4.d0 , "a.u."
write(*,*) "period : " , period , "fs"
write(*,*) "period/4.d0 : " , period/4.d0 , "fs"

E=8.77d13
E=wattcm22au(E)
write(*,*) "E : ", E
end program test
