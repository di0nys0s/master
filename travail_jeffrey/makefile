#FC=ifort
FC=mpif90 
#FOPTS=-g 
#MKL="/opt/intel/mkl/10.0.3.020/"
#MKL="/software/intel/composer_xe_2013.0.079/mkl/"
MKL="/opt/intel/composer_xe_2013.0.079/mkl/"
#FOPTS=-L$(MKL)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -g -CB -traceback -implicitnone -check all -debug all -gen-interfaces -warn interfaces -check uninit
FOPTS=-L$(MKL)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread 

scanE: scan_E.f90 travail.o morse.o potspec_maj.o spo_maj.o cmplxmath.o pulse.o complex_1d_double_ex1.o
	$(FC) $(FOPTS) morse.o potspec_maj.o spo_maj.o cmplxmath.o pulse.o complex_1d_double_ex1.o  travail.o scan_E.f90 -o scanE
travail.o: travail.f90 morse.o potspec_maj.o spo_maj.o cmplxmath.o pulse.o complex_1d_double_ex1.o
	$(FC) $(FOPTS) -c complex_1d_double_ex1.o pulse.o cmplxmath.o morse.o potspec_maj.o spo_maj.o travail.f90 	

test: morse.o test.f90
	$(FC) morse.o test.f90 -o test

morse.o: morse.f90
	$(FC) $(FOPTS) -c morse.f90

potspec_maj.o: potspec_maj.f
	$(FC) $(FOPTS) -c potspec_maj.f

spo_maj.o : spo_maj.f90 cmplxmath.o  complex_1d_double_ex1.o
	$(FC) $(FOPTS) complex_1d_double_ex1.o cmplxmath.o -c spo_maj.f90

cmplxmath.o: cmplxmath.f
	$(FC) $(FOPTS) -c cmplxmath.f

pulse.o: pulse.f90 spo_maj.o complex_1d_double_ex1.o
	$(FC) $(FOPTS) complex_1d_double_ex1.o spo_maj.o -c pulse.f90

mkl_dfti.o : mkl_dfti.f90
	$(FC) $(FOPTS) -c mkl_dfti.f90

complex_1d_double_ex1.o: mkl/10.0.3.020/examples/dftf/source/complex_1d_double_ex1.f90
	cd mkl/10.0.3.020/examples/dftf/ && make libem64t  && cd - && cp mkl/10.0.3.020/examples/dftf/complex_1d_double_ex1.o . && cp mkl/10.0.3.020/examples/dftf/fft.mod .
run: travail
	./travail

fort.10: travail
	./travail

fort.11: travail
	./travail
clean:
	rm travail *.o

cscanE : scanE.c
	icc -openmp scanE.c -o cscanE
scan_Empi: scan_Empi.f90 travail.o morse.o potspec_maj.o spo_maj.o cmplxmath.o pulse.o complex_1d_double_ex1.o
	$(FC) $(FOPTS) morse.o potspec_maj.o spo_maj.o cmplxmath.o pulse.o complex_1d_double_ex1.o  travail.o scan_Empi.f90 -o scan_Empi
