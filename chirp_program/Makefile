wpac_pulse_tot_impul.out:      wpackha_pulse_tot_fft_corr_menage.o eval.o spo_corr.o aire_pulse_tot_impul.o potspec.o \
                morse.o cmplxmath.o 
# 
	ifort    -o wpac_pulse_tot_impul.out wpackha_pulse_tot_fft_corr_menage.o eval.o spo_corr.o aire_pulse_tot_impul.o potspec.o \
                  morse.o \
               cmplxmath.o dqage_dep.o  
#
aire_pulse_tot_impul.o:		aire_pulse_tot_impul.f
		ifort   -c -O aire_pulse_tot_impul.f 
spo_corr.o:	spo_corr.f
		ifort   -c -O spo_corr.f 
wpackha_pulse_tot_fft_corr_menage.o:	wpackha_pulse_tot_fft_corr_menage.f
		ifort -c -O wpackha_pulse_tot_fft_corr_menage.f 
eval.o:	eval.f
		ifort  -c -O eval.f 
potspec.o:	potspec.f
		ifort  -c -O potspec.f 
morse.o:	morse.f
		ifort   -c -O morse.f 
cmplxmath.o:	cmplxmath.f
		ifort   -c -O cmplxmath.f 
