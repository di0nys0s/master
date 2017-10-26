       program testfreq

       double precision t,to,freq,freqt,beta
       double precision calcfreq
       integer ntemps

       write(*,*)calcfreq(0.0d0,0.001d0,1.d0,0.0d0)
       end
       double precision function calcfreq (freq0,beta,t,t0)
       double precision freqt,freq0,beta,t,t0

       calcfreq=freq0+beta*(t-t0)
       return


      end
