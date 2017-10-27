        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 25 15:19:31 2012
        MODULE TRAVAIL__genmod
          INTERFACE 
            SUBROUTINE TRAVAIL(E0IN,TITLE,PULSETYPE,WIR,PHASE,LE0WATTCM2&
     &,TC,TE,TF,IE0,LOGFILE,T,NTPS,EP,NPOS,V,X,ID,DT,NT,XMU12,POT,DELR, &
     &CHI1IN,CHI2IN,MASSREDUITE,PBFIN)
              INTEGER(KIND=4), INTENT(IN) :: NT
              INTEGER(KIND=4), INTENT(IN) :: V
              INTEGER(KIND=4), INTENT(IN) :: NPOS
              REAL(KIND=8), INTENT(IN) :: E0IN
              CHARACTER(LEN=50), INTENT(IN) :: TITLE
              INTEGER(KIND=4), INTENT(IN) :: PULSETYPE
              REAL(KIND=8), INTENT(IN) :: WIR
              REAL(KIND=8), INTENT(IN) :: PHASE
              INTEGER(KIND=4), INTENT(IN) :: LE0WATTCM2
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: TE
              REAL(KIND=8), INTENT(IN) :: TF
              INTEGER(KIND=4), INTENT(IN) :: IE0
              INTEGER(KIND=4), INTENT(IN) :: LOGFILE
              REAL(KIND=8) :: T(NT)
              INTEGER(KIND=4), INTENT(IN) :: NTPS
              REAL(KIND=8), INTENT(IN) :: EP(V,NPOS)
              REAL(KIND=8), INTENT(IN) :: X(NPOS)
              INTEGER(KIND=4), INTENT(IN) :: ID
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: XMU12(NPOS)
              REAL(KIND=8), INTENT(IN) :: POT(2,NPOS)
              REAL(KIND=8), INTENT(IN) :: DELR
              COMPLEX(KIND=8), INTENT(IN) :: CHI1IN(NPOS)
              COMPLEX(KIND=8), INTENT(IN) :: CHI2IN(NPOS)
              REAL(KIND=8), INTENT(IN) :: MASSREDUITE
              REAL(KIND=8), INTENT(OUT) :: PBFIN
            END SUBROUTINE TRAVAIL
          END INTERFACE 
        END MODULE TRAVAIL__genmod
