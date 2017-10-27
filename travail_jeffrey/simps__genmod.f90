        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 25 15:19:32 2012
        MODULE SIMPS__genmod
          INTERFACE 
            SUBROUTINE SIMPS(FUNC,VINT,DELTI,NPL)
              INTEGER(KIND=4) :: NPL
              COMPLEX(KIND=8) :: FUNC(NPL)
              REAL(KIND=8) :: VINT
              REAL(KIND=8) :: DELTI
            END SUBROUTINE SIMPS
          END INTERFACE 
        END MODULE SIMPS__genmod
