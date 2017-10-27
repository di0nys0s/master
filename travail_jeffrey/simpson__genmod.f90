        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 25 15:19:32 2012
        MODULE SIMPSON__genmod
          INTERFACE 
            SUBROUTINE SIMPSON(N,H,FI,S)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8), INTENT(IN) :: FI(N)
              REAL(KIND=8), INTENT(OUT) :: S
            END SUBROUTINE SIMPSON
          END INTERFACE 
        END MODULE SIMPSON__genmod
