        !COMPILER-GENERATED INTERFACE MODULE: Wed Mar 20 19:11:02 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BOUNDRYINTERP__genmod
          INTERFACE 
            SUBROUTINE BOUNDRYINTERP(THOURS,ZQSTEMP1,NZQSTEMP,ZQSTIME,  &
     &ZQSTEMP)
              INTEGER(KIND=4) :: NZQSTEMP
              REAL(KIND=4) :: THOURS
              REAL(KIND=4) :: ZQSTEMP1
              REAL(KIND=4) :: ZQSTIME(NZQSTEMP)
              REAL(KIND=4) :: ZQSTEMP(NZQSTEMP)
            END SUBROUTINE BOUNDRYINTERP
          END INTERFACE 
        END MODULE BOUNDRYINTERP__genmod
