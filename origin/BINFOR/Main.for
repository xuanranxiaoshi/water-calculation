C#################################################################
C       MAIN PROGRAM
C#################################################################    
        PROGRAM LAKE
        INCLUDE 'PAR.INC'     
        INCLUDE 'BASE.INC' 

       OPEN(600,FILE='../TIME.DAT')
       READ(600,*) MDT
       READ(600,*) NDAYS
       CLOSE(600)
 
       OPEN(600,FILE='../GIRD.DAT')
       READ(600,*) NOD
       READ(600,*) CEL
       CLOSE(600)      
        
       OPEN(600,FILE='../DEPTH.DAT')
       READ(600,*) HM1
       READ(600,*) HM2
       CLOSE(600)        
        
       OPEN(600,FILE='../BOUNDARY.DAT')
       READ(600,*) NZ
       READ(600,*) NQ
       READ(600,*) NZQ
       READ(600,*) NHQ       
       READ(600,*) NWE
       READ(600,*) NDI                     
       CLOSE(600)        

        OPEN(800,FILE='../CALTIME.DAT')
        READ(800,*)         
        READ(800,*) DT
        CLOSE(800)
        STIME=0.0
        JTN=NDAYS
        
        OPEN(700,FILE='../MODEL.DAT')
        READ(700,*)   
        READ(700,*) NZQM        
        READ(700,*) NWQ
        READ(700,*) NSED
        CLOSE(700)        
        
         JT0=1
        ICON=0        
        NDT=INT(DT)
        K0=MDT/DT
        
        CALL PRE2                                   !�������ǰ�����ļ�       
        CALL Take_Boundary_for_TwoD                 !����߽�����

	  IF (ICON.EQ.0) THEN                            !�����ͬ�������
        OPEN(20,FILE='../OUTPUT/ZUV.OUT')
        OPEN(24,FILE='../OUTPUT/SHA.OUT')
        OPEN(28,FILE='../OUTPUT/H2U2V2.OUT')
        END IF
        OPEN(901,FILE='../OUTPUT/XY-TEC.DAT')       !���tecplot
        OPEN(902,FILE='../OUTPUT/SHA-TEC.DAT')       !���tecplot
        CALL OUTPUT2
        
        DO 900 JT=JT0,JTN     
*******************************************************�߽�Ĳ�ֵ    
C        WRITE(*,*) 'JT=====',JT
         DO 30 L=1,NZ
          IF(JT.NE.JTN) DZT(L)=(ZT(JT+1,L)-ZT(JT,L))/K0
30       CONTINUE
         DO 31 L=1,NQ
          IF(JT.NE.JTN) DQT(L)=(QT(JT+1,L)-QT(JT,L))/K0
31       CONTINUE
*******************************************************�߽�Ĳ�ֵ            
****************************************************************         
          IF (JT.EQ.JT0) THEN            !���������
          TAL=0
          DO 400 II=1,CEL
          TAL=TAL+AREA(II)
400       CONTINUE
          END IF
****************************************************************  
C         QLUA=QIT(JT)*0.5/TAL*DT
          QLUA=0
C        
          DO 800 KT=1,K0        
            IF(NZQM==1)CALL TWOD                              !����ˮ��
            IF(NWQ==1)CALL TWODWQ                             !����ˮ�� 
            IF(NSED==1)CALL SEDIMENTS                         !������ɳ           
            IF(MOD(INT(DT*KT),MDT).EQ.0) CALL OUTPUT2
800       ENDDO
900     CONTINUE
        CLOSE(4)
        END
*************************************************************************��ʼ��������        
       BLOCK DATA                                                  
        INCLUDE 'PAR.INC'    
        INCLUDE 'BASE.INC'
        
        DATA TG,THETA/9.81,0.6/           
        DATA M1SORT,M2SORT,M3SORT/0,2,15/                           !���ֳ�й�ʺ�������
               
       END
C#####################################################################        
C#################################################################
C       SUBROUTINE TAKING IN BOUNDARY CONDITIONS FOR              �߽���������
C       2-D CALCULATION
C#################################################################

C#####################################################################
        SUBROUTINE PRE2                                     !ǰ�����ļ�
!!      PREPROCESSING OF 2D WATER BODY COMPUTATION
        INCLUDE 'PAR.INC'        
        INCLUDE 'BASE.INC' 
        
        DIMENSION NW(4)
        DIMENSION PBC0(MSORT),HDC(CEL)
        DIMENSION XP0(NOD),YP0(NOD)
        
        DATA  AMU/1.02E-06/
        DATA DB/0.018E-03,0.06E-03/
        
        OPEN(12,FILE='../SOURCES/PNAC.DAT')
        OPEN(13,FILE='../SOURCES/PNAP.DAT')
        OPEN(14,FILE='../SOURCES/PKLAS.DAT')
        OPEN(15,FILE='../SOURCES/PZBC.DAT')
        OPEN(16,FILE='../SOURCES/MBZ.DAT')
        OPEN(17,FILE='../SOURCES/MBQ.DAT')
        
        READ(12,*)
        DO I=1,CEL
          READ(12,*) NO,(NAC(J,I),J=1,4)        
        END DO
        READ(13,*)
        DO I=1,CEL
          READ(13,*) NO,(NAP(J,I),J=1,4)
        END DO
        READ(14,*)
        DO I=1,CEL
          READ(14,*) NO,(KLAS(J,I),J=1,4)
        END DO
        READ(15,*)
        DO I=1,CEL
        READ(15,*) ZBC(I)
        ENDDO
        
***********************************����߽��ļ���������        
        READ(16,*)NNZ0                  !ˮλ�߽�����
        DO I=1,NZ
        READ(16,*) NO,MBZ(I),NNZ(I)
        ENDDO
        
        DO I=1,NZ         
        DO J=1,NZ-I
            IF (MBZ(J+1)<MBZ(J))THEN
                MBZMAX=MBZ(J)
                MBZ(J)=MBZ(J+1)
                MBZ(J+1)=MBZMAX
                NNZMAX=NNZ(J)
                NNZ(J)=NNZ(J+1)
                NNZ(J+1)=NNZMAX    
            ENDIF
        ENDDO
        ENDDO   
        
        READ(17,*)NNQ0                !�����߽�����
        DO I=1,NQ
        READ(17,*) NO,MBQ(I),NNQ(I)
        ENDDO

        DO I=1,NQ         
        DO J=1,NQ-I
            IF (MBQ(J+1)<MBQ(J))THEN
                MBQMAX=MBQ(J)
                MBQ(J)=MBQ(J+1)
                MBQ(J+1)=MBQMAX
                NNQMAX=NNQ(J)
                NNQ(J)=NNQ(J+1)
                NNQ(J+1)=NNQMAX    
            ENDIF
        ENDDO
        ENDDO        
        
***********************************����߽��ļ���������         
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)
        CLOSE(15)
        CLOSE(16)
        CLOSE(17)
***********************************��������꣬��ȡ��Сֵ
        OPEN(1,FILE='../SOURCES/PXY.DAT')
        READ(1,*)           
        DO I=1,NOD
         READ(1,*) NO,XP(I),YP(I)  
         XP(I)= XP(I)   
         YP(I)= YP(I)          
        END DO
        CLOSE(1)

        DO I=1,NOD
         XP0(I)= XP(I)
         YP0(I)= YP(I)       
        END DO       
        XIMIN=MINVAL(XP0)
        YIMIN=MINVAL(YP0)        
        DO I=1,NOD
         XP(I)= XP(I)-XIMIN
         YP(I)= YP(I)-YIMIN         
        END DO        
***********************************��������꣬��ȡ��Сֵ����ת��               
******************************************************   ��ʼˮλ    
        OPEN(801,FILE='../INITIALLEVEL.DAT')        
        READ(801,*)
        DO I=1,CEL
        READ(801,*) Z1(I)
        ENDDO
        CLOSE(801)    
******************************************************   ��ʼˮλ�޸� 
******************************************************   ��ʼ����u�޸�    
        OPEN(802,FILE='../INITIALU1.DAT')        
        READ(802,*)
        DO I=1,CEL
        READ(802,*) U1(I)
        ENDDO
        CLOSE(802)    
******************************************************   ��ʼ����u�޸� 
******************************************************   ��ʼ����v�޸�    
        OPEN(803,FILE='../INITIALV1.DAT')        
        READ(803,*)
        DO I=1,CEL
        READ(803,*) V1(I)
        ENDDO
        CLOSE(803)    
******************************************************   ��ʼ����v�޸�         
******************************************************   �����޸�     
        OPEN(804,FILE='../CV.DAT')
        READ(804,*)        
        DO I=1,CEL
        READ(804,*) FNC0(I)
        ENDDO
        CLOSE(804)    
******************************************************   �����޸�         
******************************************************   ��Ԫ������߳��Լ��Ƕȼ���        
        DO 40 I=1,CEL
        !IF (NAP(1,I).EQ.0) GOTO 40
        ZB1(I)=ZBC(I)+HM1
        XC(I)=0
        YC(I)=0
        NV(I)=4
        NA=NAP(4,I)
        !IF (NA.EQ.0.OR.NA.EQ.NAP(1,I)) NV(I)=3
        DO 30 J=1,NV(I)
        NW(J)=NAP(J,I)
        XC(I)=XC(I)+XP(NW(J))
        YC(I)=YC(I)+YP(NW(J))
30      CONTINUE
         XC(I)=XC(I)/NV(I)
         YC(I)=YC(I)/NV(I)
         XP1=XP(NW(1))
         XP2=XP(NW(2))
         XP3=XP(NW(3))
         YP1=YP(NW(1))
         YP2=YP(NW(2))
         YP3=YP(NW(3))
        AREA(I)=(XP1*YP2-YP1*XP2-XP1*YP3+YP1*XP3+XP2*YP3-YP2*XP3)/2.
        IF (NV(I).EQ.4) THEN
          XP4=XP(NW(4))
          YP4=YP(NW(4))
          AREA(I)=(XP1*YP3-YP1*XP3-XP1*YP4+YP1*XP4+XP3*YP4-YP3*XP4)/2.
     #          +AREA(I)
        END IF
        DO 35 J=1,NV(I)
            N1=NW(J)
            N2=NW(MOD(J,NV(I))+1)
            DX=XP(N1)-XP(N2)
            DY=YP(N2)-YP(N1)
            SIDE(J,I)=SQRT(DX*DX+DY*DY)
             IF(SIDE(J,I).GT.0.) THEN
                SINF(J,I)=DX/SIDE(J,I)
                COSF(J,I)=DY/SIDE(J,I)
             END IF
35      CONTINUE
40    CONTINUE
******************************************************   ��Ԫ������߳��Լ��Ƕȼ���    
C
******************************************************   ��ʼˮ�����  
        DO I=1,CEL
        ZB1(I)=ZBC(I)+HM1
        ENDDO       
        IF (JT0.GT.1) THEN
        OPEN (3, FILE='HUVZ.DAT', FORM='UNFORMATTED', STATUS='OLD')
        READ (3) H1,U1,V1,Z1,KLAS
        CLOSE (3)
        ELSE
        DO 190 I=1,CEL
        !IF (NAP(1,I).EQ.0) GOTO 190
        IF (Z1(I).LE.ZBC(I)) THEN
        H1(I)=HM1
        Z1(I)=ZB1(I)
        ELSE
        H1(I)=Z1(I)-ZBC(I)
        END IF
190     CONTINUE
        END IF
******************************************************   ��ʼˮ�����
        
        DO 210 I=1,CEL
        !IF (NAP(1,I).EQ.0) GOTO 210
        DTA(I)=DT/AREA(I)
        FNC(I)=TG*FNC0(I)*FNC0(I)
*****************************************������ˮ�����ģ��
C        FNC(I)=TG*(FNC0(I)+0.01/H1(I))*(FNC0(I)+0.01/H1(I))
*****************************************������ˮ�����ģ��
        Z2(I)=Z1(I)
        H2(I)=H1(I)
        U2(I)=U1(I)
        V2(I)=V1(I)
        DO  J=1,NV(I)
        SLCOS(J,I)=SIDE(J,I)*COSF(J,I)
        SLSIN(J,I)=SIDE(J,I)*SINF(J,I)
        ENDDO
210     CONTINUE

        END
C#####################################################################
C#####################################################################
C#################################################################
C       SUBROUTINE TAKING IN BOUNDARY CONDITIONS FOR              �߽���������
C       2-D CALCULATION
C#################################################################
       SUBROUTINE Take_Boundary_for_TwoD
        INCLUDE 'PAR.INC'        
        INCLUDE 'BASE.INC'     
       
       CHARACTER pointname*80
       REAL STIME1        
C#################################################################ˮλ�߽�           
        DO K=1,NNZ0
!            WRITE(*,*) NNZ0
          write(unit=pointname,fmt="(i4.4)") K            
	    OPEN(2,FILE='../BOUNDE/NZ/'//'NZ'//trim(pointname)
     1//'.dat')            
            READ(2,*)NZTEMP            
            DO I=1,NZTEMP
            READ(2,*) QZSTIME1(I),QZSTEMP1(I)
            ENDDO
            
            DO I=1,NDAYS
            STIME1=STIME+(I-1)/(24.*3600./MDT)
              CALL BOUNDRYinterp(STIME1,ZTTEMP,NZTEMP,QZSTIME1,QZSTEMP1)
              DO J=1,NZ
                  IF(NNZ(J)==K)THEN
              ZT(I,J)=ZTTEMP
c              WRITE(*,*) I,NNZ(J),ZT(I,NNZ(J))
                  ENDIF
              ENDDO
            END DO
          CLOSE(2)          
       ENDDO
C#################################################################ˮλ�߽� 
C#################################################################�����߽� 
        DO K=1,NNQ0
          write(unit=pointname,fmt="(i4.4)") K            
	    OPEN(2,FILE='../BOUNDE/NQ/'//'NQ'//trim(pointname)
     1//'.dat')            
            READ(2,*)NQTEMP            
            DO I=1,NQTEMP
            READ(2,*) QZSTIME2(I),QZSTEMP2(I)
            ENDDO
            
            DO I=1,NDAYS
            STIME1=STIME+(I-1)/(24.*3600./MDT)
              CALL BOUNDRYinterp(STIME1,QTTEMP,NQTEMP,QZSTIME2,QZSTEMP2)
              DO J=1,NQ
                  IF(NNQ(J)==K)THEN
              QT(I,J)=QTTEMP
                  ENDIF
              ENDDO
            END DO
          CLOSE(2)          
        ENDDO            
C#################################################################�����߽�             
        END
C#################################################################�߽��ֵ
      SUBROUTINE BOUNDRYinterp(THOURS,ZQSTEMP1,NZQSTEMP,ZQSTIME,ZQSTEMP)
      REAL ZQSTIME(NZQSTEMP),ZQSTEMP(NZQSTEMP)
      REAL THOURS,ZQSTEMP1
      
	DO I=1,NZQSTEMP-1
         if(THOURS>=ZQSTIME(I) .and. THOURS<=ZQSTIME(I+1))then
         ZQSTEMP1=ZQSTEMP(I)+(ZQSTEMP(I+1)-ZQSTEMP(I))/
	1   (ZQSTIME(I+1)-ZQSTIME(I))*(THOURS-ZQSTIME(I))
         endif
      ENDDO       
      END
C#################################################################�߽��ֵ
C#####################################################################
        SUBROUTINE TWOD                                        !���������
C        MAIN CONTROL OF 2D WATER BODY COMPUTATION
        INCLUDE 'PAR.INC'       
        INCLUDE 'BASE.INC'
C
        VMIN=0.001

        CALL PRECOR                 !ͨ��������һʱ����FLUX(K,J,I)������������һʱ����WH\WU\WV

C#################################################################ͨ��H1��WH,���H2
        DO 10 I=1,CEL
           HI=H1(I)
           BI=ZBC(I)
           UI=U1(I)
           VI=V1(I)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   ����䲽��
          IF(NV(I).EQ.4) THEN
           SIDEX=AMIN1(0.5*(SIDE(1,I)+SIDE(3,I)),
     #              0.5*(SIDE(2,I)+SIDE(4,I)))
          ELSE
           SIDES=0.5*(SIDE(1,I)+SIDE(2,I)+SIDE(3,I))
           SIDEX=SQRT((SIDES-SIDE(1,I))*(SIDES-SIDE(2,I))*
     #       (SIDES-SIDE(3,I))/SIDES)
          END IF
          HSIDE=AMAX1(H1(I),hm1)
          DT2(I)=SIDEX/(U1(I)+SQRT(9.81*HSIDE))
          DT2(I)=AMIN1(DT,DT2(I))
          DT2(I)=AMAX1(DT2(I),DT/10.0)
          !IF(I==13527)THEN
          !ENDIF    
          DTA(I)=1.0*DT2(I)/(1.0*AREA(I))          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   ����䲽��
          WDTA=1.00*DTA(I) 
           H2(I)=AMAX1(HI-WDTA*WH(I)+QLUA,HM1)    !
           Z2(I)=H2(I)+BI
           IF (H2(I).LE.HM1) THEN
             U2(I)=0.
             V2(I)=0.
           !  GOTO 10
           !END IF
             ELSE
           IF (H2(I).LE.HM2) THEN
            U2(I)=SIGN(AMIN1(VMIN,ABS(UI)),UI)
            V2(I)=SIGN(AMIN1(VMIN,ABS(VI)),VI)
             !U2(I)=0.02
             !V2(I)=0.02
           ! GOTO 10
           !END IF
             ELSE
           QX1=HI*UI
           QY1=HI*VI
           DTAU=WDTA*WU(I)
           DTAV=WDTA*WV(I)
           
            FNCC=FNC(I)
            WSF=FNCC*SQRT(UI*UI+VI*VI)/HI**0.33333
             U2(I)=(QX1-DTAU-DT*WSF*UI)/H2(I)
             V2(I)=(QY1-DTAV-DT*WSF*VI)/H2(I)
        !IF (H2(I).GT.HM2) THEN
          U2(I)=SIGN(AMIN1(ABS(U2(I)),5.0),U2(I))
          V2(I)=SIGN(AMIN1(ABS(V2(I)),5.0),V2(I))
        !ELSE
        !  U2(I)=SIGN(AMIN1(ABS(U2(I)),H2(I)),U2(I))
        !  V2(I)=SIGN(AMIN1(ABS(U2(I)),H2(I)),V2(I))
        !END IF
             ENDIF
             ENDIF
             
!!!
10      CONTINUE
C
C##################################################################ǰ�������ɺ�ֵ��ǰһʱ��        
        DO 60 I=1,CEL              
        H1(I)=H2(I)
        Z1(I)=Z2(I)
        U1(I)=U2(I)
        V1(I)=V2(I)
        W1(I)=SQRT(U1(I)*U1(I)+V1(I)*V1(I))
        W2(I)=SQRT(U2(I)*U2(I)+V2(I)*V2(I))
60      CONTINUE
C
        END

C#####################################################################
C#####################################################################
C#####################################################################  ����4���ߵ�ͨ��
        SUBROUTINE PRECOR
C        FLUX BALANCE FOR 2D CELLS
        INCLUDE 'PAR.INC'       
        INCLUDE 'BASE.INC'       
C
       
        DO 110 I=1,CEL
        !!IF (NAP(1,I).EQ.0) GOTO 110
        !WRITE(*,*) JT,KT,I,NV(I) 
        WH(I)=0
        WU(I)=0
        WV(I)=0
        DO L=1,4
        FLR(L)=0
        ENDDO
        HI=AMAX1(H1(I),HM1)
        
        IF (HI.GT.HM2) THEN
        UI=U1(I)
        VI=V1(I)
        ELSE
        UI=SIGN(VMIN,UI)
        VI=SIGN(VMIN,VI)
        END IF

        BI=ZBC(I)
        ZI=AMAX1(Z1(I),ZB1(I))
C
        !IF (HI.GT.HM1) GOTO 30
        IF (HI.LE.HM1) THEN
        DO J=1,4
        DO K=1,4
        FLUX(K,J,I)=0
        ENDDO
        ENDDO       
        DO 20 K=1,NV(I)
        NC=NAC(K,I)
        IF (NC.EQ.0.AND.KLAS(K,I).NE.0) CALL PRECOR1(I)
        IF (NC.NE.0.AND.H1(NC).GT.HM1) CALL PRECOR1(I)
20      CONTINUE
        ELSE
        CALL PRECOR1(I)            
        ENDIF  
        
110	CONTINUE
C
        DO 200 I=1,CEL
        !IF (NAP(1,I).EQ.0) GOTO 200
c        WRITE(*,*) I,NV(I)
        DO 190 J=1,NV(I)
        SL=SIDE(J,I)
        SLCA=SLCOS(J,I)
        SLSA=SLSIN(J,I)
        FLR(2)=FLUX(2,J,I)+FLUX(4,J,I)
        FLR(3)=FLUX(3,J,I)
        WH(I)=WH(I)+SL*FLUX(1,J,I)
        WU(I)=WU(I)+SLCA*FLR(2)-SLSA*FLR(3)
        WV(I)=WV(I)+SLSA*FLR(2)+SLCA*FLR(3)
190     CONTINUE
200     CONTINUE
        END

        SUBROUTINE PRECOR1(I)
        INCLUDE 'PAR.INC'        
        INCLUDE 'BASE.INC'
        DIMENSION KKNV(4)
        
        DO 100 J=1,NV(I)
          DO 35 L=1,4
            FLR(L)=0            
35        CONTINUE
        KKNV(J)=0
        
        SL=SIDE(J,I)
	  KP=KLAS(J,I)
C
        COSA=COSF(J,I)
        SINA=SINF(J,I)
        QL(1)=HI
        QL(2)=UI*COSA+VI*SINA
        QL(3)=VI*COSA-UI*SINA
        CL=SQRT(9.81*HI)
        FIL=QL(2)+2*CL
        NC=NAC(J,I)
        IF (NC.EQ.0) THEN
          HC=0
          BC=0
          ZC=0
          UC=0
          VC=0
        ELSE
          HC=AMAX1(H1(NC),HM1)
          BC=ZBC(NC)
          ZC=AMAX1(ZBC(NC),Z1(NC))
          UC=U1(NC)
          VC=V1(NC)
        END IF
C
        IF(KKNV(J).EQ.0)THEN         
36      IF (KP.GE.1.AND.KP.LE.8.OR.KP.GE.10) THEN
        CALL BOUNDA(J,I)
         DO K=1,4
         FLUX(K,J,I)=FLR(K)
         ENDDO
        KKNV(J)=1
        END IF
        ENDIF
C
!40    IF (HI.LE.HM1.AND.HC.LE.HM1) GOTO 100
        IF(KKNV(J).EQ.0)THEN  
        IF (HI.LE.HM1.AND.HC.LE.HM1)THEN
         DO K=1,4
         FLUX(K,J,I)=0
         ENDDO
        KKNV(J)=1
        ENDIF 
      ENDIF

        IF(KKNV(J).EQ.0)THEN      
        IF (ZI.LE.BC) THEN
        FLR(1)=-C1*HC**1.5
        FLR(2)=HI*QL(2)*ABS(QL(2))
        FLR(3)=0
        FLR(4)=4.905*HI*HI
         DO K=1,4
         FLUX(K,J,I)=FLR(K)
         ENDDO
        KKNV(J)=1
        END IF
      ENDIF
C
        IF(KKNV(J).EQ.0)THEN      
        IF (ZC.LE.BI) THEN
        FLR(1)=C1*HI**1.5
        FLR(2)=FLR(1)*QL(2)
        FLR(3)=FLR(1)*QL(3)
        FLR(4)=0
         DO K=1,4
         FLUX(K,J,I)=FLR(K)
         ENDDO
        KKNV(J)=1
        END IF
      ENDIF
C
        IF(KKNV(J).EQ.0)THEN
        IF (HI.LE.HM2) THEN
        IF (ZC.GT.ZI)  THEN
        DH=AMAX1(ZC-BI,HM1)
        UN=-C1*SQRT(DH)
        FLR(1)=DH*UN
        FLR(2)=FLR(1)*UN
        FLR(3)=FLR(1)*(VC*COSA-UC*SINA)
        FLR(4)=4.905*HI*HI
        ELSE
        FLR(1)=HI*C1*SQRT(HI)
        FLR(2)=0
        FLR(3)=0
        FLR(4)=4.905*HI*HI
        END IF
         DO K=1,4
         FLUX(K,J,I)=FLR(K)
         ENDDO
        KKNV(J)=1
        END IF
      ENDIF
C
        IF(KKNV(J).EQ.0)THEN      
        IF (HC.LE.HM2) THEN
        IF (ZI.GT.ZC) THEN
        DH=AMAX1(ZI-BC,HM1)
        UN=C1*SQRT(DH)
        FLR(1)=DH*UN
        FLR(2)=FLR(1)*UN
        HC1=ZC-BI
        FLR(3)=FLR(1)*QL(3)
        FLR(4)=4.905*HC1*HC1
        ELSE
        FLR(1)=-HC*C1*SQRT(HC)
        FLR(2)=HI*QL(2)*QL(2)
        FLR(3)=0
        FLR(4)=4.905*HI*HI
        END IF
         DO K=1,4
         FLUX(K,J,I)=FLR(K)
         ENDDO
        KKNV(J)=1
        END IF
        ENDIF
C
        IF(KKNV(J).EQ.0)THEN      
60      QR(1)=AMAX1(ZC-BI,HM1)
        DIS=AMIN1(HC/QR(1),1.5)
        UR=UC*COSA+VC*SINA
        QR(2)=UR*DIS
        IF (HC.LE.HM2.OR.QR(1).LE.HM2) QR(2)=SIGN(VMIN,UR)
        QR(3)=VC*COSA-UC*SINA
        CALL OSHER
        FLR(2)=FLR(2)+(1-DIS)*HC*UR*UR/2
C
CC      IF (KP.EQ.9.OR.KP.EQ.20) GOTO 90
!!!        DO 70 K=1,4
!!!70      IF (NAC(K,NC).EQ.I) GOTO 80
!!!80      IF (NC.EQ.0.OR.K.GT.4) GOTO 90
!!!        FLUX(1,K,NC)=-FLR(1)
!!!        FLUX(2,K,NC)=FLR(2)
!!!        FLUX(3,K,NC)=FLR(3)
!!!        ZA=SQRT(FLR(4)/4.905)+BI
!!!        HC=0
!!!        IF (ZA.GT.BC) HC=ZA-BC
!!!        FLUX(4,K,NC)=4.905*HC*HC
C
	DO 95 K=1,4
        FLUX(K,J,I)=FLR(K)
95    CONTINUE
      ENDIF
      
100     CONTINUE      
        END
C#####################################################################  ����4���ߵ�ͨ��        
C#####################################################################  ����߽絥Ԫ�ߵ�ͨ��
        SUBROUTINE BOUNDA(J,I)
C        BOUNDARY TREATMENT OF 2D BOUNDARY CELLS
        INCLUDE 'PAR.INC'
        !INCLUDE 'TIME.INC'        
        !INCLUDE 'GIRD.INC'
        !INCLUDE 'BOUNDARY.INC'
        !INCLUDE 'DEPTH.INC'        
        INCLUDE 'BASE.INC' 
        
        DIMENSION WZ(NHQ), WQ(NHQ), QB(NZ)        
        DATA S0/0.0002/,DX2/5000/,BRDTH/100/

        
        !SF(U,H,FNI)=FNI*U*ABS(U)/H**1.33333
C        FOR SUPERCRITICAL FLOW
        IF (QL(2).GT.CL) THEN
        FLR(1)=HI*QL(2)
        FLR(2)=FLR(1)*QL(2)
        FLR(4)=4.905*HI*HI
        FLR(3)=FLR(1)*QL(3)
        RETURN
        END IF
C
        FLR(3)=0
        IF (QL(2).GT.0.) FLR(3)=HI*QL(2)*QL(3)
C
C       GIVEN DISCHARGE BOUNDARY CONDITION KP=10
        IF (KP.EQ.10) THEN
          CALL CHOICE (NQ,MBQ,I,II)
          FLR(1)=-(QT(JT,II)+DQT(II)*KT)
          FLR(1)=FLR(1)/SIDE(J,I)
          QB2=FLR(1)*FLR(1)
          HB0=HI
        DO  K=1,20
         W=FIL-FLR(1)/HB0
         HB=W*W/39.24
         IF (ABS(HB0-HB).LE.0.005) EXIT
         HB0=HB0*0.5+HB*0.5
        ENDDO
20       IF (HB.LE.1) THEN
          FLR(2)=0
           ELSE
          FLR(2)=QB2/HB
         END IF
        FLR(3)=0
        FLR(4)=4.905*HB*HB
        RETURN
        END IF
C
C        GIVEN RATING CURVE AS BOUNDARY CONDITION
        IF (KP.EQ.3) THEN
        HR0=HI
        CALL CHOICE (NZQ,MBZQ,I,II)
        DO 30 I1=1,NHQ
        WZ(I1)=ZW(I1,II)
30      WQ(I1)=QW(I1,II)
        DO 40 I1=1,20
        ZR0=HR0+BI
        CALL LAQP (ZR0,CQ,WZ,WQ,NHQ)
        W=FIL-CQ/HR0
        HR=W*W/39.24
        IF (ABS(HR-HR0).LE.0.001) EXIT
        HR0=HR
40      CONTINUE
50      FLR(1)=CQ
        FLR(2)=CQ*CQ/HR
        HB=(HI+HR)/2
        FLR(4)=4.905*HB*HB
        RETURN
        END IF

CC      GIVEN WATER LEVEL BOUNDARY CONDITION (INFLOW)
**************************************************************  ԭˮλ�߽�DBZ������ʽ
!        IF (KP.EQ.1.AND.QL(2).LT.0.) THEN                       
!        CALL CHOICE (NZ,MBZ,I,II)
!        FIA=QL(2)+6.264*SQRT(HI)
!c        if(I==4789)write(*,*)Z1(I),ZT(JT,II)       
!        FNI=FNC(I)
!c        FNI=TG*(FNC0(I)+0.01/H1(I))*(FNC0(I)+0.01/H1(I))
!C        DZB=S0*DX2
!        TGS0=TG*S0
!        IF (JT.EQ.JT0.AND.KT.EQ.1) QB(II)=HI*QL(2)
!c        HB=ZT(JT,II)+DZT(II)*(KT-0.5)-BI+DZB
!        HB=ZT(JT,II)+DZT(II)*KT-BI             !+DZB
!C        HB=ZT(JT,II)+DZT(II)*KT-BI
!        HB1=HB+DZT(II)
!        UB=QB(II)/HB
!        DX1=DT*(UB+SQRT(9.81*HB))
!        ALP=DX1/DX2
!        HA=AMAX1(ALP*HI+(1-ALP)*HB,HM1)
!        UA=(ALP*HI*QL(2)+(1-ALP)*QB(II))/HA
!        SFA=SF(UA,HA,FNI)
!        UB1=UA-6.264*(SQRT(HB1)-SQRT(HA))+TGS0-SFA
!        SFB=SF(UB1,HB1,FNI)
!        UB1=UB1+(SFA-SFB)/2
!        CQ=HB1*UB1
!        FLR(1)=(QB(II)+CQ)/2
!        QB(II)=CQ
!c        FLR(1)=HB*QL(2) 
!        FLR(2)=FLR(1)*(UB+UB1)/2
!c        FLR(2)=FLR(1)*QL(2)
!        HB=(HB+HB1)/2                            !-DZB
!        FLR(4)=4.905*HB*HB
!        RETURN
!        END IF
!!C
!        IF (KP.EQ.1.AND.QL(2).GE.0.) THEN
!        CALL CHOICE(NZ,MBZ,I,II)
!        HB1=ZT(JT,II)+DZT(II)*KT-BI-DZB
!        URB=HI*QL(2)/HB1
!c        URB=6.264*SQRT(HI)*QL(2)/HB1
!c        URB=QL(2)
!        FLR(1)=HB1*URB
!        FLR(2)=FLR(1)*URB
!        FLR(4)=4.905*HB1*HB1
!c        if(I==5061)write(*,*)I,Z1(I),ZT(JT,II),HI,HB1!,AREA(I)
!c        if(I==4789)write(*,*)I,Z1(I),ZT(JT,II),HI,HB1!,AREA(I)
!        RETURN
!        END IF
**************************************************************  ԭˮλ�߽�DBZ������ʽ

**************************************************************  ��ˮλ�߽�ѭ��������ʽ
        IF (KP.EQ.1) THEN
        CALL CHOICE (NZ,MBZ,I,II)
        HB1=ZT(JT,II)+DZT(II)*KT-BI
        FIAL=QL(2)+6.264*SQRT(HI)
        UR0 = QL(2)
        URB = UR0
        DO IURB=1,30
        FIAR= URB-6.264*SQRT(HB1)
        URB = (FIAL+FIAR)*(FIAL-FIAR)*(FIAL-FIAR)/HB1/313.92       
        IF (ABS(URB-UR0).LE.0.0001) THEN
        EXIT 
        ENDIF
        UR0 = URB
        ENDDO
        FLR(1)=HB1*URB
        FLR(2)=FLR(1)*URB
        FLR(4)=4.905*HB1*HB1
        RETURN
        ENDIF
**************************************************************  ��ˮλ�߽�ѭ��������ʽ       

C       LAND BOUNDARY CONDITION
60      IF (KP.EQ.4) THEN
        FLR(1)=0
        FLR(2)=0
        FLR(3)=0
        FLR(4)=4.905*HI*HI
        RETURN
        END IF
C
C       FREE UNIFORM OUTFLOW BOUNDARY CONDITION
        IF (KP.EQ.5) THEN
        QL(2)=AMAX1(QL(2),0.)
        FLR(1)=HI*QL(2)
        FLR(2)=FLR(1)*QL(2)
        FLR(4)=4.905*HI*HI
        RETURN
        END IF
C
C       WEIR BOUNDARY CONDITION
        IF (KP.EQ.6) THEN
        NE=I
        IF (NC.NE.0) NE=MIN0(I,NC)
        CALL CHOICE (NWE,MBW,NE,I1)
        TOP=TOPW(I1)
        IF (ZI.LT.TOP.OR.ZC.LT.TOP) THEN
        KP=4
        !GOTO 60
          FLR(1)=0
          FLR(2)=0
          FLR(3)=0
          FLR(4)=4.905*HI*HI
        END IF
        IF (ZI.GT.TOP.AND.ZC.LT.TOP) THEN
        FLR(1)=C0*(ZI-TOP)**1.5
        FLR(2)=FLR(1)*QL(2)
        FLR(3)=FLR(1)*QL(3)
        FLR(4)=4.905*(TOP-BI)*(TOP-BI)
        RETURN
        END IF
        IF (ZI.LT.TOP.AND.ZC.GT.TOP) THEN
        FLR(1)=-C0*(ZC-TOP)**1.5
        FLR(2)=FLR(1)*AMIN1(UC*COSA+VC*SINA,0.)
        FLR(3)=FLR(1)*(VC*COSA-UC*SINA)
        FLR(4)=4.905*(ZI-BI)*(ZI-BI)
        RETURN
        END IF
        DZ=ABS(ZI-ZC)
        IF (ZI.LE.ZC) THEN
        HD=ZI-TOP
        UN=AMIN1(UC*COSA+VC*SINA,0.)
        VT=VC*COSA-UC*SINA
        ELSE
        HD=ZC-TOP
        UN=AMAX1(QL(2),0.)
        VT=QL(3)
        END IF
        SH=HD+DZ
        CE=AMIN1(1.0,1.05*(DZ/SH)**0.33333)
        IF (ZI.LT.ZC.AND.UN.GT.0.) UN=0
        FLR(1)=SIGN(CE*C1*SH**1.5,ZI-ZC)
        FLR(2)=FLR(1)*ABS(UN)
        FLR(3)=FLR(1)*VT
        FLR(4)=4.905*(TOP-BI)*(TOP-BI)
        RETURN
        END IF
C
C       POSSIBLY BREACHED DYKE
        IF (KP.EQ.7) THEN
          CALL CHOICE (NDI,MDI,I,I1)
          TOP=TOPD(I1)
         IF (ZI.GT.TOP.OR.ZC.GT.TOP) THEN
            KP=0
            KLAS(J,I)=0
            CQ=QD(ZI,ZC,TOP)
            CB=BRDTH/SL
            FLR(1)=CQ*CB
            FLR(2)=CB*SIGN(CQ*CQ/HB,CQ)
            FLR(4)=4.905*HB*HB
            RETURN
          ELSE
            KP=4
          FLR(1)=0
          FLR(2)=0
          FLR(3)=0
          FLR(4)=4.905*HI*HI
          END IF
        END IF
        RETURN
        END


        SUBROUTINE LAQP(X,Y,A,B,MS)
C       ONE-VARIABLE QUADRATIC INTERPOLATION
        PARAMETER (NHQ=5)
        DIMENSION A(NHQ),B(NHQ)
         ILAQ=0
        DO 10 I=1,MS-3
        !IF (X.LT.A(I+1)) GOTO 20
         IF (X.LT.A(I+1)) THEN
             ILAQ=1
             EXIT
         ENDIF
10      CONTINUE
        IF (ILAQ==1) THEN        
        IF (I.GT.1.AND.X-A(I).LT.A(I+1)-X) I=I-1
        X0=A(I)
        X1=A(I+1)
        X2=A(I+2)
        IF (ABS(X0-X1).LT.0.01.OR.ABS(X1-X2).LT.0.01) THEN
        Y=B(I+1)
        ELSE
        !GOTO 30
        !END IF
        U=(X-X1)*(X-X2)/(X0-X1)/(X0-X2)
        V=(X-X0)*(X-X2)/(X1-X0)/(X1-X2)
        W=(X-X0)*(X-X1)/(X2-X0)/(X2-X1)
        Y=U*B(I)+V*B(I+1)+W*B(I+2)
        END IF
        RETURN
        ELSE
        I=MS-2
20      IF (I.GT.1.AND.X-A(I).LT.A(I+1)-X) I=I-1
        X0=A(I)
        X1=A(I+1)
        X2=A(I+2)
        IF (ABS(X0-X1).LT.0.01.OR.ABS(X1-X2).LT.0.01) THEN
        Y=B(I+1)
        ELSE
        !GOTO 30
        !END IF
        U=(X-X1)*(X-X2)/(X0-X1)/(X0-X2)
        V=(X-X0)*(X-X2)/(X1-X0)/(X1-X2)
        W=(X-X0)*(X-X1)/(X2-X0)/(X2-X1)
        Y=U*B(I)+V*B(I+1)+W*B(I+2)
        END IF
30      RETURN
        ENDIF
        END
C#####################################################################  ����߽絥Ԫ�ߵ�ͨ��
C#####################################################################  ���ұ߽絥Ԫ
        SUBROUTINE CHOICE(N,MB,NE,NJ)
C        LOOK UP AN ITEM IN A LIST
C       PARAMETER (NMAX=24)
        DIMENSION MB(N)
        NA=1
        NB=N
        IF (MB(NA).EQ.NE) THEN
        NJ=NA
        RETURN
        END IF
        IF (MB(NB).EQ.NE) THEN
        NJ=NB
        RETURN
        END IF
        DO WHILE (.TRUE.)
10      NJ=(NA+NB)/2
        ME=MB(NJ)
        IF (NE.EQ.ME) RETURN
        IF (NE.GT.ME) THEN
        NA=NJ+1
        IF (MB(NA).EQ.NE) THEN
        NJ=NA
        RETURN
        END IF
        ELSE
        NB=NJ-1
        IF (MB(NB).EQ.NE) THEN
        NJ=NB
        RETURN
        END IF
        END IF
        !GOTO 10
         END DO
        END

!        FUNCTION MOD(I,J)
!C        MODULO-J COUNTER
!        MOD=I-J*(I/J)
!        END

        FUNCTION QD(ZL,ZR,ZB)
C        COMPUTING DISCHARGE THROUGH BREACHED DIKE (FREE OR SUBMERGED)
C        ZL--LEFT-SIDE WATER LEVEL, ZR--RIGHT-SIDE, ZB--BOTTOM ELEVATION
C        SIGMA=(2*CM*CM)**0.33333, FI=CM*SQRT(2*G)
	PARAMETER (CM=0.384, SIGMA=0.667, FI=4.43)
	ZU=AMAX1(ZL,ZR)
	ZD=AMIN1(ZL,ZR)
	H0=ZU-ZB
	HS=ZD-ZB
	DELTA=HS/H0
	IF (DELTA.LE.SIGMA) THEN
	QD=SIGN(CM*H0**1.5, ZL-ZR)
	ELSE
	DH=ZU-ZD
	IF (DH.GT.0.09) THEN
	QD=SIGN(FI*HS*SQRT(DH),ZL-ZR)
	ELSE
	QD=SIGN(FI*HS*0.3*DH/0.1,ZL-ZR)
	END IF
	END IF
        END
C#####################################################################  ���ұ߽絥Ԫ
C#####################################################################  ������ˮ����Ԫͨ��
C#####################################################################
        SUBROUTINE OSHER
C        FLUX COMPUTATION BASED ON OSHER SCHEME FOR 16 CASES
        INCLUDE 'PAR.INC'
        !INCLUDE 'TIME.INC'        
        !INCLUDE 'GIRD.INC'
        !INCLUDE 'BOUNDARY.INC'
        !INCLUDE 'DEPTH.INC'        
        INCLUDE 'BASE.INC'       

        HL=QL(1)
        HR=QR(1)
        UL=QL(2)
        UR=QR(2)
        VL=QL(3)
        VR=QR(3)
        CR=SQRT(9.81*HR)
        FIR=UR-2*CR
        UA=(FIL+FIR)/2
        CA=ABS((FIL-FIR)/4)
        DO  I=1,4
        FLR(I)=0
        ENDDO
C
        IF (UL.LT.CL.AND.UR.GE.-CR) THEN
        K1=1
        !GOTO 10
        ELSE
        IF (CA.LT.UA) THEN
        K2=10
        !GOTO 20
        END IF
        IF (UA.GE.0.0.AND.UA.LT.CA) THEN
        K2=20
        !GOTO 20
        END IF
        IF (UA.GE.-CA.AND.UA.LT.0.0) THEN
        K2=30
        !GOTO 20
        END IF
        IF (UA.LT.-CA) K2=40                               
        END IF       
        
        
        IF (UL.GE.CL.AND.UR.GE.-CR) THEN
        K1=2
        !GOTO 10
        ELSE
        IF (CA.LT.UA) THEN
        K2=10
        !GOTO 20
        END IF
        IF (UA.GE.0.0.AND.UA.LT.CA) THEN
        K2=20
        !GOTO 20
        END IF
        IF (UA.GE.-CA.AND.UA.LT.0.0) THEN
        K2=30
        !GOTO 20
        END IF
        IF (UA.LT.-CA) K2=40        
        END IF
        
        IF (UL.LT.CL.AND.UR.LT.-CR) THEN
        K1=3
        !GOTO 10
        ELSE
        IF (CA.LT.UA) THEN
        K2=10
        !GOTO 20
        END IF
        IF (UA.GE.0.0.AND.UA.LT.CA) THEN
        K2=20
        !GOTO 20
        END IF
        IF (UA.GE.-CA.AND.UA.LT.0.0) THEN
        K2=30
        !GOTO 20
        END IF
        IF (UA.LT.-CA) K2=40 
        END IF
        
        IF (UL.GE.CL.AND.UR.LT.-CR) THEN
        K1=4
        ELSE
        IF (CA.LT.UA) THEN
        K2=10
        !GOTO 20
        END IF
        IF (UA.GE.0.0.AND.UA.LT.CA) THEN
        K2=20
        !GOTO 20
        END IF
        IF (UA.GE.-CA.AND.UA.LT.0.0) THEN
        K2=30
        !GOTO 20
        END IF
        IF (UA.LT.-CA) K2=40 
        ENDIF
        
20      K12=K1+K2        

C
        IF (K12.EQ.11) THEN
        CALL QS(2,1)
        RETURN
        END IF
C
        IF (K12.EQ.21) THEN
        CALL QS(3,1)
        RETURN
        END IF
C
        IF (K12.EQ.31) THEN
        CALL QS(5,1)
        RETURN
        END IF
C
        IF (K12.EQ.41) THEN
        CALL QS(6,1)
        RETURN
        END IF
C
        IF (K12.EQ.12) THEN
        CALL QS(1,1)
        RETURN
        END IF
C
        IF (K12.EQ.22) THEN
        CALL QS(1,1)
        CALL QS(2,-1)
        CALL QS(3,1)
        RETURN
        END IF
C
        IF (K12.EQ.32) THEN
        CALL QS(1,1)
        CALL QS(2,-1)
        CALL QS(5,1)
        RETURN
        END IF
C
        IF (K12.EQ.42) THEN
        CALL QS(1,1)
        CALL QS(2,-1)
        CALL QS(6,1)
        RETURN
        END IF
C
        IF (K12.EQ.13) THEN
        CALL QS(2,1)
        CALL QS(6,-1)
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.23) THEN
        CALL QS(3,1)
        CALL QS(6,-1)
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.33) THEN
        CALL QS(5,1)
        CALL QS(6,-1)
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.43) THEN
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.14) THEN
        CALL QS(1,1)
        CALL QS(6,-1)
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.24) THEN
        CALL QS(1,1)
        CALL QS(2,-1)
        CALL QS(3,1)
        CALL QS(6,-1)
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.34) THEN
        CALL QS(1,1)
        CALL QS(2,-1)
        CALL QS(5,1)
        CALL QS(6,-1)
        CALL QS(7,1)
        RETURN
        END IF
C
        IF (K12.EQ.44) THEN
        CALL QS(1,1)
        CALL QS(2,-1)
        CALL QS(7,1)
        END IF
        END

C#####################################################################
C#####################################################################  ������ˮ����Ԫͨ��
C#####################################################################
C##################################################################### ������ˮ����Ԫͨ��
        SUBROUTINE QS(K,J)
C        COMPUTATION OF FLUX BASED ON LEFT, RIGHT OR INTERMEDIATE STATE
        INCLUDE 'PAR.INC'
        !INCLUDE 'TIME.INC'        
        !INCLUDE 'GIRD.INC'
        !INCLUDE 'BOUNDARY.INC'
        !INCLUDE 'DEPTH.INC'        
        INCLUDE 'BASE.INC'         
        
        DIMENSION F(4)
C
        IF (K.EQ.1) THEN
        CALL QF(HL,UL,VL,F)
         DO I=1,4
         FLR(I)=FLR(I)+F(I)*J
         ENDDO
        END IF
C
        IF (K.EQ.2) THEN
        US=FIL/3
        HS=US*US/9.81
        CALL QF(HS,US,VL,F)
         DO I=1,4
         FLR(I)=FLR(I)+F(I)*J
         ENDDO
        END IF
C
        IF (K.EQ.3) THEN
        FIL=FIL-UA
        HA=FIL*FIL/39.24
        CALL QF(HA,UA,VL,F)
         DO I=1,4
         FLR(I)=FLR(I)+F(I)*J
         ENDDO
        END IF
C
        IF (K.EQ.5) THEN
        FIR=FIR-UA
        HA=FIR*FIR/39.24
        CALL QF(HA,UA,VR,F)
         DO I=1,4
         FLR(I)=FLR(I)+F(I)*J
         ENDDO
        END IF
C
        IF (K.EQ.6) THEN
        US=FIR/3
        HS=US*US/9.81
        CALL QF(HS,US,VR,F)
         DO I=1,4
         FLR(I)=FLR(I)+F(I)*J
         ENDDO
        END IF
C
        IF (K.EQ.7) THEN
        CALL QF(HR,UR,VR,F)
         DO I=1,4
         FLR(I)=FLR(I)+F(I)*J
         ENDDO
        END IF
C
!10      DO 20 I=1,4
!20      FLR(I)=FLR(I)+F(I)*J
        END

C#####################################################################
C#####################################################################
        SUBROUTINE QF(H,U,V,F)
C        COMPUTATION OF FLUX COMPONENTS
        DIMENSION F(4)
        F(1)=H*U
        F(2)=F(1)*U
        F(3)=F(1)*V
        F(4)=4.905*H*H
        END

C##################################################################### ������ˮ����Ԫͨ��
C#####################################################################
C##########################################################################
C                                                                         #
C##########################################################################  ˮ��ģ��
        SUBROUTINE TWODWQ 
C        MAIN CONTROL OF 2D WATER BODY COMPUTATION
        INCLUDE 'PAR.INC'
        !INCLUDE 'TIME.INC'        
        !INCLUDE 'GIRD.INC'
        !INCLUDE 'BOUNDARY.INC'
        !INCLUDE 'DEPTH.INC'        
        INCLUDE 'BASE.INC' 
C


C
        END

C#####################################################################
C#####################################################################

C##########################################################################
C                                                                         #
C##########################################################################   ��ɳģ��
        SUBROUTINE SEDIMENTS
C       FLUX BALANCE FOR 2D CELLS
        INCLUDE 'PAR.INC'
        !INCLUDE 'TIME.INC'        
        !INCLUDE 'GIRD.INC'
        !INCLUDE 'BOUNDARY.INC'
        !INCLUDE 'DEPTH.INC'        
        INCLUDE 'BASE.INC'       
        
        DIMENSION SOURCE(CEL,MSORT),WHS(CEL,MSORT),ALF(CEL,MSORT)
        DIMENSION DIFFUSE(CEL,MSORT),SXCYC(4)
        DATA GMAR/1690/
        DATA SK,SM/0.08,0.5/   !SK=0.01                   SK Юɳ��ϵ��
        DATA FKN/0.010/
C

        END
C##########################################################################  ���ģ��
        SUBROUTINE OUTPUT2
C       STORAGE AND OUTPUT OF COMPUTED RESULTS
        INCLUDE 'PAR.INC'       
        INCLUDE 'BASE.INC' 
        
C        OUTPUT FOR 2D WATER BODIES
        DO 10 I=1,CEL
        IF (H2(I).LE.HM1) THEN
        H2(I)=0
        Z2(I)=0
C       Z2(I)=ZBC(I)
        W2(I)=0
        FI2(I)=0
        !GOTO 10
        ELSE
        W1(I)=SQRT(U1(I)*U1(I)+V1(I)*V1(I))
        W2(I)=SQRT(U2(I)*U2(I)+V2(I)*V2(I))
        FI2(I)=FI(U2(I),V2(I))*57.298
        END IF        
10      CONTINUE
        WRITE (28,100) JT,KT,INT(DT),NDT*KT/3600,NSF,NSF2,ALPHA,CQL,
     1	 INE
        WRITE (20,100) JT,KT,INT(DT),NDT*KT/3600,NSF,NSF2,ALPHA,CQL,
     1	 INE
C        WRITE (*,*) CEL      
        WRITE (28,200) (H2(I),I=1,CEL)         
        WRITE (20,300) (Z2(I),I=1,CEL)
        WRITE (28,400) (U2(I),I=1,CEL)
        WRITE (28,500) (V2(I),I=1,CEL)
        WRITE (20,600) (W2(I),I=1,CEL)
        WRITE (20,700) (FI2(I),I=1,CEL)
100     FORMAT (' ',//1X,3HJT=,I5, 2X,3HKT=,I5, 2X,3HDT=,I3,1X,3HSEC,
     1  2X,2HT=,I2,1X,1HH, 2X,4HNSF=,I2,1H/,I2, 2X,4HWEC=,F4.2,1H/,
     2	2X,4HCQL=,F4.2, 2X,4HINE=,I1)
200     FORMAT (' ',/5X, 3HH2=, 10(/5X, 10F7.2))
300     FORMAT (' ',/5X, 3HZ2=, 10(/5X, 10F7.2))
400     FORMAT (' ',/5X, 3HU2=, 10(/5X, 10F7.2))
500     FORMAT (' ',/5X, 3HV2=, 10(/5X, 10F7.2))
600     FORMAT (' ',/5X, 3HW2=, 10(/5X, 10F7.2))
700     FORMAT (' ',/5X, 3HFI=, 10(/5X, 10F7.2))
        WRITE (*,*) 'JT=',JT,' KT=',KT
********************************************************���tecplot��ʽ 
        !OPEN(20,FILE='..\OUTPUT\HUV.DAT')  Z2W2FI2
        !OPEN(24,FILE='..\OUTPUT\SHA.DAT')
        !OPEN(28,FILE='..\OUTPUT\H2U2V2.DAT')    H2U2V2    
       write(901,*) ' VARIABLES = "X", "Y", "H2", "Z2","U2","V2","W2"'       
       write(901,*) 'ZONE N=24889, E=24020, DATAPACKING=BLOCK, 
	1ZONETYPE=FEQUADRILATERAL'
       write(901,*) 'VARLOCATION=([3-7]=CELLCENTERED)'
     	 write(901,*) (XP(I)+XIMIN,I=1,NOD)
       write(901,*) (YP(I)+YIMIN,I=1,NOD)
       write(901,*) (H2(I),I=1,CEL)
       write(901,*) (Z2(I),I=1,CEL)
       write(901,*) (U2(I),I=1,CEL)
       write(901,*) (V2(I),I=1,CEL)
       write(901,*) (W2(I),I=1,CEL) 
	DO I = 1, CEL
      write (901,*) (NAP(K,I),K=1,4)
      ENDDO     
********************************************************���tecplot��ʽ        
      !WRITE(*,*) 'CEL=7657  ','U2=',U2(7657),'  V2=',V2(7657)    
      !WRITE(*,*) 'CEL=4093  ','U2=',U2(4093),'  V2=',V2(4093)          
        
        END

!        SUBROUTINE OUTPUTSQ
!C        STORAGE AND OUTPUT OF COMPUTED RESULTS
!        INCLUDE 'CJ.INC'
!        INCLUDE 'CJ2.INC'
!        
!        WRITE(24,*) '######## SHA ##########','JT=',JT,'KT=',KT
!        WRITE(24,1000) (SHASUM(I),I=1,CEL)
!        WRITE(24,*) '######## ZBC ##########','JT=',JT,'KT=',KT
!        WRITE(24,1010) (ZBC2(I),I=1,CEL)
!1000    FORMAT(1X,10F10.3) 
!1010    FORMAT(1X,10F10.5)
!********************************************************���tecplot��ʽ        
!       write(902,*) ' VARIABLES = "X", "Y", "SHASUM", "ZBC2", "SXTWO"'       
!       write(902,*) 'ZONE N=24889, E=24020, DATAPACKING=BLOCK, 
!	1ZONETYPE=FEQUADRILATERAL'
!       write(902,*) 'VARLOCATION=([3-5]=CELLCENTERED)'
!     	 write(902,*) (XP(I)+523000.,I=1,NOD)
!       write(902,*) (YP(I)+3268000.,I=1,NOD)
!       write(902,*) (SHASUM(I),I=1,CEL)
!       write(902,*) (ZBC2(I),I=1,CEL)
!       write(902,*) (SXTWO(I,1),I=1,CEL)
!	DO I = 1, CEL
!      write (902,*) (NAP(K,I),K=1,4)
!      ENDDO     
********************************************************���tecplot��ʽ   
        !END

        FUNCTION FI(X,Y)
C        COMPUTATION OF DIRECTIONAL ANGLE OF VELOCITY VECTOR FOR 2D CELLS
        !IF (X*Y.NE.0.0) GOTO 10
        IF (X*Y.NE.0.0) THEN
10      W=ATAN2(ABS(Y),ABS(X))
        !IF (X*Y.GT.0.0) GOTO 30
        !IF (Y.GT.0.0) GOTO 50
        IF (X*Y.GT.0.0) THEN
          IF (X.GT.0.0) THEN
          FI=W
          ELSE
          FI=3.14159+W
          ENDIF
        ELSE
          IF (Y.GT.0.0)THEN
          FI=3.14159-W
          ELSE
          FI=6.2832-W    
          ENDIF
        ENDIF               
        ELSE            
        IF (X.EQ.0.0.AND.Y.GE.0.0) FI=1.5708
        IF (X.EQ.0.0.AND.Y.LT.0.0) FI=4.7124
        IF (Y.EQ.0.0.AND.X.GE.0.0) FI=0.0
        IF (Y.EQ.0.0.AND.X.LT.0.0) FI=3.1416
        RETURN
        ENDIF     
        END
C#####################################################################
C#####################################################################
