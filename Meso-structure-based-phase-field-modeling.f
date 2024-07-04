C ================================================================================================ C
C User Subroutine UEL for Abaqus: two elements for the split scheme operator phase and displacement problem:
C     Type 2: C2D8 displacement rectangular element
C     Type 4: C2P4 phase-field rectangular element
C ================================================================================================ C
C Material properties to be given through the input file (*.inp):
C
C For Type 2 element (stress-strain):
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = Yield stress (sig_y)
C PROPS(4) = Hardening modulus (H)
C PROPS(5) = Critical plastic strain (eps_pl_crit)
C PROPS(6) = Thickness of the element (t)
C PROPS(7) = Density (rho)
C PROPS(8) = Aniso energy degradation switch (if 1 - yes, 0 - no)
C PROPS(9) = Plasticuty switch (if 1 - yes, 0 - no)
C PROPS(10) = Length scale parameter (lc)
C PROPS(11) = Crack surface energy (gc)
C
C For Type 4 element (phase field):
C PROPS(1) = Length scale parameter (lc)
C PROPS(2) = Crack surface energy (gc)
C PROPS(3) = Thickness of the element (t)
C PROPS(4) = Elastic switch  (if 1 - yes, 0 - no)
C
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C            by 3 - (N_phase+N_stress+N_UMAT)/3 (to be changed for each model)
C
C NSTVTT - solution dependent variables for the displacement element
C            (displacements, strains, stresses,energies, phase, etc.)
C NSTVTO - solution dependent variables for the phase-field element
C            (phase, energy history)
C NSTV - overall solution dependent variables (NSTVTO+NSTVTT+4), where
C           the additional 4 variables are the: time and iteration number
C
C ================================================================================================ C
C Comments on solution dependent variables
C Stress/strain element
C SVARS(1-112): SDV(1-28)x4(or 1)
C               SDV(1) - X translation
C               SDV(2) - Y translation
C               SDV(3) - X normal strain
C               SDV(4) - Y normal strain
C               SDV(5) - XY engineering shear strain
C               SDV(6-9) - Elastic strains (x, y, z, xy)
C               SDV(10-13) - Plastic strains
C               SDV(14) - Eq. plastic strain
C               SDV(15-18) - Stresses 
C               SDV(19) - Hydrostatic stress
C               SDV(20) - von Mises stress
C               SDV(21) - plastic energy
C               SDV(22) - tensile elastic energy
C               SDV(23) - potential strain energy
C               SDV(24) - phase-field
C               SDV(25-28) - El. strain at the beginning of the step
C               SVARS(113-120): RHS(t_{n-1}) - 8(or 6) components, previous internal
C                                    force vector for HHT (dynamic) for each DOF
C
C Phase-field element
C SVARS(1-8):   SDV(1-2)x4(or 1)
C               SDV(1) - phase-field
C               SDV(2) - history energy
C               SVARS(9-12): RHS(t_{n-1}), 4(or 3) components previous internal
C                                    force vector for HHT (dynamic) for each DOF
C
C ================================================================================================ C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     -------------------------------------------------------------------------------------------- C
      INCLUDE 'ABA_PARAM.INC'
C     -------------------------------------------------------------------------------------------- C
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
     2 TEN=10.D0,PHCALCMAX=0.95D0,DEPSCR=0.1D0,DENGMAX=1,DTMIN=1.0D-9,
     3 N_ELEM=100000,NSTVTT=28,NSTVTO=2,NSTV=34)

      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
     
      INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY,NALL,CNT

      REAL*8 AINTW(4),XII(4,2),XI(2),dNdxi(NNODE,2),
     1    VJACOB(2,2),dNdx(NNODE,2),VJABOBINV(2,2),AN(4),BP(2,NDOFEL),
     2    DP(2),SDV(NSTV),BB(3,NDOFEL),CMAT(3,3),EPS(3),STRESS(3),
     3    VNI(2,NDOFEL),ULOC(2),PHASENOD(NNODE),AMASS(NDOFEL,NDOFEL),
     4    EIGV(3),ALPHAI(3),CMATG(6,6),VECTI(3),EPSZ(4),
     5    EPSC(6),ASTIFF(NDOFEL,NDOFEL),RHSINI(NDOFEL),
     6    RHSK(NDOFEL),DEPS(3),EELAS(4),EPLAS(4),SFULL(4),FLOW(4),
     7    CMATP(4,4),EPSP(4),DEPLAS(4),DDSDDEEQ(4,4),
     8    FLOWG(4),PHMAX(4),PNEWDTIP(4)

      REAL*8 DTM,THCK,HIST,CLPAR,GCPAR,EMOD,ENU,PARK,ENG,ENGK,ENGEG,
     1 ENGKG,ENGDG,ENGPG,ENGD,PLSWT,ANISOSWT,EQPLAS,YSHARD,YIELDS,
     2 DLAMB,VNEVEZ,ALPHA,ENGPL,REITER
C
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
C     ============================================================================================ C
C     History variables
      ENGKG = ZERO
      ENGEG = ZERO
      ENGDG = ZERO
      ENGPG = ZERO
      NALL  = NSTVTO+NSTVTT  
C     ******************************************************************************************** C
C     ******************************************************************************************** C
C     Constructing elemet TYPE 2 (stress/strain - displacement)
C     ******************************************************************************************** C
C     ******************************************************************************************** C
      IF ((JTYPE.EQ.TWO)) THEN
C     ============================================================================================ C
C     Time an iteration variables
C     ============================================================================================ C
      IF (TIME(2).EQ.ZERO) THEN
          TIMEZ    = -999.D0
          TIMEZLOC = -99.D0
          DO K2=1,NSTV
             DO K3=1,4
                USRVAR(JELEM,K2,K3) = ZERO         ! 34 variables, 4 IPs
             END DO
          END DO
      ELSE
          TIMEZ    =  USRVAR(JELEM,NALL+1,1)       ! 31
          TIMEZLOC =  TIME(2)-DTIME                ! 
      ENDIF
C
      DTZERO = USRVAR(JELEM,NALL+2,1)
C
      IF (TIMEZ.LT.TIMEZLOC) THEN                  ! New iteration
          USRVAR(JELEM,NALL+1,1) = TIMEZLOC        ! 31 - TIMEZLOC   
          USRVAR(JELEM,NALL+2,1) = DTIME           ! 32 - DTIME 
          USRVAR(JELEM,NALL+3,1) = ZERO            ! 33 - 0
          USRVAR(JELEM,NALL+4,1) = ZERO            ! 34 - 0
      ELSE
          IF (DTZERO.GT.DTIME*(ONE+TOLER)) THEN
C         -----   New correcting iteration   -----
              USRVAR(JELEM,NALL+2,1) = DTIME
              USRVAR(JELEM,NALL+3,1) = USRVAR(JELEM,NALL+3,1)+ONE
              USRVAR(JELEM,NALL+4,1) = ZERO
          ELSE
C         -----   New local step   -----
              USRVAR(JELEM,NALL+4,1) = USRVAR(JELEM,NALL+4,1)+ONE
          ENDIF
      ENDIF      
      REITER   = USRVAR(JELEM,NALL+3,1)
      STEPITER = USRVAR(JELEM,NALL+4,1)
C     ============================================================================================ C
C     Additional plasticity control
C     ============================================================================================ C
      IF ((REITER.EQ.ZERO).AND.(STEPITER.EQ.ZERO)) THEN
          PLSWTGLOB = ONE
          USRVAR(JELEM,NALL+1,2) = PLSWTGLOB
      ELSE
          PLSWTGLOB = USRVAR(JELEM,NALL+1,2)
      ENDIF
C       
C     ============================================================================================ C
C     Material parameters
C     ============================================================================================ C
      EMOD     = PROPS(1)
      ENU      = PROPS(2)
      YIELDS   = PROPS(3)
      VHMOD    = PROPS(4)
      EQPLASCR = PROPS(5)
      THCK     = PROPS(6)
      DENS     = PROPS(7)
      ANISOSWT = PROPS(8)
      PLSWT    = PROPS(9)
      CLPAR    = PROPS(10)
      GCPAR    = PROPS(11)
C      
      PARK     = TOLER
      ELAMEL   = EMOD*ENU/((ONE+ENU)*(ONE-TWO*ENU))
      ELAMEG   = EMOD/(TWO*(ONE+ENU))
C     ============================================================================================ C
C     Initial preparations
C     ============================================================================================ C
      DO K1 = 1, NDOFEL                      
         DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
         END DO
         RHSK(K1) = ZERO
         RHSINI(K1) = ZERO
         DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
            AMASS (K2,K1) = ZERO
            ASTIFF(K2,K1) = ZERO
         END DO
      END DO
C     ============================================================================================ C
C     Local coordinates and weights
C     ============================================================================================ C
      XII(1,1) = -ONE/THREE**HALF
      XII(1,2) = -ONE/THREE**HALF
      XII(2,1) =  ONE/THREE**HALF
      XII(2,2) = -ONE/THREE**HALF
      XII(3,1) =  ONE/THREE**HALF
      XII(3,2) =  ONE/THREE**HALF
      XII(4,1) = -ONE/THREE**HALF
      XII(4,2) =  ONE/THREE**HALF
      INNODE   =  FOUR
      DO I=1,INNODE
         AINTW(I) = ONE
      END DO
C
C     ============================================================================================ C
C     Determining maximum phase value in the element
C     ============================================================================================ C
      PHELEMAX = ZERO
      DO K1=1,4
         PHMAX(K1)    = ZERO
         PNEWDTIP(K1) = TEN
      END DO
      DO K1=1,INNODE
         IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
              PHMAX(K1) = USRVAR(JELEM,NSTVTT+1,K1)
         ELSE
              PHMAX(K1) = USRVAR(JELEM,24,K1)
         ENDIF        
      END DO
      PHELEMAX = MAXVAL(PHMAX)
C
C     ============================================================================================ C
C     Calculating properties at each integration point
C     ============================================================================================ C
      DO INPT = 1,INNODE
C     -------------------------------------------------------------------------------------------- C
C     Variable initilization
      DO I=1,NSTVTT
         SDV(I) = SVARS(NSTVTT*(INPT-1)+I)
      END DO
C     -------------------------------------------------------------------------------------------- C
C     Local coordinates of the integration point
      XI(1) = XII(INPT,1)
      XI(2) = XII(INPT,2) 
C     Shape functions and local derivatives
      CALL SHAPEFUN(AN,dNdxi,XI)
C     Shape functions
      IY = ZERO
      DO I = 1,NNODE
         IX=IY+1
         IY=IX+1
         VNI(1,IX)=AN(I)
         VNI(1,IY)=ZERO
         VNI(2,IX)=ZERO
         VNI(2,IY)=AN(I)
      END DO
C     Jacobian
      DO I = 1,2
         DO J = 1,2
            VJACOB(I,J) = ZERO
            DO K = 1,NNODE
               VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
            END DO
         END DO
      END DO
C        
      DTM = ZERO
      DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
      IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT
      ENDIF
C     Inverse of Jacobian
      VJABOBINV(1,1) =  VJACOB(2,2)/DTM
      VJABOBINV(1,2) = -VJACOB(1,2)/DTM
      VJABOBINV(2,1) = -VJACOB(2,1)/DTM
      VJABOBINV(2,2) =  VJACOB(1,1)/DTM
C        
C     Derivatives of shape functions respect to global ccordinates
      DO K = 1,NNODE
         DO I = 1,2
           dNdx(K,I) = ZERO
           DO J = 1,2
              dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
           END DO
         END DO
      END DO
C
C     Calculating B matrix (B=LN)
      IY=0
      DO INODE=1,NNODE
          IX=IY+1
          IY=IX+1
          BB(1,IX)= dNdx(INODE,1)
          BB(1,IY)= ZERO
          BB(2,IX)= ZERO
          BB(2,IY)= dNdx(INODE,2)
          BB(3,IX)= dNdx(INODE,2)
          BB(3,IY)= dNdx(INODE,1)
      END DO
C
C     -------------------------------------------------------------------------------------------- C
C     Nodal displacements
C     -------------------------------------------------------------------------------------------- C
      DO J=1,2
         ULOC(J) = ZERO
      END DO
      DO J=1,2
         DO I=1,NDOFEL
            ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
         END DO
      END DO  
      DO J=1,2
         SDV(J) = ULOC(J)
      END DO
C     -------------------------------------------------------------------------------------------- C
C     Nodal phase-field
C     -------------------------------------------------------------------------------------------- C
      IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         PHASE = USRVAR(JELEM,NSTVTT+1,INPT)
      ELSE
         PHASE = USRVAR(JELEM,24,INPT)
      ENDIF
C
      IF (PHASE.GT.ONE) THEN
         PHASE = ONE
      ELSEIF (PHASE.LT.ZERO) THEN
         PHASE = ZERO
      ENDIF
C
      SDV(24) = PHASE
C     -------------------------------------------------------------------------------------------- C
C     Calculating strain
C     -------------------------------------------------------------------------------------------- C
      DO J=1,3
         EPS(J)  = ZERO
         DEPS(J) = ZERO
      END DO
      DO I=1,3
         DO J=1,NDOFEL
            EPS(I) = EPS(I)+BB(I,J)*U(J)    
         END DO
      END DO

      DO J=1,3
         DEPS(J) = EPS(J)-SDV(J+2)
      END DO
      DO J=1,3
         SDV(J+2) = EPS(J)
      END DO
C     -------------------------------------------------------------------------------------------- C
C     Recovering elastic and plastic strains from previous step
C     -------------------------------------------------------------------------------------------- C
      DO K1=1,4 
         EELAS(K1) = SDV(K1+5)            ! elastic strain
         EPLAS(K1) = SDV(K1+9)            ! plastic strain
         EPSZ(K1)  = SDV(K1+5)
      END DO
C
      EQPLAS = SDV(14)                    ! equivalent plastic strain
C
      EELAS(1) = EELAS(1) + DEPS(1)
      EELAS(2) = EELAS(2) + DEPS(2)
      EELAS(4) = EELAS(4) + DEPS(3)
C
C     -------------------------------------------------------------------------------------------- C
C     Calculating eigenvalue decomposition
C     -------------------------------------------------------------------------------------------- C
C     Only updating the stiffness matrix in the fist iteration step
      DO K1=1,6
         EPSC(K1) = ZERO
      END DO
C
      IF ((STEPITER.LE.(3-REITER)).AND.(REITER.LT.4)) THEN
C ------ Normal case: the stiffness is refreshed in every iteration
C        based on the actual strain state
         DO K1=1,4
            EPSC(K1)   = EELAS(K1)
            SDV(24+K1) = EPSC(K1)
         END DO
C     
      ELSEIF ((REITER.GE.4).AND.(STEPITER.EQ.ZERO)) THEN
C ------ If maximum trials are exceeded, the strain state from
C        the original (converged) step is taken
         DO K1=1,4
            EPSC(K1)=EPSZ(K1)
            SDV(24+K1)=EPSC(K1)
         END DO
      ELSE
C ------ If the stiffness is not refreshed, the strain is recovered
C        from an older NR iteration
         DO K1=1,4
            EPSC(K1)=USRVAR(JELEM,K1+24,INPT)
            SDV(24+K1)=EPSC(K1)
         END DO
      ENDIF
C       
      VALMDEPS = MAXVAL(ABS(DEPS))
C
      IF (VALMDEPS.GT.(0.1)) THEN
          USRVAR(JELEM,NALL+1,2)=ZERO
          PLSWTGLOB=ZERO
      ENDIF
C
C     -------------------------------------------------------------------------------------------- C
C     Calculating degrated elastic stiffness matrix
C     -------------------------------------------------------------------------------------------- C
      DO I=1,6
         DO J=1,6
            CMATG(I,J)=ZERO
         END DO
      END DO
C     
      IF (PHASE.LT.TOLER) THEN
          CMATG(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMATG(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMATG(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
          CMATG(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMATG(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMATG(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMATG(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMATG(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMATG(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
          CMATG(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
          CMATG(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
          CMATG(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
          DO I=1,6
              DO J=1,6
                  CMATG(I,J)=CMATG(I,J)*((ONE-PHASE)**TWO+PARK)
              END DO
          END DO
      ELSE
          CALL HMAT(EPSC,CMATG,ENU,EMOD,ANISOSWT,PLSWT,PHASE,PARK)
      ENDIF
C        
      DO I=1,4
         DO J=1,4
            CMATP(I,J)=CMATG(I,J)
         END DO
      END DO
C
C     -------------------------------------------------------------------------------------------- C
C
C     ------------------ Testing the yield criterion -------------------
C
C     -------------------------------------------------------------------------------------------- C
C     Expending the 2D stress to 3D (as in plane-strain sig_z~=0)
C     Strain vectors
      DO K1=1,4
         EPSP(K1)   = EELAS(K1)
         DEPLAS(K1) = ZERO
      END DO
C
C     Initial stress prediction
      DO K1=1,4
         SFULL(K1)=ZERO
      END DO
      DO I=1,4
         DO J=1,4
            SFULL(I) = SFULL(I)+CMATP(I,J)*EPSP(J)
         END DO 
      END DO     
C     Hardening yield strength
C
      YSHARD = (YIELDS+EQPLAS*VHMOD)*((ONE-PHASE)**TWO+PARK)
       
C     Stress invariants
C       
      SMISES = (SFULL(1)-SFULL(2))*(SFULL(1)-SFULL(2)) + 
     1         (SFULL(2)-SFULL(3))*(SFULL(2)-SFULL(3)) + 
     2         (SFULL(3)-SFULL(1))*(SFULL(3)-SFULL(1)) 
      SMISES = SMISES+SIX*SFULL(4)*SFULL(4)
      SMISES = SQRT(SMISES/TWO)
      FFUN   = SMISES-YSHARD
C
      DLAMB=ZERO
      IF ((PLSWT.GT.TOLER).AND.(PLSWTGLOB.GT.TOLER)) THEN
C ==========	Testing yield criterion =============
C
      IF ((FFUN.GT.TOLER).AND.(PHELEMAX.LT.PHCALCMAX)) THEN
          FLOW(1) = TWO*SFULL(1)-SFULL(2)-SFULL(3)
          FLOW(2) = TWO*SFULL(2)-SFULL(1)-SFULL(3)
          FLOW(3) = TWO*SFULL(3)-SFULL(2)-SFULL(1)
          FLOW(4) = SIX*SFULL(4)
          DO K1=1,4
             FLOW(K1) = FLOW(K1)*HALF/SMISES
          END DO
          DO K1=1,4
             FLOWG(K1) = FLOW(K1)
          END DO
C
C   Calculating the relationship between EQPLAS and DLAMB
          VNEVEZ=ZERO
          DO K1=1,4
            DO K2=1,4
               VNEVEZ=VNEVEZ+FLOW(K2)*CMATP(K2,K1)*FLOWG(K1)
            END DO
          END DO
          VNEVEZ=VNEVEZ+(VHMOD*((ONE-PHASE)**TWO+TOLER))
          DLAMB=FFUN/VNEVEZ
C
          DO K1=1,4
            DEPLAS(K1)=DLAMB*FLOWG(K1)
          END DO
          EQPLAS=EQPLAS+DLAMB
C
C   Updating stresses and strains
          DO K2=1,4
             DO K1=1,4
                SFULL(K2)=SFULL(K2)-CMATP(K2,K1)* DEPLAS(K1)
             END DO
             EELAS(K2)=EELAS(K2)-DEPLAS(K2)
             EPLAS(K2)=EPLAS(K2)+DEPLAS(K2)
          END DO
C       
C   Tangent matrix      
          DO K2=1,4
             DO K4=1,4
                DDSDDEEQ(K2,K4)= ZERO
             END DO
          END DO
C
          DO K1=1,4
             DO K2=1,4
                DO K3=1,4
                   DO K4=1,4
            DDSDDEEQ(K2,K4)= DDSDDEEQ(K2,K4)+CMATP(K2,K1)*FLOWG(K1)*
     1        FLOW(K3)*CMATP(K3, K4)
                   END DO
                END DO
             END DO
          END DO
C       
          SMISES = (SFULL(1)-SFULL(2))*(SFULL(1)-SFULL(2)) + 
     1             (SFULL(2)-SFULL(3))*(SFULL(2)-SFULL(3)) + 
     2             (SFULL(3)-SFULL(1))*(SFULL(3)-SFULL(1)) 
          SMISES = SMISES+SIX*SFULL(4)*SFULL(4)
          SMISES = SQRT(SMISES/TWO)        
C             
          DO K2=1,4
             DO K1=1,4
                CMATP(K2,K1)=CMATP(K2,K1)-DDSDDEEQ(K2,K1)/VNEVEZ
             END DO
          END DO
         
       ENDIF
       ENDIF
C
C     Converting every result back to 2D
C
      DO K1=1,4
         SDV(K1+5) = EELAS(K1)
         SDV(K1+9) = EPLAS(K1)
      END DO
      SDV(14) = EQPLAS
C
C     Stresses
      STRESS(1) = SFULL(1)
      STRESS(2) = SFULL(2)
      STRESS(3) = SFULL(4)
      DO J=1,4
         SDV(J+14) = SFULL(J)
      END DO
      HYDRO = (SFULL(1)+SFULL(2)+SFULL(3))/THREE
      SDV(19) = HYDRO
      SDV(20) = SMISES

C
C     Materials stiffness matrix
      CMAT(1,1)=CMATP(1,1)
      CMAT(1,2)=CMATP(1,2)
      CMAT(2,1)=CMATP(2,1)
      CMAT(2,2)=CMATP(2,2)
      CMAT(1,3)=CMATP(1,4)
      CMAT(3,1)=CMATP(4,1)
      CMAT(2,3)=CMATP(2,4)
      CMAT(3,2)=CMATP(4,2)
      CMAT(3,3)=CMATP(4,4)
C
C     -------------------------------------------------------------------------------------------- C
C     Calculating elastic ENERGY
C     -------------------------------------------------------------------------------------------- C
      ENG  = ZERO
      ENGP = ZERO
      ENGK = ZERO
C
C     Kinetic enery        
      ENGK = ZERO
C
C     Elastic strain enery        
      DO I=1,4
         EPSC(I)=EELAS(I)
      END DO
      EPSC(5)=ZERO
      EPSC(6)=ZERO
      CALL EIGOWN(EPSC,EIGV,ALPHA,ALPHAI,VECTI)
      IF (ANISOSWT.LT.TOLER) THEN
         DO K1=1,3
            ALPHAI(K1)=ONE
         END DO
         ALPHA=ONE
      ENDIF
      IF (PLSWT.GT.TOLER) THEN
          DO K1=1,3
             ALPHAI(K1)=ONE
          END DO
      ENDIF
C      
      ENGP=(ELAMEL*(ALPHA*(EIGV(1)+EIGV(2)+EIGV(3)))**TWO)/
     1  TWO+ELAMEG*((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*
     2  ALPHAI(2))**TWO+(EIGV(3)*ALPHAI(3))**TWO)
      ENGN=(ELAMEL*((ONE-ALPHA)*(EIGV(1)+EIGV(2)+EIGV(3)))**
     1  TWO)/TWO+ELAMEG*((EIGV(1)*(ONE-ALPHAI(1)))**TWO+(EIGV(2)*
     2  (ONE-ALPHAI(2)))**TWO+(EIGV(3)*(ONE-ALPHAI(3)))**TWO)
C
C     Plastic energy
C        
      PLPEN=ZERO
      IF (EQPLASCR.EQ.ZERO) THEN
          PLPEN=ZERO        
      ELSE
          ENGMAXEL = (YIELDS+VHMOD*(DEPSCR+ONE)*EQPLASCR)**TWO/
     1                SIX/ELAMEG
          ENGMAXPL = YIELDS*(DEPSCR+ONE)*EQPLASCR+VHMOD*
     1                ((DEPSCR+ONE)*EQPLASCR)**TWO/TWO
          ENGMAXCR = GCPAR/TWO/CLPAR
C         
          PLPEN = ENGMAXCR-ENGMAXEL-ENGMAXPL
          PLPEN = (PLPEN+ABS(PLPEN))/TWO
          PLPEN = TWO*PLPEN/(DEPSCR*EQPLASCR)**TWO
      ENDIF
C
      ENGPL=ZERO
      IF (EQPLAS.GT.EQPLASCR) THEN
         ENGPL=HALF*(EQPLAS-EQPLASCR)**TWO*PLPEN
      ELSE
         ENGPL=ZERO
      ENDIF
      ENGPL = ENGPL+EQPLAS*(YIELDS+HALF*EQPLAS*VHMOD)
C        
C
C     ---------   Energy increment and time step control -------------
C
      DENG=ENGPL+ENGP-SDV(21)-SDV(22)
C        
      DENGMAXV=GCPAR/TWO/CLPAR*DENGMAX
      IF ((LFLAGS(1).EQ.1).OR.(LFLAGS(1).EQ.11)) THEN
C --------- If automatic time integration is used
          IF ((DENG.GT.DENGMAXV).AND.(REITER.LT.4)) THEN
C --------- If energy criterion is violated (dE>dE_max)
              IF (DTIME.GT.DTMIN*(ONE+TOLER)) THEN
C --------- If time-step still can be reduced
                  PNEWDTIP(INPT)=DENGMAXV/DENG/TEN
                  DTIMENEXT=PNEWDTIP(INPT)*DTIME
                  IF (DTIMENEXT.LT.(DTMIN*(ONE+TOLER))) THEN
                      PNEWDTIP(INPT)=PNEWDTIP(INPT)*DTMIN/DTIMENEXT
                  ENDIF
              ELSE
C --------- If time-step is already too small
              PNEWDTIP(INPT)=ONE 
              ENDIF
          ENDIF
      ENDIF
C
      IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
      ELSE
         SDV(21)=ENGPL
         SDV(22)=ENGP
      ENDIF
      ENGKG=ENGKG+ENGK*DTM*THCK
      ENGEG=ENGEG+(ENGP*((ONE-PHASE)**TWO+PARK)+ENGN)*DTM*THCK
      ENGPG=ENGPG+ENGPL*((ONE-PHASE)**TWO+PARK)*DTM*THCK
      SDV(23)=((ENGPL+ENGP)*((ONE-PHASE)**TWO+PARK)+ENGN)
        
        
      

      
C
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
      DO K=1,NDOFEL
         DO L=1,NDOFEL
            DO I=1,3
               DO J=1,3
            ASTIFF(K,L)=ASTIFF(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1       BB(J,L)*DTM*THCK
               END DO
            END DO
         END DO
      END DO
C
C     ==================================================================
C     Calculating element mass matrix
C     ==================================================================
      DO K=1,NDOFEL
         DO L=1,NDOFEL
            DO I=1,2
            AMASS(K,K)=AMASS(K,K)+AINTW(INPT)*VNI(I,K)*VNI(I,L)
     1       *DTM*THCK*DENS
            END DO
         END DO
      END DO
C       
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
      IF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
         DO K1=1,NDOFEL
           DO K2=1,2
              RHSK(K1) = RHSK(K1)-AINTW(INPT)*VNI(K2,K1)*A(K1)*
     1                   DTM*THCK*DENS
          END DO
         END DO
        ENDIF
C
        DO K1=1,NDOFEL
         DO K4=1,3
           RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1      THCK
         END DO
        END DO
C       
C     -------------------------------------------------------------------------------------------- C
C     Uploading solution dep. variables
C     -------------------------------------------------------------------------------------------- C
        DO I=1,NSTVTT
         SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
         IF (LFLAGS(3).EQ.5) THEN
         ELSE
          USRVAR(JELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
         ENDIF
        END DO
       END DO
C
C     -------------------------------------------------------------------------------------------- C
C     Jacobien for the element
C     -------------------------------------------------------------------------------------------- C
      IF ((LFLAGS(1).EQ.1).OR.(LFLAGS(1).EQ.2)) THEN
C     Simple static calculation
         DO K=1,NDOFEL
          DO L=1,NDOFEL
             AMATRX(K,L)=AMATRX(K,L)+ASTIFF(K,L)
          END DO
         END DO
C
         ELSEIF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
C         Dynamic calculation
C
          PARALPHA=PARAMS(1)
          PARBETA=PARAMS(2)
          DO I=1,NDOFEL
           RHSINI(I)=SVARS(I+NSTVTT*INNODE)
           SVARS(I+NSTVTT*INNODE)=RHS(I,1)
          END DO
C
          IF (LFLAGS(3).EQ.4) THEN
           DO I=1,NDOFEL
            RHS(I,1)=ZERO
           END DO
          ELSE
           DO I=1,NDOFEL
            RHS(I,1)=RHS(I,1)*(ONE+PARALPHA)-RHSINI(I)*PARALPHA+
     1        RHSK(I)
           END DO
          ENDIF
C          
          IF ((LFLAGS(3).EQ.4).OR.(LFLAGS(3).EQ.6)) THEN
          DO K=1,NDOFEL
            DO L=1,NDOFEL
              AMATRX(K,L)=AMATRX(K,L)+AMASS(K,L)
            END DO
           END DO         
C           
          ELSEIF (LFLAGS(3).EQ.1) THEN
           DADU=ONE/(PARBETA*DTIME**TWO)
           DO K=1,NDOFEL
            DO L=1,NDOFEL
              AMATRX(K,L)=AMATRX(K,L)+AMASS(K,L)*DADU+(ONE+PARALPHA)*
     1         ASTIFF(K,L)
            END DO
           END DO 
          ENDIF
        ENDIF       
C     New time increment
       PNEWDTE=MINVAL(PNEWDTIP)
       IF (PNEWDTE.LT.(ONE+TOLER)) THEN
        PNEWDT=PNEWDTE
      ENDIF
C
C     ******************************************************************************************** C
C     ******************************************************************
C     Constructing elemet TYPE 4 (damage phase-field)
C     ******************************************************************
C     ******************************************************************************************** C
      ELSEIF ((JTYPE.EQ.FOUR)) THEN
C     ============================================================================================ C
      REITER   = USRVAR(JELEM-N_ELEM,NALL+3,1)
      STEPITER = USRVAR(JELEM-N_ELEM,NALL+4,1)
C     ============================================================================================ C
C     Material parameters
C     ============================================================================================ C
      CLPAR = PROPS(1)
      GCPAR = PROPS(2)
      THCK  = PROPS(3)
      ELSWT = PROPS(4)
C     ============================================================================================ C
C     Initial preparations
C     ============================================================================================ C
      DO K1 = 1, NDOFEL                      
         DO KRHS = 1, NRHS
            RHS(K1,KRHS) = ZERO
         END DO
         DO K2 = 1, NDOFEL
            AMATRX(K2,K1) = ZERO
         END DO
      END DO
C     ============================================================================================ C
C     Local coordinates and weights
C     ============================================================================================ C
      XII(1,1) = -ONE/THREE**HALF
      XII(1,2) = -ONE/THREE**HALF
      XII(2,1) =  ONE/THREE**HALF
      XII(2,2) = -ONE/THREE**HALF
      XII(3,1) =  ONE/THREE**HALF
      XII(3,2) =  ONE/THREE**HALF
      XII(4,1) = -ONE/THREE**HALF
      XII(4,2) =  ONE/THREE**HALF
      INNODE = FOUR
      DO I=1,INNODE
         AINTW(I) = ONE
      END DO
C     ============================================================================================ C
C     Calculating properties at each integration point
C     ============================================================================================ C
      DO INPT=1,INNODE
C     -------------------------------------------------------------------------------------------- C
C     -------------------------------------------------------------------------------------------- C
C     Initializing solution dependent variables (phase,history)
      DO I=1,NSTVTO
         SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
      END DO
C
C     Local coordinates of the integration point
      XI(1) = XII(INPT,1)
      XI(2) = XII(INPT,2) 
C     Shape functions and local derivatives
      CALL SHAPEFUN(AN,dNdxi,XI)
C     Jacobian
      DO I = 1,2
         DO J = 1,2
            VJACOB(I,J) = ZERO
            DO K = 1,NNODE
               VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
            END DO
         END DO
      END DO
C        
      DTM = ZERO
      DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
      IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
      ENDIF
C     Inverse of Jacobian
      VJABOBINV(1,1) =  VJACOB(2,2)/DTM
      VJABOBINV(1,2) = -VJACOB(1,2)/DTM
      VJABOBINV(2,1) = -VJACOB(2,1)/DTM
      VJABOBINV(2,2) =  VJACOB(1,1)/DTM
C        
C     Derivatives of shape functions respect to global ccordinates
      DO K = 1,NNODE
         DO I = 1,2
            dNdx(K,I) = ZERO
            DO J = 1,2
               dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
            END DO
         END DO
      END DO
C
C     Calculating B matrix (B=LN)
      DO INODE=1,NNODE
         BP(1,INODE)=dNdx(INODE,1)
         BP(2,INODE)=dNdx(INODE,2)
      END DO
C
C     -------------------------------------------------------------------------------------------- C
C     Nodal phase-field
C     -------------------------------------------------------------------------------------------- C
      PHASE  = ZERO
      DPHASE = ZERO
      DO I=1,NDOFEL
         PHASE=PHASE+AN(I)*U(I)
      END DO
      DO I=1,NDOFEL
         DPHASE=DPHASE+AN(I)*DU(I,1)
      END DO
      SDV(1)=PHASE
C
C     Gradient
      DO I=1,2
         DP(I)=ZERO
      END DO
      DO I=1,2
         DO J=1,NNODE
            DP(I)=DP(I)+BP(I,J)*U(J)
         END DO
      END DO
C
C     -------------------------------------------------------------------------------------------- C
C     Calculating elastic ENERGY history
C     -------------------------------------------------------------------------------------------- C
      IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
C          ENGN=USRVAR(JELEM-N_ELEM,21,INPT)+USRVAR(JELEM-N_ELEM,22,INPT)                 ! 21--plastic energy, 22-elastic energy
C          IF (ELSWT.GT.ZERO) THEN
C              ENGN=ENGN-GCPAR/TWO/CLPAR
C          ENDIF
          
          ENGN_ELA_THRES = 0.00144D0
          ENGN_PLA_THRES = 0.08D0
          
          ENGN_ELA = max(USRVAR(JELEM-N_ELEM,22,INPT) - ENGN_ELA_THRES, 0.0d0)
          ENGN_PLA = max(USRVAR(JELEM-N_ELEM,21,INPT) - ENGN_PLA_THRES, 0.0d0)
          
          ENGN = ENGN_ELA + ENGN_PLA
           
      ELSE
          ENGN=USRVAR(JELEM-N_ELEM,NSTVTT+2,INPT)
      ENDIF
C       
      HISTN = USRVAR(JELEM-N_ELEM,NSTVTT+2,INPT)
      IF (ENGN.GT.HISTN) THEN
          HIST=ENGN
      ELSE
          HIST=HISTN
      ENDIF
      SDV(2)=HIST
      
      



C     -------------------------------------------------------------------------------------------- C
C     Calculating fracture energy for history output
C     -------------------------------------------------------------------------------------------- C
C
      ENGD=ZERO
      ENGD=PHASE**TWO/TWO/CLPAR*DTM*THCK*GCPAR
C
      DO J=1,2
          ENGD=ENGD+DP(J)*DP(J)*CLPAR*DTM*THCK/TWO*GCPAR
      END DO
      ENGDG=ENGDG+ENGD
C
C     -------------------------------------------------------------------------------------------- C
C     Calculating element stiffness matrix
C     -------------------------------------------------------------------------------------------- C
      DO I=1,NNODE
         DO K=1,NNODE
            DO J=1,2
               AMATRX(I,K) = AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1                       THCK*GCPAR*CLPAR*AINTW(INPT)
            END DO
            AMATRX(I,K) = AMATRX(I,K)+AN(I)*AN(K)*DTM*THCK*
     1                    AINTW(INPT)*(GCPAR/CLPAR+TWO*HIST)
         END DO
      END DO
C        
C     -------------------------------------------------------------------------------------------- C
C     Internal forces (residual vector)
C     -------------------------------------------------------------------------------------------- C
      DO I=1,NDOFEL
         DO J=1,2
            RHS(I,1) = RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR*
     1                 AINTW(INPT)*DTM*THCK
         END DO
         RHS(I,1) = RHS(I,1)-AN(I)*AINTW(INPT)*DTM*THCK*
     1              ((GCPAR/CLPAR+TWO*HIST)*PHASE-TWO*HIST)
      END DO
C
C     -------------------------------------------------------------------------------------------- C
C     Uploading solution dep. variables
C     -------------------------------------------------------------------------------------------- C
      DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
         IF (LFLAGS(3).EQ.5) THEN
         ELSE
             USRVAR(JELEM-N_ELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
         ENDIF
      END DO
C     -------------------------------------------------------------------------------------------- C
C     -------------------------------------------------------------------------------------------- C
      END DO
C     -------------------------------------------------------------------------------------------- C
      DO I=1,NDOFEL
         RHSINI(I)=SVARS(I+NSTVTO*INNODE)
         SVARS(I+NSTVTO*INNODE)=RHS(I,1)
      END DO
      IF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
          PARALPHA=PARAMS(1)
          DO I=1,NDOFEL
             RHS(I,1)=RHS(I,1)*(ONE+PARALPHA)+RHSINI(I)*PARALPHA
          END DO
          DO I=1,NNODE
             DO K=1,NNODE
                AMATRX(I,K)=AMATRX(I,K)*(ONE+PARALPHA)
             END DO
          END DO
      ENDIF
C     ******************************************************************************************** C
C     ******************************************************************************************** C
      ENDIF
C     ******************************************************************************************** C

      ENERGY(1)=ENGKG
      ENERGY(2)=ENGEG
      ENERGY(4)=ENGPG
      ENERGY(7)=ENGDG
      RETURN
      END

C     ============================================================================================ C
C     ============================================================================================ C
C     Shape functions for square elements
C     ============================================================================================ C
C     ============================================================================================ C
      SUBROUTINE SHAPEFUN(AN,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
C
C     Values of shape functions as a function of local coord.
      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
C
C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END
C
C     ============================================================================================ C
C     ============================================================================================ C
C     Eigenstrains from Voigt notation
C     ============================================================================================ C
C     ============================================================================================ C

      SUBROUTINE EIGOWN(EPS,EIGV,ALPHA,ALPHAI,VECTI)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,
     1  TS=27.D0,THREE=3.D0,HALF=0.5D0,TOLER=1.0D-12,FOUR=4.D0,
     2  CNTN=100,TOLERE=1.0D-12)
       INTEGER I, J, K
       REAL*8 EPS(6), EIGV(3), ALPHAI(3), VECTI(3)
       REAL*8 PC, QC, ALPHA, DISC, PI, CNT
C
       PI=FOUR*ATAN(ONE)
C Scaling the strain vector
       VMAXE=MAXVAL(ABS(EPS))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)/VMAXE
        END DO
       ENDIF
C    
C   Calculating eigenvalues
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
C
C   Depressed coefficients    
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
C
       DO I=1,3
        EIGV(I)=ZERO
       END DO
       CNT=ZERO
       IF (ABS(DISC).LT.TOLER) THEN
        IF ((ABS(QC).LT.TOLER).AND.(ABS(PC).LT.TOLER)) THEN
         EIGV(1)=VECTI(1)/THREE
         EIGV(2)=VECTI(1)/THREE
         EIGV(3)=VECTI(1)/THREE
        ELSE
         EIGV(1)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(2)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(3)=THREE*QC/PC+VECTI(1)/THREE
         IF (EIGV(1).GT.EIGV(3)) THEN
          EONE=EIGV(1)
          EIGV(1)=EIGV(3)
          EIGV(3)=EONE
         ENDIF
        ENDIF
       ELSE
        DO I=1,3
         EIGV(I)=VECTI(1)/THREE+TWO*(MONE*PC/THREE)**HALF*
     1   COS(ONE/THREE*ACOS(MONE*QC/TWO*(TS/(MONE*PC**THREE))**
     2   HALF)+TWO*I*PI/THREE)
        END DO
       ENDIF
C       
       ALPHA=ZERO
       IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
        ALPHA=ONE
       ENDIF
       DO K1=1,3
        ALPHAI(K1)=ZERO
        IF (EIGV(K1).GT.TOLER) THEN
         ALPHAI(K1)=ONE
        ENDIF
       END DO
C
C    Rescaling eigenvalues       
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)*VMAXE
        END DO
        DO K1=1,3
         EIGV(K1)=EIGV(K1)*VMAXE
        END DO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
       ENDIF
C    
       RETURN
      END
C     ============================================================================================ C
C     ============================================================================================ C
C
C     ============================================================================================ C
C     ============================================================================================ C
       SUBROUTINE HMAT(EPSI,CMAT,ENU,EMOD,ANISOSWT,PLSWT,PHASE,PARK)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.D0,
     1  TOLER=1.0D-12,SIX=6.D0,FT=50.D0,THREE=3.D0,HALF=0.5D0,
     2  TS=27.D0,CNTM=1000,TEN=10.D0,TOLERE=1.0D-12,TOLD=1.0D-7)
       INTEGER I, J, K, NDIM
       REAL*8 EPS(6),VECTI(3),CLMAT(3,3),EIGV(3),ALPHAI(3),
     1  DLDI(3,3),DLDIDI(3,3,3),DIDE(3,6),DIDEDE(3,6,6),
     2  DLDEDE(3,6,6),CMAT(6,6),DLDE(3,6),EPSI(6),
     3  EIGVEC(3,3),EIGVAL(3,3),EPSNEW(3,3)
       REAL*8 VMAXE, ALPHA, ANISOSWT, DENOM, ENU, EMOD,
     1 PHASE, PC, QC, CNT, DISC, PLSWT
C
C Rescaling the strain vector
       VMAXE=MAXVAL(ABS(EPSI))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPSI(K1)/VMAXE
        END DO
       ELSE
        DO K1=1,6
         EPS(K1)=EPSI(K1)
        END DO
       ENDIF
C       
C   Calculating invariants
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
C
C   Calculating eigenvalues
       CALL EIGOWN(EPS,EIGV,ALPHA,ALPHAI,VECTI)
       IF (PLSWT.GT.TOLER) THEN
        DO K1=1,3
         ALPHAI(K1)=ONE
        END DO
       ENDIF
C
C	Initialising CMAT
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=ZERO
         END DO
        END DO
C
C ************ Starting spectral decomposition ************************        
       IF ((MINVAL(EIGV).GE.-TOLER).OR.(ANISOSWT.LT.TOLER)) THEN
C All eigenvalues are in tension, therefore the stiffness matrix can be
C     degraded without the decomposition
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=CMAT(I,J)*((ONE-PHASE)**TWO+PARK)
         END DO
        END DO
C
       ELSEIF ((MAXVAL(EIGV).LT.-TOLER).AND.(PLSWT.LT.HALF)) THEN
C All eigenvalues are in compression, therefore the stiffness
C     matrix does not need to be degraded 
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
C
       ELSE
C   Calculating materials stiffness matrix
        DO I=1,3
         DO J=1,3
          CLMAT(I,J)=ZERO
         END DO
        END DO
        GAMMA=ENU/(ONE-TWO*ENU)
        CLMAT(1,1)=((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,2)=((ONE-ALPHAI(2)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(3,3)=((ONE-ALPHAI(3)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,2)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,3)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,3)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,1)=CLMAT(1,2)
        CLMAT(3,1)=CLMAT(1,3)
        CLMAT(3,2)=CLMAT(2,3)
        DO I=1,3
         DO J=1,3
          CLMAT(I,J)=CLMAT(I,J)*EMOD/(ONE+ENU)
         END DO
        END DO
C
C Adding small permutation if two eigenvalues are close or equal 
        IF (ABS(DISC).LT.TOLD) THEN
         CALL JACOBYEIG(EPS,EIGVEC,EIGVAL)
         IF (ABS(EIGV(2)).LT.TOLER) THEN
          EIGV(2)=-0.1
         ELSE
          EIGV(2)=EIGV(2)*1.1D0 
         ENDIF
C
C Reconstructing strain matrix from the modified eigenstrains         
         DO I=1,3
          DO J=1,3
           EPSNEW(I,J)=ZERO
           DO K=1,3
            EPSNEW(I,J)=EPSNEW(I,J)+EIGVEC(I,K)*EIGV(K)*EIGVEC(J,K)
           END DO
          END DO
         END DO
         EPS(1)=EPSNEW(1,1)
         EPS(2)=EPSNEW(2,2)
         EPS(3)=EPSNEW(3,3)
         EPS(4)=EPSNEW(1,2)*TWO
         EPS(5)=EPSNEW(1,3)*TWO
         EPS(6)=EPSNEW(2,3)*TWO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
         PC=VECTI(2)-VECTI(1)**TWO/THREE
         QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
         DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)       
        ENDIF
C        
C   Calculating derivatives of eigenvalues respect to invariants
        DO I=1,3
         DO J=1,3
          DLDI(I,J)=ZERO
         END DO
        END DO
        DO K=1,3
         DENOM=THREE*EIGV(K)**TWO-VECTI(1)*TWO*EIGV(K)+
     1   VECTI(2)
         IF (ABS(DENOM).LT.TOLER) THEN
          WRITE(7,*) 'EPS: ',EPS
          WRITE(7,*) 'EIGV: ',EIGV
          WRITE(7,*) 'VECTI: ',VECTI
          WRITE(7,*) 'DENOM: ',DENOM
          WRITE(7,*) 'Denominator is close to 0.'
C          CALL XIT
         ENDIF
         DLDI(K,1)=EIGV(K)**TWO/DENOM
         DLDI(K,2)=MONE*EIGV(K)/DENOM
         DLDI(K,3)=ONE/DENOM
        END DO 
C
        DO I=1,3
         DO J=1,3
          DO K=1,3
           DLDIDI(K,I,J)=ZERO
          END DO
         END DO
        END DO
        DO K=1,3
         DENOM=THREE*EIGV(K)**TWO-VECTI(1)*TWO*EIGV(K)+
     1   VECTI(2)
         DLDIDI(K,1,1)=TWO*EIGV(K)*DLDI(K,1)/DENOM-EIGV(K)**
     1   TWO*(SIX*EIGV(K)*DLDI(K,1)-TWO*VECTI(1)*DLDI(K,1)-
     2   TWO*EIGV(K))/DENOM**TWO
         DLDIDI(K,1,2)=TWO*EIGV(K)*DLDI(K,2)/DENOM-EIGV(K)**
     1   TWO*(SIX*EIGV(K)*DLDI(K,2)-TWO*VECTI(1)*DLDI(K,2)+
     2   ONE)/DENOM**TWO
         DLDIDI(K,1,3)=TWO*EIGV(K)*DLDI(K,3)/DENOM-EIGV(K)**
     1   TWO*(SIX*EIGV(K)*DLDI(K,3)-TWO*VECTI(1)*DLDI(K,3))/
     2   DENOM**TWO
         DLDIDI(K,2,1)=EIGV(K)*(SIX*EIGV(K)*DLDI(K,1)-TWO*
     1   VECTI(1)*DLDI(K,1)-TWO*EIGV(K))/DENOM**TWO-
     2   DLDI(K,1)/DENOM
         DLDIDI(K,2,2)=EIGV(K)*(SIX*EIGV(K)*DLDI(K,2)-TWO*
     1   VECTI(1)*DLDI(K,2)+ONE)/DENOM**TWO-DLDI(K,2)/DENOM
         DLDIDI(K,2,3)=EIGV(K)*(SIX*EIGV(K)*DLDI(K,3)-TWO*
     1   VECTI(1)*DLDI(K,3))/DENOM**TWO-DLDI(K,3)/DENOM
         DLDIDI(K,3,1)=MONE*(SIX*EIGV(K)*DLDI(K,1)-TWO*
     1   VECTI(1)*DLDI(K,1)-TWO*EIGV(K))/DENOM**TWO
         DLDIDI(K,3,2)=MONE*(SIX*EIGV(K)*DLDI(K,2)-TWO*
     1   VECTI(1)*DLDI(K,2)+ONE)/DENOM**TWO
         DLDIDI(K,3,3)=MONE*(SIX*EIGV(K)*DLDI(K,3)-TWO*
     1   VECTI(1)*DLDI(K,3))/DENOM**TWO
        END DO
C
C   Calculating derivatives of invariants respect to matrix components
        DO I=1,3
         DO J=1,6
          DIDE(I,J)=ZERO
         END DO
        END DO
        DIDE(1,1)=ONE
        DIDE(1,2)=ONE
        DIDE(1,3)=ONE
        DIDE(2,1)=EPS(2)+EPS(3)
        DIDE(2,2)=EPS(1)+EPS(3)
        DIDE(2,3)=EPS(2)+EPS(1)
        DIDE(2,4)=MONE*EPS(4)/TWO
        DIDE(2,5)=MONE*EPS(5)/TWO
        DIDE(2,6)=MONE*EPS(6)/TWO
        DIDE(3,1)=EPS(2)*EPS(3)-EPS(6)**TWO/FOUR
        DIDE(3,2)=EPS(1)*EPS(3)-EPS(5)**TWO/FOUR
        DIDE(3,3)=EPS(2)*EPS(1)-EPS(4)**TWO/FOUR
        DIDE(3,4)=EPS(5)*EPS(6)/FOUR-EPS(4)*EPS(3)/TWO
        DIDE(3,5)=EPS(4)*EPS(6)/FOUR-EPS(5)*EPS(2)/TWO
        DIDE(3,6)=EPS(5)*EPS(4)/FOUR-EPS(6)*EPS(1)/TWO
C
        DO I=1,6
         DO J=1,6
          DO K=1,3
           DIDEDE(K,I,J)=ZERO
          END DO
         END DO
        END DO
        DIDEDE(2,4,4)=MONE*HALF
        DIDEDE(2,5,5)=MONE*HALF
        DIDEDE(2,6,6)=MONE*HALF
        DIDEDE(2,1,2)=ONE
        DIDEDE(2,2,1)=ONE
        DIDEDE(2,1,3)=ONE
        DIDEDE(2,3,1)=ONE
        DIDEDE(2,2,3)=ONE
        DIDEDE(2,3,2)=ONE
        DIDEDE(3,1,2)=EPS(3)
        DIDEDE(3,1,3)=EPS(2)
        DIDEDE(3,1,6)=MONE*EPS(6)/TWO
        DIDEDE(3,2,3)=EPS(1)
        DIDEDE(3,2,5)=MONE*EPS(5)/TWO
        DIDEDE(3,3,4)=MONE*EPS(4)/TWO
        DIDEDE(3,4,4)=MONE*EPS(3)/TWO
        DIDEDE(3,4,5)=EPS(6)/FOUR
        DIDEDE(3,4,6)=EPS(5)/FOUR
        DIDEDE(3,5,5)=MONE*EPS(2)/TWO
        DIDEDE(3,5,6)=EPS(4)/FOUR
        DIDEDE(3,6,6)=MONE*EPS(1)/TWO
        DO I=1,6
         DO J=1,I-1
          DIDEDE(3,I,J)=DIDEDE(3,J,I)
         END DO
        END DO
C
C    Calculating derivatives of eigenvalues respect to matrix components
        DO I=1,3
         DO J=1,6
          DLDE(I,J)=ZERO
         END DO
        END DO
        DO K=1,3
         DO J=1,6
          DO I=1,3
           DLDE(K,J)=DLDE(K,J)+DLDI(K,I)*DIDE(I,J)
          END DO
         END DO
        END DO
C
        DO I=1,6
         DO J=1,6
          DO K=1,3
           DLDEDE(K,I,J)=ZERO
          END DO
         END DO
        END DO
        DO N=1,3
         DO I=1,6
          DO J=1,6
           DO M=1,3
            DLDEDE(N,I,J)=DLDEDE(N,I,J)+DLDI(N,M)*DIDEDE(M,I,J)
             DO P=1,3
              DLDEDE(N,I,J)=DLDEDE(N,I,J)+DLDIDI(N,M,P)*
     1         DIDE(P,J)*DIDE(M,I)
             END DO
            END DO
           END DO
          END DO
         END DO
C
C   Stiffness matrix
        DO I=1,6
         DO J=1,6
          DO N=1,3
           DO M=1,3
            CMAT(I,J)=CMAT(I,J)+DLDE(M,J)*CLMAT(N,M)*DLDE(N,I)
     1       +EIGV(M)*CLMAT(M,N)*DLDEDE(N,I,J)
           END DO
          END DO
         END DO
        END DO
       ENDIF
       RETURN
       END
C     ============================================================================================ C
C     ============================================================================================ C
C
C     ============================================================================================ C
C     ============================================================================================ C
      SUBROUTINE JACOBYEIG(EPS,EIGVEC,EIGVAL)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.D0,
     1  TOLER=1.0D-12,SIX=6.D0,FT=50.D0,THREE=3.D0,HALF=0.5D0,
     2  TS=27.D0,CNTM=1000,TEN=10.D0,TOLERE=1.0D-12)
      INTEGER I, J, K, N, CNT
      REAL*8 EPS(6),EIGVEC(3,3),EIGVAL(3,3)
      REAL*8 B2, BAR, BETA, COEFF, S, C, CS, SC, VAG
C
      DO I=1,3
         DO J=1,3
            EIGVEC(I,J) = ZERO
            EIGVAL(I,J) = ZERO
         END DO
      END DO
C       
      DO I=1,3
         EIGVEC(I,I) = ONE
      END DO
C       
      EIGVAL(1,1) = EPS(1)
      EIGVAL(2,2) = EPS(2)
      EIGVAL(3,3) = EPS(3)
      EIGVAL(1,2) = EPS(4)/TWO
      EIGVAL(1,3) = EPS(5)/TWO
      EIGVAL(2,3) = EPS(6)/TWO
      EIGVAL(2,1) = EIGVAL(1,2)
      EIGVAL(3,1) = EIGVAL(1,3)
      EIGVAL(3,2) = EIGVAL(2,3)
C     
      B2 = ZERO
      DO I=1,3
         DO J=1,3
            IF (I.NE.J) THEN
               B2 = B2 + EIGVAL(I,J)**TWO
            ENDIF
         END DO
      END DO
      BAR = B2/(THREE*THREE)/TWO
C      
      CNT=ONE
      DO WHILE ((B2.GT.TOLER).AND.(CNT.LT.CNTM))
         DO I=1,2
            DO J=I+1,3
               IF (EIGVAL(J,I)**TWO.LE.BAR) THEN
               ELSE
                  B2 = B2 - TWO*EIGVAL(J,I)**TWO
                  BAR = HALF*B2/(THREE*THREE)
                  BETA = (EIGVAL(J,J)-EIGVAL(I,I))/(TWO*EIGVAL(J,I))
                  COEFF = HALF*BETA/SQRT(ONE+BETA**TWO)
                  S = SQRT(MAX(HALF+COEFF,ZERO))
                  C = SQRT(MAX(HALF-COEFF,ZERO))
                  DO K=1,3
                     CS =  C*EIGVAL(I,K)+S*EIGVAL(J,K)
                     SC = -S*EIGVAL(I,K)+C*EIGVAL(J,K)
                     EIGVAL(I,K) = CS
                     EIGVAL(J,K) = SC
                  END DO
                  DO K=1,3
                     CS =  C*EIGVAL(K,I)+S*EIGVAL(K,J)
                     SC = -S*EIGVAL(K,I)+C*EIGVAL(K,J)
                     EIGVAL(K,I) = CS
                     EIGVAL(K,J) = SC
                     CS =  C*EIGVEC(K,I)+S*EIGVEC(K,J)
                     SC = -S*EIGVEC(K,I)+C*EIGVEC(K,J)
                     EIGVEC(K,I) = CS
                     EIGVEC(K,J) = SC
                  END DO
              ENDIF
          END DO
       END DO
       CNT=CNT+ONE
       END DO
C
C ---------- Sorting eigenvalues and vectors ---------------
        VAG=ZERO
        IF (EIGVAL(1,1).GT.EIGVAL(2,2)) THEN
            VAG=EIGVAL(1,1)
            EIGVAL(1,1)=EIGVAL(2,2)
            EIGVAL(2,2)=VAG
            DO I=1,3
                VAG=EIGVEC(I,1)
                EIGVEC(I,1)=EIGVEC(I,2)
                EIGVEC(I,2)=VAG
            END DO
        ENDIF
        IF (EIGVAL(1,1).GT.EIGVAL(3,3)) THEN
            VAG=EIGVAL(1,1)
            EIGVAL(1,1)=EIGVAL(3,3)
            EIGVAL(3,3)=VAG
            DO I=1,3
                VAG=EIGVEC(I,1)
                EIGVEC(I,1)=EIGVEC(I,3)
                EIGVEC(I,3)=VAG
            END DO
        ENDIF 
        IF (EIGVAL(2,2).GT.EIGVAL(3,3)) THEN
            VAG=EIGVAL(2,2)
            EIGVAL(2,2)=EIGVAL(3,3)
            EIGVAL(3,3)=VAG
            DO I=1,3
                VAG=EIGVEC(I,2)
                EIGVEC(I,2)=EIGVEC(I,3)
                EIGVEC(I,3)=VAG
            END DO
        ENDIF        
       RETURN
      END            
C
C     ============================================================================================ C
C     ============================================================================================ C
C     UMAT for visulization
C     NOTE: N_ELEM has to be changed according to the UEL !!!!!
C     ============================================================================================ C
C     ============================================================================================ C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
      PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=100000,NSTV=34) 
      DATA NEWTON,TOLER/40,1.D-6/ 
C       
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,4)
C 
C     -------------------------------------------------------------------------------------------- C
C     PROPS(1) - Young's modulus 
C     PROPS(2) - Poisson ratio 
C     -------------------------------------------------------------------------------------------- C
C
C	Elastic properties
C
      EMOD = PROPS(1)
      ENU  = PROPS(2)
      EG   = EMOD/(TWO*(ONE+ENU))
      EG2  = EG*TWO
      ELAM = EG2*ENU/(ONE-TWO*ENU)
C
C	Stiffness tensor
C
      DO K1=1, NTENS
         DO K2=1, NTENS
            DDSDDE(K2, K1)=0.0
         END DO
      END DO
C
      DO K1=1, NDI
         DO K2=1, NDI
            DDSDDE(K2, K1)=ELAM
         END DO
         DDSDDE(K1, K1)=EG2+ELAM
      END DO 
C
      DO K1=NDI+1, NTENS
         DDSDDE(K1, K1)=EG
      END DO
C
C	Calculate Stresses
C
      DO K1=1, NTENS
         DO K2=1, NTENS
            STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
         END DO
      END DO 
C
      NELEMAN=NOEL-TWO*N_ELEM
      
      IF (NPT.EQ.3) THEN
         NPT=4
      ELSEIF (NPT.EQ.4) THEN
         NPT=3
      ENDIF
C       
      DO I=1,NSTATV
         STATEV(I)=USRVAR(NELEMAN,I,NPT)
      END DO
C       
      RETURN
      END      
      
