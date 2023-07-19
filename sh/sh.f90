!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   cir :    Homoclinic Bifurcation in an Electronic Circuit
!                (the same equations as in demo tor)
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION R, D, E

       R=PAR(1)
       D=PAR(2)
       E=PAR(3)

       F(1)= U(2)
       F(2)= U(3)
       F(3)= U(4)
       F(4)= (R-1)*U(1)-2*U(3)+D*U(1)**3-E*U(1)**5

       IF(IJAC.EQ.0)RETURN

       DFDU(1,1)=0
       DFDU(1,2)=1.0
       DFDU(1,3)=0
       DFDU(1,4)=0

       DFDU(2,1)=0
       DFDU(2,2)=0
       DFDU(2,3)=1.0
       DFDU(2,4)=0

       DFDU(3,1)=0
       DFDU(3,2)=0
       DFDU(3,3)=0
       DFDU(3,4)=1.0

       DFDU(4,1)=(R-1)+3*D*U(1)**2-5*E*U(1)**4
       DFDU(4,2)=0
       DFDU(4,3)=-2.0
       DFDU(4,4)=0

      IF(IJAC.EQ.1)RETURN

!      *Parameter derivatives
       DFDP(1,1)=0
       DFDP(2,1)=0
       DFDP(3,1)=0
       DFDP(4,1)=U(1)

       DFDP(1,2)=0
       DFDP(2,2)=0
       DFDP(3,2)=0
       DFDP(4,2)=U(1)**3

       DFDP(1,3)=0
       DFDP(2,3)=0
       DFDP(3,3)=0
       DFDP(4,3)=-U(1)**5


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

!----------------------------------------------------------------------
! Problem parameters (only PAR(1-9) are available to the user) :

       U(1)=0.0
       U(2)=0.0
       U(3)=0.0
       U(4)=0.0
       PAR(1)=0.0

!----------------------------------------------------------------------
! If IEQUIB >0 put initial equilibrium in PAR(11+i), i=1,...,NDIM :

      PAR(12) = 0.0
      PAR(13) = 0.0
      PAR(14) = 0.0
      PAR(15) = 0.0


      END SUBROUTINE STPNT

      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
! Homoclinic bifurcations COMMON block needed here :
      COMMON /BLHOM/ ITWIST,ISTART,IEQUIB,NFIXED,NPSI,NUNSTAB,NSTAB,NREV
      INTEGER ITWIST,ISTART,IEQUIB,NFIXED,NPSI,NUNSTAB,NSTAB,NREV
      INTEGER I


! If IEQUIB =0 put analytic equilibrium in PAR(11+i), i=1..NDIM

      IF(IEQUIB.EQ.0)THEN
        DO I=1,NDIM
          PAR(11+I)= 0.0
        ENDDO
      ENDIF

      END SUBROUTINE PVLS

      !----------------------------------------------------------------------
            SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
      !     ---------- ----

            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NDIM,ICP(*),NBC,IJAC
            DOUBLE PRECISION, INTENT(IN) :: PAR(*),U0(NDIM),U1(NDIM)
            DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
            DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

             FB(1)=U0(1)
             FB(2)=U1(1)
             FB(3)=U0(3)
             FB(4)=U1(3)

            END SUBROUTINE BCND
      !----------------------------------------------------------------------

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT
