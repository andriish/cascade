      DOUBLE PRECISION FUNCTION dcasrn(i)

      DOUBLE PRECISION DR(1)
      DOUBLE PRECISION TWOM24
      PARAMETER (TWOM24 = 2.D0**(-24))
      Real Rvec(2)
	Integer i
      CALL RANLUX(RVEC,2)
      DR(1) = DBLE (RVEC(1)) + DBLE (RVEC(2)) * TWOM24
c      CALL RM48 (DR,1)
      dcasrn = DR(1)


      RETURN
      END


      SUBROUTINE dcasrnv(RVEC,LENV)

      Integer lenv,I
      Double Precision RVEC(LENV)
      Double Precision dcasrn,t
      Do I=1,Lenv
         t=dcasrn(i)
         if(t.ne.t) then
            write(6,*) ' DCASRNv: serious problem with Rn = ',t
         elseif(t.gt.1d0.or.t.lt.0.d0) then
            write(6,*) ' DCASRNv: out of range Rn = ',t
         endif
         Rvec(i)=t
      enddo

      RETURN
      END
