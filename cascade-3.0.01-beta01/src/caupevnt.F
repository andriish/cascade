      SUBROUTINE CAUPEVNT
 
      Implicit None
 
#include "caslujets.inc"
#include "caspara.inc"
#include "caspara1.inc"
#include "casherup.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "casludat2.inc"

C...CASCADE process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
      
      Integer NCT,MCT(4000,2),IMO,IDA

      Integer kb1,kb2,kb3,kdq1,kdq2,kdq3,kd1,kd2,kd3,kq

      Integer nupp
	Integer I,J,ncall,iorg
	Data ncall/0/
	Integer Icol
      integer k3(100)
      Integer Ktest,Ktest1
      
      Integer Iup,KUPPO(100)
      
	ncall = ncall+1
      IF(ILHA.eq.-1) then
c        write(6,*) ' write LHA partonic event record to file '
c        write(6,*) ' starting caupevnt ', ncall,nf1,nf2,nia1,nia2,n
c        call pylist(1)
	  IDPRUP = IPRO
        XWGTUP = 1
	  SCALUP = 0.1
	  AQEDUP = -1.
	  AQCDUP = -1.
	  NUP=1
	  ISTUP(NUP) = -1
	  IDUP(NUP) = K(2,2)
	  ICOLUP(1,NUP) = 0
	  ICOLUP(2,NUP) = 0
	  Do J = 1,5
	    PUP(J,NUP) = P(2,J)
	  Enddo
	  NUP=NUP+1
	  ISTUP(NUP) = -1
	  IDUP(NUP) = K(1,2)
	  ICOLUP(1,NUP) = 0
	  ICOLUP(2,NUP) = 0
	  Do J = 1,5
	   PUP(J,NUP) = P(1,J)
	  Enddo

        DO 150 I=21,27
        MCT(I,1)=0
        MCT(I,2)=0
  150   CONTINUE
 
c ala PYTHIA
        NCT=100
        Do i=21,27

          KQ=KCHG(PYCOMP(K(I,2)),2)*ISIGN(1,K(I,2))
          IF(K(I,1).EQ.3.OR.K(I,1).EQ.13.OR.K(I,1).EQ.14) THEN
            IF(K(I,4).NE.0.AND.(KQ.EQ.1.OR.KQ.EQ.2).AND.MCT(I,1).EQ.0)
     &      THEN
              IMO=MOD(K(I,4)/MSTU(5),MSTU(5))
              IDA=MOD(K(I,4),MSTU(5))
c            write(6,*) K(I,4),MSTU(5),IMO,IDA
              IF(IMO.NE.0.AND.MOD(K(IMO,5)/MSTU(5),MSTU(5)).EQ.I.AND.
     &        MCT(IMO,2).NE.0) THEN
                MCT(I,1)=MCT(IMO,2)
              ELSEIF(IMO.NE.0.AND.MOD(K(IMO,4),MSTU(5)).EQ.I.AND.
     &        MCT(IMO,1).NE.0) THEN
                MCT(I,1)=MCT(IMO,1)
              ELSEIF(IDA.NE.0.AND.MOD(K(IDA,5),MSTU(5)).EQ.I.AND.
     &        MCT(IDA,2).NE.0) THEN
                MCT(I,1)=MCT(IDA,2) 
              ELSE
                NCT=NCT+1
                MCT(I,1)=NCT
              ENDIF
            ENDIF
            IF(K(I,5).NE.0.AND.(KQ.EQ.-1.OR.KQ.EQ.2).AND.MCT(I,2).EQ.0)
     &      THEN
              IMO=MOD(K(I,5)/MSTU(5),MSTU(5))
              IDA=MOD(K(I,5),MSTU(5))
              IF(IMO.NE.0.AND.MOD(K(IMO,4)/MSTU(5),MSTU(5)).EQ.I.AND.
     &        MCT(IMO,1).NE.0) THEN
                MCT(I,2)=MCT(IMO,1)
              ELSEIF(IMO.NE.0.AND.MOD(K(IMO,5),MSTU(5)).EQ.I.AND.
     &        MCT(IMO,2).NE.0) THEN
                MCT(I,2)=MCT(IMO,2)
              ELSEIF(IDA.NE.0.AND.MOD(K(IDA,4),MSTU(5)).EQ.I.AND.
     &        MCT(IDA,1).NE.0) THEN
                MCT(I,2)=MCT(IDA,1)
              ELSE
                NCT=NCT+1
                MCT(I,2)=NCT
              ENDIF
            ENDIF
          ENDIF
        Enddo
        
        Do I=21,26
          if(I.eq.21.or.I.eq.23.or.I.eq.25.or.I.eq.26) then
	      NUP=NUP+1
            If(I.eq.21.or.I.eq.23) then 
              ICOLUP(1,NUP)=MCT(I,1)
              ICOLUP(2,NUP)=MCT(I,2)
            Elseif(I.eq.25) then
              ICOLUP(1,NUP)=MCT(27,1)
            Elseif(I.eq.26) then
              ICOLUP(2,NUP)=MCT(27,2)
            Endif
            MOTHUP(1,NUP)=K(I,3)
            MOTHUP(2,I)=0

	      VTIMUP(NUP) = V(I,5)
            SPINUP(NUP) = 9
	      IDUP(NUP) = K(I,2)

	      Do J = 1,5
	        PUP(J,NUP) = P(I,J)
	      End do
          endif
        End do 
      else
      if(IPRO.eq.501) then
         write(6,*) ' caupevnt: IPRO =',IPRO,' not implemented for IHLA=10 '
         write(6,*) ' caupevnt:program stop '
         stop
      endif
      
c      write(6,*) ' starting caupevnt ', ncall,nf1,nf2,nia1,nia2,n
c      call pylist(1)
      If(k(nf1,1).ne.1) k(nf1,1)=-1
      if(nf2.ge.1) then 
         If(k(nf2,1).ne.1) k(nf2,1)=-1
      endif
      if(Iabs(K(2,2)).gt.1000) then
      k(nia2,1)=-1
      k(nia2+1,1)=-1
      endif
      if(iabs(K(1,2)).gt.1000) then
      k(nia1,1)=-1
      k(nia1+1,1)=-1
      endif
c      call pyedit(11)
c      write(6,*) ' after pyedit ', ncall,nf1,nf2,nia1,nia2,n
c      call pylist(1)
c checking diquark and replace it with incoming quark
c for finding it back, incoming quark is 5000+ID      
      kq=-9999
      Ktest = -999
      Ktest1 = -999
      Do I=3,N
        if(K(I,2).eq.10443.or.K(I,2).eq.20443.or.
     &     K(I,2).ge.200000) then
           Ktest=K(I,2)
           Ktest1 = 123
           K(I,2) = 123
           endif  
        if(iabs(K(I,2)).gt.1000.and.iabs(K(I,2)).lt.5000) then
c          write(6,*) ' diquark found ... ',I,K(I,2),K(I,3)
          iorg = K(I,3)
          KB1=MOD(K(iorg,2)/1000,10)
          KB2=MOD(K(iorg,2)/100,10)
          KB3=MOD(K(iorg,2)/10,10)
          
          KDQ1=MOD(K(I,2)/1000,10)
          KDQ2=MOD(K(I,2)/100,10)
          KDQ3=MOD(K(I,2)/10,10)
          kd1=kb1-kdq1
          kd2=kb2-kdq2
          kd3=kb3-kdq3
c          write(6,*) ' origin ',iorg,K(iorg,2)
c          write(6,*) ' kb1 ',kb1,kb2,kb3
c          write(6,*) ' kdq1 ',kdq1,kdq2,kdq3
c          write(6,*) ' kd1 ',kd1,kd2,kd3
          if(kd1.ne.0) write(6,*) 'caupevnt: kd1 ',kd1 
          if(kd2.eq.0.and.kd3.ne.0) then 
            kq = kd3
          elseif(kd2.ne.0.and.kd3.ne.0) then 
c            write(6,*) 'caupevnt: kd2 ',kd2 
            kq = kd2+kd3
          endif
c          write(6,*) ' kq = ',kq
          if(kq.eq.-9999)  write(6,*) ' ERROR kq = ',kq
          Do J=1,3
          P(I,J) = P(iorg,J) - P(I,J)
          Enddo
          P(I,5) = PYMASS(kq)
          P(I,5) = 0
          P(I,4) = Dsqrt(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2)
          K(I,2) = kq+isign(1,K(iorg,2))*5000
c          write(6,*) 'kq ... ',I,kq,K(I,2)
        endif        
      End do
c      call pylist(1)
      NUP=0
	NUP=NUP+1
	ISTUP(NUP) = -1
	IDUP(NUP) = K(2,2)
	ICOLUP(1,NUP) = 0
	ICOLUP(2,NUP) = 0
	Do J = 1,5
	 PUP(J,NUP) = P(2,J)
	Enddo
	NUP=NUP+1
	ISTUP(NUP) = -1
	IDUP(NUP) = K(1,2)
	ICOLUP(1,NUP) = 0
	ICOLUP(2,NUP) = 0
	Do J = 1,5
	 PUP(J,NUP) = P(1,J)
	Enddo
30    continue	
	Icol = 0
      Do I=3,N
       if(K(I,2).eq.123) then 
         K(I,2)=Ktest
         endif
	 IF(K(I,1).eq.1.or.K(I,1).eq.2) then
c          write(6,*) ' start loop ',i,k(i,1),k(i,2),N
	    NUP=NUP+1
	    ICOLUP(1,NUP) = 0
	    ICOLUP(2,NUP) = 0
          MOTHUP(1,nup)=3	
          MOTHUP(2,nup)=1	
	    ISTUP(NUP) = 1
	    IDUP(NUP) = K(I,2)
	    if(K(i,1).eq.2) then 
	      Icol = Icol+1
		if(K(i,2).gt.0.and.iabs(K(i,2)).lt.10) then
 	         ICOLUP(1,NUP) = 500+icol
		   ICOLUP(2,NUP) = 0
		   elseif(K(i,2).gt.0.and.iabs(K(i,2)).gt.1000.and.Ktest1.ne.123) then
               write(6,*) K(I,2),Ktest1,Ktest
		   ICOLUP(1,NUP) = 500+icol
  	         ICOLUP(2,NUP) = 0
c               write(6,*)  ' parton A K',K(I,2)
               K(I,2) = K(I,2)-isign(1,K(I,2))*5000
	         IDUP(NUP) = K(I,2)
     	         ISTUP(NUP) = -1 
               if(K(I,3).eq.1) then
                 MOTHUP(1,nup) = 2
                 elseif(K(I,3).eq.2) then
                 MOTHUP(1,nup) = 1
                 endif
c               write(6,*) ' parton  A ',IDUP(NUP),K(i,2)
		   elseif(K(i,2).lt.0.and.iabs(K(i,2)).lt.10) then
  	         ICOLUP(1,NUP) = 0
		   ICOLUP(2,NUP) = 500+icol
		   elseif(K(i,2).lt.0.and.iabs(K(i,2)).gt.1000.and.Ktest1.ne.123) then
		   ICOLUP(1,NUP) = 0
 	         ICOLUP(2,NUP) = 500+icol
c               write(6,*)  ' parton B K',K(I,2)
               K(I,2) = K(I,2)-isign(1,K(I,2))*5000
	         IDUP(NUP) = K(I,2)
	         ISTUP(NUP) = -1
               if(K(I,3).eq.1) then
                 MOTHUP(1,nup) = 2
                 elseif(K(I,3).eq.2) then
                 MOTHUP(1,nup) = 1
                 endif
c               write(6,*) ' parton B ',IDUP(NUP),K(I,2)
		   elseif(K(i,2).eq.21) then
c look for gluons with KF=21
		      if(ICOLUP(1,NUP-1).gt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 500+icol
			ICOLUP(2,NUP) = 500+icol-1
                  if(MOTHUP(1,NUP-1).eq.1.or.MOTHUP(1,NUP-1).eq.2)then
			ICOLUP(2,NUP) = 500+icol
			ICOLUP(1,NUP) = 500+icol-1
                  endif
		      elseif(ICOLUP(1,NUP-1).lt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 500+icol-1
			ICOLUP(2,NUP) = 500+icol
                  if(MOTHUP(1,NUP-1).eq.1.or.MOTHUP(1,NUP-1).eq.2)then
			ICOLUP(2,NUP) = 500+icol-1
			ICOLUP(1,NUP) = 500+icol
                  endif
			endif
		   elseif(K(i,2).gt.100) then
		   else
  	         ICOLUP(1,NUP) = 0
		   ICOLUP(2,NUP) = 0
		endif   
	    elseif(K(i,1).eq.1) then
c            write(6,*) ' here 1 ',i,k(i,2)
c            write(6,*) ' here 1 col ',ICOLUP(1,NUP-1),ICOLUP(2,NUP-1)
            if(iabs(K(i,2)).gt.1000.and.Ktest1.ne.123) then
c               write(6,*) ' incoming A K',k(i,2),nup
               K(I,2) = K(I,2)-isign(1,K(I,2))*5000
c               write(6,*) ' incoming 2nd A K',k(i,2),nup
	         IDUP(NUP) = K(I,2)
	         ISTUP(NUP) = -1
               if(K(i,2).gt.0) then 
		   ICOLUP(1,NUP) = 500+icol
  	         ICOLUP(2,NUP) = 0
               else
  	         ICOLUP(1,NUP) = 0
		   ICOLUP(2,NUP) = 500+icol
               endif 
               MOTHUP(1,nup) = K(I,3)
               if(K(I,3).eq.1) then
                 MOTHUP(1,nup) = 2
                 elseif(K(I,3).eq.2) then
                 MOTHUP(1,nup) = 1
                 endif
c               write(6,*) ' this is for incoming part A ',IDUP(NUP)
	      elseif(ICOLUP(1,NUP-1).ne.0.or.ICOLUP(2,NUP-1).ne.0) then
 		    if(ICOLUP(1,NUP-1).gt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 0
			ICOLUP(2,NUP) = 500+icol
		    elseif(ICOLUP(1,NUP-1).lt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 500+icol
			ICOLUP(2,NUP) = 0
	          endif
	      else
c            write(6,*) ' enter here ',(K(I,J),J=1,3),(P(I,J),J=1,5)
            endif
	    endif
       elseif(K(I,2).LT.1000.or.Ktest1.eq.123) then
	    NUP=NUP+1
	    ISTUP(NUP) = 3
	    ICOLUP(1,NUP) = 0
	    ICOLUP(2,NUP) = 0
	    IDUP(NUP) = K(I,2)
	    MOTHUP(2,NUP) = 0  ! special
	    if(K(I,3).eq.1) then 
	      MOTHUP(1,NUP) = 2
	    elseif(K(I,3).eq.2) then 
	      MOTHUP(1,NUP) = 1
	    else
	      MOTHUP(1,NUP) = 3
	      MOTHUP(2,NUP) = 1
	    endif
100       continue
	 Else
          goto 111
       Endif
	 VTIMUP(NUP) = V(I,5)
	 SPINUP(NUP) = 9

	 Do J = 1,5
	   PUP(J,NUP) = P(I,J)
	 Enddo
	
111   Continue
	Enddo 
      
c      write(6,*) 'before reshuffling'
c      call pylist(7)
      do I=3,NUP
c      write(6,*) ' checking here ',ISTUP(I)
      if(ISTUP(I).lt.0) then 
        iorg=mothup(1,i)
	  Do J = 1,5
	   PUP(J,iorg) = PUP(J,I) 
         PUP(J,I) = 0
	  Enddo
	  ISTUP(iorg) = ISTUP(I)
	  IDUP(iorg) =  IDUP(I)
	  IDUP(iorg) =  IDUP(I)
	  ICOLUP(1,iorg) = ICOLUP(1,I) 
	  ICOLUP(2,iorg) = ICOLUP(2,I) 
        
	  ISTUP(I) = 3
	  IDUP(I) =  0
	  ICOLUP(1,I) = 0
	  ICOLUP(2,I) = 0
	  MOTHUP(1,1) =  0
 	  MOTHUP(2,1) =  0
	  MOTHUP(1,2) =  0
 	  MOTHUP(2,2) =  0
      Endif
      enddo
      
      nupp =0
      do i=1,nup
      if(idup(i).eq.0) then
      else
         nupp=nupp+1
         istup(nupp)=istup(i)
         icolup(1,nupp) = icolup(1,i)	
         icolup(2,nupp) = icolup(2,i)	
	   IDUP(nupp) =  IDUP(I)
	   MOTHUP(1,nupp) = MOTHUP(1,i)  
 	   MOTHUP(2,nupp) = MOTHUP(2,i)  
	   Do J = 1,5
	     PUP(J,nupp) = PUP(J,I) 
	   Enddo
      endif
      enddo
            
      
      nup=nupp
	IDPRUP = IPRO
      XWGTUP = 1
	SCALUP = 0.1
	AQEDUP = -1.
	AQCDUP = -1.
      
	if(ncall.lt.5) then 
c         write(6,*) ' before print'
c         call pylist(7)
      endif
c         call pylist(7)
c         write(6,*) ' before print'
c         call pylist(7)
      endif
      
      IF(ILHA.ge.10) then        
c         write(6,*) ' starting pyevnt: ncall=',ncall,' nup =',nup
c         call pylist(1)
c         call pylist(7)
c	   if(ncall.lt.5) call pylist(7)
	   do i=1,4000
	   do j=1,5
	   p(i,j)=0
	   k(i,j)=0
	   enddo
	   enddo
         IF(NUP.LT.100) then 
	      CALL PYEVNT
         ELSE
            write(6,*) ' too many partons to be handled by PYTHIA ', NUP
         ENDIF
c         write(6,*) ' after pyevnt'
c         call pylist(1)
c	   if(ncall.lt.50) call pylist(1)
      elseif(ILHA.le.-2) then
c read in LHA event recors
        write(6,*) ' caupevnt dummy version of read in record '      
	  Read(21) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &   (IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),
     &   ICOLUP(2,I),(PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I),I=1,NUP)
        write(6,*) 'caupevnt nup=',NUP 
        I=0
        DO 470 IUP=1,NUP
          KUPPO(IUP)=0
          I=I+1
          KUPPO(IUP)=I
          K(I,1)=21
          K(I,2)=IDUP(IUP)
          IF(IDUP(IUP).EQ.0) K(I,2)=90
          K(I,3)=0
          IF(MOTHUP(1,IUP).GE.3) K(I,3)=KUPPO(MOTHUP(1,IUP))
          K(I,4)=0
          K(I,5)=0
          DO 460 J=1,5
            P(I,J)=PUP(J,IUP)
  460     CONTINUE
          V(I,5)=VTIMUP(IUP)
  470   CONTINUE
        N=I
 
        call pylist(2)
     
	else
c      write(6,*) ' cascade RGupevnt '
	   Write(21) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &   (IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),
     &   ICOLUP(2,I),(PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I),I=1,NUP)
	
c         write(6,*) ' cascade caupevnt ',ncall
	   if(ncall.lt.5) call pylist(7)
c	   if(ncall.lt.5) call pylist(1)
	endif
c	call pylist(1)
 
      RETURN
      END 
 
