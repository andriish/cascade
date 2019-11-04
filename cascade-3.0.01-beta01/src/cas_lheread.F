      Subroutine cas_lheinit
      Implicit None
#include "caslujets.inc"
#include "caseffic.inc"
#include "caslhe.inc"

      character lhefile1*132
       LOGICAL ex
      Integer lun,ipr, ievent
      character* 200 string
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
      
      lun = 25
      ievent=0
      lhefile1=trim(lhefile)
      inquire(FILE=lhefile,EXIST=ex)
      if(ex) then
	   Write(6,*) '*     File exists ',trim(lhefile1)
         OPEN(lun,FILE=lhefile1,STATUS='OLD',ERR=100)
         else
         write(6,*) '*           input lhe file ',lhefile1,' not found    *'
         write(6,*) '*           Program stopped                         *'
         write(6,*) '*****************************************************'
         stop
      endif
      avgi = 0.
      sd = 0.
 1    continue
      string=' '
      read(lun,fmt='(a)',err=777,end=666) string
c     powheg files start in 1st column, aMCatNLO files at 3rd column
c      write(6,*) string(3:9)
      if(string(1:6).eq.'<init>'.or.string(3:8).eq.'<init>') then
         read(lun,*) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         write(6,101) idbmup(1),idbmup(2),ebmup(1),ebmup(2)
  101    Format(' *     Beam Id = ',I6,1X,I6,'                       *',/, 
     &          ' *     Beam energies =  ',E10.4,1X,E10.4, '        *')    
         write(6,102) pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2)
  102    Format(' *     collinear pdf Id = ',I4,1X,I4,'                  *'/
     &     ' *                  ',I7,1X,I7,'                  *')    
         write(6,103) idwtup,nprup
  103    Format(' *     Idwtup = ',I6,' nprup = ',I6,'                *') 
         do ipr=1,nprup
            read(lun,*) xsecup(ipr),xerrup(ipr),xmaxup(ipr),lprup(ipr)                     
            avgi = avgi + xsecup(ipr)
            sd = sd + xerrup(ipr)**2
         
         enddo
         goto 1
      elseif(string(1:6).eq.'<event'.or.string(3:8).eq.'<event') then
         ievent = ievent + 1
         goto 1
      else 
         goto 1
      endif
 
999   return
         
100	Write(6,*) ' Error! File open failed for ',lhefile1,'  STOP '
      return
777   Write(6,*) ' Error! reading file   '
      return
666   Write(6,*) '*     End of file. Events read: = ',ievent,'     *'
      if(nevent.le.0) nevent = ievent
      
      if(idwtup.eq.1) then      
         normhepmc = 1d0/dabs(avgi)
      elseif(idwtup.eq.3) then
         normhepmc = 1d0
      elseif(idwtup.eq.-4) then      
         normhepmc = 1d0/dabs(avgi)
      else
         normhepmc = 1d0/dabs(avgi)
      endif
C xsection in LHE is given in pb 
C CASCADE wants xsection in nb:
      avgi =avgi/1000.
      sd = sqrt(sd)/1000.
      
      write(6,6001) avgi, sd
6001  format(' *     LHE xsection [nb]:',G10.3,'+-',G10.3,'      *')     
      rewind lun
      return     
      End
      subroutine lhrwgt_id_value_ind(string,iid,iendid,value)
      implicit none
      character *(*) string
      integer iid,iendid
      real * 8 value
      integer j,k,l,ios
      j=1
      l=len(string)
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j+4).ne.'<wgt ') goto 777
      j=j+5
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j+1).ne.'id') goto 777
      j=j+2
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j).ne.'=') goto 777
      j = j + 1
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      select case (string(j:j))
      case ('"',"'")
         iid = j
         iendid = index(string(iid+1:),string(iid:iid))+iid
         j=iendid+1
      case default
         goto 777
      end select
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      if(string(j:j).ne.'>') goto 777
      j = j+1
      do while(string(j:j).eq.' '.and.j.lt.l)
         j = j + 1
      enddo
      k = index(string,'</wgt>')
      if(k.eq.0) goto 777
      read(unit=string(j:k-1),fmt=*,iostat=ios) value
      if(ios.ne.0) goto 777
      return
c Something went wrong:
 777  iid = -1
      end
      
      Subroutine cas_lheread
c this is a copy of lheread from powheg      
      Implicit None
#include "caslujets.inc"
#include "casludat1.inc"
#include "casludat2.inc"
#include "casludat3.inc"
#include "caslhe.inc"
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
      character* 200 string
      Integer lun,i,j,iup,i1,i2,i3,i4,ILBL,ISDE1,ISDE2,IUP1,IUP2,NLBL,NMAT
      Integer KUPPO(100),ILAB(100),IUPDAU
      
      Character *4 name
      character *4 cdum
      character* 200 namepdf     
      
      double Precision scalepdf
      
      Integer ncall
      data ncall/0/
      double precision x1,x2
      Common/myxval/x1,x2
      
      Double Precision  pi, phi1,phi2,kt2min,kt2max,kt21,kt22,rn1,rn2,rn3
      Double Precision wt1,wt2,wt1_pdf,wt2_pdf,wtmax,wt_tmd,wt_ratio
      double Precision xfa(-6:6),xfb(-6:6),xfa_max(-6:6),xfb_max(-6:6)
      Double Precision xfamax,xfbmax,scalfac
      Integer Id
	Double Precision DETOT,DBETA(2,3),STHETA(2),SPHI(2),cmpart(4),cm(4),sqrts,e1,e2,p1,p2,shat,shat1
      Double Precision shatr,shats,pms(2)
      Double Precision Psave(20,5)
      Double Precision dcasrn,dot1,dot
      
      Double Precision db,dbp,dga,dgabp,ptest
      
      Double Precision x1_bar, x2_bar, sss, mu_pdf1, mu_pdf2, xsi_hard(2)
      Double Precision weightnew,value
      integer iid,iendid,i_ev,nlheWid
            
      Integer itest,Iopt
      Integer Ntry, ntry_tmd
      External dcasrn,dot1,dot
      pi= 4.D0*datan(1.D0)
c initialise tmp array
      do i=1,100
        kuppo(i) = 0
        ilab(i) = 0
      end do

      lun = 25
      
 1    continue
      string=' '
      read(lun,fmt='(a)',err=777,end=666) string
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      i_ev=INDEX(STRING,'<event')
c      write(6,*) 'cas_lheread: ',i_ev, string(i_ev:i_ev+5)
      if(string(i_ev:i_ev+5).eq.'<event') then
c      if(string(1:6).eq.'<event'.or.string(3:8).eq.'<event') then
c on error try next event. The error may be caused by merging
c truncated event files. Thus we are tolerant about it.
c Only on EOF return with no event found
         read(lun,*,end=998,err=778)nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(lun,*,end=998,err=1) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         end do
c         write(6,*) ' cas_lheread: ',hepmcout,lhefile,Ilhe,iTMDw,Collfile
c check for scale of process:
         scale_process = scalup
cccc endnew for powheg           
         xidP1 = pup(4,1)/ebmup(1)
         xidP2 = pup(4,2)/ebmup(2)
         idP1 = idup(1)
         idP2 = idup(2)
         x1=xidP1
         x2=xidP2   
         x1_bar=0
         x2_bar=0
         xsi_hard(1) = (x2+x1_bar)/(x1+x2_bar)
         xsi_hard(2) = (x1+x2_bar)/(x2+x1_bar)
         sss=4*ebmup(1)*ebmup(2)
         if(Iscale.eq.0) then 
            scale_process = scalup
         elseif(Iscale.eq.1) then	
            scale_process = max(x1*sqrt(xsi_hard(1)*SSS),x2*sqrt(xsi_hard(2)*SSS))
c         elseif(Iscale.eq.2) then 
c            if(scalepdf.le.0) then
c                 write(6,*) ' Error cas_lheread: Iscale = ',Iscale,' has scalepdf = ',scalepdf
c                 stop
c            endif
c            scale_process = scalepdf
         else
            write(6,*) ' Error cas_lheread: Iscale = ',Iscale,' not defined '
            stop
         endif
         scalup_cas = scalup
c Radek       
         if(scale_process .le. 1.5d0 ) scale_process = 1.5d0
Cend Radek           
         mu_pdf1 = scale_process
         mu_pdf2 = scale_process
         weight = xwgtup
         string=' '
         read(lun,fmt='(a)',err=777,end=666) string
c         write(6,*) ' check for additional params ', string
         i_ev=INDEX(STRING,'</event')
567         if(i_ev.le.0) then 
c             write(6,*) ' ananlysing additonal paramters ',string
             if(INDEX(STRING,'<rwgt>').gt.0) then
                nlheWid=0
5667            read(lun,fmt='(a)',err=777,end=666) string 
c                write (6,*) ' after <rwgt> ', string
                if (INDEX(STRING,'<wgt id=').gt.0 ) then
                   nlheWid=nlheWid+1
c call subroutine from POWHEG
                   call lhrwgt_id_value_ind(string,iid,iendid,value)
c                   write(6,*) ' wgt id = : ',string(iid:iendid),value
                   if(nlheWid.eq.lheWid) weight = value
c                write(6,*) ' string:',string(iid:iendid),'  CWGT = ',Collfile
                   if(string(iid+1:iendid-1).eq.Collfile) then
                      weightnew = value
                      weight = weightnew
                   endif
                   goto 5667
                endif 
            endif
c            write(6,*) ' test weight: lheWid = ',lheWid,weight
            if(INDEX(STRING,'#pdf  ').gt.0) then
               read(string,*,end=998,err=1) name, idP1,idP2,xidP1,xidP2,scalepdf,pdf1,pdf2
c              write(6,*) ' read infos ',name, idP1,idP2,xidP1,xidP2,pdf1,pdf2,scalepdf
            elseif(INDEX(STRING,'   #pdf1pdf2').gt.0) then
c            write(6,*) ' pdf string found ',string(1:12)
               idP1 = -9999
               idP2 = -9999
               xidP1 = -9999.
               xidP2 = -9999.
               kt2idP1 = -9999.
               kt2idP2 = -9999.
               read(string,*,end=998,err=1) name, pdf1pdf2 
               idP1 = idup(1)
               idP2 = idup(2)
               xidP1 = pup(4,1)/ebmup(1)
               xidP2 = pup(4,2)/ebmup(2)
               if(iTMDw.ge.1) then
                  idP1 = idup(1)
                  idP2 = idup(2)
                  xidP1 = pup(4,1)/ebmup(1)
                  xidP2 = pup(4,2)/ebmup(2)
                  sss=4*ebmup(1)*ebmup(2)
                  kt2idP1 = pup(1,1)**2 + pup(2,1)**2
                  kt2idP2 = pup(1,2)**2 + pup(2,2)**2   
                  mu_pdf1 = scalup
                  mu_pdf2 = scalup

                  if(iTMDw.eq.2) then
               
                     x1=xidP1
                     x2=xidP2   
                     x2_bar=0
                     x1_bar=0
                     xsi_hard(1) = (x2+x1_bar)/(x1+x2_bar)
                     xsi_hard(2) = (x1+x2_bar)/(x2+x1_bar)	
                     mu_pdf1 = x1*sqrt(xsi_hard(1)*SSS)
                     mu_pdf2 = x2*sqrt(xsi_hard(2)*SSS)
               
                  endif
                  CALL cauniglu(idbmup(1),xidP1,kt2idP1,mu_pdf1,xfa) 
                  CALL cauniglu(idbmup(2),xidP2,kt2idP2,mu_pdf2,xfb) 
                  if(idP1.eq.21) idP1=0
                  if(idP2.eq.21) idP2=0
                  pdf_weight=1.d0/pdf1pdf2
                  if(kt2idP1.ge.1E-5 ) pdf_weight = pdf_weight*xfa(idP1)
                  if(kt2idP1.ge.1E-5 ) pdf_weight = pdf_weight*xfb(idP2)
                  weight = weight * pdf_weight
                  if(pdf_weight.le.0.or.pdf1pdf2.le.0) then
                     write(6,*) ' TMD reweighting problem: new TMD weight ', pdf_weight
                     write(6,*) ' parton 1 (id,x,kt2,mu) ',idP1,xidP1,kt2idP1,mu_pdf1,xfa(idP1)
                     write(6,*) ' parton 2 (id,x,kt2,mu) ',idP2,xidP2,kt2idP2,mu_pdf2,xfb(idP2)
                     write(6,*) ' pdf weights ',idP1,idP2, pdf_weight,xfa(idP1)*xfb(idP2),pdf1pdf2
                  endif 
               endif
            endif
c            write(6,*) ' end ', string
            read(lun,fmt='(a)',err=777,end=666) string
            i_ev=INDEX(STRING,'</event')
            if(i_ev.le.0) goto 567 
            
         else
            write(6,*) ' end of event reached, read new event' 
         endif
         
c         write(6,*) ' after pdfweights '
         ncall = ncall + 1 
         Nf_part = nup-2
         do  I=1,MSTU(4)
           do  J=1,5
           K(I,J)=0
           P(I,J)=0.D0
           V(I,J)=0.d0
           end do
         end do
         
c         write(6,*) ' after read ',ncall
         if (ncall .le.5) call pylist(7)
         K(1,1)=21
         K(1,2)=idbmup(1)
         K(2,2)=idbmup(2)
         K(2,1)=21
         P(1,4)=ebmup(1)
         P(1,3)=ebmup(1)
         P(1,5) = 0
         P(2,4)=ebmup(2)
         P(2,3)=-ebmup(2)
         P(2,5) = 0
         I=2
c         write(6,*) ' sign of 1st aprticel ',sign(1d0,PUP(3,1)),sign(1d0,PUP(3,2))

C this is taken form pythia (pyrand)
        KUPPO(1)=1
        KUPPO(2)=2

         N=I
C...Store final state partons for user-defined processes
        DO 490 IUP=1,NUP
          N=N+1
          K(N,1)=1
          IF(ISTUP(IUP).EQ.2.OR.ISTUP(IUP).EQ.3) K(N,1)=11
          K(N,2)=IDUP(IUP)
          IF(IDUP(IUP).EQ.0) K(N,2)=90
          IF(MOTHUP(1,IUP).EQ.0) THEN
            K(N,3)=KUPPO(IUP)
          ELSE
            K(N,3)=MOTHUP(1,IUP)
          ENDIF
          K(N,4)=0
          K(N,5)=0
C...Search for daughters of intermediate colourless particles.
          IF(K(N,1).EQ.11.AND.KCHG(PYCOMP(K(N,2)),2).EQ.0) THEN
            DO 475 IUPDAU=IUP+1,NUP
              IF(MOTHUP(1,IUPDAU).EQ.IUP.AND.K(N,4).EQ.0) K(N,4)=
     &        N+IUPDAU-IUP
              IF(MOTHUP(1,IUPDAU).EQ.IUP) K(N,5)=N+IUPDAU-IUP
  475       CONTINUE
          ENDIF
          DO 480 J=1,5
            P(N,J)=PUP(J,IUP)
  480     CONTINUE
          V(N,5)=VTIMUP(IUP)
  490   CONTINUE
c         write(6,*) ' before color flow '
c         call pylist(2)
C...Arrange colour flow for user-defined processes
        NLBL=0
        DO 540 IUP1=1,NUP
          I1=I+IUP1
          IF(KCHG(PYCOMP(K(I1,2)),2).EQ.0) GOTO 540
          IF(K(I1,1).EQ.1) K(I1,1)=3
          IF(K(I1,1).EQ.11) K(I1,1)=14
C...Find a not yet considered colour/anticolour line.
          DO 530 ISDE1=1,2
            IF(ICOLUP(ISDE1,IUP1).EQ.0) GOTO 530
            NMAT=0
            DO 500 ILBL=1,NLBL
              IF(ICOLUP(ISDE1,IUP1).EQ.ILAB(ILBL)) NMAT=1
  500       CONTINUE
            IF(NMAT.EQ.0) THEN
              NLBL=NLBL+1
              ILAB(NLBL)=ICOLUP(ISDE1,IUP1)
C...Find all others belonging to same line.
              I3=I1
              I4=0
              DO 520 IUP2=IUP1+1,NUP
                I2=2+IUP2
                DO 510 ISDE2=1,2
                  IF(ICOLUP(ISDE2,IUP2).EQ.ICOLUP(ISDE1,IUP1)) THEN
                    IF(ISDE2.EQ.ISDE1) THEN
                      K(I3,3+ISDE2)=K(I3,3+ISDE2)+I2
                      K(I2,3+ISDE2)=K(I2,3+ISDE2)+MSTU(5)*I3
                      I3=I2
                    ELSEIF(I4.NE.0) THEN
                      K(I4,3+ISDE2)=K(I4,3+ISDE2)+I2
                      K(I2,3+ISDE2)=K(I2,3+ISDE2)+MSTU(5)*I4
                      I4=I2
                    ELSEIF(IUP2.LE.2) THEN
                      K(I1,3+ISDE1)=K(I1,3+ISDE1)+I2
                      K(I2,3+ISDE2)=K(I2,3+ISDE2)+I1
                      I4=I2
                    ELSE
                      K(I1,3+ISDE1)=K(I1,3+ISDE1)+MSTU(5)*I2
                      K(I2,3+ISDE2)=K(I2,3+ISDE2)+MSTU(5)*I1
                      I4=I2
                    ENDIF
                  ENDIF
  510           CONTINUE
  520         CONTINUE
            ENDIF
  530     CONTINUE
  540   CONTINUE

c         write(6,*) ' after color flow '       
c         call pylist(2)
        if(Ilhe.eq.1) then
c generate kt according to TMD density
           ntry = 0
           do i=1,4        
             cmpart(i) = P(3,I) + P(4,i)
           end do 
           shat=dot(cmpart,cmpart)
           do  I=1,4
             CM(I) =P(1,I)+P(2,I)
           end do
           sqrts=sqrt(dot(cm,cm))
           do j=1,19
             do i=1,5
               Psave(J,I)=P(J,I)
             end do
           end do
455        Continue
           do j=1,19
             do i=1,5
               P(J,I)=Psave(J,I)
             end do
           end do
c Radek       
           if(scale_process .le. 1.5d0 ) scale_process = 1.5d0
Cend Radek           
           kt2min = 0.1
ccc           kt2max = scale_process**2
           kt2max = scalup**2
           ntry_tmd = 0
           Wtmax = 1.
           scalfac = 8.
544        RN1=dcasrn() 
           RN2=dcasrn()
c           write(6,*) RN1,RN2
           phi1=dcasrn()*2.*pi
           phi2=dcasrn()*2.*pi
           kt21 = 0.001
           ntry_tmd = ntry_tmd + 1
           if(ntry_tmd.gt.200000) then
           	  write(6,*) ' cas_lheread: tmd ntry  '
              goto 545
           end if
           wt1 = 1.
           wt1_pdf = 1.
           xfamax = 0.
           id=idP1
           if(idP1.eq.21) id=0
           xfa(id) = 0 .
           if(xidP1.le.0.95.and.abs(idbmup(1)).eq.2212) then 
             kt21 = kt2min*(kt2max/kt2min)**RN1   
             wt1 = kt21                 
             id=idP1
             if(idP1.eq.21) id=0
             CALL cauniglu(idbmup(1),xidP1,kt2min,mu_pdf1,xfa_max) 
             xfamax=xfa_max(id)        
             CALL cauniglu(idbmup(1),xidP1,kt21,mu_pdf1,xfa)
c             write(6,*) ' 1st ',id,xidP1,kt2min,mu_pdf1,xfa_max,kt21
             wt1=wt1*xfa(id)
             wt1_pdf=xfa(id)/pdf1 
           endif      
           kt22 = 0.001 
           wt2 = 1.
           wt2_pdf = 1.
           xfbmax = 0.
           id=idP2
           if(idP2.eq.21) id=0
           xfb(id) = 0 .
           if(xidP2.le.0.95.and.abs(idbmup(2)).eq.2212) then
             kt22 = kt2min*(kt2max/kt2min)**RN2 
             wt2 = kt22 
             CALL cauniglu(idbmup(2),xidP2,kt2min,mu_pdf2,xfb_max)
             xfbmax=xfb_max(id)
             CALL cauniglu(idbmup(2),xidP2,kt22,mu_pdf2,xfb)
             wt2=wt2*xfb(id)
             wt2_pdf=xfb(id)/pdf2
c             write(6,*) 'wt2 ',xfb(id),idP2,pdf2,kt22,scale
           endif
           wt_tmd = wt1*wt2
           RN3=dcasrn()
c           wtmax=max(xfamax,xfbmax)*scalfac 
           wtmax=max(1.,xfamax)*max(1.,xfbmax)*1.5
c             write(6,*) ' xfa_max = ',xfamax,' xfb_max = ',xfbmax,' kt2min = ',kt2min,xidP2
c             write(6,*) ' xidP1 = ',xidP1,' xidP2 = ',xidP2 
           if(wtmax.le.0) then
             write(6,*) ' xfa_max = ',xfamax,' xfb_max = ',xfbmax,' kt2min = ',kt2min,xidP2
             write(6,*) ' xidP1 = ',xidP1,' xidP2 = ',xidP2 
             goto 545
           endif
           wt_ratio = wt_tmd/wtmax
c              write(6,*) ' cas_lheread: wt_tmd,  wtmax ',wt_tmd/wtmax,wt_tmd,wtmax,kt22,mu_pdf1,mu_pdf2
           if(wt_ratio.gt.1.) then
              write(6,*) ' cas_lheread: wt_tmd > wtmax ',wt_tmd/wtmax,wt_tmd,wtmax,kt22
           end if  
           if(wt_ratio.le.RN3) goto 544
545        continue           


c           write(6,*) ' cas_lheread: sclae process ',scale_process,' scalup = ',scalup
           scalup_cas = scalup
           x1 = (P(3,4)+P(3,3))/sqrts
           x2 = (P(4,4)-P(4,3))/sqrts

c              call pylist(1)
           P(3,1) = sqrt(kt21)*cos(phi1)
           P(3,2) = sqrt(kt21)*sin(phi1)
           P(3,5) = -sqrt(abs(dot1(3,3)))
           P(4,1) = sqrt(kt22)*cos(phi2)
           P(4,2) = sqrt(kt22)*sin(phi2)   
           P(4,5) = -sqrt(abs(dot1(4,4)))
           PMS(1) = P(3,1)**2+P(3,2)**2 
           PMS(2) = P(4,1)**2+P(4,2)**2
c              write(6,*) ' shat ',shat   
           if(shat.le.0) goto 455
           if((shat-PMS(1)-PMS(2))**2-4D0*PMS(1)*PMS(2).LE.0D0) then
c                write(6,*) ' problem with shat and masses ',sqrt(shat),sqrt(pms(1)),sqrt(pms(2))
                goto 455
           endif
           if((shat-P(3,5)**2 - P(4,5)**2 - P(3,5)**2*P(4,5)**2/shat - PMS(1) -PMS(2)).le.0) then 
c                write(6,*) ' 2nd problem with shat and masses '
                goto 455
           endif 
           e1 = (shat - p(3,5)**2 + p(4,5)**2)/2./sqrt(shat)
           e2 = (shat + p(3,5)**2 - p(4,5)**2)/2./sqrt(shat)
           p1 = ((shat + p(3,5)**2 + p(4,5)**2)**2 - 4.*p(3,5)**2*p(4,5)**2)/4./shat
           p1 = sqrt(p1) 
           p2=-p1
           DETOT=P(3,4)+P(4,4)
           DBETA(2,1)=(P(3,1)+P(4,1))/DETOT
           DBETA(2,2)=(P(3,2)+P(4,2))/DETOT
           DBETA(2,3)=(P(3,3)+P(4,3))/DETOT

            
              
           ntry = ntry + 1
           if(Ntry.gt.100) then
              write(6,*) ' cas_lhread ntry limit reached = ',ntry
              goto 456
           endif
           if(sqrt(DBETA(2,1)**2+DBETA(2,2)**2+DBETA(2,3)**2).ge.1) then
                write(6,*) ' boost too large ',sqrt(DBETA(2,1)**2+DBETA(2,2)**2+DBETA(2,3)**2)
                write(6,*) ' beta ',DBETA(2,1),DBETA(2,2),DBETA(2,3)
              goto 455
           endif

           CALL PYROBO(3,4,0d0,0d0,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3)) 
           SPHI(2)=PYANGL(P(3,1),P(3,2))
           CALL PYROBO(3,4,0d0,-SPHI(2),0d0,0d0,0d0)
c              write(6,*) ' after phi '
c              call pylist(1)
           STHETA(2)=PYANGL(P(3,3),P(3,1))
           CALL PYROBO(3,4,-STHETA(2),0d0,0d0,0d0,0d0)
c              write(6,*) ' after theta '
           P(3,3)=p1
           P(4,3)=p2
           P(3,4)=e1
           P(4,4)=e2             

           CALL PYROBO(5,N,0d0,0d0,0d0,0d0,-DBETA(2,3)) 
c              write(6,*) ' in cm frame: q2 = ',sqrt(abs(dot1(3,3))),sqrt(abs(dot1(4,4))),sqrt(abs(dot1(5,5)))    
c              call pylist(1)
           DETOT=P(3,4)+P(4,4)
              CALL PYROBO(3,0,STHETA(2),SPHI(2),DBETA(2,1),DBETA(2,2),DBETA(2,3))

c              call pylist(1)
           P(3,5) = -sqrt(abs(dot1(3,3)))
           P(4,5) = -sqrt(abs(dot1(4,4)))
           x1 = (P(3,4)+P(3,3))/sqrts
           x2  = (P(4,4)-P(4,3))/sqrts
c              write(6,*) ' before check:  x1, x2 ',x1,x2,sqrt(kt21),sqrt(kt22)
           if(x1.ge.1.or.x2.ge.1) then
c              write(6,*) ' rejected end:  x1, x2 c',x1,xidP1,x2,xiDP2,sqrt(kt21),sqrt(kt22),sqrt(shat)
c              write(6,*) ' x1*x2 ',x1*x2,xidP1*xidP2,xidP1*xidP2-x1*x2
c                 write(6,*) ' rejected boost ', sqrt(DBETA(2,1)**2+DBETA(2,2)**2+DBETA(2,3)**2),DBETA(2,3)
c                 call pylist(7)
c                 call pylist(1)
              goto 455
              else
c                 write(6,*) ' accepted end:  x1, x2 ',x1,xidP1,x2,xiDP2,sqrt(kt21),sqrt(kt22),sqrt(shat)
c                 write(6,*) ' accepted x1*x2 ',x1*x2,xidP1*xidP2,xidP1*xidP2-x1*x2
c                 write(6,*) ' accepted boost ', sqrt(DBETA(2,1)**2+DBETA(2,2)**2+DBETA(2,3)**2),DBETA(2,3)
c                 call pylist(7)
c                 call pylist(1)
           end if
456        continue  
        
       End If      
c       write(6,*) ' cascade weight ',weight
       else
         goto 1
      endif
c no event found:
      return
 777  continue
      print *,"Error in reading 1 time events"
      print *,string
      nup=0
      return
 778  continue
      print *,"Error in reading events", ncall
      print *,string
      nup=0
      return
 666  continue
      print *,"reached EOF after ",ncall,' events '
      print *,string
      nup=0
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      nup=0      
 999  Continue
      return
      end
 
      Subroutine cas_lheprocess
      Implicit None
#include "caslujets.inc"
#include "caspart.inc"
#include "casludat1.inc"
#include "caspar.inc"
#include "caspara5.inc"
#include "casluco.inc"
#include "casgki.inc"
#include "casscal.inc"
#include "casshwr.inc"
#include "caslhe.inc"
	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar

      DOUBLE PRECISION DETOT,DBETA,STHETA,SPHI
      COMMON /MYBOOST/ DETOT,DBETA(2,3),STHETA(2),SPHI(2)
      Double Precision XPY,SH
      Integer KFL
      COMMON /MYPROC/KFL(3,2),XPY(2),SH
      Double Precision dboost
      Integer NVMAX
	Parameter (NVMAX=5000)
      Double Precision PS(nvmax,5),PSAVE(2,5)
      Integer KS(nvmax,5)
      Double Precision dot1,dot
      Integer I,J,NS,ii,nfin,Nn,Jn
      Integer IPU1,IPU2
     
      Double Precision sqrts
      Double Precision yph,yph_bar
      
      Double Precision PFIN(4),QMAX, PTfin(100),Mfin(100)
      DOUBLE PRECISION QG21,QG22
      COMMON/SEMIH/ QG21,QG22      
	Integer kf
	Common /cbran1/ kf
      double precision dx,dkt2,dp,xfa(-6:6)
      Integer Iremn
      common/casflags/Iremn	

      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

      Double Precision chec,p1,p2,p3,p4
      Double Precision pc1,pc2,pc3,pc4
      
      double precision x1,x2
      Common/myxval/x1,x2
	Integer MSTJO41,MSTJO45
      
      Integer nfsh1,nfsh2,nfin1,nfin2
      
      Integer IPU1O,IPU2O
      
      Integer ktest4,ktest5,nold

C... CALCULATE CMS ENERGY
      DO 20 I=1,4
         CM(I) =P(1,I)+P(2,I)
   20 CONTINUE
C BOOST TO CMS
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4))
      DO  I=1,2
         DO  J=1,5
            PBEAM(I,J)=P(I,J)
            KBEAM(I,J)=K(I,J)  
         End do
      End do    

      KINT(1,2)= K(1,2)
      KINT(2,2)= K(2,2)     
      SSS = DOT(CM,CM)
      sqrts=dsqrt(sss)
      ISEMIH = 1
      KFL(2,1)=K(3,2)
      KFL(2,2)=K(4,2)
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      PYVAR(2)=2.D0*DOT1(1,2)
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(5)=PYVAR(1)/2.

      if(Ilhe.ne.1) then
        xg1     = (P(3,4)+P(3,3))/sqrts
        xg     = (P(4,4)-P(4,3))/sqrts
       else
        xg1     = (P(3,4)+P(3,3))/sqrts
        xg     = (P(4,4)-P(4,3))/sqrts

      endif
      PYVAR(31) = xg1
      PYVAR(32) = xg 
c      write(6,*) ' cas_lheread ', xg1,xg
      if(xg1.gt.1.or.xg.gt.1) then 
        write(6,*) ' cas_lheprocess: x > 1 :', xg1,xg,sqrts
        call pylist(7)
        call pylist(1)
        LST(21)=201
        Return
        
      endif
c variables for initial shower
      QG21 = P(3,1)**2+P(3,2)**2
      QG22 = P(4,1)**2+P(4,2)**2
c check here, if we want to have DGLAP type cut on kts:
       if(ICCFM.eq.0) then
        if(sqrt(QG21).gt.scalup/2..or.sqrt(qg22).gt.scalup/2.) then
            LST(21)=200
            Return
         endif
      Endif
      KF = K(1,2)
c  call once cauniglu with dummy aruments to set up grid 
      dx=0.01
      dkt2= 1.
      dp= 100.
      CALL cauniglu(2212,dx,dkt2,dp,xfa)
c
      
      
      do I=1,5
            PSAVE(1,I)=P(1,I)
            PSAVE(2,I)=P(2,I)
      End do
      DBETA(2,1)=0d0
      DBETA(2,2)=0d0
      SPHI(2)=0d0
      STHETA(2)=0d0
      dboost = (PYVAR(31)-PYVAR(32))/
     &            (PYVAR(31)+PYVAR(32))
      DBETA(2,3)=dboost
c      call pylist(1)
C...Transform to parton-parton cms, boost parameters in double precision.
      CALL PYROBO(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
c      write(6,*) ' before event record rearrangement '
c      call pylist(1)
C...Save event record in  cms
      do I=1,N
         do J=1,5
            KS(I,J)=K(I,J)
            PS(I,J)=P(I,J)
         end do
      end do
C...Partons with colour information in hadronic cms frame.
      DO 60  I=20,MSTU(4)
         DO 60 J=1,5
            K(I,J)=0
            P(I,J)=0.D0
   60 V(I,J)=0.d0
      NS=20
c set pt to a small value to avid problems when calling uPDFs      
      if(abs(P(3,1)).le.0.0001) P(3,1)=0.0001
      if(abs(P(3,2)).le.0.0001) P(3,2)=0.0001
      if(abs(P(4,1)).le.0.0001) P(4,1)=0.0001
      if(abs(P(4,2)).le.0.0001) P(4,2)=0.0001
      do  J=1,5
c         N=NS+7
         K(NS+1,J)=K(3,J)
         P(NS+1,J)=P(3,J)
         K(NS+3,J)=K(4,J)
         P(NS+3,J)=P(4,J)
      end do
ccccc      
      nfin = 0
      Do I=5,N
        nfin = nfin + 1
        Nn = NS+4+nfin
        Do J=1,5        
         K(Nn,J)=K(I,J)
         P(Nn,J)=P(I,J)
        End do

      End do

      do i=5,N
         K(I,1) = 13
cc      write(6,*) ' K... ',i,K(i,4),K(i,5)
      enddo

      N=Nn

C...Origin and Colour info for incoming parton
      K(NS+1,1)=13
      K(NS+3,1)=13
c      write(6,*) ' after copy all partons ',N,NS,nfin
c      MSTU(15)=2
c      call pylist(7)
c      call pylist(2)

      do i=21,23,2
c      write(6,*) ' K... ',i,K(i,4),K(i,5)
         ktest4=k(i,4)
         ktest5=k(i,5)
         k(i,4)=0
         k(i,5)=0
         if(mod(ktest4/10000,10000).eq.3) k(i,4)=21*MSTU(5)
         if(mod(ktest5/10000,10000).eq.3) k(i,5)=21*MSTU(5)
         if(mod(ktest4/10000,10000).eq.4) k(i,4)=23*MSTU(5)
         if(mod(ktest5/10000,10000).eq.4) k(i,5)=23*MSTU(5)
         if(mod(ktest4,10000).eq.3) k(i,4)=k(i,4)+21
         if(mod(ktest5,10000).eq.3) k(i,5)=k(i,5)+21
         if(mod(ktest4,10000).eq.4) k(i,4)=k(i,4)+23
         if(mod(ktest5,10000).eq.4) k(i,5)=k(i,5)+23
         if(mod(ktest4/10000,10000).gt.4) k(i,4)=(20+mod(ktest4/10000,10000))*MSTU(5)
         if(mod(ktest5/10000,10000).gt.4) k(i,5)=(20+mod(ktest5/10000,10000))*MSTU(5)
         if(mod(ktest4,10000).gt.4) k(i,4)=k(i,4)+(20+mod(ktest4,10000))
         if(mod(ktest5,10000).gt.4) k(i,5)=k(i,5)+(20+mod(ktest5,10000))
c      write(6,*) ' K... ',i,K(i,4),K(i,5)
      enddo
      do i=25,N
c      write(6,*) ' before K... ',i,K(i,4),K(i,5)
         ktest4 = k(i,4)
         ktest5 = k(i,5)
         k(i,4)=0
         k(i,5)=0
         if(mod(ktest4/10000,10000).eq.3) k(i,4)=21*MSTU(5)
         if(mod(ktest5/10000,10000).eq.3) k(i,5)=21*MSTU(5)
         if(mod(ktest4/10000,10000).eq.4) k(i,4)=23*MSTU(5)
         if(mod(ktest5/10000,10000).eq.4) k(i,5)=23*MSTU(5)
         if(mod(ktest4,10000).eq.3) k(i,4)=k(i,4)+21
         if(mod(ktest5,10000).eq.3) k(i,5)=k(i,5)+21
         if(mod(ktest4,10000).eq.4) k(i,4)=k(i,4)+23
         if(mod(ktest5,10000).eq.4) k(i,5)=k(i,5)+23
         if(mod(ktest4/10000,10000).gt.4) k(i,4)=(20+mod(ktest4/10000,10000))*MSTU(5)
         if(mod(ktest5/10000,10000).gt.4) k(i,5)=(20+mod(ktest5/10000,10000))*MSTU(5)
         if(mod(ktest4,10000).gt.4) k(i,4)=k(i,4)+(20+mod(ktest4,10000))
         if(mod(ktest5,10000).gt.4) k(i,5)=k(i,5)+(20+mod(ktest5,10000))
c      write(6,*) ' after K... ',i,K(i,4),K(i,5)
      enddo
      
c      write(6,*) ' before rearraning '
c      call pylist(2)

        Nold = N
c trick when lepton are between colored particles
        Do I=25,Nold
        if(K(i,1).eq.11.or.(iabs(K(I,2)).ge.11.and.Iabs(K(I,2)).le.16).or.
     &       (Iabs(K(I,2)).ge.22.and.Iabs(K(I,2)).le.25)) then 
c          write(6,*) ' adding parton ',I
          if(iabs(K(I,2)).le.16) then
            do j=1,5
              k(n+1,j) = k(i,j)
              p(n+1,j) = p(i,j)
              k(i,j) = 0
            end do
            n=n+1
          endif 
c fake here an entry with negleglible energy to not screw up pyshow
          do j=1,5
            p(i,j) = 0.0001 
          end do
        endif
        End do
c       write(6,*) ' after rearraning Nold, N, nf1',Nold,N,nf1
c       call pylist(2)
     
      
       IPY(47) = N
c 
c branching type: MSTJ(41) = 1 only QCD branching allowed
      MSTJO41=MSTJ(41)

C...Generate timelike parton shower (if required)
      if(IPY(13).EQ.1) then
         MSTJ(41)=1
c check for heavy flavor production in PS: max flavor for g->qq
         MSTJO45 = MSTJ(45)
cc         MSTJ(45) = 3
        do j=1,4
          Pfin(j)=0.
        end do
c final state particles
        nf1=25
        if(nfin.gt.100) then
          write(6,*) ' cas_lheprocess:  nfin > 100. stop '
          stop
        endif
        nfin1 = 9999
        nfin2 = 0
        do i=nf1,N
          if(iabs(K(I,1)).eq.3.or.iabs(K(I,1)).eq.1) then
            if(i.lt.nfin1) nfin1=i
            if(i.gt.nfin2) nfin2=i
            PTfin(I) = P(I,1)**2 + P(I,2)**2
            Mfin(I) = P(I,5)**2  
c          write(6,*) ' K = ',I, K(I,1), K(I,2) , nfin
          endif
          do J=1,4
            PFIN(J) = PFIN(J) + P(I,J)  
          End do
c          write(6,*) ' adding ',i,(P(i,j),j=1,4)
        End do
        QMAX=dabs(dot(PFIN,PFIN)) 
c        write(6,*) ' before final state shower ',nf1,nfin,Qmax
c        write(6,*) ' final state PS ',IFinal
        IF(IFINAL.EQ.1) THEN
           Qmax = 0
           Do i=nfin1,nfin2
             QMAX = QMAX+ PTfin(I) + Mfin(I)
           end do 
           Qmax = Qmax/nfin
           QMAX = 4d0*QMAX
        ELSEIF(IFINAL.EQ.2) THEN
           QMAX=dot(PFIN,PFIN)       
        ELSEIF(IFINAL.EQ.3) THEN
           Qmax = 0
           Do i=nfin1,nfin2
             QMAX = QMAX+ Mfin(I)
           end do 
           Qmax = Qmax/nfin
           QMAX = 4d0*QMAX
        ELSEIF(IFINAL.EQ.4) THEN
           Qmax = SCALUP**2
           QMAX = 4d0*QMAX
        ENDIF
        
        QMAX = QMAX * SCALFAF
        QMAX=DSQRT(QMAX)
c        write(6,*) ' before final state shower ',nf1,nfin,Qmax
        nfsh1 = nfin1
        nfsh2 = nfin2-nfin1   
c        write(6,*) ' before final state shower ',nf1,nfsh1,nfsh2,nfin1,nfin2
c shower only hadronic parts, no leptons in a color single system 
        CALL PYSHOW(nfsh1,-nfsh2,Qmax)       

c         call pylist(2)
        MSTJ(45) = MSTJO45
        MSTJ(41) = MSTJO41
      endif
      IPY(47) = N

      
      IPU1=21
      IPU2=23
      IPY(41)=K(1,2)
      IPY(42)=K(2,2)
      IRES(1) = 1
      IRES(2) = 1
      ILEPTO =0 
      IPU1O=IPU1
      IPU2O=IPU2
      IF(IPY(14).GE.1) THEN
c         write(6,*) ' before initial state shower '
c         call pylist(1)

         CALL cascps(IPU1,IPU2)
         IF(LST(21).NE.0) THEN
            IF(LST(21).EQ.55.OR.LST(21).EQ.58) THEN
            ELSE
               RETURN
            ENDIF
         ENDIF
      ENDIF

      IPY(47) = N
c      write(6,*) ' before cremn, N = ',N
c      call pylist(1)
c      call pylist(2) 
      
      IRES(1) = 1
      IRES(2) = 1
      ILEPTO =0 
      CALL CAREMN(IPU1,IPU2)
c      write(6,*) ' cas_lheprocess after caremn',IPY(48)
      IF(IPY(48).GE.1) THEN
         LST(21)=47+IPY(48)
c         write(6,*) ' cas_lheread: LST(21) = ',LST(21),IPY(48)
         RETURN
      ENDIF
     
c      write(6,*) ' cas_lheprocess after caremn'
c      call pylist(1)
c rearrange event record
      ii = 0
      Do i=3,4
        ii = ii+1
        DO  J=1,5
         K(N+ii,J)=K(i,J)
         P(N+ii,J)=P(i,J)
        End do
      End do
      Nia1=3
      Nia2=5

	K(21,3)=3
	K(21,1)=21
	K(23,3)=4
	K(23,1)=21

C...Zero irrelevant lines with K(I,1)<0
      do i=5,19
c      K(I,4)=0
      enddo

             
c      write(6,*) ' before pyprep '
c      call pylist(2) 
      MSTU(28)=0
      if (Iremn.eq.1.and.p6frag.eq.1) CALL PYPREP(0)
      if(MSTU(28).ne.0) then
         write(6,*) ' pythia error MSTU(28) ',MSTU(28)
         call pylist(1)
         call pylist(7)
      endif
      if(MSTU(24).ne.0) then
         write(6,*) ' pythia error MSTU(24) ',MSTU(24)
         call pylist(1)
         call pylist(7)
      endif
c      write(6,*) ' after pyprep '
c      call pylist(1)

      DO 190 I=1,4
         P(1,I)=PBEAM(1,I)
         P(2,I)=PBEAM(2,I)
  190 CONTINUE


	if(Iremn.eq.1) then 
c check for enegy momentum conservation
         DO 170 I=1,N
            IF(P(I,5).LT.0.0D0.OR.K(I,1).GT.2) GOTO 170
            DO 160 J=1,4
               IF(ABS(P(I,J)).LE.1.E-6) THEN
                  P(I,J)=0.0D0
               ENDIF
  160       CONTINUE
  170    CONTINUE
         CHEC = 5D-4
c         CHEC = 5D-6
         CHEC = 1D-3
         CHEC = 1D-2
         P1 =ABS(P(1,1)+P(2,1))
         P2 =ABS(P(1,2)+P(2,2))
         P3 =ABS(P(1,3)+P(2,3))
         P4 =ABS(P(1,4)+P(2,4))
         PC1 = (ABS(PYP(0,1)) - P1)/P4
         PC2 = (ABS(PYP(0,2)) - P2)/P4
         PC3 = (ABS(PYP(0,3)) - P3)/P4
         PC4 = (ABS(PYP(0,4)) - P4)/P4
         IF(DABS(PC1).GT.CHEC.
     +      OR.DABS(PC2).GT.CHEC.
     +      OR.DABS(PC3).GT.CHEC.
     +      OR.DABS(PC4).GT.CHEC) THEN
            write(6,*) ' cas_lheread: energy of final particles not correct'
     +      ,chec,LST(21)
c            write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
c            write(6,*) ' x1 = ',pyvar(31),' x2 = ',pyvar(32)
c            write(6,*) ' xg1 = ',xg1,'xg2_bar',xg2_bar
c            write(6,*) ' xg2 = ',xg2,'xg1_bar',xg1_bar
            call pylist(1)
            LST(21) = 100
            return
         ENDIF
      ENDIF
      
c check whether PS was generated 
      if(IPU1O.lt.IPU1.or.IPU2O.lt.IPU2) then
         if((P(IPU1,1).eq.P(IPU1O,1)).and.(P(IPU2,1).eq.P(IPU2O,1))) then
c            write(6,*) ' no PS generated ',IPU1,IPU1O,IPU2,IPU2O
c            call pylist(1)
cccc            LST(21) = 100
         endif
      endif       
c      call pylist(1)

      return
      end
      
