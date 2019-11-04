      subroutine input_steer

      implicit none
      Integer Isteer_new
      character *512 TXT 
	Logical ex
      Integer Irr	
      External Pydata
c initialize default parameters      
      call casini
      call dummysub(pydata)
      Isteer_new = 0 
c check if old or new steering file is used      
   20 Read(5,101,END=100) TXT
      If(TXT(1:1).EQ.'*') then
c         WRITE(6,103) '*',TXT
         GOTO 20
      elseif(TXT(1:1).EQ.'&') then
         Isteer_new = 1
      Endif
      Rewind (5)
      If(Isteer_new.eq.0 ) then
         call steer
         call cascha
         call pytcha
	   inquire(FILE='pythia64-update-particle.dat',EXIST=ex)
         if(ex) then
            open(3,FILE='pythia64-update-particle.dat', FORM='formatted',
     +      STATUS='OLD',IOSTAT=IRR,ERR=220)
            write(6,*) ' updating PYTHIA particle data with file ',
     +      'pythia64-update-particle.dat'
	      call pyupda(3,3)
         else 
            write(6,*) '  PYTHIA particle data file not existing ' 
            open(3,FILE='pythia64-update-particle.dat',
     +      FORM='formatted',STATUS='NEW',
     +      IOSTAT=IRR,ERR=220)
	      call pyupda(1,3)
            write(6,*) ' pythia64-update-particle.dat written '
            write(6,*) ' NOTE: since no update of file ',
     +      'no 2 or 3S resonances can be produced'
	   Endif
         close(3)
      elseif(Isteer_new.eq.1) then
         call steer_new
      else
         write(6,*) ' no valid steering file found; Isteer_new = ',Isteer_new
      endif
      
      
  100 Return

  101 Format(A512)

      Return
      
  220 write(6,*) ' error opening  file for PYTHIA update '
      STOP
      END
      subroutine dummysub(a)
      return
      end

