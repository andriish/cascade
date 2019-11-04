c get dcasrn for bases/spring
      FUNCTION DRN(IRN)
	Implicit None
	Integer IRN
      double precision drn1,DRN
      Double Precision dcasrn
	EXTERNAL dcasrn
      DRN1= dcasrn(IRN)
      DRN = DRN1
	RETURN
	END
	
