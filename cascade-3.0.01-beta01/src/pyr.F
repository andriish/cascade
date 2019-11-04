      FUNCTION PYR(IDUMMY)

      double precision pyr,dcasrn
	external dcasrn

      PYR= dcasrn(idummy)

      RETURN
      END
