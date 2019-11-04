      Integer Nf_part
      Double precision normhepmc,weight
      Common/lheinfo/normhepmc,weight,Nf_part

      Integer idP1,idP2
      double Precision xidP1,xidP2,pdf1,pdf2,scale_process,scalup_cas,pdf1pdf2,kt2idP1,kt2idP2,pdf_weight
      common/lhepdfinfo/idP1,idP2,xidP1,xidP2,pdf1,pdf2,scale_process,scalup_cas,pdf1pdf2,kt2idP1,kt2idP2,pdf_weight

      character*512 hepmcout,lhefile,Collfile
      Integer Ilhe,iTMDw,Iscale,lheWid
      common/CALHE/lhefile,Ilhe,iTMDw,Collfile,Iscale,lheWid
      Common/CAHEPMCOUT/hepmcout

c replace nr of event with nr of event found in file
	Integer Nevent
	Common/steer1/Nevent
