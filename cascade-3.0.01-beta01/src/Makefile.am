## Copyright (C) 2006 Free Software Foundation, Inc.
## This Makefile.am is free software; the Free Software Foundation
## gives unlimited permission to copy, distribute and modify it.

AM_FFLAGS = -I$(srcdir)/../include 
#AM_FFLAGS += -I$(PYTHIA_SOURCE) 
AM_FFLAGS += $(MYFFLAGS)
AM_CPPFLAGS = -I$(srcdir)/../include 
AM_CPPFLAGS+= -I$(HEPMC)/include  
#AM_CPPFLAGS+= -I$(PYTHIA8_DIR)/include
AM_CPPFLAGS+= -I$(BOOST_CFLAGS)
#EXTRA_DIST = $(top_srcdir)/src/pythia-bkdata/pydata.f

AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = cascade
#cascade_SOURCES =  casmain.F  drn.F  upevnt.F upinit.F $(top_srcdir)/src/pythia-bkdata/pydata.f
cascade_SOURCES =  casmain.F  drn.F  upevnt.F upinit.F pydata.f
#cascade_FFLAGS = -I$(top_srcdir)/include -I$(PYTHIA_SOURCE) $(MYFFLAGS) 
#cascade_FFLAGS = -I$(top_srcdir)/include -I$(PYTHIA_SOURCE) $(MYFFLAGS) 
conv_source=convhepmc-helpers.F
if NO_HEPMC
   conv_source+=convhepmc-dummy.cc
else
   conv_source+=convhepmc.cc
endif
if NO_HEPMC3
conv_source +=convhepmc3-dummy.cc 
else
conv_source +=convhepmc3.cc 
endif
if NO_TMDLIB
   tmdlib_source=tmdlib-dummy.F
endif
if NO_RIVET
conv_source += convrivet-dummy.cc 
else
conv_source += convrivet.cc 
endif
if NO_PYTHIA8
   pythia8_source=pythia8frag-dummy.cc
else
   pythia8_source=pythia8frag.cc
endif

lib_LTLIBRARIES = libcascade2.la
libcascade2_la_SOURCES = aalam.F    caps.F      cascha.F    caupevnt.F  dfun.F   \
                   event.F    mscatt.F     steer.F   xsecpsi.F   \
                   alphas.F   caps_ms.F   cascps.F    caupinit.F  dgrv_nln.F \
                   fxn1.F     multscatt.F  pyr.F ccfmstfq.F   \
                   bcasin.F   casgen.F    ccfmstfg.F  dot1.F  \
                   fxnb.F     parti.F      pytcha.F  \
                   blueml.F   caremh.F    casini.F    ccfmstfu.F  dot.F  \
                   ktstfu.F   p_gamgam.F   xglu.F \
                   caend.F    caremn.F     cutg.F     \
                   mehiggs.F  phase.F      xgx.F \
                   caprikt.F  casbran.F   caspli.F    dcasrn.F    d_xgx.F \
                   meoffsh.F  p_mscatt.F   ranums.F  xhiggs.F \
                   capsee.F   cascade.F   cauniglu.F  dfridr.F    elemnz.F \
                   minuts.F   p_semih.F    ryssha.F  xsec1.F xsecqg.F  \
                   meoffjpsi.F  dhfint.F   \
                   colorflow.F xvq.F mevq.F meqg.F  \
                   mcint.F mcgen.F \
                   xvqq.F mevqq.F ktOrtAmVqNew.F \
                   xsecchi.F meoffchi.F mev.F xv1.F \
                   xsecjj.F meoffjj.F jj_ggf.F jj_me2.F \
                   input_steer.F steer_new.F \
                   cas_lheread.F $(conv_source) $(tmdlib_source) $(pythia8_source) casbran_new.F
#                   meoffjpsi.F pdfset.F   structm.F  structp.F 


cascade_LDADD =  ./libcascade2.la ../bases51/libbases.la \
                $(HEPMCLIB) $(HEPMC3LIB) $(RIVETLIB) $(PYTHIALIB) $(TMDLIB) $(PDFLIB)  $(GSL_LIB) $(PYTHIA8LIB) ../misc/libmycern.la -lstdc++ -lm
                
	
