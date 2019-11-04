#if BUILD_DOCS
pkgdocdir = ${mandir}
#dist_
pkgdoc_DATA = cascade.pdf cascade-lhe.pdf
EXTRA_DIST = cascade.tex cascade.bbl mystyle-new.sty ladder-ccfm.eps ref.bib mcite.sty heralhc.bst \
	parton-shower-evolution-crop.pdf splitting-process-latex-crop.pdf ladder-ccfm.pdf cascade-lhe.bbl cascade-lhe.tex 
CLEANFILES =cascade.pdf cascade-lhe.pdf
DISTCLEANFILES = cascade.log cascade.idx cascade.aux cascade.dvi \
                 cascade.ilg cascade.blg cascade.ind \
                 cascade-lhe.aux cascade-lhe.idx cascade-lhe.log \
                 copy-doc 

all-local: copy-doc cascade.pdf

copy-doc:
	find ${srcdir} -type f -exec ln -s {} 2> /dev/null \;
	if test -h cascade.ind; then rm -f cascade.ind; cp ${srcdir}/cascade.ind . && chmod a+w cascade.ind; fi
	touch copy-doc

#cascade.ps: copy-doc
#	latex cascade && bibtex cascade && latex cascade && makeindex cascade && latex cascade && dvips -o cascade.ps cascade

cascade.pdf: copy-doc
	pdflatex cascade  && pdflatex cascade && makeindex cascade && pdflatex cascade 

cascade-lhe.pdf: copy-doc
	pdflatex cascade-lhe  && pdflatex cascade-lhe &&  pdflatex cascade-lhe 

#endif
