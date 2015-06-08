
#define the suffixes
.SUFFIXES:	.tex

#macro definitions
SHELL = /bin/sh
LATEX = pdflatex

.tex.aux:	
	$(LATEX) $<


all:	Presentation1.pdf

Presentation1.pdf: Presentation1.tex *.tex
	$(LATEX) Presentation1

force:
	$(LATEX) Presentation1

view:	Presentation1.pdf
	/usr/bin/okular Presentation1.pdf &


clean:
	rm -f *.aux *.nav *.snm *.toc *.log *.dvi *.out *~

