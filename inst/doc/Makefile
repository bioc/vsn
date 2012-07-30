UNAME := $(shell uname)
all: vsn likelihoodcomputations convergence2 clean

vsn: vsn.tex
ifeq ("$(UNAME)", "CYGWIN_NT-6.1-WOW64")
		cmd /C texi2dvi --pdf vsn.tex
else
		texi2dvi --pdf vsn.tex
endif

likelihoodcomputations: likelihoodcomputations.tex
ifeq ("$(UNAME)", "CYGWIN_NT-6.1-WOW64")
		cmd /C texi2dvi --pdf likelihoodcomputations.tex
else
		texi2dvi --pdf likelihoodcomputations.tex
endif

convergence2: convergence2.tex
		cp -p ../scripts/convergence2.pdf .

clean:
		-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
		-$(RM) -f vsn-* likelihoodcomputations-* Rplots.pdf