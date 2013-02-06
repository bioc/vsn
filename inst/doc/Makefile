all: vsn likelihoodcomputations convergence2 clean

vsn: vsn.tex
	pdflatex vsn.tex

likelihoodcomputations: likelihoodcomputations.tex
	pdflatex likelihoodcomputations.tex

convergence2: convergence2.tex
	cp -p ../scripts/convergence2.pdf .

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	-$(RM) -f vsn-* likelihoodcomputations-* Rplots.pdf



