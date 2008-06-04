all: vsn likelihoodcomputations convergence2 clean

vsn: vsn.tex
	pdflatex vsn
	pdflatex vsn

likelihoodcomputations: likelihoodcomputations.tex
	pdflatex likelihoodcomputations
	pdflatex likelihoodcomputations

convergence2: convergence2.tex
	cp -p ../scripts/convergence2.pdf .

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	rm -f *.dvi *.ps vsn-* likelihoodcomputations-* Rplots.pdf


