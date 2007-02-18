all: vsn incremental convergence2 clean

vsn: vsn.tex
	pdflatex vsn
	pdflatex vsn

incremental: incremental.tex
	pdflatex incremental	
	pdflatex incremental

convergence2: convergence2.tex
	cp -p ../scripts/convergence2.pdf .

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	rm -f *.dvi *.ps vsn-* incremental-* 


