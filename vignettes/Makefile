# vignettes actually being built here
BUILD = 0-vsn 1-vsn 3-likelihoodcomputations

# prebuilt vignettes that are copied from the SRCS_DIR directory
STUBS = 4-convergence
SRCS_DIR = ../inst/vignettes

all: pdf copy

pdf: $(addsuffix .pdf,$(BUILD))

copy: $(STUBS)

%.pdf: %.tex
	pdflatex $*
	pdflatex $*

%: %.tex
	cp -p $(SRCS_DIR)/$*.pdf .

clean:
	$(RM) *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	$(RM) $(addsuffix -*,$(BUILD))
