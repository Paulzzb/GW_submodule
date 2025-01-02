name=kssolvug

all: $(name).pdf

$(name).pdf: $(name).tex
	pdflatex $(name)
	bibtex   $(name).aux
	pdflatex $(name)
	pdflatex $(name)

clean:
	rm -f $~ *.aux *.idx *.log *.blg *.bbl *.toc *.dvi *.ps *.out *.gz \
			 *.fdb_latexmk *.fls
cleanless:
	rm -f $~ *.aux *.idx *.log *.blg *.toc *.dvi *.ps *.out *.gz \
			 *.fdb_latexmk *.fls
cleanall:
	rm -f $~ *.aux *.idx *.log *.blg *.bbl *.toc *.dvi *.ps *.out *.gz \
			 *.fdb_latexmk *.fls \
	         figure/*-eps-converted-to.pdf $(name).pdf
