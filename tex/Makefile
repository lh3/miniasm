.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:miniasm.pdf

miniasm.pdf:miniasm.tex miniasm.bib
		pdflatex miniasm; bibtex miniasm; pdflatex miniasm; pdflatex miniasm;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ miniasm.pdf
