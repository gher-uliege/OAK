BIBDIR=$(HOME)/Biblio
SRC=$(wildcard *.tex)
NAME=$(patsubst %.tex,%,$(SRC))
WEBSITE=modb.oce.ulg.ac.be:/var/lib/mediawiki/upload/Alex/OAK/doc

all: pdf html all.bib

all.bib: $(BIBDIR)/Journals/full.bib $(wildcard $(BIBDIR)/*.bib) 
	cat $(BIBDIR)/Journals/full.bib $(BIBDIR)/*bib > all.bib

pdf: 
	pdflatex oak.tex	
	bibtex oak
	pdflatex oak.tex	
	pdflatex oak.tex	

html:
	htlatex oak.tex "oak,xhtml,jsmath" " -cmozhtf"

clean:
	rm -f oak.4ct oak.log oak.xref oak.4tc oak.idv oak.aux oak.dvi oak.lg oak.tmp oak.out

upload:
	scp oak.pdf oak*.html oak.css $(WEBSITE)