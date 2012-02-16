WEBSITE=modb.oce.ulg.ac.be:/var/lib/mediawiki/upload/Alex/OAK/doc

all: pdf html

bib:
	cat /home/abarth/Biblio/Journals/full.bib /home/abarth/Biblio/*bib > all.bib

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