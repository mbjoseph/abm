# Makefile for the scientific article in project-template

all: output/article_text.docx output/article_text.pdf output/full_article.pdf

output/article_text.docx: article_text.md library.bib
	pandoc -H format.sty -V fontsize=12pt --bibliography library.bib --csl=ecology.csl --latex-engine=xelatex article_text.md -o output/article_text.docx

output/article_text.pdf: article_text.md library.bib
	# reset the starting page to 1 
	./rwstyle 1 > format.sty

	pandoc -H format.sty -V fontsize=12pt --bibliography library.bib --csl=ecology.csl --latex-engine=xelatex article_text.md -o output/article_text.pdf

output/full_article.pdf: article_text.md library.bib
	# merge pdfs
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=output/full_article.pdf output/article_text.pdf
	
	# reset the starting page to 1 
	./rwstyle 1 > format.sty
