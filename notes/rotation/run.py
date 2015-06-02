import os

#run latex and bibtex on my thesis
os.system('pdflatex rotation.tex')
os.system('bibtex rotation')
os.system('pdflatex rotation.tex')
os.system('pdflatex rotation.tex')
os.system('open rotation.pdf')
