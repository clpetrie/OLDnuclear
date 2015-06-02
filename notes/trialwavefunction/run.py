import os

#run latex and bibtex on my thesis
os.system('pdflatex trial.tex')
os.system('bibtex trial')
os.system('pdflatex trial.tex')
os.system('pdflatex trial.tex')
os.system('open trial.pdf')
