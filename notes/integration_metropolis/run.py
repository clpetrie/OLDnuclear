import os

os.system('pdflatex integmet.tex')
os.system('bibtex integmet')
os.system('pdflatex integmet.tex')
os.system('pdflatex integmet.tex')
os.system('open integmet.pdf')
