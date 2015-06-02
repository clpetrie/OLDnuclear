import os

#run latex and bibtex on my thesis
os.system('pdflatex quantummc.tex')
os.system('bibtex quantummc')
os.system('pdflatex quantummc.tex')
os.system('pdflatex quantummc.tex')
os.system('open quantummc.pdf')
