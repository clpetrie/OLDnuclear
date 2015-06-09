import os

os.system('pdflatex vars.tex')
os.system('pdflatex subr.tex')
os.system('pdftk vars.pdf subr.pdf cat output program.pdf')
os.system('gnome-open program.pdf')
