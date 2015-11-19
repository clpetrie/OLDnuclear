#!/bin/bash

#get filename
filename="$1"
if [ "$filename" = "" ]
then
 filename="data.dat"
fi

grep 'total energy-propagated' $filename > mytemp.dat

gnuplot plot.gp -p

rm mytemp.dat
