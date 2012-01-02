#!/bin/sh

X0=$1
Y0=$2
X=$3
Y=$4

#echo "($X - $X0)^2 + ($Y - $Y0)^2" | bc
echo "($X - $X0)^2 + ($Y - $Y0)^2" | gp -f -q
