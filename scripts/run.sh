#!/bin/bash

output=../examples/output
if [ $# -eq 1 ]; then
    output=$1
fi

for input in ../examples/input/*
do
     echo executing sampler-pnt on $input and $output

     echo `../bin/sampler-pnt $input 3 $output`

     #
     # If you have a PowerPC, replace "lpdata50NT.dat with lpdata50.dat"
     #

     echo executing sampler-loop on $input and $output

     echo `../bin/sampler-loop $input 3 ../dat/lpdata50NT.dat $output`
done
