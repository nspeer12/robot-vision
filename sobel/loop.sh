#!/bin/bash

for i in {1..300}
do
    fout=$"output/output_${i}.pgm"
    thresh=$"${i}"
    ./a.out garb34.pgm fout thresh
    
    #echo "output/output_${i}.pgm" "${i}"
done

