#!/bin/tcsh

foreach i(`seq 0 1 23`)
    echo $i
    ./random_matrix_generator 500 $i > matrix_$i.dat
end
