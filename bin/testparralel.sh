#/bin/sh

./faststi -f parms/parms_parallel_tests.txt > tmp2.txt &
./faststi -f parms/parms_parallel_tests.txt > tmp3.txt &
wait
echo finished
