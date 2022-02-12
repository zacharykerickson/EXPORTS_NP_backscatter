#!/bin/bash

instruments=(FLBBRR FLBBSR BBFL2WW BB2FLSG BB9RR BB9SR HS6RR FLBBLF MCOMSBGC)
len=${#instruments[@]}

echo "Running comparisons with the following variables:"
echo "DIST_THRES = 5 # km"
echo "TIME_THRES = 6/24 # days"
echo "DTEMP_THRES = 0.5 # deg C"
echo "DSAL_THRES = 0.1 # psu"
echo ""

for ((i=0; i<$len; i++)); do
for ((j=i+1; j<$len; j++)); do

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "running: compare_inst.py ${instruments[$j]} ${instruments[$i]}"
python compare_inst.py ${instruments[$j]} ${instruments[$i]}
done
done
