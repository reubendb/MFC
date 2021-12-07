#!/bin/bash

# Color ANSI Escape Sequences
FG_RED='\033[0;31m'
FG_GREEN='\033[0;32m'
FG_ORANGE='\033[0;33m'
FG_NONE='\033[0m'

### Note: Designed to be run from MFC root directory, not tests/

rm -rf ./tests/*/D ./tests/*/*.inp ./tests/*/p_all ./tests/*/*.out

mytests=( 1d_bubscreen kapila_shocktube sod_shocktube vacuum_generation )

ntest=${#mytests[@]}

npass=0
nfail=0
i=1

echo -----------------------------------------------

cd tests
for mytest in "${mytests[@]}"; do
    cd $mytest
        
        #Run test case
        /usr/bin/env python3 ./input.py pre_process > pre_process.out
        /usr/bin/env python3 ./input.py  simulation > simulation.out
        
        cd check
            check_file=$(echo *)
        cd ..
        
        #Check that the files are the same
        rm -f diff.out
        diff check/$check_file D/$check_file > diff.out

        mytest="Test $i of $ntest: $mytest"
        #Print if not
        if [ -s diff.txt ]; then
            echo -e $mytest": "$FG_RED"Test failed! Output files are different."$FG_NONE
            ((++nfail))
        elif [ ! -f "D/$check_file" ]; then
            echo -e $mytest": "$FG_RED"Test failed! Output file was not found."$FG_NONE
            ((++nfail))
        else
            echo -e $mytest": "$FG_GREEN"Test passed!"$FG_NONE
            ((++npass))
        fi

        ((++i))

    cd ..
done

echo -----------------------------------------------
echo -e ---- $FG_RED$nfail Tests failed$FG_NONE
echo -e ---- $FG_GREEN$npass Tests passed$FG_NONE

# Proper Exit Code
#      0 - Success
# Not(0) - Fail (1 or more tests)
exit $nfail
