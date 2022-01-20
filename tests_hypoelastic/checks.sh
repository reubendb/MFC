#!/bin/bash

### Note: Designed to be run from MFC root directory

rm -rf ./tests_hypoelastic/*/D ./tests_hypoelastic/*/*.inp ./tests_hypoelastic/*/p_all ./tests_hypoelastic/*/*.out

mytests=( 1d_sod 1d_impact 1d_2materials 2d_5wave )

ntest=${#mytests[@]}

npass=0
nfail=0
i=1

echo -----------------------------------------------

cd tests_hypoelastic
for mytest in "${mytests[@]}"; do
    cd $mytest
        
        #Run test case
        ./input.py pre_process > pre_process.out
        ./input.py  simulation > simulation.out
        
        cd check
            check_file=$(echo *)
        cd ..
        
        #Check that the files are the same
        rm -f diff.out
        diff check/$check_file D/$check_file > diff.out

        mytest="Test $i of $ntest: $mytest"
        #Print if not
        if [ -s diff.out ]; then
            echo $mytest: Test failed! Output files are different.
            ((++nfail))
        elif [ ! -f "D/$check_file" ]; then
            echo $mytest: Test failed! Output file was not found.
            ((++nfail))
        else
            echo $mytest: Test passed!
            ((++npass))
        fi

        ((++i))

    cd ..
done

echo -----------------------------------------------
echo ---- $nfail Tests failed
echo ---- $npass Tests passed
