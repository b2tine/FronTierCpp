#!/bin/bash
 
for infile in ReferenceInputs/*;
do
    testcase=$(basename $infile)
    echo $testcase
    testcase=${testcase#in-}
    echo $testcase
    exit 0

    echo "Testing: $(basename $testcase)"
    ../iFluid -d 2 -i $testcase -o "ReferenceSolutions/out-$(basename $testcase)"

    echo "Comparing Test Simulation to Correct"
    echo "Comparing CorrectSimulations/ReferenceOutput/"$(basename $test)-tests"/ outputs to CorrectSimulations/ReferenceOutput/"$(basename $test)"/ outputs"
    for output in CorrectSimulations/ReferenceOutput/"$(basename $test)"/*;
    do
        if [[ $output =~ "ts" ]]
        then
            echo "Comparing $output to CorrectSimulations/ReferenceOutput/"$(basename $test)-tests"/"$(basename $output)""
            diff $output CorrectSimulations/ReferenceOutput/"$(basename $test)-tests"/"$(basename $output)" 
            if [ $? -ne 0 ]; then
                echo "$(basename $test) test FAILED, generated files NOT identical"
                exit -1
            fi
            echo "$($output) PASSED"
        fi
    done
    rm -rf CorrectSimulations/ReferenceOutput/"$(basename $test)-tests"
    echo "$(basename $test) test PASSED, generated files ARE identical"
    exit 0 
done

