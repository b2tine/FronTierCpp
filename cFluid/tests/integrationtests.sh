#!/bin/bash

rm -rf ReferenceSolutions/test-*

for infile in ReferenceInputs/in-*; do

    infile=$(basename $infile)
    testcase=${infile#in-}
    
    ./runtestcase.sh $testcase
    if [ "$?" -ne 0 ]; then
        exit 1
    fi

    ./testcompare.sh $testcase
    if [ "$?" -ne 0 ]; then
        exit 1
    fi

done

exit 0
