#!/bin/bash

casename=$1

testdir="ReferenceSolutions/test-$casename"
refdir="ReferenceSolutions/refsol-$casename"

for refsol in $refdir/{intfc,state}*; do
    filename=$(basename $refsol)
    checksol="$testdir/$filename"
    if [ -e "$checksol" ]; then
        if ! diff $checksol $refsol; then
            echo "FAIL: $checksol NOT EQUAL TO $refsol"
            exit 1
        fi
    else
        echo "FAIL: $checksol NOT FOUND"
        exit 1
    fi
done

echo "PASS: $casename"
exit 0
