#!/bin/bash

testcase=$1

testdir="test-$testcase"
refdir="refsol-$testcase"

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

echo "PASS: $testcase"
exit 0
