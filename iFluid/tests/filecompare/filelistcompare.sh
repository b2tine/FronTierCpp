#!/bin/bash

testcase="rt2d"

testdir="test-$testcase"
refdir="refsol-$testcase"

for refsol in $refdir/{intfc,state}*; do
    filename=$(basename $refsol)
    checksol="$testdir/$filename"
    if [ -e "$checksol" ]; then
        if ! diff $checksol $refsol; then
            echo "FAIL: $checksol NOT equal to $refsol"
            exit 1
        fi
    else
        echo "FAIL: $checksol not found"
        exit 1
    fi
done

echo "PASS: $testcase"
exit 0
