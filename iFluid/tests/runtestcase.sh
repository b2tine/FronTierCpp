#!/bin/bash

casename=$1

if [ -d "ReferenceSolutions/test-$casename" ]; then
    rm -rf "ReferenceSolutions/test-$casename"
fi

echo "Running testcase: $casename"

../iFluid -d 2 -i ReferenceInputs/in-$casename -o ReferenceSolutions/test-$casename

if [ "$?" -ne 0 ]; then
    echo "RUNTIME ERROR: $casename"
    exit 1
fi

exit 0
