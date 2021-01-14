#!/bin/bash

./iFluid -d 2 -i in-cyl2d-VREMAN -o out-cyl2d-slip &
./iFluid -d 2 -i in-cyl2d-VREMAN-noslip -o out-cyl2d-noslip &
