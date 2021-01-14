#!/bin/bash

./parachute -d 3 -i input-C9-VREMAN/in-C9-nVnG-ball-VREMAN-v15 -o out-nvng-ball-slip-v15 &
./parachute -d 3 -i input-C9-VREMAN/in-C9-nVnG-ball-VREMAN-noslip-v15 -o out-nvng-ball-noslip-v15 &

./parachute -d 3 -i input-C9-VREMAN/in-C9-nVnG-ball-VREMAN-v50 -o out-nvng-ball-slip-v50 &
./parachute -d 3 -i input-C9-VREMAN/in-C9-nVnG-ball-VREMAN-noslip-v50 -o out-nvng-ball-noslip-v50 &
