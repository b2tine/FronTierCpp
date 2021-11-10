#!/bin/bash

INPUTDIR="input-brandon-FWE-v50"

nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-pointmass-VREMAN-v50 -o out-C9-nVnG-pointmass-VREMAN-v50 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-ball-VREMAN-v50 -o out-C9-nVnG-ball-VREMAN-v50 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-box-VREMAN-v50 -o out-C9-nVnG-box-VREMAN-v50 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-human-VREMAN-v50 -o out-C9-nVnG-human-VREMAN-v50 &

nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedpointmass-VREMAN-v50 -o out-C9-nVnG-fixedpointmass-VREMAN-v50 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedball-VREMAN-v50 -o out-C9-nVnG-fixedball-VREMAN-v50 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedbox-VREMAN-v50 -o out-C9-nVnG-fixedbox-VREMAN-v50 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedhuman-VREMAN-v50 -o out-C9-nVnG-fixedhuman-VREMAN-v50 &

