#!/bin/bash

INPUTDIR="input-brandon-FWE-v12"

nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-pointmass-VREMAN-v12 -o out-C9-nVnG-pointmass-VREMAN-v12 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-ball-VREMAN-v12 -o out-C9-nVnG-ball-VREMAN-v12 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-box-VREMAN-v12 -o out-C9-nVnG-box-VREMAN-v12 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-human-VREMAN-v12 -o out-C9-nVnG-human-VREMAN-v12 &

nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedpointmass-VREMAN-v12 -o out-C9-nVnG-fixedpointmass-VREMAN-v12 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedball-VREMAN-v12 -o out-C9-nVnG-fixedball-VREMAN-v12 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedbox-VREMAN-v12 -o out-C9-nVnG-fixedbox-VREMAN-v12 &
nohup ./parachute -d 3 -i ${INPUTDIR}/in-C9-nVnG-fixedhuman-VREMAN-v12 -o out-C9-nVnG-fixedhuman-VREMAN-v12 &

