#!/bin/bash


./iFluid -d 3 -i in-box3d-noturb -o out-ball3d-VREMAN &
./iFluid -d 3 -i in-box3d-VREMAN -o out-ball3d-VREMAN &
./iFluid -d 3 -i in-ball3d-CGAL-noturb -o out-ball3d-CGAL-noturb &
./iFluid -d 3 -i in-ball3d-CGAL-VREMAN -o out-ball3d-CGAL-VREMAN &
./iFluid -d 3 -i in-cyl3d-CGAL-noturb -o out-cyl3d-CGAL-noturb &
./iFluid -d 3 -i in-cyl3d-CGAL-VREMAN -o out-cyl3d-CGAL-VREMAN &
