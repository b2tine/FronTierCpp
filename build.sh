#!/bin/bash

function config_defiant {
    CONF="--with-mpich-dir=/usr/local/mpich-3.2.1 --with-petsc-dir=/usr/local/petsc --with-hdf4-dir=/usr/local/hdf4-2.10 --with-gd-inc=/usr/include --with-gd-lib=/usr/lib/x86_64-linux-gnu --with-cgal-inc=/usr/include --with-cgal-lib=/usr/lib/x86_64-linux-gnu"
}


HOST=`uname -n`

if [[ "$HOST" == defiant ]]; then
    echo "Machine is recognized as defiant."
    config_defiant
fi


autoreconf -vif

echo "configure options: ${CONF}"

./configure ${CONF}

make
