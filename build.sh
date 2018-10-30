#!/bin/bash


function config_alpha {
    CONF="CFLAGS='-Wformat-overflow=0' CXXFLAGS='-Wformat-overflow=0' --with-mpich-dir=/usr/local/pkg/mpich-3.2.1 --with-petsc-dir=/usr/local/pkg/petsc --with-hdf4-dir=/usr/local/pkg/hdf4-2.10 --with-gd-inc=/usr/include --with-gd-lib=/usr/lib/x86_64-linux-gnu --with-cgal-inc=/usr/include --with-cgal-lib=/usr/lib/x86_64-linux-gnu"
}

function config_defiant {
    CONF="--with-mpich-dir=/usr/local/mpich-3.2.1 --with-petsc-dir=/usr/local/petsc --with-hdf4-dir=/usr/local/hdf4-2.10 --with-gd-inc=/usr/include --with-gd-lib=/usr/lib/x86_64-linux-gnu --with-cgal-inc=/usr/include --with-cgal-lib=/usr/lib/x86_64-linux-gnu"
}


HOST=`uname -n`

if [[ "$HOST" == alpha ]]; then
    echo "Machine recognized as alpha."
    config_alpha
elif [[ "$HOST" == defiant ]]; then
    echo "Machine recognized as defiant."
    config_defiant
fi


autoreconf -vif

echo "configure options: ${CONF}"

./configure ${CONF}

make
