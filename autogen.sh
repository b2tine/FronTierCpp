#!/bin/sh

export LD_LIBRARY_PATH="/usr/local/hdf4-2.10/lib"

autoreconf --verbose --install --force
./configure --with-hdf4-dir=/usr/local/hdf4-2.10
make clean
make
