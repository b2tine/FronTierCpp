#!/bin/sh

autoreconf -vif
./configure --with-hdf4-dir=/usr/local/hdf4-2.10
make
