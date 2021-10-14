#!/bin/bash

cd ../
homedir=$(pwd)
cd lib

library=lapack-3.4.2
echo 'opening tarball'
tar -zxf $library'.tgz'
echo 'finished opening tarball'
cd $library
    #cp make.inc.example make.inc
    cat ../make_lapack.inc | sed -e s+XXMAKEDIR+$homedir+g > make.inc
    make blaslib 
    make lapacklib 
    make install
cd ..
