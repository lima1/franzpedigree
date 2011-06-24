#!/bin/sh
cp configure.ac.WIN32 configure.ac
./configure --host=i386-mingw32 
make
make clean
make
man -t man/FRANz.1 | ps2pdf - > man/FRANz.pdf
cp extras/installer/win/FRANz.nsi .
cp extras/installer/win/EnvVarUpdate.nsh .
makensis FRANz.nsi
rm FRANz.nsi
rm EnvVarUpdate.nsh
cp configure.ac.unix configure.ac
