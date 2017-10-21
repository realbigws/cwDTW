#!/bin/bash

if [ -f "wletdtw" ] || [ -f "mer2seq" ]
then
	echo " executable files 'wletdtw' and 'mer2seq' already been compiled. "
	exit 1
fi


mkdir -p Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
mv bin/wletdtw ../
mv bin/mer2seq ../
cd ../
rm -rf Release


