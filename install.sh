#!/bin/bash

if [ -f "cwDTW" ] && [ -f "seq2sig" ] && [ -f "sig2peak" ]
then
	echo " executable files 'cwDTW', 'seq2sig' and 'sig2peak' already been compiled. "
	exit 1
fi


mkdir -p Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
mv bin/wletdtw ../cwDTW
mv bin/mer2seq ../seq2sig
mv bin/sig2peak ../sig2peak
cd ../
rm -rf Release


