#!/bin/bash

if [ -f "cwDTW" ] && [ -f "seq2sig" ]
then
	echo " executable files 'cwDTW' and 'seq2sig' already been compiled. "
	exit 1
fi


mkdir -p Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
mv bin/wletdtw ../cwDTW
mv bin/mer2seq ../seq2sig
cd ../
rm -rf Release


