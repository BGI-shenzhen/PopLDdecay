#!/bin/sh
#$ -S /bin/sh
echo Start Time : 
date
chmod 755  ../bin/*
echo g++	-g -O2 LD_Decay.cpp		-lz  -L/usr/lib/ -L./include/zlib/	-o	../bin/PopLDdecay
g++	-g -O2 LD_Decay.cpp		-lz -L/usr/lib64/   -L/usr/lib/	-L./include/zlib/	-o	../bin/PopLDdecay 
echo done see the [ ../bin/PopLDdecay ]
echo End Time : 
date
