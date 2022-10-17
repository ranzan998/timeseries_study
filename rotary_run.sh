#!/bin/bash

f77 -c -I/usr/include wavetest_ureal.f
f77 -o test wavetest_ureal.o -L/usr/lib/x86_64-linux-gnu -lnetcdff
./test

f77 -c -I/usr/include wavetest_uimag.f
f77 -o test wavetest_uimag.o -L/usr/lib/x86_64-linux-gnu -lnetcdff
./test


f77 -c -I/usr/include wavetest_vreal.f
f77 -o test wavetest_vreal.o -L/usr/lib/x86_64-linux-gnu -lnetcdff
./test

f77 -c -I/usr/include wavetest_vimag.f
f77 -o test wavetest_vimag.o -L/usr/lib/x86_64-linux-gnu -lnetcdff
./test

