#!/bin/bash
reset
file=ising
gcc -O3 -W $file.c -o $file.e -lm
./$file.e 32 1000
#ipython3 -i problema2b.py
#ipython3 -i problema2d.py
ipython3 -i problema2e.py
#ipython3 -i correlacion2b.py
