#!/bin/bash
export OMP_NUM_THREADS=48
./gapmis -a ./data/100/queries100K.fa -b ./data/100/targets100.fa -t 48
./gapmis -a ./data/150/queries100K.fa -b ./data/150/targets100.fa -t 48
./gapmis -a ./data/200/queries100K.fa -b ./data/200/targets100.fa -t 48

./gapmis -a ./data/100/queries100K.fa -b ./data/100/targets100.fa -t 32
./gapmis -a ./data/150/queries100K.fa -b ./data/150/targets100.fa -t 32
./gapmis -a ./data/200/queries100K.fa -b ./data/200/targets100.fa -t 32

./gapmis -a ./data/100/queries100K.fa -b ./data/100/targets100.fa -t 16
./gapmis -a ./data/150/queries100K.fa -b ./data/150/targets100.fa -t 16
./gapmis -a ./data/200/queries100K.fa -b ./data/200/targets100.fa -t 16

./gapmis -a ./data/100/queries100K.fa -b ./data/100/targets100.fa -t 8
./gapmis -a ./data/150/queries100K.fa -b ./data/150/targets100.fa -t 8
./gapmis -a ./data/200/queries100K.fa -b ./data/200/targets100.fa -t 8

./gapmis -a ./data/100/queries100K.fa -b ./data/100/targets100.fa -t 4
./gapmis -a ./data/150/queries100K.fa -b ./data/150/targets100.fa -t 4
./gapmis -a ./data/200/queries100K.fa -b ./data/200/targets100.fa -t 4

./gapmis -a ./data/100/queries100K.fa -b ./data/100/targets100.fa -t 1
./gapmis -a ./data/150/queries100K.fa -b ./data/150/targets100.fa -t 1
./gapmis -a ./data/200/queries100K.fa -b ./data/200/targets100.fa -t 1
