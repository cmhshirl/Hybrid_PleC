#!/bin/sh

for i in $(seq 1 30)
  do 
    time ./Stochastic $i &
    sleep 1
  done
wait

