#!/bin/sh

for i in $(seq 1 2)
  do 
    time ./cell_lsodar $i &
    sleep 1
  done
  wait 
