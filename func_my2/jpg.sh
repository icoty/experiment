#!/bin/bash

list=`ls *json`

for i in $list
do
  jpg=`echo $i | cut -d '.' -f 1`".jpg"
  if [ -f $jpg ]; then
    echo $jpg" exist."
    continue
  fi
  python3 graph.py $i
done
