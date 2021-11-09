#!/bin/bash

if [ $# -lt 1 ]; then
  echo "Usage: $0 <dirs>"
  exit 1
fi

for var in "$@" ; do
  list=`ls -d ~/4+6f/$var`
  for i in $list ; do
    if [ -d $i ] ; then # Ignore non-directory arguments
      echo $i
      echo -ne "Conf: "
      ls $i/Configs/ckpoint_lat.* | wc -l
      echo -ne "Out:  "
      ls $i/Out/pipi.* | wc -l
      echo -ne "Done: "
      grep "Total measurement time: " $i/Out/pipi.* | wc -l
      echo -ne "HDF5: "
      for h in $i/pipi/wall_llll.* ; do
        h5dump $h | wc -l | grep -w 808
      done | wc -l
      echo
    fi
  done
done
