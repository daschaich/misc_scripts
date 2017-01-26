#!/bin/bash
# ------------------------------------------------------------------
# Estimate total amount of time the given Wilson flow measurement will take
# based on the progress it makes in a minute
if [ $# != 1 ]; then
  echo "Usage: $0 <file>"
  exit 1
fi

file=$1

START=`tail -n 1 $file | awk '{print $2}'`
sleep 60
tail -n 1 $file | awk -v s=$START '{print s, "-->", $2, "==>", 128.0 / (60.0 * ($2 - s)), "hours"}'
# ------------------------------------------------------------------
