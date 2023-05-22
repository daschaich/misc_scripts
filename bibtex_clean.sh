#!/bin/bash

if [ $# -lt 1 ]; then
  echo "Usage: $0 <file>"
  exit 1
fi

file=$1

sed -i -E "s/author =/author        =/" $file
sed -i -E "s/title =/title         =/" $file
sed -i -E "s/journal =/journal       =/" $file
sed -i -E "s/volume =/volume        =/" $file
sed -i -E "s/year =/year          =/" $file
sed -i -E "s/pages =/pages         =/" $file
sed -i -E "s/doi =/doi           =/" $file
sed -i -E "s/eprint =/eprint        =/" $file
