#!/bin/bash

# Save the current version just in case something goes wrong
cp ~/public_html/Tables.html ~/public_html/Tables.bak

cd /nfs/beowulf02/anna2/KSNHYP/FA-ks481216/
echo "<html><head><title>Errors and Wilson flow measurements</title></head>" > ~/public_html/Tables.html
echo "<body>" >> ~/public_html/Tables.html

date=`date`
echo "<p>Last update: $date</p>" >> ~/public_html/Tables.html

# First pass for table of contents
echo "<ul>" >> ~/public_html/Tables.html
for dir in Run_APBC8_[123]* Run_2HYPAPBC8_[123]* ; do
  echo "<li><a href=#$dir>$dir</a></li>" >> ~/public_html/Tables.html
done
echo "</ul>" >> ~/public_html/Tables.html

for dir in Run_APBC8_[123]* Run_2HYPAPBC8_[123]* ; do
  echo $dir
  cd $dir
  lastUp=`grep Last *_*.html | tail -n 1 | awk '{sub(/<p>/,"(");print $2, $3, $4, $5, $6, $7, $8")"}'`
  echo "<p><b><a name=$dir></a>$dir</b> $lastUp</p>" >> ~/public_html/Tables.html
  echo "<pre>" >> ~/public_html/Tables.html
  python ~/scripts/check_ensembles.py >> ~/public_html/Tables.html
  echo "</pre>" >> ~/public_html/Tables.html
  cd ../
done

echo "</body></html>" >> ~/public_html/Tables.html
