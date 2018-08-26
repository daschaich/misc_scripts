#!/bin/bash

# Save the current version just in case something goes wrong
cp ~/public_html/KS_nHYP_FA/modeTables.html ~/public_html/KS_nHYP_FA/modeTables.bak

cd /nfs/beowulf03/beowulf02/anna2/KSNHYP/FA-ks481216/
echo "<html><head><title>Errors and Wilson flow measurements</title></head>" > ~/public_html/KS_nHYP_FA/modeTables.html
echo "<body>" >> ~/public_html/KS_nHYP_FA/modeTables.html

date=`date`
echo "<p>Last update: $date</p>" >> ~/public_html/KS_nHYP_FA/modeTables.html

# First pass for table of contents
echo "<ul>" >> ~/public_html/KS_nHYP_FA/modeTables.html
for dir in Run_APBC16_[123]* Run_unsAPBC12_[234]* Run_APBC12_[234]* Run_2HYPAPBC12_[234]* Run_new8_2448 Run_new8_3264 ; do
  echo "<li><a href=#$dir>$dir</a></li>" >> ~/public_html/KS_nHYP_FA/modeTables.html
done
echo "</ul>" >> ~/public_html/KS_nHYP_FA/modeTables.html

for dir in Run_APBC16_[123]* Run_unsAPBC12_[234]* Run_APBC12_[234]* Run_2HYPAPBC12_[234]* Run_new8_2448 Run_new8_3264 ; do
  echo $dir
  cd $dir
  lastUp=`grep Last *_*.html | tail -n 1 | awk '{sub(/<p>/,"(");print $2, $3, $4, $5, $6, $7, $8")"}'`
  echo "<p><b><a name=$dir></a>$dir</b> $lastUp</p>" >> ~/public_html/KS_nHYP_FA/modeTables.html
  echo "<pre>" >> ~/public_html/KS_nHYP_FA/modeTables.html
  python ~/scripts/check_ensembles.py >> ~/public_html/KS_nHYP_FA/modeTables.html
  echo "</pre>" >> ~/public_html/KS_nHYP_FA/modeTables.html
  cd ../
done

echo "</body></html>" >> ~/public_html/KS_nHYP_FA/modeTables.html
