#!/bin/bash

# Save the current version just in case something goes wrong
cp ~/public_html/KS_nHYP_FA/12fTables.html ~/public_html/KS_nHYP_FA/12fTables.bak

cd /nfs/beowulf03/beowulf02/anna2/KSNHYP/FA-ks481216/
echo "<html><head><title>Errors and Wilson flow measurements</title></head>" > ~/public_html/KS_nHYP_FA/12fTables.html
echo "<body>" >> ~/public_html/KS_nHYP_FA/12fTables.html

date=`date`
echo "<p>Last update: $date</p>" >> ~/public_html/KS_nHYP_FA/12fTables.html

# First pass for table of contents
echo "<ul>" >> ~/public_html/KS_nHYP_FA/12fTables.html
for dir in Run_APBC12_[1234]* ; do
  echo "<li><a href=#$dir>$dir</a></li>" >> ~/public_html/KS_nHYP_FA/12fTables.html
done
echo "</ul>" >> ~/public_html/KS_nHYP_FA/12fTables.html

for dir in Run_APBC12_[1234]* ; do
  echo $dir
  cd $dir
  latest=`ls -tr *_*.html | tail -n 1`
  lastUp=`grep Last $latest | awk '{sub(/<p>/,"("); sub(/<\/p>/,")"); print}'`
  echo "<p><b><a name=$dir></a>$dir</b> $lastUp</p>" >> ~/public_html/KS_nHYP_FA/12fTables.html
  echo "<pre>" >> ~/public_html/KS_nHYP_FA/12fTables.html
  python ~/scripts/check_ensembles.py >> ~/public_html/KS_nHYP_FA/12fTables.html
  echo "</pre>" >> ~/public_html/KS_nHYP_FA/12fTables.html
  cd ../
done

echo "</body></html>" >> ~/public_html/KS_nHYP_FA/12fTables.html
