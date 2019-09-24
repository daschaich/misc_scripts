echo "  <p>Last update: `date`</p>"
echo "  <table width=\"4000px\">"
echo "    <tr>"
for m in 4 2 1 05 ; do
  echo "      <th width=\"500px\">m=0."$m" Unflowed</th>"
  echo "      <th width=\"500px\">m=0."$m" Wilson-flowed</th>"
done
echo "    </tr>"
echo "    <tr>"
for m in 4 2 1 05 ; do
  echo "      <td width=\"500px\"><a href=\"poly_hist_m"$m".gif\"><img src=\"poly_hist_m"$m".gif\" width=\"450px\"></a></td>"
  echo "      <td width=\"500px\"><a href=\"Wpoly_hist_m"$m".gif\"><img src=\"Wpoly_hist_m"$m".gif\" width=\"450px\"></a></td>"
done
echo "    </tr>"
echo "    <tr>"
# Explicit heights are a bit of a hack to improve alignment...
# TODO!!! Need to correct bash's sorting whenever couplings have different numbers of digits after the decimal point...
for m in 4 2 1 05 ; do
  echo "      <td width=\"500px\">"
  for i in poly_hist_b*_m0.$m.pdf ; do
    echo "        <a href=\""$i"\"><img src=\""${i/pdf/png}"\" width=\"450px\" height=\"300px\"></a><br>"
  done
  echo "      </td>"
  echo "      <td width=\"500px\">"
  for i in poly_hist_b*_m0.$m.pdf ; do
    echo "        <a href=\"W"$i"\"><img src=\"W"${i/pdf/png}"\" width=\"450px\" height=\"300px\"></a><br>"
  done
  echo "      </td>"
done
echo "    </tr>"
echo "  </table>"
echo "</body></html>"
