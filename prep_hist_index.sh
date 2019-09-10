echo "  <p>Last update: `date`</p>"
echo "  <ul>"
echo "    <li><a href="poly_hist_m05.gif">poly_hist_m05.gif</a></li>"
echo "    <li><a href="Wpoly_hist_m05.gif">Wpoly_hist_m05.gif</a></li>"
echo "    <br>"
echo "    <li><a href="poly_hist_m1.gif">poly_hist_m1.gif</a></li>"
echo "    <li><a href="Wpoly_hist_m1.gif">Wpoly_hist_m1.gif</a></li>"
echo "    <br>"
echo "    <li><a href="poly_hist_m2.gif">poly_hist_m2.gif</a></li>"
echo "    <li><a href="Wpoly_hist_m2.gif">Wpoly_hist_m2.gif</a></li>"
echo "    <br>"
echo "    <li><a href="poly_hist_m4.gif">poly_hist_m4.gif</a></li>"
echo "    <li><a href="Wpoly_hist_m4.gif">Wpoly_hist_m4.gif</a></li>"
for m in 05 1 2 4 ; do
  for i in poly_hist_b*_m0.$m.pdf ; do
    echo "    <br>"
    echo "    <li><a href="$i">$i</a></li>"
    echo "    <li><a href="W$i">W$i</a></li>"
  done
done
echo "  </ul>"
echo "</body></html>"
