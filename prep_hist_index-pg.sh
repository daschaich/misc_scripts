echo "  <p>Last update: `date`</p>"
echo "  <ul>"
echo "    <li><a href="poly_hist.gif">poly_hist.gif</a></li>"
echo "    <li><a href="Wpoly_hist.gif">Wpoly_hist.gif</a></li>"
for i in poly_hist_b*.pdf ; do
  echo "    <br>"
  echo "    <li><a href="$i">$i</a></li>"
  echo "    <li><a href="W$i">W$i</a></li>"
done
echo "  </ul>"
echo "</body></html>"
