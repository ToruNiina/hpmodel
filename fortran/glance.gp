set xrange[-18:18]
set yrange[-18:18]

do for [idx=0:9999] {
  plot 'trajectory.dat' i idx u 1:2 w l
}
