reset
set terminal png size 1020,1020
set encoding koi8u
set encoding utf8

unset xtics
unset xlabel
unset ytics
unset ylabel

set xrange[0:50]
set yrange[0:50]

f91(n) = sprintf("DataA_%d_1.dat",n)
f92(n) = sprintf("DataA_%d_2.dat",n)
f93(n) = sprintf("DataA_%d_3.dat",n)
f94(n) = sprintf("DataA_%d_4.dat",n)
f95(n) = sprintf("DataA_%d_5.dat",n)
f96(n) = sprintf("DataA_%d_6.dat",n)
f97(n) = sprintf("DataA_%d_7.dat",n)
f98(n) = sprintf("DataA_%d_8.dat",n)
f99(n) = sprintf("DataA_%d_9.dat",n)
f90(n) = sprintf("DataNBA_%d.dat",n)
f10(n) = sprintf("DataCL_%d.dat",n)
f11(n) = sprintf("DataMY_%d.dat",n)

fT1(n) = sprintf("DataA_%d_1T.dat",n)
fT2(n) = sprintf("DataA_%d_2T.dat",n)
fT3(n) = sprintf("DataA_%d_3T.dat",n)
fT4(n) = sprintf("DataA_%d_4T.dat",n)
fT5(n) = sprintf("DataA_%d_5T.dat",n)

set key below
do for [i=1:1000]{
   nom="AActina".i."T.png";
   set output nom
   plot fT1(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#4b03a1' title "Tensión<=20%",\
        fT2(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#06a4ca' title "Tensión<=40%",\
        fT3(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#66A61E' title "Tensión<=60%",\
        fT4(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#f8c932' title "Tensión<=80%",\
        fT5(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#FF0000' title "Tensión Rotura"
}

no

unset key
do for [i=1:1000]{
   nom="AActina".i.".png";
   set output nom
   plot f91(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#66A61E' title "Actin Bundle 1",\
        f92(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#f8c932' title "Actin Bundle 2",\
        f93(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#FF0000' title "Actin Bundle 3",\
        f94(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#4b03a1' title "Actin Bundle 4",\
        f95(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#aadc32' title "Actin Bundle 5",\
        f96(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#E6AB02' title "Actin Bundle 6",\
        f97(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#cb4679' title "Actin Bundle 7",\
        f98(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#4dbeee' title "Actin Bundle 8",\
        f99(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#7570B3' title "Actin Bundle 9",\
        f90(i) u 1:2 w p lw 1.5 pt 5 lc rgb '#666666' title "Actin NoBundle",\
        f10(i) u 1:2 w p lw 2 pt 12 lc rgb '#000004' title "C-link",\
        f11(i) u 1:2 w p lw 2 pt 13 lc rgb '#000004' title "Myo"
}
