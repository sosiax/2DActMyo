reset
set terminal png enhanced size 1024,768
set encoding koi8u
set encoding utf8
set key below
set xlabel 't [Iterations]'
set lmargin 11
set rmargin 4

set key font ",20"
set xtics font ",15"
set xlabel font ",21"
set ytics font ",15"
set ylabel font ",21"

set output 'A0NA.png'
set ylabel 'Number of Actins'
set title "NA"
plot "DataN.dat" u 1:2 w p lw 1 pt 2 lc 7 notitle,"DataN.dat" u 1:5 w p lw 1 pt 2 lc 6 notitle, "DataN.dat" u 1:($2+$5) w p lw 1 pt 2 lc 1 notitle;

set output 'A0NL.png'
set ylabel 'Number of C-Linker'
set title "NL"
plot "DataN.dat" u 1:3 w p lw 1 pt 2 lc 7 notitle,"DataN.dat" u 1:6 w p lw 1 pt 2 lc 6 notitle, "DataN.dat" u 1:($3+$6) w p lw 1 pt 2 lc 1 notitle;

set output 'A0NM.png'
set ylabel 'Number of Myosin'
set title "NM"
plot "DataN.dat" u 1:4 w p lw 1 pt 2 lc 7 notitle,"DataN.dat" u 1:7 w p lw 1 pt 2 lc 6 notitle, "DataN.dat" u 1:($4+$7) w p lw 1 pt 2 lc 1 notitle;

set output 'A0MU.png'
set ylabel 'Chemical potential'
set title "MU"
plot "DataMU.dat" u 1:2 w p ps 2 pt 2 lc rgb '#FF0000' title "Actins",\
     "DataMU.dat" u 1:3 w p ps 2 pt 2 lc rgb '#06a4ca' title "Cross-Linkers",\
     "DataMU.dat" u 1:4 w p ps 2 pt 2 lc rgb '#E6AB02' title "Myosins";

set output 'A0Fil.png'
set ylabel 'Number of Filaments'
set title "NFil"
plot "DataNF.dat" u 1:2 w p lw 1 pt 2 lc 7 notitle;

set output 'A0AeFil.png'
set ylabel 'Number of Actins in Filaments'
set title "NAeF"
plot "DataNF.dat" u 1:3 w p lw 1 pt 2 lc 7 notitle;

set key below
set output 'A0BA.png'
set ylabel 'Number of Actins in Bundle'
set title "NBA"
plot "DataABA.dat" u 1:2 w l lw 2 lc rgb '#66A61E' title "Actin Bundle 1",\
     "DataABA.dat" u 1:3 w l lw 2 lc rgb '#f8c932' title "Actin Bundle 2",\
     "DataABA.dat" u 1:4 w l lw 2 lc rgb '#FF0000' title "Actin Bundle 3",\
     "DataABA.dat" u 1:5 w l lw 2 lc rgb '#4b03a1' title "Actin Bundle 4",\
     "DataABA.dat" u 1:6 w l lw 2 lc rgb '#aadc32' title "Actin Bundle 5",\
     "DataABA.dat" u 1:7 w l lw 2 lc rgb '#E6AB02' title "Actin Bundle 6",\
     "DataABA.dat" u 1:8 w l lw 2 lc rgb '#cb4679' title "Actin Bundle 7",\
     "DataABA.dat" u 1:9 w l lw 2 lc rgb '#4dbeee' title "Actin Bundle 8",\
     "DataABA.dat" u 1:10 w l lw 2 lc rgb '#7570B3' title "Actin Bundle 9",\
     "DataABA.dat" u 1:($2+$3+$4+$5+$6+$7+$8+$9+$10) w l lw 2 lc rgb '#000004' title "All Actin Bundle";

set key below
set output 'A0TM.png'
set ylabel 'Tension of filaments'
set title "TM"
plot "DataTA_1.dat" u 1:2 w l lw 2 lc rgb '#66A61E' title "Actin Bundle 1",\
     "DataTA_2.dat" u 1:2 w l lw 2 lc rgb '#f8c932' title "Actin Bundle 2",\
     "DataTA_3.dat" u 1:2 w l lw 2 lc rgb '#FF0000' title "Actin Bundle 3",\
     "DataTA_4.dat" u 1:2 w l lw 2 lc rgb '#4b03a1' title "Actin Bundle 4",\
     "DataTA_5.dat" u 1:2 w l lw 2 lc rgb '#aadc32' title "Actin Bundle 5",\
     "DataTA_6.dat" u 1:2 w l lw 2 lc rgb '#E6AB02' title "Actin Bundle 6",\
     "DataTA_7.dat" u 1:2 w l lw 2 lc rgb '#cb4679' title "Actin Bundle 7",\
     "DataTA_8.dat" u 1:2 w l lw 2 lc rgb '#4dbeee' title "Actin Bundle 8",\
     "DataTA_9.dat" u 1:2 w l lw 2 lc rgb '#7570B3' title "Actin Bundle 9";
