#####   Prologue   #####
clear # erases the current screen or output device 
reset # all graph-related options take on their default values

######  Load data   #####

######  Set labels  #####
#set xlabel 'x'
#set ylabel 'y'
#set zlabel 'u'
#set grid

###### Plot options #####
set yrange [-0.1:1.1]

######  Plot data   #####
stats 'modele.out'
#si 2 lignes entre les blocs
do for [i=0:int(STATS_blocks-1)] {
    set term qt 0
    set key title 'free cadherins at t'.i
    plot 'modele.out' index i using 2:4 with lines notitle
    set term qt 1
    set key title 'fixed cadherins at t'.i
    plot 'modele.out' index i using 2:5 with lines notitle
    pause -1
}
