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
set surface
set dgrid3d 11,11 qnorm 2 # 11*11 est le nombre de points de la grille
#set dgrid3d 11,11
set dgrid3d splines
set style data lines
set hidden3d
#set pm3d # plot en faux 3d
set pm3d map # plot en 2d i.e. "carte" avec couleur
#set ticslevel 0.8
#set isosamples 40,40
#set view 60, 30, 1, 1 # pour voir depuis un angle particulier
#set contour base
set cbrange [0:1]
set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')

######  Plot data   #####
stats 'modele.out' using 4:5
do for [i=0:int(STATS_blocks-1)] {
    set term qt 0
    set key title 'free cadherins at t'.i
    splot 'modele.out' index i using 2:3:4 notitle
#    splot 'modele.out' index i using 2:3:4 pal notitle
    set term qt 1
    set key title 'fixed cadherins at t'.i
    splot 'modele.out' index i using 2:3:5 notitle
    pause -1
}
