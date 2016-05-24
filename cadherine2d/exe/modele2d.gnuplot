######   Prologue   #####
#clear					# n'est plus utile
reset					# all graph-related options take on their default values

######  Load data   #####

######  Set labels  #####
#set xlabel 'x'
#set ylabel 'y'
#set zlabel 'u'

###### Plot output  #####
#set terminal gif size 640, 960 animate \
#    delay 20 loop 1 optimize		# safest to give `set terminal` first
#set output "cadherins.gif"
set terminal qt size 640, 960

###### Plot options #####
set grid
set surface
set dgrid3d 99,99 qnorm 2 # 99*99 est le nombre de points de la grille
#set dgrid3d 99,99
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

###### Stat of data #####
stats 'modele2d.out' using 4:5 nooutput

######  Plot data   #####
do for [i=0:int(STATS_blocks-1)] {	# si 2 lignes entre les blocs
    ######  Multiplot   #####
    set size   1, 1
    set origin 0, 0
    set multiplot layout 2, 1 columnsfirst scale 1, 1	# 2 rows, 1 column
    #set terminal qt 0			# n'est plus utile
    set key title 'free cadherins at t'.i
    set key title 'free cadherins at t'.i
    #splot 'modele2d.out' index i using 2:3:4 notitle
    splot 'modele2d.out' index i using 2:3:4 palette notitle
    #set terminal qt 1			# n'est plus utile
    set key title 'fixed cadherins at t'.i
    splot 'modele2d.out' index i using 2:3:5 notitle
    pause -1				# n'est plus utile
    unset multiplot			# turn off automatic layout and restore the values of `size` and
					#     `origin` as they were before `set multiplot layout`
}
