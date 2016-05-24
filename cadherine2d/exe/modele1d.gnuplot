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
#set xrange [ : ]
set yrange [-0.1 : 2.1]
set grid

###### Stat of data #####
stats 'modele1d.out' nooutput

######  Plot data   #####
do for [i=0:int(STATS_blocks-1)] {	# si 2 lignes entre les blocs
    ######  Multiplot   #####
    set size   1, 1
    set origin 0, 0
    set multiplot layout 2, 1 columnsfirst scale 1, 1	# 2 rows, 1 column
    #set terminal qt 0			# n'est plus utile
    set key title 'free cadherins at t'.i
    plot 'modele1d.out' index i using 2 : 4 with lines notitle
    #set terminal qt 1			# n'est plus utile
    set key title 'fixed cadherins at t'.i
    plot 'modele1d.out' index i using 2 : 5 with lines notitle
    pause -1				# n'est plus utile
    unset multiplot			# turn off automatic layout and restore the values of `size` and
					#     `origin` as they were before `set multiplot layout`
}

######   Epilogue   #####
set output				# any previously opened output file will be closed
