#!/bin/bash

#Haiyang Kehoe
#University of Arizona
#Department of Geosciences
#19 February 2020
#Modified 18 October 2022

dir='0.20s/Aftershock_0_xcor'
cd $dir
idbp=idbp_radiators.txt
cent=idbp_radiators_cent.txt
param=../../parameters.txt
i=sum.txt

gmt gmtset FONT_LABEL 12p,Helvetica,black
gmt gmtset FONT_ANNOT_PRIMARY 12p,Helvetica,black
gmt gmtset FORMAT_GEO_MAP ddd.x
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset PS_MEDIA A2
gmt gmtset PS_PAGE_ORIENTATION landscape
gmt gmtset PS_LINE_CAP butt
gmt gmtset PS_LINE_JOIN miter

lonxc=$(awk '{print $1}' $cent)
latxc=$(awk '{print $2}' $cent)
depxc=$(awk '{print $3}' $cent)

PSFILE="idbp_bp.ps"

gmt makecpt -T0/1/0.01 -Z -Chaxby > col.cpt
gmt makecpt -Chot -T0/2.0/0.1 > myvel.cpt

north_idbp=$(awk '{print $3}' $idbp | sort | tail -1)
south_idbp=$(awk '{print $3}' $idbp | sort | head -1)
west_idbp=$(awk '{print $2}' $idbp | sort | head -1)
east_idbp=$(awk '{print $2}' $idbp | sort | tail -1)
top_idbp=$(awk '{print $4}' $idbp | sort | head -1)
bot_idbp=$(awk '{print $4}' $idbp | sort | tail -1)

north=$(awk '{print $2}' $i | sort -nr | head -n1)
south=$(awk '{print $2}' $i | sort -n | head -n1)
west=$(awk '{print $1}' $i | sort -n | head -n1)
east=$(awk '{print $1}' $i | sort -nr | head -n1)
top=$(awk '{print $3}' $i | sort -n | head -n1)
bot=$(awk '{print $3}' $i | sort -nr | head -n1)

latd=$(awk '{print $2}' $i | uniq | tail -r -2 | paste -sd- - | bc)
lond=$(awk '{print $1}' $i | uniq | tail -r -2 | paste -sd- - | bc)
depd=$(awk '{print $3}' $i | uniq | tail -r -2 | paste -sd- - | bc)

#Plot map-view
gmt psbasemap -Y5 -R$west/$east/$south/$north -JX20d/24d -B0.1NSEW -K > $PSFILE
awk -v var=$depxc '$3==var{print $1,$2,$4}' $i |\
gmt xyz2grd -R -Genv.grd -I$latd/$lond
gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >> $PSFILE
gmt grdcontour env.grd -R -J -C0.1 -L0.8/1 -O -K -W1,white >> $PSFILE
gmt psxy -R -J -A -W0.5,white,- -O -K << EOF >> $PSFILE
$west $latxc
$east $latxc
EOF
gmt psxy -R -J -A -W0.5,white,- -O -K << EOF >> $PSFILE
$lonxc $south
$lonxc $north
EOF
gmt pstext -R -J -W1.0 -G255/255/255 -N -O -K << EOF >> $PSFILE
$west $latxc B
$east $latxc BB
$lonxc $south A
$lonxc $north AA
EOF
awk '{print $2,$3,$5}' $idbp |\
gmt psxy -R -J -Ccol.cpt -W0.8,black -Sc0.4 -O -K >> $PSFILE
awk '{print $1,$2,$8,$4,$5}' $param |\
gmt psxy -R -J -Sv0.5+s+e -Cmyvel.cpt -W2.0 -K -O >> $PSFILE

#Plot Latitude-Depth cross-section
gmt psbasemap -X24 -Y13 -R$south/$north/$top/$bot -JX20d/-11 -B0.1:"Latitude":/10:"Depth (km)":nSEW -O -K >> $PSFILE
awk -v var=$lonxc '$1==var{print $2,$3,$4}' $i |\
gmt xyz2grd -R -Genv.grd -I$latd/$depd
gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >> $PSFILE
gmt grdcontour env.grd -R -J -C0.1 -L0.8/1 -O -K -W1,white >> $PSFILE
gmt psxy -R -J -A -W0.5,white,- -O -K << EOF >> $PSFILE
$south $depxc
$north $depxc
EOF
gmt pstext -R -J -W1.0 -G255/255/255 -N -O -K << EOF >> $PSFILE
$south $depxc A
$north $depxc AA
EOF
awk '{print $3,$4,$5}' $idbp |\
gmt psxy -R -J -Ccol.cpt -W0.8,black -Sc0.4 -O -K >> $PSFILE
awk '{print $2,$3,$8,$5,$6}' $param |\
gmt psxy -R -J -Sv0.5+s+e -Cmyvel.cpt -W2.0 -K -O >> $PSFILE

#Plot Longitude-Depth cross-section
gmt psbasemap -Y-13 -R$west/$east/$top/$bot -JX20d/-11 -B0.1:"Latitude":/10:"Depth (km)":nSEW -O -K >> $PSFILE
awk -v var=$latxc '$2==var{print $1,$3,$4}' $i |\
gmt xyz2grd -R -Genv.grd -I$lond/$depd
gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >> $PSFILE
gmt grdcontour env.grd -R -J -C0.1 -L0.8/1 -O -K -W1,white >> $PSFILE
gmt psxy -R -J -A -W0.5,white,- -O -K << EOF >> $PSFILE
$west $depxc
$east $depxc
EOF
gmt pstext -R -J -W1.0 -G255/255/255 -N -O -K << EOF >> $PSFILE
$west $depxc B
$east $depxc BB
EOF
awk '{print $2,$4,$5}' $idbp |\
gmt psxy -R -J -Ccol.cpt -W0.8,black -Sc0.4 -O -K >> $PSFILE
awk '{print $1,$3,$8,$4,$6}' $param |\
gmt psxy -R -J -Sv0.5+s+e -Cmyvel.cpt -W2.0 -K -O >> $PSFILE

gmt psscale -D10.0/-1.2/20.0/0.3h -B0.25+l"Normalized Amplitude" -X-24 -K -O -Ccol.cpt >> $PSFILE
gmt psscale -D10.0/-1.2/20.0/0.3h -B0.5+l"Rupture Velocity (v@-r@-/v@-S@-)" -X24 -O -Cmyvel.cpt >> $PSFILE

gmt psconvert $PSFILE -A0.2c -Tf -P

cd ..
