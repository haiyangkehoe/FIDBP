#!/bin/bash

#Haiyang Kehoe
#University of Arizona
#Department of Geosciences
#19 February 2020
#Modified 19 February 2020

gmt gmtset FONT_LABEL 12p,Helvetica,black
gmt gmtset FONT_ANNOT_PRIMARY 12p,Helvetica,black
gmt gmtset FORMAT_GEO_MAP ddd.x
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset PS_MEDIA A2
gmt gmtset PS_PAGE_ORIENTATION landscape
gmt gmtset PS_LINE_CAP butt
gmt gmtset PS_LINE_JOIN miter

dir='Pad_Aftershock_00_syn_xcor'
rm files.txt
ls $dir/3D_*.txt > files.txt
filename=$dir/hyp.txt

t_step=$(awk '{print $5}' $filename)
t_min=$(awk -v ts=$t_step '{print (1-$1)*ts}' $filename)
#Make t_min have the same sig figs as t_step:
t_min=$(echo "$(echo "$t_min + $t_step" | bc) - $t_step" | bc)
hypo_lat=$(awk '{print $3}' $filename)
hypo_lon=$(awk '{print $2}' $filename)
hypo_dep=$(awk '{print $4}' $filename)

PSFILE="results.ps"

gmt makecpt -T0/1/0.01 -Z -Chaxby > col.cpt

idx=1
while read i; do
  echo $t_min
  peak=$(awk 'BEGIN{a=0}{if ($4>0+a) {a=$4; l=$0}} END{print l}' $i)
  lonxc=$(echo $peak | awk '{print $1}')
  latxc=$(echo $peak | awk '{print $2}')
  depxc=$(echo $peak | awk '{print $3}')

  north=$(awk '{print $2}' $i | sort -nr | head -n1)
  south=$(awk '{print $2}' $i | sort -n | head -n1)
  west=$(awk '{print $1}' $i | sort -n | head -n1)
  east=$(awk '{print $1}' $i | sort -nr | head -n1)
  top=$(awk '{print $3}' $i | sort -n | head -n1)
  bot=$(awk '{print $3}' $i | sort -nr | head -n1)

  latd=$(awk '{print $2}' $i | uniq | tail -r -2 | paste -sd- - | bc)
  lond=$(awk '{print $1}' $i | uniq | tail -r -2 | paste -sd- - | bc)
  depd=$(awk '{print $3}' $i | uniq | tail -r -2 | paste -sd- - | bc)

  #echo $north $south $west $east $top $bot $latd $lond $depd

  gmt psbasemap -Y5 -R$west/$east/$south/$north -JX20d/24d -B0.1NSEW -K > $PSFILE
  if [[ -z $peak ]]; then
    awk -v var=$depxc '$3==var{print $1,$2,$4}' $i |\
    gmt xyz2grd -R -Genv.grd -I$latd/$lond -di0
    gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >>$PSFILE
    gmt psxy -R -J -G255/255/255 -W0.8,black -Sa0.4 -O -K << EOF >> $PSFILE
    $hypo_lon $hypo_lat
EOF
    pstext -R -J -F+cTL -C0.1 -D0.25/-0.25 -G255/255/255 -W1.0 -N -O -K << EOF >> $PSFILE
    $t_min s
EOF
  else
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
    gmt psxy -R -J -G255/255/255 -W0.8,black -Sa0.4 -O -K << EOF >> $PSFILE
    $hypo_lon $hypo_lat
EOF
    pstext -R -J -F+cTL -C0.1 -D0.25/-0.25 -G255/255/255 -W1.0 -N -O -K << EOF >> $PSFILE
    $t_min s
EOF
  fi

  gmt psbasemap -X24 -Y14 -R$south/$north/$top/$bot -JX20d/-10 -B0.1:"Latitude":/10:"Depth (km)":nSEW -O -K >> $PSFILE
  if [[ -z $peak ]]; then
    awk -v var=$lonxc '$1==var{print $2,$3,$4}' $i |\
    gmt xyz2grd -R -Genv.grd -I$latd/$depd -di0
    gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >>$PSFILE
    gmt psxy -R -J -G255/255/255 -W0.8,black -Sa0.4 -O -K << EOF >> $PSFILE
    $hypo_lat $hypo_dep
EOF
    pstext -R -J -F+cTL -C0.1 -D0.25/-0.25 -G255/255/255 -W1.0 -N -O -K << EOF >> $PSFILE
    $t_min s
EOF
  else
    awk -v var=$lonxc '$1==var{print $2,$3,$4}' $i |\
    gmt xyz2grd -R -Genv.grd -I$latd/$depd
    gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >> $PSFILE
    gmt grdcontour env.grd -R -J -C0.1 -L0.8/1 -O -K -W1,white >> $PSFILE
    gmt psxy -R -J -G255/255/255 -W0.8,black -Sa0.4 -O -K << EOF >> $PSFILE
    $hypo_lat $hypo_dep 
EOF
    pstext -R -J -F+cTL -C0.1 -D0.25/-0.25 -G255/255/255 -W1.0 -N -O -K << EOF >> $PSFILE
    $t_min s
EOF
    gmt psxy -R -J -A -W0.5,white,- -O -K << EOF >> $PSFILE
    $south $depxc
    $north $depxc
EOF
    gmt pstext -R -J -W1.0 -G255/255/255 -N -O -K << EOF >> $PSFILE
    $south $depxc A
    $north $depxc AA
EOF
  fi

  gmt psbasemap -Y-12 -R$west/$east/$top/$bot -JX20d/-10 -B0.1:"Latitude":/10:"Depth (km)":nSEW -O -K >> $PSFILE
  if [[ -z $peak ]]; then
    awk -v var=$latxc '$2==var{print $1,$3,$4}' $i |\
    gmt xyz2grd -R -Genv.grd -I$lond/$depd -di0
    gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >>$PSFILE 
    gmt psxy -R -J -G255/255/255 -W0.8,black -Sa0.4 -O -K << EOF >> $PSFILE
    $hypo_lon $hypo_dep 
EOF
    pstext -R -J -F+cTL -C0.1 -D0.25/-0.25 -G255/255/255 -W1.0 -N -O -K << EOF >> $PSFILE
    $t_min s
EOF
  else
    awk -v var=$latxc '$2==var{print $1,$3,$4}' $i |\
    gmt xyz2grd -R -Genv.grd -I$lond/$depd
    gmt grdimage env.grd -R -J -Ccol.cpt -n+c -O -K >> $PSFILE
    gmt grdcontour env.grd -R -J -C0.1 -L0.8/1 -O -K -W1,white >> $PSFILE
    gmt psxy -R -J -G255/255/255 -W0.8,black -Sa0.4 -O -K << EOF >> $PSFILE
    $hypo_lon $hypo_dep
EOF
    pstext -R -J -F+cTL -C0.1 -D0.25/-0.25 -G255/255/255 -W1.0 -N -O -K << EOF >> $PSFILE
    $t_min s
EOF
    gmt psxy -R -J -A -W0.5,white,- -O -K << EOF >> $PSFILE
    $west $depxc
    $east $depxc
EOF
    gmt pstext -R -J -W1.0 -G255/255/255 -N -O -K << EOF >> $PSFILE
    $west $depxc B
    $east $depxc BB
EOF
  fi

  gmt psscale -D10.0/-1.2/20.0/0.3h -B0.1+l"Normalized Amplitude" -O -Ccol.cpt >> $PSFILE
  
  gmt psconvert $PSFILE -A0.2c -Tf -P

  mv results.pdf $i.pdf

  t_min=$(echo "$t_min + $t_step" | bc)
  #t_min=$(( $t_min+$t_step ))
  idx=$(( $idx+1 ))

done < files.txt
