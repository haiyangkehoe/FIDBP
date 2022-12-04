#Haiyang Kehoe
#University of Arizona
#Department of Geosciences
#08 September 2022
#Modified 08 September 2022

MS='Aftershock_0'
ALIGN='Aftershock_1_Alignment'

#for i in {0.20,0.50,1.00}; do
for i in {0.20,0.25}; do
  if [ ! -d ${i}s ]; then
    mkdir ${i}s
    mkdir ${i}s/${MS}_xcor
    mkdir ${i}s/${MS}_syn_xcor
    cd ../${MS}_xcor
    rm 3D_*.txt
    ./plot_rf_frames_output_get_bp.x $i
    cp 3D_*.txt ../IDBP_${MS}/${i}s/${MS}_xcor
    cd ../../../Data_Hinet_Synthetic/${ALIGN}/${MS}_xcor/
    rm 3D_*.txt
    ./plot_rf_frames_output_get_bp.x $i
    cp 3D_*.txt ../../../Data_Hinet/${ALIGN}/IDBP_${MS}/${i}s/${MS}_syn_xcor/
    cd ../../../Data_Hinet/${ALIGN}/IDBP_${MS}
  fi
done
