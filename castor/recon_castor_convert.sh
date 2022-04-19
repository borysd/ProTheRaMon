#!/bin/bash

#module load castor/3.1_JPET
module load castor/3.1-tmerlin

patientID=G3P007E2

#scanner geometry name

#scanner=BARREL1
scanner=BARREL2
#scanner=BARREL3
#scanner=DUALHEAD_1x6
#scanner=DUALHEAD_2x6
#scanner=DUALHEAD_3x4
#scanner=DUALHEAD_3x4_experimental
#scanner=dualhead_ccb  #2780 ...

#for simnum in {2780..2780..20}
for simnum in {2200..2214..1}
do
    simLabel=Sim${simnum}
    pathToSim=/home/jpet/patientData/${patientID}/gate/${simLabel}
    pathToRoots=/home/jpet/patientData/${patientID}/gate/${simLabel}/Results

    listRootsTxt=/home/jpet/patientData/${patientID}/gate/${simLabel}/Results/listRoots.txt


    `ls  ${pathToRoots} | grep .root > ${listRootsTxt}` 
    
    # TOF 230
    castor-GATERootToCastor -il ${listRootsTxt} -o ${pathToSim}/patient_lm_230 -m  ${pathToSim}/patient_PET.mac -s ${scanner} -ot -TOF_reso 230
    # TOF 500
    castor-GATERootToCastor -il ${listRootsTxt} -o ${pathToSim}/patient_lm_500 -m  ${pathToSim}/patient_PET.mac -s ${scanner} -ot -TOF_reso 500
    # TOF 700
    castor-GATERootToCastor -il ${listRootsTxt} -o ${pathToSim}/patient_lm_700 -m  ${pathToSim}/patient_PET.mac -s ${scanner} -ot -TOF_reso 700

    ./recon_castor_convert_calibration.sh ${pathToSim}/patient_lm_230_df.Cdh $simnum
    ./recon_castor_convert_calibration.sh ${pathToSim}/patient_lm_500_df.Cdh $simnum
    ./recon_castor_convert_calibration.sh ${pathToSim}/patient_lm_700_df.Cdh $simnum

    #castor-recon  -df my_data.cdh -opti MLEM -it 10:16 -proj joseph  -conv gaussian,4.,4.5,3.5::psf -dim 128,128,128  -vox 3.,3.,3. -dout my_images
done