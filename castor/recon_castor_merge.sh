#!/bin/bash
#SBATCH -J mergeAttSens
#SBATCH --ntasks=20
#SBATCH -p quanta
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

#module load castor/3.1_JPET
module load castor/3.1-tmerlin

#attSource=K
attSource=J

#patientID=Exp_GTR4_SOBP_S2F5_2Gy
patientID=G3P007E1


#geometry=barrel
geometry=barrel2
#geometry=barrel3
#geometry=dualhead_1x6
#geometry=dualhead_2x6
#geometry=dualhead_3x4
#geometry=sensitivity_lm_dualhead_ccb

geomName=${geometry}_${attSource}

#attFile=umap_CT_CT_ExtROICrop_2.5x2.5x2.5_${attSource}.hdr
attFile=umap_CT_CT_2.5x2.5x2.5_${attSource}.hdr
#sensFile=barrel_df.Cdh
#sensFile=barrel2_df.Cdh
#sensFile=barrel3_df.Cdh
#sensFile=dualhead_1x6_df.Cdh
#sensFile=dualhead_2x6_df.Cdh
sensFile=${geometry}_df.Cdh

x=80
y=60
z=80

pathToPat=/home/jpet/patientData/${patientID}
pathToSens=/scratch/scratch-ssd/jpet/noTOF/

castor-recon  -df ${pathToSens}${sensFile}  -dim  ${x},${y},${z}  -vox  2.5,2.5,2.5  -th  0  -it  1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
