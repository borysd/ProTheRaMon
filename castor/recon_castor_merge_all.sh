#!/bin/bash

patientID=$1 
#patientID=G3P035E0


# READ x y z dimensions from CT file
source /home/jpet/venv38/bin/activate
xyz=$(python3 /home/jpet/scripts/castor/xyzRead.py -P ${patientID} 2>&1)
#xyz=$(python3 /home/jpet/scripts/castor/xyzReadNoCrop.py -P ${patientID} 2>&1)

lenxyz=$((${#xyz}-2))
xyz_=`echo ${xyz:1:${lenxyz}}`
IFS=',' read -ra xxyyzz <<< `echo ${xyz_}| tr -d '[:space:]'`
x=$((${xxyyzz[0]}))
y=$((${xxyyzz[1]}))
z=$((${xxyyzz[2]}))

#x=80
#y=60
#z=80

attSource=J
pathToPat=/home/jpet/patientData/${patientID}
pathToSens=/scratch/scratch-ssd/jpet/noTOF/

FILE=${pathToPat}/TPS/umap_CT_CT_ExtROICrop_2.5x2.5x2.5_${attSource}.hdr

if test -f "$FILE"; then
    echo "$FILE exists."
    attFile=umap_CT_CT_ExtROICrop_2.5x2.5x2.5_${attSource}.hdr
else
    attFile=umap_CT_CT_2.5x2.5x2.5_${attSource}.hdr
fi
#attFile=umap_CT_CT_ExtROICrop_2.5x2.5x2.5_${attSource}.hdr
#attFile=umap_CT_CT_2.5x2.5x2.5_${attSource}.hdr

# -------------------------------------------------------------------------------

geomName=barrel_${attSource}
sensFile=barrel_df.Cdh

echo "#!/bin/bash
#SBATCH -J mergeAttSens_barrel1
#SBATCH --ntasks=24
#SBATCH -p ibm
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

module load castor/3.1-tmerlin

castor-recon  -df ${pathToSens}${sensFile}  -dim ${x},${y},${z}  -vox 2.5,2.5,2.5  -th 0  -it 1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
" | sbatch

# -------------------------------------------------------------------------------

geomName=barrel2_${attSource}
sensFile=barrel2_df.Cdh

echo "#!/bin/bash
#SBATCH -J mergeAttSens_barrel2
#SBATCH --ntasks=24
#SBATCH -p ibm
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

module load castor/3.1-tmerlin

castor-recon  -df ${pathToSens}${sensFile}  -dim  ${x},${y},${z}  -vox  2.5,2.5,2.5  -th  0  -it  1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
" | sbatch

# -------------------------------------------------------------------------------

geomName=barrel3_${attSource}
sensFile=barrel3_df.Cdh

echo "#!/bin/bash
#SBATCH -J mergeAttSens_barrel3
#SBATCH --ntasks=24
#SBATCH -p ibm
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

module load castor/3.1-tmerlin

castor-recon  -df ${pathToSens}${sensFile}  -dim  ${x},${y},${z}  -vox  2.5,2.5,2.5  -th  0  -it  1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
" | sbatch

# -------------------------------------------------------------------------------

geomName=dualhead_1x6_${attSource}
sensFile=dualhead_1x6_df.Cdh

echo "#!/bin/bash
#SBATCH -J mergeAttSens_dual1
#SBATCH --ntasks=24
#SBATCH -p ibm
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

module load castor/3.1-tmerlin

castor-recon  -df ${pathToSens}${sensFile}  -dim  ${x},${y},${z}  -vox  2.5,2.5,2.5  -th  0  -it  1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
" | sbatch

# -------------------------------------------------------------------------------

geomName=dualhead_2x6_${attSource}
sensFile=dualhead_2x6_df.Cdh

echo "#!/bin/bash
#SBATCH -J mergeAttSens_dual2
#SBATCH --ntasks=24
#SBATCH -p ibm
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

module load castor/3.1-tmerlin

castor-recon  -df ${pathToSens}${sensFile}  -dim  ${x},${y},${z}  -vox  2.5,2.5,2.5  -th  0  -it  1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
" | sbatch

# -------------------------------------------------------------------------------

geomName=dualhead_3x4_${attSource}
sensFile=dualhead_3x4_df.Cdh

echo "#!/bin/bash
#SBATCH -J mergeAttSens_dual3
#SBATCH --ntasks=24
#SBATCH -p ibm
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/merge_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/merge_%j.e
##SBATCH --reservation=jpet_sim

module load castor/3.1-tmerlin

castor-recon  -df ${pathToSens}${sensFile}  -dim  ${x},${y},${z}  -vox  2.5,2.5,2.5  -th  0  -it  1:1  -fout ${pathToPat}/TPS/sens_atten_${geomName} -opti SENS -img ${pathToPat}/TPS/${attFile} -ignore-corr fdur
" | sbatch