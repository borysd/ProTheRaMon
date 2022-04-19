#!/bin/bash
#SBATCH -J castorRecon
#SBATCH --ntasks=10
#SBATCH -p quanta
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/recon_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/recon_%j.e

#module load castor/3.1_JPET
module load castor/3.1-tmerlin

#attSource=K
attSource=J

patientID=G3P007E2
#patientID=Exp_GTR4_SOBP_S1F5_2Gy

#simLabel=Sim2100

source /home/jpet/venv38/bin/activate

#va=$(python test.py 2>&1)

#xyz=$(python3 /home/jpet/scripts/castor/xyzRead.py -P ${patientID} 2>&1)

#lenxyz=$((${#xyz}-2))
#xyz_=`echo ${xyz:1:${lenxyz}}`
#IFS=',' read -ra xxyyzz <<< `echo ${xyz_}`
#x=${xxyyzz[0]}
#y=${xxyyzz[1]}
#z=${xxyyzz[2]}

#echo $x
#echo $y
#echo $z



xyz=$(python3 /home/jpet/scripts/castor/xyzRead.py -P ${patientID} 2>&1)
#xyz=$(python3 /home/jpet/scripts/castor/xyzReadNoCrop.py -P ${patientID} 2>&1)

lenxyz=$((${#xyz}-2))
xyz_=`echo ${xyz:1:${lenxyz}}`
IFS=',' read -ra xxyyzz <<< `echo ${xyz_} | tr -d '[:space:]'`
x=${xxyyzz[0]}
y=${xxyyzz[1]}
z=${xxyyzz[2]}

echo ${xyz_}

echo ${x}
echo ${y}
echo ${z}


#geomName=barrel_${attSource}
geomName=barrel2_${attSource}
#geomName=barrel3_${attSource}
#geomName=dualhead_1x6_${attSource}
#geomName=dualhead_2x6_${attSource}
#geomName=dualhead_3x4_${attSource}
##geomName=sensitivity_lm_dualhead_ccb_${attSource}

#simLabel=Sim2100

#for simnum in {2780..2780..20}
for simnum in {2200..2214..1}
do
    simLabel=Sim${simnum}

    pathToPat=/home/jpet/patientData/${patientID}
    pathToSim=/home/jpet/patientData/${patientID}/gate/${simLabel}

##    castor-recon -df ${pathToSim}/patient_lm_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -th 0 -it 20:1 -dout ${pathToPat}/recon/${simLabel} -ignore-corr fdur
    castor-recon -df ${pathToSim}/patient_lm_230_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -it 12:1 -dout ${pathToPat}/recon/${simLabel}_${attSource}_230 -ignore-corr fdur
    castor-recon -df ${pathToSim}/patient_lm_500_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -it 12:1 -dout ${pathToPat}/recon/${simLabel}_${attSource}_500 -ignore-corr fdur
    castor-recon -df ${pathToSim}/patient_lm_700_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -it 12:1 -dout ${pathToPat}/recon/${simLabel}_${attSource}_700 -ignore-corr fdur
done
