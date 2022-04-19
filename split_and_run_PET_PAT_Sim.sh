#!/bin/bash


# to arguments are needed : PATIENT ($1)  and SimLABEL ($2)  partition ($3)  ntasks ($4)

PATIENT=$1
#Sim2000 or other
simLABEL=$2
#partition name
partPARAM=$3
#number of cores
ntasksPARAM=$4

#load your virtual environment
source /home/jpet/virtplan/bin/activate

#load GATE module
module load gate/9.0


# config:
#getPlanDIR=/home/jpet/scripts/getPlan/
patientsDIR=/home/jpet/patientData/
#simLabel group (configuration)
simLABELgroup=Sim2000
#simulations configuration folder
GateSimConf=/home/jpet/GatePETSimConf/GatePETSimConf_${simLABELgroup}/
scriptsDIR=/home/jpet/scripts/

cd ${patientsDIR}

patientID=$PATIENT

#for patientID in */; do
    echo "$patientID"
    petsimDIR=${patientsDIR}${patientID}/gate/${simLABEL}
    echo $fieldsDIR

    outputBeamFile=${petsimDIR}/Beam.sh
    outputBatchFile=${petsimDIR}/batch.sh
    
    export GC_GATE_EXE_DIR=/opt/gate/v9.0/bin/
    export GC_DOT_GATE_DIR=${petsimDIR}

    if [[ -d ${petsimDIR}/.Gate ]]
    then
	`rm -r ${petsimDIR}/.Gate/*`
    fi
    
    gjs -numberofsplits ${ntasksPARAM} -clusterplatform openmosix "${petsimDIR}"/patient_PET.mac 


echo "#!/bin/bash

#Based on example from https://gooseslurm.readthedocs.io/en/latest/examples/tempdir/readme.html
#get name of the temporary directory working directory, physically on the compute-node
workdir=\"\${TMPDIR}/jpet/beamjob/\${SLURM_JOB_ID}/\${SLURM_PROCID}/data\"
if [ ! -d \"\${workdir}\" ]; then
  mkdir -p \"\${workdir}\"
fi

echo $workdir
# get submit directory
#submitdir=\"\${SLURM_SUBMIT_DIR}\"
submitdir=${petsimDIR}
echo \$submitdir
#mkdir -p /tmp/jpet-gate/gate
cd \$workdir

if [ ! -d Results ]; then
  mkdir -p Results
fi
#run gate
Gate ${petsimDIR}/.Gate/patient_PET/patient_PET\$((\$SLURM_PROCID+1)).mac

#cd Results
mv *.root \${submitdir}/Results
#mv *.bin \${submitdir}/Results
#mv *.txt \${submitdir}/Results
#mv *.out \${submitdir}/Results
#mv *.dat \${submitdir}/Results

function clean_up {
  # - cleanup
  rm -rf \"\${TMPDIR}/jpet/beamjob/\${SLURM_JOB_ID}/\${SLURM_PROCID}\"
  # - exit the script
  exit
}
# definen_up\" function when this script exits, it is run even if SLURM cancels the job
trap 'clean_up' EXIT" > $outputBeamFile
`chmod 770 ${outputBeamFile}`

if [ "$partPARAM" = "quanta" ]; then
    moduleVer="gate/9.0_Haswell"


echo "#!/bin/bash
#SBATCH -J jpet_PET_${patientID}_${simLABEL}_${counter}
#SBATCH --ntasks=${ntasksPARAM}
#SBATCH -p ${partPARAM}
#SBATCH --output=${petsimDIR}/c_%j.out
#SBATCH --error=${petsimDIR}/c_%j.err
#SBATCH --time=100:00:00
#SBATCH --hint=nomultithread
#SBATCH --reservation=jpet_sim

#module load gate/9.0_Haswell
module load ${moduleVer}
#export GC_DOT_GATE_DIR=./
export GC_DOT_GATE_DIR=${petsimDIR}
if [ ! -d \"${petsimDIR}/Results\" ]; then
  echo \"Creating ${petsimDIR}/Results directory...\"
  mkdir -p ${petsimDIR}/Results
fi
srun ${petsimDIR}/Beam.sh" > $outputBatchFile



else   #for ibm partition
    moduleVer="gate/9.0"


echo "#!/bin/bash
#SBATCH -J jpet_PET_${patientID}_${simLABEL}_${counter}
#SBATCH --ntasks=${ntasksPARAM}
#SBATCH -p ${partPARAM}
#SBATCH --output=${petsimDIR}/c_%j.out
#SBATCH --error=${petsimDIR}/c_%j.err
#SBATCH --time=100:00:00
#SBATCH --hint=nomultithread
##SBATCH --reservation=jpet_sim

#module load gate/9.0_Haswell
module load ${moduleVer}
#export GC_DOT_GATE_DIR=./
export GC_DOT_GATE_DIR=${petsimDIR}
if [ ! -d \"${petsimDIR}/Results\" ]; then
  echo \"Creating ${petsimDIR}/Results directory...\"
  mkdir -p ${petsimDIR}/Results
fi
srun ${petsimDIR}/Beam.sh" > $outputBatchFile


fi



workdir=`pwd`
echo $workdir

#RUN prepared simulations for all Fields
echo "QUEUEING SLURM JOBS"

`cd ${petsimDIR}`
`sbatch ${petsimDIR}/batch.sh`
`cd $workdir`



deactivate

