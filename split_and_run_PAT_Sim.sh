#!/bin/bash


# to arguments are needed : PATIENT ($1)  and SimLABEL ($2)  partition ($3)  ntasks ($4)

PATIENT=$1
simLABEL=$2
partPARAM=$3
ntasksPARAM=$4

source /home/jpet/virtplan/bin/activate

module load gate/9.0


# config:
getPlanDIR=/home/jpet/scripts/getPlan/
patientsDIR=/home/jpet/patientData/
CTfilename=CT_ExtROICrop_2.5x2.5x2.5.mhd
#simLabel group (configuration)
simLABELgroup=Sim1000
GateSimConf=/home/jpet/GateSimConf/GateSimConf_${simLABELgroup}/

scriptsDIR=/home/jpet/scripts/

#simLABEL=Sim1010


cd ${patientsDIR}

patientID=$PATIENT

#for patientID in */; do
    echo "$patientID"
#    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l ${simLABEL} -n 0.1 --Nunit primPr --CPUs 1000 -a prodAll -x0 -y0 -z0
    
    fieldsDIR=${patientsDIR}${patientID}/gate/${simLABEL}
    echo $fieldsDIR

    counter=1
    for fieldID in ${fieldsDIR}/*; do
	if [[ -d "${fieldID}" ]]
	then
	    echo "$fieldID"

	    
	    outputBeamFile=${fieldID}/Beam.sh
	    outputBatchFile=${fieldID}/batch.sh
	    export GC_GATE_EXE_DIR=/opt/gate/v9.0/bin/
	    export GC_DOT_GATE_DIR=${fieldID}
	    
	    if [[ -d ${fieldID}/.Gate ]]
	    then
		`rm -r ${fieldID}/.Gate/*`
	    fi
	    
	    gjs -numberofsplits ${ntasksPARAM} -clusterplatform openmosix "${fieldID}"/mainSim.mac 

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
submitdir=${fieldID}
echo \$submitdir
#mkdir -p /tmp/jpet-gate/gate
cd \$workdir

if [ ! -d Results ]; then
  mkdir -p Results
fi
#run gate
Gate ${fieldID}/.Gate/mainSim/mainSim\$((\$SLURM_PROCID+1)).mac

cd Results
mv *map*.txt \${submitdir}/Results
mv *.mhd \${submitdir}/Results
mv *.raw \${submitdir}/Results
mv *Stat*.out \${submitdir}/Results

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
else
    moduleVer="gate/9.0"
fi


echo "#!/bin/bash
#SBATCH -J jpet_${patientID}_${simLABEL}_${counter}
#SBATCH --ntasks=${ntasksPARAM}
#SBATCH -p ${partPARAM}
#SBATCH --output=${fieldID}/c_%j.out
#SBATCH --error=${fieldID}/c_%j.err
#SBATCH --time=180:00:00
#SBATCH --hint=nomultithread
#SBATCH --reservation=jpet_sim

#module load gate/9.0_Haswell
module load ${moduleVer}
#export GC_DOT_GATE_DIR=./
export GC_DOT_GATE_DIR=${fieldID}
if [ ! -d Results ]; then
  mkdir -p Results
fi
srun ${fieldID}/Beam.sh" > $outputBatchFile


	    counter=$((counter+1))
	fi

    done

    workdir=`pwd`
    echo $workdir

    #RUN prepared simulations for all Fields
    echo "QUEUEING SLURM JOBS"
    
    fieldsDIR=${patientsDIR}${patientID}/gate/${simLABEL}
    
    
    for fieldID in ${fieldsDIR}/*; do
	if [[ -d "${fieldID}" ]]
	then
	    echo "$fieldID"
	    `cd ${fieldID}`
	    `sbatch ${fieldID}/batch.sh`
	fi
    done

    `cd $workdir`

#done


deactivate

