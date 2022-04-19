#!/bin/bash

simLabelPARAM=$1
simLabelResultPARAM=$2
xx=$3
yy=$4
zz=$5
geomNamePARAM=$6

source /home/jpet/virtplan/bin/activate

# config:
petPlanDIR=/home/jpet/scripts/petPlan/
patientsDIR=/home/jpet/patientData/
#CTfilename=CT_ExtROICrop_2.5x2.5x2.5.mhd
CTfilename=CT_2.5x2.5x2.5.mhd

#simLabel group (configuration)
simLABELgroup=Sim2000
GatePETSimConf=/home/jpet/GatePETSimConf/GatePETSimConf_${simLABELgroup}/

geometryName=${geomNamePARAM}

#simLABEL=Sim1000
simLABEL=${simLabelPARAM}
#CPUs=${ntasksPARAM}

echo "Scanning input file"

input="patientsList.txt"
while IFS= read -r line
do
    patientID=${line}
    echo ${patientID}
    echo ${patientsDIR}
    simPath=${patientsDIR}${patientID}/gate/${simLabelResultPARAM}
    python3 ${petPlanDIR}petPlan_CCB.py -c${GatePETSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l${simLABEL} -s${simPath} -g${geometryName} -x${xx} -y${yy} -z${zz}
done < "$input"

deactivate

