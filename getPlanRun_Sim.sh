#!/bin/bash

simLabelPARAM=$1
ntasksPARAM=$2
xx=$3
yy=$4
zz=$5
HU=$6

source /home/jpet/virtplan/bin/activate

# config:
getPlanDIR=/home/jpet/scripts/getPlan/
patientsDIR=/home/jpet/patientData/
CTfilename=CT_ExtROICrop_2.5x2.5x2.5.mhd
#CTfilename=CT_2.5x2.5x2.5.mhd
#simLabel group (configuration)
simLABELgroup=Sim1000
GateSimConf=/home/jpet/GateSimConf/GateSimConf_${simLABELgroup}/


#simLABEL=Sim1000
simLABEL=${simLabelPARAM}
CPUs=${ntasksPARAM}


#cd ${patientsDIR}
#for patientID in */; do
#    echo "$patientID"
##    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} ${patientsDIR}${patientID} ${CTfilename} ${simLABEL} -n 1E4 -Nunit primNoPB -core 1000
##    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} ${patientsDIR}${patientID} ${CTfilename} ${simLABEL} -n 0.1 -Nunit primPr -core 1000
##    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l ${simLABEL} -n 0.1 --Nunit primPr --CPUs 1000 -a prodActor -a doseAll -a doseP -a LETdAll -a LETdP -x0 -y0 -z0
#    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l ${simLABEL} -n 0.1 --Nunit primPr --CPUs ${CPUs} -a prodAll -x${xx} -y${yy} -z${zz}
#done


echo "Scanning input file"

input="patientsListSPB.txt"
while IFS= read -r line
do
    patientID=${line}
    echo ${patientID}
#####    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l ${simLABEL} -n 0.1 --Nunit primPr --CPUs ${CPUs} -a prodAll -x${xx} -y${yy} -z${zz} -m${HU}

    # 10% protons
    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l ${simLABEL} -n 0.1 --Nunit primPr --CPUs ${CPUs} -a prodAll -a doseAll -x${xx} -y${yy} -z${zz} -m${HU}
    
    # 100% protons
#    python3 ${getPlanDIR}getPlanGate_CCB.py -c${GateSimConf} -p${patientsDIR}${patientID} -t${CTfilename} -l ${simLABEL} -n 1.0 --Nunit primPr --CPUs ${CPUs} -a prodAll -a doseAll -x${xx} -y${yy} -z${zz} -m${HU}

done < "$input"

deactivate

