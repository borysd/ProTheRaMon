#!/bin/bash


echo "Scanning input file"

#input="patientsListPET.txt"
input="patientsListSPB.txt"

while IFS= read -r line
do
    patientID=${line}
    echo ${patientID}

    ./attPlanRun_Sim.sh ${patientID}

done < "$input"
