#!/bin/bash


echo "Scanning input file"

input="patientsListPET.txt"

while IFS= read -r line
do
    patientID=${line}
#    echo ${patientID}

    ./split_and_run_PET_PAT_master.sh ${patientID}

done < "$input"
