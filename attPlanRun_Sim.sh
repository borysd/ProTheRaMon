#!/bin/bash

Patient=$1

source /home/jpet/virtplan/bin/activate

# config:
attPlanDIR=/home/jpet/scripts/attPlan/
patientsDIR=/home/jpet/patientData/
CTfilename=CT_ExtROICrop_2.5x2.5x2.5.mhd
#CTfilename=CT_2.5x2.5x2.5.mhd

patientID=${Patient}
echo ${patientID}
echo ${patientsDIR}
python3 ${attPlanDIR}attPlan.py -p${patientsDIR}${patientID} -t${CTfilename}

deactivate

