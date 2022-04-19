#!/bin/bash

simLabelPARAM=$1
partPARAM=$2
patientPARAM=$3



echo "#!/bin/bash

#SBATCH -J prodPlan_${simLabelPARAM}
#SBATCH --ntasks=10
#SBATCH -p ${partPARAM}
#SBATCH -t 08:00:00
#SBATCH --output=/home/jpet/scripts/prodPlan/slurmout/prod_%j.o
#SBATCH --error=/home/jpet/scripts/prodPlan/slurmout/prod_%j.e

source /home/jpet/virtplan/bin/activate

# config:
prodPlanDIR=/home/jpet/scripts/prodPlan/
patientsDIR=/home/jpet/patientData/
CTfilename=CT_ExtROICrop_2.5x2.5x2.5.mhd
#CTfilename=CT_2.5x2.5x2.5.mhd
#simLabel group (configuration)
simLABELgroup=Sim1000
GateSimConf=/home/jpet/GateSimConf/GateSimConf_\${simLABELgroup}/


simLABEL=${simLabelPARAM}

# PatientList=


patientID=${patientPARAM}

cd \${patientsDIR}

#for patientID in */; do
    echo \"\$patientID\"

#    python3 \${prodPlanDIR}prodPlan.py -p\${patientsDIR}\${patientID} -t\${CTfilename} -l\${simLABEL}
#    python3 \${prodPlanDIR}prodPlan_29032021.py -p\${patientsDIR}\${patientID} -t\${CTfilename} -l\${simLABEL}
    python3 \${prodPlanDIR}prodPlan.py -p\${patientsDIR}\${patientID} -t\${CTfilename} -l\${simLABEL}
    
#done

deactivate" | sbatch
