#!/bin/bash

patID=$1
#patID=G4P133E0


#./split_and_run_PET_PAT_Sim.sh G3P052E0 Sim2111 quanta 200

partition=quanta
#partition=ibm
ntasks=200

for simnum in {2600..2600..1}
#for simnum in {2100..2600..100}
do
    simLabel=Sim${simnum}

#    ./split_and_run_PET_PAT_Sim.sh ${patID} ${simLabel} ${partition} ${ntasks}
    ./split_and_run_PET_PAT_Sim_quanta.sh ${patID} ${simLabel} ${partition} ${ntasks}
done
