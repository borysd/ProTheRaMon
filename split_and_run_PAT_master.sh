#!/bin/bash

patID=$1
#patID=G4P024E0
#patID=G3P120E0
#patID=G3P013E1
#patID=Exp_GTR4_SPB_R115_8Gy


#./split_and_run_PAT_Sim.sh G3P052E0 Sim2111 quanta 200

partition=quanta
#partition=ibm

ntasks=400
#ntasks=600

for simnum in {1000..1000}
do
    simLabel=Sim${simnum}

    ./split_and_run_PAT_Sim.sh ${patID} ${simLabel} ${partition} ${ntasks}
done
