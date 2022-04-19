#!/bin/bash


patientPARAM=$1
#queuePARAM

#patientPARAM=G3P013E1
#patientPARAM=Exp_Sphere_2F
#slurmQueue=ibm
slurmQueue=quanta

simPath=/home/jpet/patientData/


for i in {1000..1000}
do
    simuLabel=Sim${i}
    echo $simuLabel
    ./prodPlanRun_Sim.sh ${simuLabel} ${slurmQueue} ${patientPARAM}

    finalPath=${simPath}${patientPARAM}/gate/${simuLabel}

    #read subdirectories where Results subdir should exists 
#    for subdir in ${finalPath}/*
#    do
#	echo ${subdir}
#	gzip ${subdir}/*Stop.txt
#	gzip ${subdir}/*Prod.txt
#	echo "done"
#    done

done