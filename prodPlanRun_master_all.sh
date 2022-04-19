#!/bin/bash


input="/home/jpet/scripts/patientsListSPB.txt"

while IFS= read -r line;
do 
    PATID=$line
    ./prodPlanRun_master.sh $PATID
done <"$input"