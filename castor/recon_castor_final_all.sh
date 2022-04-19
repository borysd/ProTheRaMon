#!/bin/bash
#SBATCH -J castorRecon
#SBATCH --ntasks=10
#SBATCH -p quanta
#SBATCH -t 24:00:00
#SBATCH --output=/home/jpet/scripts/castor/slurmOut/recon_%j.o
#SBATCH --error=/home/jpet/scripts/castor/slurmOut/recon_%j.e

#module load castor/3.1_JPET
module load castor/3.1-tmerlin

#attSource=K
attSource=J

#patientID=Exp_GTR4_SOBP_S2F5_2Gy
patientID=$1

source /home/jpet/venv38/bin/activate

xyz=$(python3 /home/jpet/scripts/castor/xyzRead.py -P ${patientID} 2>&1)
#xyz=$(python3 /home/jpet/scripts/castor/xyzReadNoCrop.py -P ${patientID} 2>&1)

lenxyz=$((${#xyz}-2))
xyz_=`echo ${xyz:1:${lenxyz}}`
IFS=',' read -ra xxyyzz <<< `echo ${xyz_} | tr -d '[:space:]'`
x=${xxyyzz[0]}
y=${xxyyzz[1]}
z=${xxyyzz[2]}

echo ${xyz_}

echo ${x}
echo ${y}
echo ${z}


#geomName=barrel_${attSource}
#geomName=barrel2_${attSource}
##geomName=barrel3_${attSource}
#geomName=dualhead_1x6_${attSource}
#geomName=dualhead_2x6_${attSource}
#geomName=dualhead_3x4_${attSource}

#simLabel=Sim2100

#for simnum in {2600..2614..1}
#for simnum in {2720..2780..20}
for simnum in {2200..2214..1}
do
    echo $simnum
    case $simnum in
	2100|2101|2102|2103|2104|2105|2106|2107|2108|2109|2110|2111|2112|2113|2114) geomName=barrel_${attSource};;
	2200|2201|2202|2203|2204|2205|2206|2207|2208|2209|2210|2211|2212|2213|2214) geomName=barrel2_${attSource};;
	2300|2301|2302|2303|2304|2305|2306|2307|2308|2309|2310|2311|2312|2313|2314) geomName=barrel3_${attSource};;
	2400|2401|2402|2403|2404|2405|2406|2407|2408|2409|2410|2411|2412|2413|2414) geomName=dualhead_1x6_${attSource};;
	2500|2501|2502|2503|2504|2505|2506|2507|2508|2509|2510|2511|2512|2513|2514) geomName=dualhead_2x6_${attSource};;
#	2600|2601|2602|2603|2604|2605|2606|2607|2608|2609|2610|2611|2612|2613|2614) geomName=dualhead_3x4_Field1_${attSource};;
	2600|2601|2602|2603|2604|2605|2606|2607|2608|2609|2610|2611|2612|2613|2614) geomName=dualhead_3x4_${attSource};;
	2700|2720|2740|2760|2780) geomName=sensitivity_lm_dualhead_ccb_${attSource};;
	*) geomName=dualhead_3x4_${attSource}
    esac
    echo $geomName

    simLabel=Sim${simnum}

    pathToPat=/home/jpet/patientData/${patientID}
    pathToSim=/home/jpet/patientData/${patientID}/gate/${simLabel}

    #castor-recon -df acquisition_lm_df.Cdh -dim x,y,z -vox x,y,z -sens sens_atten_calib_it1.hdr -th 0 -it 20:1 -dout recon
##    castor-recon -df ${pathToSim}/patient_lm_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -th 0 -it 20:1 -dout ${pathToPat}/recon/${simLabel} -ignore-corr fdur

    castor-recon -df ${pathToSim}/patient_lm_230_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -it 12:1 -dout ${pathToPat}/recon/${simLabel}_${attSource}_230 -ignore-corr fdur
    castor-recon -df ${pathToSim}/patient_lm_500_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -it 12:1 -dout ${pathToPat}/recon/${simLabel}_${attSource}_500 -ignore-corr fdur
    castor-recon -df ${pathToSim}/patient_lm_700_df.Cdh -dim ${x},${y},${z} -vox 2.5,2.5,2.5 -sens ${pathToPat}/TPS/sens_atten_${geomName}_it1.hdr -it 12:1 -dout ${pathToPat}/recon/${simLabel}_${attSource}_700 -ignore-corr fdur

done
