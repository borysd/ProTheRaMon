#!/bin/bash

#process converted by Castor file and include calibration factor


input=$1
output=$1".tmp"
simLabel=$2

# PUT your own values !!!!
case $simLabel in
    2100|2101|2102|2103|2104|2105|2106|2107|2108|2109|2110|2111|2112|2113|2114) calibrationF="1.31748E-08";;
    2200|2201|2202|2203|2204|2205|2206|2207|2208|2209|2210|2211|2212|2213|2214) calibrationF="1.476252E-08";;
    2300|2301|2302|2303|2304|2305|2306|2307|2308|2309|2310|2311|2312|2313|2314) calibrationF="1.497636E-08";;
    2400|2401|2402|2403|2404|2405|2406|2407|2408|2409|2410|2411|2412|2413|2414) calibrationF="8.395512E-09";;
    2500|2501|2502|2503|2504|2505|2506|2507|2508|2509|2510|2511|2512|2513|2514) calibrationF="1.0235124E-08";;
    2600|2601|2602|2603|2604|2605|2606|2607|2608|2609|2610|2611|2612|2613|2614) calibrationF="6.998856E-09";;
    2700|2701|2702|2703|2704|2705|2706|2707|2708|2709|2710|2711|2712|2713|2714|2720|2740|2760|2780) calibrationF="6.998856E-09";;
    *) calibrationF="6.998856E-09"
esac

#read line by line and create temporary file .tmp (replace line with Calibration factor)
while IFS= read -r line
do
    new_line="Calibration factor: ""$calibrationF"

    if [[ "$line" == *"Calibration factor:"* ]]; then
    	echo "$new_line" >> "$output"
    else
	echo "$line" >> "$output"
    fi
done < "$input"

#move tmp file to old input file with the new Calibration factor
mv "$output" "$input"