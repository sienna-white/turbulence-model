#!/bin/bash
#SBATCH --nodes 1
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J pmax_v_ws
#SBATCH -t 00:40:00

##### SBATCH -account m1266

module load PrgEnv-gnu
module load python

out="pressure_vs_ws.csv"
rm $out


#######################################

linspace() { # 1=start / 2=end / 3=points
    step=$(bc <<< "scale=10; ($2 - $1) / ($3 - 1)")
    echo $step
    seq "$1" "$step" "$2"
}

points=30
############# PRESSURE #############
start=0.0000001 
end=0.001 

pressure=$(linspace "$start" "$end" "$points")

############# WS ##################
start=-0.0000001 
end=0.0000001  

ws0=$(linspace "$start" "$end" "$points")

#######################################

echo "ws ="
echo $ws0

echo "pressure = "
echo $pressure

cd .. 

# Loop through each value in the vector and pass it to the function
for ws00 in $ws0; do
    for px0 in $pressure; do
    python watercolumn-sw.py $out $ws00 $px0
    sleep 5
    done
done