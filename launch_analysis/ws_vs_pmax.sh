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

out="ws_vs_pmax.csv"
rm $out


#######################################

linspace() { # 1=start / 2=end / 3=points
    step=$(bc <<< "scale=10; ($2 - $1) / ($3 - 1)")
    echo $step
    seq "$1" "$step" "$2"
}

points=30
############# PMAX #############
start=0.0005 
end=0.5
pmax=$(linspace "$start" "$end" "$points")
echo $pmax


############# WS ##################
start=-0.0000001 
end=0.0000001  
ws=$(linspace "$start" "$end" "$points")
echo $ws

#######################################

echo "ws ="
echo $ws

echo "pmax = "
echo $pmax

cd .. 

# Loop through each value in the vector and pass it to the function
for ws0 in $ws; do
    for pmax0 in $pmax; do
    python watercolumn-sw.py $out $ws0 $pmax0
    sleep 5
    done
done