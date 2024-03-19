out="pressure_vs_ws.csv"
rm $out
############# PRESSURE #############
start=0.0000001 
end=0.001 
points=10

# Calculate the step size
step=$(bc <<< "scale=10; ($end - $start) / ($points - 1)")

# Generate the vector using seq command
pressure=$(seq "$start" "$step" "$end")
#######################################

############# WS ##################
start=-0.0000001 
end=0.0000001  
points=10

# Calculate the step size
step=$(bc <<< "scale=10; ($end - $start) / ($points - 1)")
echo $step
ws=$(seq "$start" "$step" "$end")
#######################################

echo "ws ="
echo $ws

echo "pressure = "
echo $pressure


# Loop through each value in the vector and pass it to the function
for ws0 in $ws; do
    for px0 in $pressure; do
    python watercolumn-sw.py $out $ws0 $px0
    sleep 10
    done
done