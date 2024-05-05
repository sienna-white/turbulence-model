


file="pressure_vs_ws.csv" 
export pmax=0.1
output_csv="pressure_vs_ws_pmax=${pmax}.csv"
export output_csv

tail -n +2 "$file" | while IFS=, read -r arg1 arg2 || [ -n "$arg1" ]; do
    # Print the two values to the script
    echo "pressure: $arg1"
    echo "ws: $arg2"
    export pressure=$arg1
    export ws=$arg2
    python ../watercolumn-sw.py 
    wait 

done 