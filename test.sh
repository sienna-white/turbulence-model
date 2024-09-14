

file="launch_analysis/phasing_options.csv" 
export pressure=2e-6
output_csv="tidalphase_vs_windphase_pressure=${pressure}_wind=3.5.csv"
export output_csv

tail -n +2 "$file" | while IFS=, read -r arg1 arg2 || [ -n "$arg1" ]; do
    # Print the two values to the script
    echo "tidal_phase: $arg1"
    echo "wind_phase: $arg2"
    export tidal_phase=$arg1
    export wind_phase=$arg2
    python watercolumn-sw_2species.py & 
    # sleep 0.5 

done 

# &>scriptname.out