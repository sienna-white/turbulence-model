
# cd figures/phase/
# rm *.png 
# cd ../../

file="launch_analysis/phasing_options.csv" 
# file="launch_analysis/temp_vs_wind_options.csv"
export pressure=2e-6
output_csv="tidalphase_vs_windphase_pressure=${pressure}_wind=3.5.csv"
# output_csv="temp_vs_windphase_pressure=${pressure}_wind=0.csv"

rm $output_csv
export output_csv

echo "tidal_phase: $arg1"
echo "wind_phase: $arg2"
export tidal_phase="0"
# export del_temp=$arg1
export wind_phase="0"
python watercolumn-sw_2species.py & 


# tail -n +2 "$file" | while IFS=, read -r arg1 arg2 || [ -n "$arg1" ]; do
#     # Print the two values to the script
#     echo "tidal_phase: $arg1"
#     echo "wind_phase: $arg2"
#     export tidal_phase=$arg1
#     # export del_temp=$arg1
#     export wind_phase=$arg2
#     python watercolumn-sw_2species.py & 

#     # sleep 0.5 

# done 

# &>scriptname.out