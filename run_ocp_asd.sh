#!/bin/bash

# -------------------------------------------------------------
# Configuration
# -------------------------------------------------------------

# Define the 10 switching times as "n switching_time_in_ns"
# (n is just a label; the switching time is used in the input file)
switching_times=(
  "1 1.18281"
  "2 2.36561"
  "3 3.54842"
  "4 4.73122"
  "5 5.91403"
  "6 7.09683"
  "7 8.27964"
  "8 9.46244"
  "9 10.6453"
  "10 11.8281"
)

# Path to the input file and the executable
input_file="/mnt/c/Users/abhin/Desktop/fortran_code_Abhinav2/Input_File.in"
executable="/mnt/c/Users/abhin/Desktop/fortran_code_Abhinav2/OCP"

# Define the final output folder for the copied and renamed outputs.
# (This will hold 20 files: 10 for pulse.in and 10 for path.in)
out_dir="/mnt/c/Users/abhin/Desktop/ASD_ocp_pulses_paths"
mkdir -p "$out_dir"

# -------------------------------------------------------------
# Loop over each switching time and run the executable
# -------------------------------------------------------------
for pair in "${switching_times[@]}"; do
    # Split the pair into the label (n) and the switching time value
    read -r n t_value <<< "$pair"
    
    # Format the switching time in Fortran double precision (e.g., "1.18281d0")
    time_fortran=$(printf "%.5fd0" "$t_value")
    
    # Create a temporary folder for this run.
    # The folder is labeled as "n.0tau_0" (e.g., "1.0tau_0", "2.0tau_0", etc.)
    folder=$(printf "%d.0tau_0" "$n")
    mkdir -p "$folder"
    
    echo "Running for n = $n with switching_time = $time_fortran in folder '$folder'..."
    
    # Copy the input file into the run folder
    cp "$input_file" "$folder/Input_File.in"
    
    # Modify the input file to update the switching_time.
    # (This assumes that the input file has a line starting with "switching_time")
    sed -i "s/^switching_time .*/switching_time $time_fortran/" "$folder/Input_File.in"
    
    # Run the executable within the run folder.
    # The output (stdout) is captured in output.log.
    (cd "$folder" && "$executable" > output.log)
    
    # After the run, check for the expected output files and copy them to the final output folder.
    # Rename External_P.out to "n.0tau_0_pulse.in" and Path_P.out to "n.0tau_0_path.in"
    if [ -f "$folder/External_P.out" ]; then
        cp "$folder/External_P.out" "$out_dir/$(printf "%d.0tau_0_pulse.in" "$n")"
    else
        echo "Warning: External_P.out not found in folder '$folder'."
    fi
    
    if [ -f "$folder/Path_P.out" ]; then
        cp "$folder/Path_P.out" "$out_dir/$(printf "%d.0tau_0_path.in" "$n")"
    else
        echo "Warning: Path_P.out not found in folder '$folder'."
    fi
    
    echo "Completed run for n = $n."
    echo "------------------------------------------------------"
done

echo "All runs completed! Final output files are located in: $out_dir"
