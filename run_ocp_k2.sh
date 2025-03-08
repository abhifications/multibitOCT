#!/bin/bash

# Define the range of switching times
start=0.1
end=3.0
increment=0.1

# Path to the input file and executable
input_file="/mnt/c/Users/abhin/Desktop/fortran_code_Abhinav2/Input_File.in"
executable="/mnt/c/Users/abhin/Desktop/fortran_code_Abhinav2/OCP"

# Loop over switching times
for time in $(seq -f "%.1f" $start $increment $end); do
    # Convert to Fortran-style double precision (e.g., "1.0d0")
    time_fortran=$(printf "%.1fd0" $time)

    # Create a folder for this switching time
    folder=$(printf "%.1fd0" $time)
    mkdir -p "$folder"

    # Copy the input file to the folder
    cp "$input_file" "$folder/Input_File.in"

    # Modify the input file to update switching_time
    sed -i "s/^switching_time .*/switching_time $time_fortran/" "$folder/Input_File.in"

    # Run the executable
    (cd "$folder" && "$executable" > output.log)

    echo "Completed for switching_time = $time_fortran"
done

echo "All runs completed!"
