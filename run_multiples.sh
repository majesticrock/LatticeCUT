#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 name value1 [value2 ...]"
    exit 1
fi

name=$1
shift
values=("$@")

config_file="params/tc.config"

if [ ! -f "$config_file" ]; then
    echo "Error: Configuration file '$config_file' not found!"
    exit 1
fi

# Loop over each value in the values array
for value in "${values[@]}"; do
    # Replace the value corresponding to the key == name in the params.config file
    # Use sed to find the line containing the key and replace its value
    sed -i "s/^$name.*/$name $value/" "$config_file"

    # Execute the script ./exec.sh
    ./build/latticecut $config_file
done
