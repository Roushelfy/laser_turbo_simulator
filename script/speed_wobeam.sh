#!/bin/bash

# Run simulator with different speed

speed="10 50 100 150 200 250 300 350 400 450 500"

distance="50 100 150 200 250 300 350 400 450 500 550 600"

target_string="lost"

for j in $distance
    do
       for i in $speed
        do
            echo "Distance: $j"
            echo "Speed: $i"
            #../build/simulate -object_speed=$i -object_distance=$j >> ../data/speed.txt
            attempt=0
            while ((attempt < 10)); do
                radius=$(echo "0.01 * $j" | bc)
                output=$(../build/simulate -object_speed=$i -object_distance=$j -object_moving_radius=$radius)
                echo "$output" | grep -q "$target_string"
                if [ $? -ne 0 ]; then
                    #echo "$output" >> ../data/speed_wobeam_$j.txt
                    break
                fi
                ((attempt++))
            done
        done
    done
