#!/bin/bash

#test pid

kp=$(seq 0.4 0.01 1.3)
kd=$(seq 0.7 0.01 1.4)

target_string="lost"

try_number=5
max_background_processes=8  # 控制同时运行的后台进程数量
current_background_processes=0

for j in $kd; do
    for i in $kp; do
        (
            attempt=0
            success=0
            start_time=$(date +%s.%N)  # 开始时间
            while ((attempt < try_number)); do
                output=$(../build/simulate -kp=$i -kd=$j)
                echo "$output" | grep -q "$target_string"
                if [ $? -ne 0 ]; then
                    ((success++))
                fi
                attempt=$((attempt + 1))
            done
            end_time=$(date +%s.%N)   # 结束时间
            duration=$(echo "$end_time - $start_time" | bc)  # 计算持续时间
            echo "Kp: $i Kd: $j Success: $success Time: $duration" >> pid_result.txt
        ) &
        ((current_background_processes++))

        if ((current_background_processes >= max_background_processes)); then
            wait -n
            ((current_background_processes--))
        fi
    done
done

wait

sort -k6,6nr -k8,8nr pid_result.txt -o sorted_pid_result.txt
