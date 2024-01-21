import matplotlib.pyplot as plt
import numpy as np
import os

# 假设的速度列表
speeds = [10, 50, 100 ,150, 200, 250 ,300, 350, 400, 450 ,500]

# 假设的距离列表
distances = [50 ,100, 150, 200 ,250, 300, 350 ,400 ,450 ,500 ,550 ,600]

for speed in speeds:
    snr_data = []
    for dis in distances:
        file_path = f'../data/dis_{dis}/speed_{speed}/snr.txt'
        if os.path.exists(file_path):
            with open(file_path, 'r') as file:
                data = [float(line.strip()) for line in file if float(line.strip())>0]
                #only append positive data
                snr_data.append(data)
    # 绘制小提琴图
    plt.figure()
    plt.violinplot(snr_data)
    plt.title(f'Speed {speed} - SNR Distribution')
    plt.xlabel('Distance(m)')
    plt.ylabel('SNR')
    plt.xticks(range(1, len(distances) + 1), distances)

    # 保存图表
    plt.savefig(f'../data/speed_{speed}_snr_violin_withbeam.png')
    plt.close()
