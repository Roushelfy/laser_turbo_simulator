import matplotlib.pyplot as plt
import numpy as np
import os

# 假设的速度列表
speeds = [10]

# 假设的距离列表
distances = [100,200,300,400]

for speed in speeds:
    snr_data = []
    for dis in distances:
        file_path = f'../data/dis_{dis}/speed_{speed}/snr.txt'
        if os.path.exists(file_path):
            with open(file_path, 'r') as file:
                data = [float(line.strip()) for line in file]
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
    plt.savefig(f'../data/speed_{speed}_snr_violin_withbeam.pdf')
    plt.close()

    #绘制箱型图
    plt.figure()
    bplot=plt.boxplot(snr_data, patch_artist=True, showfliers=False)
    plt.title(f'Speed {speed} - SNR Distribution')
    plt.xlabel('Distance(m)')
    plt.ylabel('SNR')
    #设置为虚线
    plt.grid(axis='y', linestyle='--')
    plt.xticks(range(1, len(distances) + 1), distances)

        # 设置所有箱体的颜色
    box_color = '#FF9999'  # 箱体颜色
    for patch in bplot['boxes']:
        patch.set_facecolor(box_color)

    # 设置中位数线的颜色
    median_color = 'blue'
    for median in bplot['medians']:
        median.set_color(median_color)

    # 保存图表
    plt.savefig(f'../data/speed_{speed}_snr_box_withbeam.pdf')
    plt.close()

   