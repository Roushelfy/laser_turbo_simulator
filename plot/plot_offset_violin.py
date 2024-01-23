import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
from cycler import cycler
import matplotlib.colors as mcolors
## Global configuration

# constrained_layout automatically adjusts subplots and decorations like legends and colorbars
# so that they fit in the figure window while still preserving, as best they can, the logical
# layout requested by the user. constrained_layout is similar to tight_layout, but uses a constraint
# solver to determine the size of axes that allows them to fit.

linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as (0, (1, 1)) or ':'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')]  # Same as '-.'

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),
     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

# ['#3853a4', '#146533', '#ed1f24', '#708191', '#faa51a', '#b9519f']
# ['#F27970', '#BB9727',  '#54B345', '#32B897', '#05B9E2', '#8983BF']
palette = ['#EE6677', '#4477AA', '#8ECFC9', '#FFBE7A', '#BEB8DC', '#E7DAD2']  # Genshin by Yuhan
patterns = ["", "/" , "\\" , "x", ".", "o"] #  "|" , "-" , "+" , "x", "o", "O", ".", "*" ]
linestyles = ['-', '--', ':', '-.', (0, (3, 1, 1, 1, 1, 1)), (5, (10, 3))]
markers = ['v', '*', '.', 's', '1', 'x']

plt.rcParams["figure.constrained_layout.use"] = True

plt.rcParams['figure.figsize'] = [4.0, 3.0]
plt.rcParams['figure.dpi'] = 80
plt.rcParams['savefig.dpi'] = 100

plt.rcParams['font.size'] = 18
plt.rcParams["font.family"] = "Arial"
plt.rcParams['legend.fontsize'] = 'medium'
plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['legend.edgecolor'] = 'white'
plt.rcParams['legend.framealpha'] = 0.9
plt.rcParams['figure.titlesize'] = 'medium'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['axes.prop_cycle'] = cycler(color=palette) + cycler(linestyle=linestyles)


transparency = 0.5

def get_transparent_color(color):
    c = mcolors.hex2color(color)
    c = tuple(map(lambda x: x*transparency + (1.0-transparency), mcolors.hex2color(c)))
    hex_color = '#{:02X}{:02X}{:02X}'.format(int(c[0] * 255), int(c[1] * 255), int(c[2] * 255))
    return hex_color

# `colors` for border lines and `colors_fill` for filled areas
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colors_fill = list(map(get_transparent_color, colors))
# 假设的速度列表
speeds = [10]

# 假设的距离列表
distances = [100,200,300,400]

for speed in speeds:
    offset_data = []
    for dis in distances:
        file_path = f'../data/dis_{dis}/speed_{speed}/offset.txt'
        if os.path.exists(file_path):
            with open(file_path, 'r') as file:
                data = [float(line.strip())*100 for line in file]
                #only append positive data
                offset_data.append(data)
    # 绘制小提琴图
    plt.figure()
    plt.violinplot(offset_data)
    plt.xlabel('Distance(m)')
    plt.ylabel('offset(cm)')
    plt.xticks(range(1, len(distances) + 1), distances)

 # 保存图表
    plt.savefig(f'../data/speed_{speed}_offset_violin_withbeam.png')
    plt.close()

    #绘制箱型图
    plt.figure()
    bplot=plt.boxplot(offset_data,patch_artist=True, boxprops={
                              'edgecolor': colors[0],
                              'facecolor': colors_fill[0],
                              'linewidth': 2,
                              'hatch': patterns[0],
                          },
                          medianprops={ 'color': colors[0], 'linewidth': 2 },
                          whiskerprops={ 'color': colors[0], 'linewidth': 2 },
                          capprops={ 'color': colors[0], 'linewidth': 2 },showfliers=False)
    plt.xlabel('Distance(m)')
    plt.ylabel('offset(cm)')
    #设置为虚线
    plt.grid(axis='y', linestyle='--')
    plt.xticks(range(1, len(distances) + 1), distances)

    # 保存图表
    plt.savefig(f'../data/speed_{speed}_offset_box_withbeam.pdf')
    plt.close()

   