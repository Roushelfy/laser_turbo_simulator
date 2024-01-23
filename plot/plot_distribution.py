import matplotlib.pyplot as plt

def read_section(file, start_line,end_line):
    data = {}
    for line in file:
        if line.strip() == start_line:
            break
    for line in file:
        if line.strip() == end_line:
            break
        key, value = line.split(':')
        data[float(key)] = float(value)
        print(key, value)
    return data

def plot_data(data, title, filename):
    keys = list(data.keys())
    values = [data[k] for k in keys]
    if title == 'Laser Offset Distribution':
        keys = [k * 100 for k in keys]
    plt.figure(figsize=(10, 10))
    plt.bar(keys, values, width=0.8 * (keys[1] - keys[0]))
    plt.title(title)
    if title == 'Laser Offset Distribution':
        plt.xlabel('Distance (cm)')
    else:
        plt.xlabel('Intensity(normalized)')
    plt.ylabel('Frequncy')
    plt.savefig(filename)  # 保存图表为图片文件

def main():
    with open('../data/record.txt', 'r') as file: 
        distance_data = read_section(file, 'laser_offset:', 'tag_intensity:')
        file.seek(0)
        intensity_data = read_section(file, 'tag_intensity:','')

    plot_data(distance_data, 'Laser Offset Distribution', 'laser_offset_distribution.pdf')
    plot_data(intensity_data, 'Tag Intensity Distribution', 'tag_intensity_distribution.pdf')

if __name__ == '__main__':
    main()
