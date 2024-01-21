import os
import matplotlib.pyplot as plt
import numpy as np
# Function to parse data from a given file path
def parse_data_from_file(file_path):
    object_speeds = []
    average_offsets = []
    average_intensities = []

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            for part in parts:
                key, value = part.split(':')
                if key == 'object_speed':
                    object_speeds.append(int(value))
                elif key == 'average_offset':
                    average_offsets.append(float(value))
                elif key == 'average_intensity':
                    average_intensities.append(float(value))

    return np.array(object_speeds), np.array(average_offsets), np.array(average_intensities)

# Function to plot and save data
def plot_and_save_data(speed, directory, suffixes, labels, colors, markers):
    # Dictionary to store parsed data
    parsed_data = {}

    # Parse and store data from both 'wbeam' and 'wobeam' files
    for suffix in suffixes:
        file_name = f"speed_{suffix}_{speed}.txt"
        file_path = os.path.join(directory, file_name)
        parsed_data[suffix] = parse_data_from_file(file_path)

    # Determine global maximum intensity for normalization across both datasets
    max_intensity = max(
        max(parsed_data['wbeam'][2]), 
        max(parsed_data['wobeam'][2])
    )

    # Plotting
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    for suffix, label, color, marker in zip(suffixes, labels, colors, markers):
        object_speeds = parsed_data[suffix][0]
        average_offsets = parsed_data[suffix][1] * 100  # Assuming the offset is given in meters, convert to cm
        average_intensities = parsed_data[suffix][2] / max_intensity  # Normalize the intensities

        # First plot (Normalized power)
        axs[0].plot(object_speeds, average_intensities, color=color, marker=marker, linestyle='-', label=label)

        # Second plot (Offset in cm)
        axs[1].plot(object_speeds, average_offsets, color=color, marker=marker, linestyle='-', label=label)

    # Setting labels, titles, and legends
    axs[0].set_xlabel('Velocity (m/s)')
    axs[0].set_ylabel('Normalized power')
    axs[0].set_title('Power vs. Velocity')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].set_xlabel('Velocity (m/s)')
    axs[1].set_ylabel('Offset (cm)')
    axs[1].set_title('Offset vs. Velocity')
    axs[1].grid(True)
    axs[1].legend()

    plt.tight_layout()

    # Save the figure
    plt.savefig(os.path.join(directory, f"combined_plot_speed_{speed}.png"))
    plt.close()

# Assuming the directory "../data" contains the text files.
data_directory = "../data"  # This is a placeholder path
speeds = [100, 200, 300, 400, 500, 600, 700, 800]  # Replace with actual speed values
suffixes = ['wbeam', 'wobeam']  # Suffixes for files with and without beamwidth adjustment
labels = ['Beamwidth adjustment', 'w/o Beamwidth adjustment']  # Labels for the plot lines
colors = ['blue', 'red']  # Colors for the plot lines
markers = ['o', 's']  # Markers for the plot lines

# Iterate over each speed and generate plots
for speed in speeds:
    plot_and_save_data(speed, data_directory, suffixes, labels, colors, markers)
