import numpy as np
import matplotlib.pyplot as plt
from processLVIS import lvisGround
import argparse
import os
import glob

def get_file_bounds(filename):
    """ Function for reading the bounds of a file"""
    lvis = lvisGround(filename, onlyBounds=True)
    return lvis.bounds

def is_within_bounds(bounds, minX, minY, maxX, maxY):
    """ Function to check if the bounds intersect with the specific region"""
    return not(bounds[2] < minX or bounds[0] > maxX or bounds[3] < minY or bounds[1] > maxY)

def select_files(file_list, minX, minY, maxX, maxY):
    selected_files = []
    for filename in file_list:
        bounds = get_file_bounds(filename)
        print(f"File: {filename}, Bounds: {bounds}")  # Debugging statement
        if is_within_bounds(bounds, minX, minY, maxX, maxY):
            selected_files.append(filename)
    return selected_files

def plot_waveform(filename, index):
    """ Function to plot waveform at given index """
    # Read the LVIS data
    lvis = lvisGround(filename)
    
    # Set elevations to ensure the 'z' attribute is available
    lvis.setElevations()
    
    # Get the waveform at the specified index
    z, waveform = lvis.getOneWave(index)
    
    # Plot the waveform
    plt.figure()
    plt.plot(z, waveform)
    plt.title(f'Waveform at Index {index}')
    plt.xlabel('Elevation (m)')
    plt.ylabel('Amplitude')
    plt.show()

if __name__ == "__main__":
    # Parser set up 
    parser = argparse.ArgumentParser(description='Select LVIS files covering Pine Island Glacier and plot a waveform')
    parser.add_argument('folder_path', type=str, help='Path to the folder containing LVIS files')
    parser.add_argument('index', type=int, help='Index of the waveform to plot')

    # Parse arguments
    args = parser.parse_args()

    # Get list of all files in the specified folder
    file_list = glob.glob(os.path.join(args.folder_path, '*.h5'))

    
    x0 = -100.0
    y0 = -75.0
    x1 = -90.0
    y1 = -70.0

    # Select files covering Pine Island Glacier
    selected_files = select_files(file_list, x0, y0, x1, y1)

    # Print selected files
    for file in selected_files:
        print(f"Selected file: {file}")
    
    # Plot the waveform from the first selected file
    plot_waveform(selected_files[0], args.index)
    