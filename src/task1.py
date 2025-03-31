import numpy as np
import matplotlib.pyplot as plt
from processLVIS import lvisGround
import argparse
import os

def plot_waveform(filename, index):
    """
    Plots the waveform from an LVIS file at a specified index.

    Parameters:
    filename (str): Path to the LVIS file.
    index (int): Index of the waveform to plot.

    Raises:
    FileNotFoundError: If the LVIS file does not exist.
    IndexError: If the index is out of range.
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file {filename} does not exist.")
    
    # Read the LVIS data
    lvis = lvisGround(filename)
    
    # Set elevations to ensure the 'z' attribute is available
    lvis.setElevations()
    
    # Check if the index is within the valid range
    if index < 0 or index >= len(lvis.z):
        raise IndexError(f"Index {index} is out of range. Valid range is 0 to {len(lvis.z) - 1}.")
    
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
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Plot LVIS waveform')
    parser.add_argument('filename', type=str, help='Path to the LVIS file')
    parser.add_argument('index', type=int, help='Index of the waveform to plot')
    
    # Parse arguments
    args = parser.parse_args()
    
    try:
        # Plot the waveform
        plot_waveform(args.filename, args.index)
    except (FileNotFoundError, IndexError) as e:
        print(f"Error: {e}")