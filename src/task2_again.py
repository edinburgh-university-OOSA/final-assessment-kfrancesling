import os
import numpy as np 
import argparse # cmd line args
import matplotlib.pyplot as plt # for plotting 
from handleTiff import tiffHandle 
from processLVIS import lvisGround
import rasterio # raster data handling 
from rasterio.merge import merge # 
from glob import glob
from pyproj import Transformer # cordinate tranformation 
import tracemalloc # memory tracking

def norm_lon(lon):
    """
    function to normalise longitude to 0-360 degree range.
    
    arguments:
        lon (float): Longitude in degrees.
        
    returns:
        float: Normalised longitude.
    """
    return lon % 360

def writeDEM(lvis, res, outName):
    """
    function to write LVIS elevation data to a GeoTIFF file
    
    arguments:
        lvis (lvisGround): LVIS data object.
        res (float): Resolution of the DEM.
        outName (str): Output filename for the GeoTIFF.
    """

    tiff = tiffHandle(outName) # creating tiff handle object for writing the geotiff 
    tiff.minX = np.min(lvis.x) # min longitude for the geotiff from reproj xcoords 
    tiff.maxY = np.max(lvis.y) # min latitude for the geotiff from reproj ycoords
    tiff.res = res # defining spatial res
    tiff.nX = int((np.max(lvis.x) - np.min(lvis.x)) / res) + 1 # calc number of pixels in x dir (+1 to cover edges)
    tiff.nY = int((np.max(lvis.y) - np.min(lvis.y)) / res) + 1 # calc number of pixels in y dir (+1 to cover edges)

    # initiaise DEM grid -- with no-data value
    dem_grid = np.full((tiff.nY, tiff.nX), -999.0) # numpy array with nodata val
    x_indices = ((lvis.x - tiff.minX) / res).astype(int) # converting xcoords to column indicies 
    y_indices = ((tiff.maxY - lvis.y) / res).astype(int) # converting ycoords to row indicies

    # ensure that indicies are within grid bounds to prevent errors
    valid_x = (x_indices >= 0) & (x_indices < tiff.nX)
    valid_y = (y_indices >= 0) & (y_indices < tiff.nY)
    valid_indices = valid_x & valid_y # combine x + y to ensure flll row-colum validity
    
    # valid data into DEM grid
    dem_grid[y_indices[valid_indices], x_indices[valid_indices]] = lvis.zG[valid_indices]

    # write this to GeoTiff
    tiff.writeTiff(dem_grid, lvis.x, lvis.y, res, filename=outName, epsg=3031)
    print("GeoTIFF written to", outName)

def merge_tifs(output_dir, merged_output):
    """
    function to merge individual GeoTIFF files into one.
    
    arguments:
        output_dir (str): Directory containing individual GeoTIFFs.
        merged_output (str): Filename for the merged GeoTIFF.
    """
    input_tifs = glob(os.path.join(output_dir, "*.tif")) # using glob to create list of TIFFs

    # check for files
    if not input_tifs:
        print("No files found to merge.")
        return

    print(f"Merging {len(input_tifs)} files...")
    src_files_to_mosaic = [rasterio.open(tif) for tif in input_tifs] # opening files creating rasterio objs
    mosaic, out_trans = merge(src_files_to_mosaic) # merge the files

    out_meta = src_files_to_mosaic[0].meta.copy() # prepare metadata
    out_meta.update({ # update metadata
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "compress": "LZW"
    })
    # write the merged file
    with rasterio.open(merged_output, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Merged file saved as {merged_output}")

def plot_dem(file_path):
    """
    function to plot the DEM using matplotlib.
    
    arguments:
        file_path (str): path to the GeoTIFF file to plot.
    """
    with rasterio.open(file_path) as src:
        dem = src.read(1)  # read the first band
        extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)

    plt.figure(figsize=(10, 8))
    plt.imshow(dem, cmap='terrain', extent=extent)
    plt.colorbar(label='Elevation (m)')
    plt.title('DEM')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid(False)
    plt.show()

def get_bounds(lvis_object):
    """
    function to get spatial bounds of the LVIS data.
    
    arguments:
        lvis_object (lvisGround): LVIS data object.
        
    returns:
        tuple: min and max longitude and latitude.
    """
    return np.min(lvis_object.lon), np.min(lvis_object.lat), np.max(lvis_object.lon), np.max(lvis_object.lat)

def calculate_overlap(file_bounds, pig_bounds):
    """
    function to calculate overlap area between file bounds and PIG bounds to find the best file.
    
    arguments:
        file_bounds (tuple): min/max longitude and latitude of the file.
        pig_bounds (tuple): min/max longitude and latitude of the PIG.
        
    return:
        float: the overlapping area
    """
    overlap_lon = max(0, min(file_bounds[2], pig_bounds[2]) - max(file_bounds[0], pig_bounds[0])) # overlap in longitiude 
    overlap_lat = max(0, min(file_bounds[3], pig_bounds[3]) - max(file_bounds[1], pig_bounds[1])) # overlap in latitude
    return overlap_lon * overlap_lat # calculate the total area of overlap

def process_pig_files(directory, resolution, output_dir):
    """
    function to process LVIS files to create a DEM covering the PIG.
    
    arguments:
        directory (str): Directory containing LVIS files.
        resolution (float): Resolution of the DEM.
        output_dir (str): Directory to save the GeoTIFF files.
    """
    pig_bounds = (norm_lon(-102.0), -75.4, norm_lon(-99.0), -74.6) # PIG bounds
    best_file = None
    best_overlap = 0 # preparing to track largest overlap with PIG

    tracemalloc.start()

    # file iteration
    for filepath in glob(os.path.join(directory, "*.h5")):
        b = lvisGround(filepath, onlyBounds=True)
        file_bounds = get_bounds(b)

        overlap = calculate_overlap(file_bounds, pig_bounds) # calculate geo overlap
        if overlap > best_overlap: 
            best_file = filepath
            best_overlap = overlap

    # process best file
    if best_file:
        print("Best file selected:", best_file)
        b = lvisGround(best_file, onlyBounds=True)
        b_min_lon, b_min_lat, b_max_lon, b_max_lat = get_bounds(b) # retrive bounds 

        # divide file into chunks for processing (based on range)
        chunk_size_x = (b_max_lon - b_min_lon) / 20
        chunk_size_y = (b_max_lat - b_min_lat) / 20 

    # process each chunk
        for x0 in np.arange(b_min_lon, b_max_lon, chunk_size_x):
            x1 = x0 + chunk_size_x # loop over long chunks
            for y0 in np.arange(b_min_lat, b_max_lat, chunk_size_y):
                y1 = y0 + chunk_size_y # loop over lat chunks 
                # initialise lvis for chunk
                lvis = lvisGround(best_file, minX=x0, minY=y0, maxX=x1, maxY=y1, setElev=True)
                # check for data within chunk 
                if lvis.nWaves > 0:
                    lvis.setElevations()
                    lvis.estimateGround()
                    # reproject coords to sys for antartica 
                    transformer = Transformer.from_crs("epsg:4326", "epsg:3031", always_xy=True)
                    lvis.x, lvis.y = transformer.transform(lvis.lon, lvis.lat)
                    # write tiff
                    outName = f"{output_dir}/DEM_x_{x0}_y_{y0}.tif"
                    writeDEM(lvis, resolution, outName)
        # memory usage reporting - reporting most used lines 
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')

        print("[ Top 10 memory usage ]")
        for stat in top_stats[:10]:
            print(stat)
        # merge the geotiffs
        merged_output = f"{output_dir}/Merged.tif"
        merge_tifs(output_dir, merged_output)

        # Plot the merged DEM
        plot_dem(merged_output)

if __name__ == "__main__":
    """
    Main execution block to parse arguments and process files.
    """
    parser = argparse.ArgumentParser(description=("Create DEM from LVIS data"))
    parser.add_argument("directory", type=str, help="Directory containing LVIS files")
    parser.add_argument("resolution", type=float, help="Resolution of the DEM in meters")
    parser.add_argument("output_dir", type=str, help="Directory to save the merged GeoTIFF file")
    args = parser.parse_args()

    process_pig_files(args.directory, args.resolution, args.output_dir)