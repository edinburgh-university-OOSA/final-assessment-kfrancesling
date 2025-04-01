import tracemalloc  # Tracing memory allocations
import numpy as np
import argparse  # Handle command-line arguments
from pyproj import Transformer  # Coordinate transformations
from handleTiff import tiffHandle
from processLVIS import lvisGround
import rasterio
from rasterio.merge import merge
from glob import glob  # File handling
import os  # OS interactions

tracemalloc.start()

class plotLVIS(lvisGround):
    '''A class to handle LVIS data and write DEMs'''

    def reprojectLVIS(self, outEPSG):
        '''Reproject geolocation data'''
        inProj = "epsg:4326"
        outProj = f"epsg:{outEPSG}"
        transformer = Transformer.from_crs(inProj, outProj, always_xy=True)
        self.x, self.y = transformer.transform(self.lon, self.lat)

    def writeDEM(self, res, outName):
        '''Write LVIS ground elevation data to a GeoTIFF using the tiffHandle class'''
        if self.x is None or self.y is None:
            raise ValueError("ReprojectLVIS must be called before writeDEM to set x and y coordinates.")
        print("Writing GeoTIFF...")

        # Create an instance of tiffHandle
        tiff = tiffHandle(outName)

        # Set the necessary attributes for geotransform
        tiff.minX = np.min(self.x)
        tiff.maxY = np.max(self.y)
        tiff.res = res
        tiff.nX = int((np.max(self.x) - np.min(self.x)) / res) + 1
        tiff.nY = int((np.max(self.y) - np.min(self.y)) / res) + 1

        # Initialize the DEM grid
        dem_grid = np.full((tiff.nY, tiff.nX), -999.0)

        # Fill the grid
        x_indices = ((self.x - tiff.minX) / res).astype(int)
        y_indices = ((tiff.maxY - self.y) / res).astype(int)

        valid_x = (x_indices >= 0) & (x_indices < tiff.nX)
        valid_y = (y_indices >= 0) & (y_indices < tiff.nY)
        valid_indices = valid_x & valid_y

        dem_grid[y_indices[valid_indices], x_indices[valid_indices]] = self.zG[valid_indices]

        # Write the DEM to a GeoTIFF file
        tiff.writeTiff(dem_grid, self.x, self.y, res, filename=outName, epsg=3031)
        print("GeoTIFF written to", outName)

def merge_tifs(input_dir, output_tif):
    '''Merge multiple GeoTIFF files into one'''
    # Use glob to find tif files in the directory
    input_tifs = glob(os.path.join(input_dir, "*.tif"))

    if not input_tifs:
        print("No files found to merge.")
        return

    src_files_to_mosaic = [rasterio.open(tif) for tif in input_tifs]

    mosaic, out_trans = merge(src_files_to_mosaic)

    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "compress": "LZW",
    })

    with rasterio.open(output_tif, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Merged file saved as {output_tif}")

def get_pine_island_bounds():
    '''Define bounds for the Pine Island Glacier region'''
    x0 = -102.00  # Min longitude
    y0 = -75.4    # Min latitude
    x1 = -99.00   # Max longitude
    y1 = -74.6    # Max latitude
    return x0, x1, y0, y1

def getCmdArgs():
    '''Get command line arguments'''
    p = argparse.ArgumentParser(description=("Create DEM from LVIS data"))
    p.add_argument("filename", type=str, help="Path to the LVIS file")
    p.add_argument("resolution", type=float, help="Resolution of the DEM in meters")
    p.add_argument("output_dir", type=str, help="Directory to save the DEM files and merged GeoTIFF")
    return p.parse_args()

if __name__ == "__main__":
    '''Main block'''

    args = getCmdArgs()
    filename = args.filename
    resolution = args.resolution
    output_dir = args.output_dir

    # Use Pine Island Glacier bounds
    x0, x1, y0, y1 = get_pine_island_bounds()

    # Define smaller chunks for more efficient processing
    chunk_size_x = (x1 - x0) / 20
    chunk_size_y = (y1 - y0) / 20

    # Loop over spatial subsets
    for x_start in np.arange(x0, x1, chunk_size_x):
        x_end = x_start + chunk_size_x
        for y_start in np.arange(y0, y1, chunk_size_y):
            y_end = y_start + chunk_size_y

            try:
                lvis = plotLVIS(filename, minX=x_start, minY=y_start, maxX=x_end, maxY=y_end, setElev=True)

                if lvis.nWaves > 0:
                    print(f"Data found in region: {x_start}, {y_start}, {x_end}, {y_end}")
                    lvis.setElevations()
                    lvis.estimateGround()
                    
                    # Reproject data
                    lvis.reprojectLVIS(3031)

                    # Write DEM
                    outName = f"{output_dir}/DEM_x_{x_start}_y_{y_start}.tif"
                    lvis.writeDEM(resolution, outName)
                    print(f"DEM saved to {outName}")
            except AttributeError:
                print(f"Tile skipped: {x_start}, {y_start}")

    # Merge all output files into a single file
    merged_output = f"{output_dir}/Merged.tif"
    merge_tifs(output_dir, merged_output)

    # Stop tracing memory allocations
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')

    print("[ Top 10 memory usage ]")
    for stat in top_stats[:10]:
        print(stat)

    tracemalloc.stop()
