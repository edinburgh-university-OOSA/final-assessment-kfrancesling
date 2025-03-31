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
        "compress": "LZW"
    })

    with rasterio.open(output_tif, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Merged file saved as {output_tif}")

def getCmdArgs():
    '''Get command line arguments'''
    p = argparse.ArgumentParser(description=("Create DEM from LVIS data"))
    p.add_argument("filename", type=str, help="Path to the LVIS file")
    p.add_argument("resolution", type=float, help="Resolution of the DEM in meters")
    p.add_argument("input_dir", type=str, help="Directory containing tiled GeoTIFF files")
    p.add_argument("output_dir", type=str, help="Directory to save the merged GeoTIFF file")
    return p.parse_args()

if __name__ == "__main__":
    '''Main block'''

    # Start tracing memory allocations
    tracemalloc.start()

    # Read the command line arguments
    args = getCmdArgs()
    filename = args.filename
    resolution = args.resolution
    input_dir = args.input_dir
    output_dir = args.output_dir

    # Create an instance of the class with "onlyBounds" flag
    b = plotLVIS(filename, onlyBounds=True)

    # Print the bounds to check the range
    print("Bounds of the data:", b.bounds)

    # Define smaller chunks for more efficient processing
    chunk_size_x = (b.bounds[2] - b.bounds[0]) / 20
    chunk_size_y = (b.bounds[3] - b.bounds[1]) / 20

    # Loop over spatial subsets
    for x0 in np.arange(b.bounds[0], b.bounds[2], chunk_size_x):
        x1 = x0 + chunk_size_x
        for y0 in np.arange(b.bounds[1], b.bounds[3], chunk_size_y):
            y1 = y0 + chunk_size_y

            try:
                lvis = plotLVIS(filename, minX=x0, minY=y0, maxX=x1, maxY=y1, setElev=True)

                if lvis.nWaves > 0:
                    print(f"Data found in region: {x0}, {y0}, {x1}, {y1}")
                    lvis.setElevations()
                    lvis.estimateGround()
                    
                    # Reproject data
                    lvis.reprojectLVIS(3031)

                    # Write DEM
                    outName = f"{output_dir}/DEM_x_{x0}_y_{y0}.tif"
                    lvis.writeDEM(resolution, outName)
                    print(f"DEM saved to {outName}")
            except AttributeError:
                print(f"Tile skipped: {x0}, {y0}")

    # Merge all output files into a single file
    merged_output = f"{output_dir}/Merged.tif"
    merge_tifs(output_dir, merged_output)

    # Stop tracing memory allocations
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')

    print("[ Top 10 memory usage ]")
    for stat in top_stats[:10]:
        print(stat)
