
from pathlib import Path
import geopandas as gpd
import os
import zipfile
import numpy as np
import richdem as rd
from os import walk
from os.path import join
from pathlib import Path
import rasterio
from rasterio.merge import merge

# download UK grid data https://github.com/OrdnanceSurvey/OS-British-National-Grids
# download UK administrative borders, including coast line http://www.diva-gis.org/gdata
# download OS DEM data from http://www.ukso.org/static-maps/ordnance-survey-terrain-50.html

# I have exploded the 10k UK grid to be able to have access to each cell
# In QGIS this is done by going to Vector -> Geometry Tools -> Multipart to Singleparts
# The coast line vector data needs to be reprojected to EPSG 27700 and convert to line
# so that the sjoin method can look only to the cells that intersect the coast line


# input data
os_dem_data = '/Users/corina/3keel/data/terr50_gagg_gb/data' # path to OS DEM data folder
uk_10_km_grid = "/Users/corina/3keel/data/os_bng_grids_10km_explode.gpkg" # path to exploaded UK 10 km grid
uk_coast_line = "/Users/corina/3keel/data/GBR_adm0_27700_line.gpkg" # path to UK coast line

# out data
uk_coast_data = '/Users/corina/3keel/data/filetred_data' # path to filtered OS DEM tiles that interesct the coast line
slope_data = '/Users/corina/3keel/data/slope/' # path to output slope data
merged_slope_tiles = '/Users/corina/3keel/data/merged_slope.tif' # path to slope mosaic

def read_vector_data():
    """Read the grid and the boundaries vector files as geopandas
    """
    grid = gpd.read_file(uk_10_km_grid)
    coast_line = gpd.read_file(uk_coast_line)
    coast_line = coast_line[['NAME_ISO', 'geometry']]
    return grid, coast_line


def get_zipped_files():
    """Get all the zipped files that have been downloaded from the OS
    """
    for subdir, _, files in os.walk(os_dem_data):
        for file in files:
            if file.endswith(".zip"):
                yield Path(subdir, file)


def filter_10k_grid():
    """Filter the 10k UK grid to get only the tiles that overlap with the UK coast line
    """
    grid, coast_line = read_vector_data()
    # gpd.sjoin returns a new GeoDataFrame with the geometries for each object on the left dataframe
    # repeated for each geometry they intersect in the right, with the index of the object in the right
    intersect_grid_with_line_coast = gpd.sjoin(coast_line, grid, predicate='intersects')
    remove_duplicates = intersect_grid_with_line_coast.drop_duplicates(subset=['tile_name']) # remove duplicated tile names
    return remove_duplicates[['NAME_ISO', 'tile_name']]


def unzip_files_along_uk_coast():
    """Match the filtered coast cells with the name of the zipped files
    """
    zipped = get_zipped_files()
    filtered_grid = filter_10k_grid()
    # get only the file name to the zipped files
    for zip_files in zipped:
        filename = Path(zip_files).stem

    for tile in filtered_grid['tile_name']:
        to_lower = tile.lower() #convert the tile name to lowercase
        spl = filename.split('_')[1] #get the second element in the file name after the underscore
        spl2 = filename.split('_')[2] #get the third element in the file name after the underscore
        rm_digits = ''.join(filter(lambda x: not x.isdigit(), to_lower)) #remove the numbers from the tiles to be able to pass it as a folder in the path to the zipped files
        filtered_zipped_files = Path(os_dem_data, rm_digits, f"{to_lower}_{spl}_{spl2}.zip") #reconstruct the path to the zip files but get only the files that are on the coast
        print(f"Unzipping: {filtered_zipped_files}")
        if filtered_zipped_files.is_file():
            # if the files exist, unzip it
            zip_ref = zipfile.ZipFile(filtered_zipped_files) # create zipfile object
            zip_ref.extractall(uk_coast_data) # extract file to dir
            zip_ref.close() # close file


def create_slope_from_dem():
    """Genereting slope from DEM
    """
    allfiles = [join(uk_coast_data, f) for _, _, files in walk(uk_coast_data) for f in files if f.endswith(".asc")]
    for asc_files in allfiles:
        print(f"Generate slope for: {asc_files}")
        dem = rd.LoadGDAL(asc_files, no_data=np.nan)
        slope = rd.TerrainAttribute(dem, attrib='slope_degrees') # generate slope from DEM
        with rasterio.open(asc_files) as src:
            kwargs = src.meta.copy()
            kwargs.update({
                'driver': "GTiff",
                })
            filename = Path(asc_files).stem # keep only the file name from the whole path
            rm_ext = filename.removesuffix('.asc') # remove the extension from the filename to
            slope_files = Path(slope_data, f'{rm_ext}.tif') # dynamically generate path for slope data

            with rasterio.open(slope_files, 'w', **kwargs) as dst:
                dst.write(slope, 1)


def create_mosaic():
    """Creating a mosaic out of all the tiles
    """
    allfiles = [join(slope_data, f) for _, dirs, files in walk(slope_data) for f in files if f.endswith(".tif")]
    src_files_to_mosaic = []
    for fp in allfiles:
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)
    print("Creating the slope mosaic")
    mosaic, out_trans = merge(src_files_to_mosaic)
    return mosaic, out_trans


def save_mosaic():
    """Saves the mosaic to disk
    """
    mosaic, out_trans = create_mosaic()
    out_meta={'driver': 'GTiff',
            'dtype': 'float32',
            'nodata': np.nan,
            'width': mosaic.shape[2],
            'height': mosaic.shape[1],
            'count': 1,
            'crs': rasterio.crs.CRS.from_epsg(27700),
            'transform': out_trans,
            "compress": "lzw"}
    print(f'Saving the slope mosaic to: {merged_slope_tiles}')
    with rasterio.open(merged_slope_tiles, "w", **out_meta) as dest:
        dest.write(mosaic)


def main():
    unzip_files_along_uk_coast()
    create_slope_from_dem()
    save_mosaic()

if __name__ == "__main__":
    main()