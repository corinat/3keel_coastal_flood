import geopandas as gpd
import os
import zipfile
import richdem as rd
from os import walk
from os.path import join
from pathlib import Path
import rasterio
from rasterio.merge import merge
from rasterstats import zonal_stats
import pandas as pd
import numpy as np
from datetime import datetime
import glob

pd.options.mode.chained_assignment = None

# download UK grid data https://github.com/OrdnanceSurvey/OS-British-National-Grids
# download UK administrative borders, including coast line http://www.diva-gis.org/gdata
# download OS DEM data from http://www.ukso.org/static-maps/ordnance-survey-terrain-50.html


# input data
os_dem_data = "/Users/corina/3keel_coastal_flood/data/terr50_gagg_gb/data"  # path to OS DEM data folder
uk_coast_line = "/Users/corina/3keel_coastal_flood/data/GBR_adm/GBR_adm0.shp"  # path to UK coast line
voronoi_data = (
    "/Users/corina/3keel_coastal_flood/data/2100_SeaLevelRise&Surge_Voronoi.gpkg"
)
os_bng_grids = "/Users/corina/3keel_coastal_flood/data/os_bng_grids.gpkg"

# out data
uk_coast_data = "/Users/corina/3keel_coastal_flood/data/filetred_data/"  # path to filtered OS DEM tiles that interesct the coast line
slope_data = (
    "/Users/corina/3keel_coastal_flood/data/slope/"  # path to filtered slope data
)
merged_slope_tiles = (
    "/Users/corina/3keel_coastal_flood/data/merged_slope.tif"  # path to slope mosaic
)
centroid_stats = "/Users/corina/3keel_coastal_flood/data/centroid_stats.gpkg"  # path to the location where the centroids filel is saved
slope_on_1km_grid = "/Users/corina/3keel_coastal_flood/data/uk_grid_voronoi_slope.gpkg"  # path to final product


def polygon_to_line():
    """Convert polygon to line for coast data and transform it from 4326 to 27700"""
    print("Convert polygon to line for coast data and transform it from 4326 to 27700")
    coast_line = gpd.read_file(uk_coast_line)
    coast = coast_line.copy()  # copy GeoDataFrame
    coast.geometry = coast_line.geometry.boundary
    coast = coast[["NAME_ISO", "geometry"]]
    return coast.to_crs(27700)


def read_uk_grid():
    """Extracting the 1km and 10km grids from os_bng_grids.gpkg"""
    grid_10k = gpd.read_file(os_bng_grids, layer="10km_grid")
    grid_1k = gpd.read_file(os_bng_grids, layer="1km_grid")
    return grid_10k, grid_1k


def get_zipped_files():
    """Get all the zipped files that have been downloaded from the OS"""
    for subdir, _, files in os.walk(os_dem_data):
        for file in files:
            if file.endswith(".zip"):
                yield Path(subdir, file)


def filter_10k_grid():
    """Filter the 10k UK grid to get only the tiles that overlap with the UK coast line"""
    coast_line = polygon_to_line()
    _, exploded_grid_10k = explode_grid()
    # gpd.sjoin returns a new GeoDataFrame with the geometries for each object on the left dataframe
    # repeated for each geometry they intersect in the right, with the index of the object in the right
    intersect_grid_with_line_coast = gpd.sjoin(
        coast_line, exploded_grid_10k, predicate="intersects"
    )
    remove_duplicates = intersect_grid_with_line_coast.drop_duplicates(
        subset=["tile_name"]
    )  # remove duplicated tile names
    return remove_duplicates[["NAME_ISO", "tile_name"]]


def unzip_files_along_uk_coast():
    """Match the filtered coast cells with the name of the zipped files"""
    zipped = get_zipped_files()
    filtered_grid = filter_10k_grid()
    # get only the file name to the zipped files
    for zip_files in zipped:
        filename = Path(zip_files).stem

    for tile in filtered_grid["tile_name"]:
        to_lower = tile.lower()  # convert the tile name to lowercase
        spl = filename.split("_")[
            1
        ]  # get the second element in the file name after the underscore
        spl2 = filename.split("_")[
            2
        ]  # get the third element in the file name after the underscore
        rm_digits = "".join(
            filter(lambda x: not x.isdigit(), to_lower)
        )  # remove the numbers from the tiles to be able to pass it as a folder in the path to the zipped files
        filtered_zipped_files = Path(
            os_dem_data, rm_digits, f"{to_lower}_{spl}_{spl2}.zip"
        )  # reconstruct the path to the zip files but get only the files that are on the coast
        print(f"Unzipping: {filtered_zipped_files}")
        if filtered_zipped_files.is_file():
            # if the files exist, unzip it
            zip_ref = zipfile.ZipFile(filtered_zipped_files)  # create zipfile object
            zip_ref.extractall(uk_coast_data)  # extract file to dir
            zip_ref.close()  # close file


def create_slope_from_dem():
    """Genereting slope from DEM"""
    allfiles = [
        join(uk_coast_data, f)
        for _, _, files in walk(uk_coast_data)
        for f in files
        if f.endswith(".asc")
    ]
    for asc_files in allfiles:
        print(f"Generate slope for: {asc_files}")
        dem = rd.LoadGDAL(asc_files, no_data=np.nan)
        slope = rd.TerrainAttribute(
            dem, attrib="slope_degrees"
        )  # generate slope from DEM
        with rasterio.open(asc_files) as src:
            kwargs = src.meta.copy()
            kwargs.update(
                {
                    "driver": "GTiff",
                }
            )
            filename = Path(
                asc_files
            ).stem  # keep only the file name from the whole path
            rm_ext = filename.removesuffix(
                ".asc"
            )  # remove the extension from the filename to
            slope_files = Path(
                slope_data, f"{rm_ext}.tif"
            )  # dynamically generate path for slope data

            with rasterio.open(slope_files, "w", **kwargs) as dst:
                dst.write(slope, 1)


def create_mosaic():
    """Creating a mosaic out of all the tiles"""
    allfiles = [
        join(slope_data, f)
        for _, _, files in walk(slope_data)
        for f in files
        if f.endswith(".tif")
    ]  # looping in the directory where all the slopw tiles have been save and create a list to these files
    src_files_to_mosaic = [
        rasterio.open(fp) for fp in allfiles
    ]  # opening all the slope files in the list created above with rasterio
    print("Creating the slope mosaic")
    mosaic, out_trans = merge(src_files_to_mosaic)  # creating the mosaic
    return mosaic, out_trans


def save_mosaic():
    """Saves the mosaic to disk"""
    mosaic, out_trans = create_mosaic()
    out_meta = {
        "driver": "GTiff",
        "dtype": "float32",
        "nodata": np.nan,
        "width": mosaic.shape[2],
        "height": mosaic.shape[1],
        "count": 1,
        "crs": rasterio.crs.CRS.from_epsg(27700),
        "transform": out_trans,
        "compress": "lzw",
    }
    print(f"Saving the slope mosaic to: {merged_slope_tiles}")
    with rasterio.open(merged_slope_tiles, "w", **out_meta) as dest:
        dest.write(mosaic)


def filter_1km_grid():
    """Getting only the cells from 1 km grid that overlap with the coast line"""
    exploded_grid_1k, _ = explode_grid()
    coast_line = polygon_to_line()
    intersect_grid_with_line_coast = gpd.sjoin(
        coast_line, exploded_grid_1k, predicate="intersects"
    )  # create a spatial join to get the ikm grid cells that intersect the coast line
    get_index = intersect_grid_with_line_coast[
        "index1"
    ].tolist()  # get the index of the cells that intersect the coast line and save it as a list
    filt_1km_grid = exploded_grid_1k[
        exploded_grid_1k["index1"].isin(get_index)
    ]  # based on the index list created above get the cells from the exploded 1km grid
    print("Creating filtered 1km grid along UK coast")
    return filt_1km_grid


def explode_grid():
    """Explode grids to be able to reach each cell inside the grid"""
    grid_10k, grid_1k = read_uk_grid()
    grid_1k["index1"] = grid_1k.index
    exploded_grid_1k = grid_1k.explode(
        index_parts=False
    )  # explode the 1km grid to be able to reach each cell
    exploded_grid_10k = grid_10k.explode(
        index_parts=False
    )  # explode the 1km grid to be able to reach each cell
    return exploded_grid_1k, exploded_grid_10k


def get_centroid():
    """Getting the centroid for multiple polygons"""
    filtered_1km_grid = filter_1km_grid()
    print("Creating centroids for each 1km grid cell along the UK coast")
    cent = filtered_1km_grid.centroid
    cent.to_file(centroid_stats, driver="GPKG")


def produce_zonal_stats(vector, raster, affine, metrics, nodata):
    """Producing zonal stats"""
    gdf = gpd.read_file(vector)

    stats = zonal_stats(
        gdf,
        raster,
        affine=affine,
        stats=metrics,
        nodata=nodata,
    )
    stats = pd.DataFrame(stats)
    stats = gdf.join(stats, how="left")
    print("Finished calculating zonal stats")
    return stats


def get_stats():
    """Getting the value of the slope mosaic for each centroid"""
    dem = rasterio.open(merged_slope_tiles)
    array = dem.read(1)
    print("Getting the value of the slope where the centroid overlaps the slope mosaic")
    print("This step takes around 8 minutes to complete")
    stats = produce_zonal_stats(
        vector=centroid_stats,
        raster=array,
        affine=dem.profile["transform"],
        metrics="mean",
        nodata=np.nan,
    )
    return stats


def attach_slope_value_to_1km_grid():
    """Getting the values from voronoi file  for those cells that intersect with the voronoi
    polygons and the slope from the centroids and add them to the 1km grid
    """
    voronoi = gpd.read_file(voronoi_data)
    voronoi = voronoi.to_crs(27700)
    point = get_stats()
    poly = filter_1km_grid()
    join_poly_point = poly.sjoin(point, how="left")
    join_poly_point = join_poly_point.drop(["index_right"], axis=1)
    join_voronoi_uk_grid = join_poly_point.sjoin(voronoi, how="inner").drop_duplicates(
        subset=["tile_name"]
    )
    join_voronoi_uk_grid.rename(columns={"mean": "slope_value"}, inplace=True)
    slope_1km_uk_grid = join_voronoi_uk_grid[
        [
            "tile_name",
            "slope_value",
            "LT_05_local",
            "LT_50_local",
            "LT_95_local",
            "SWRL_MT",
            "slrbLT_50",
            "swrl-slr50",
            "Add-Swrl5",
            "geometry",
        ]
    ]
    print("Create 1km grid with slope values")
    slope_1km_uk_grid.to_file(slope_on_1km_grid, driver="GPKG")


def remove_intermediary_data():
    """Remove intermediary data"""
    rm_data = [merged_slope_tiles, centroid_stats]
    for files in glob.glob(f"{uk_coast_data}*", recursive=True):
        os.remove(files)
    for files in glob.glob(f"{slope_data}*", recursive=True):
        os.remove(files)
    for data in rm_data:
        os.remove(data)


def main():
    start_time = datetime.now()
    unzip_files_along_uk_coast()
    create_slope_from_dem()
    save_mosaic()
    get_centroid()
    attach_slope_value_to_1km_grid()
    remove_intermediary_data()
    end_time = datetime.now()
    print("Duration: {}".format(end_time - start_time))


if __name__ == "__main__":
    main()
