"""
This script is used to convert netcdf to geotiff
"""
import os

import rioxarray
import xarray
import numpy as np
import glob
import pandas as pd

#%%-------------------- NetCDF to Tiff (reproject to UTM-CA EPSG:32611) ---------------------------------------------------------
def manage_crs(xds, input_crs, output_crs="EPSG:32611"):
    # modify GeoTransform string to reflect 1000m spatial resolution
    # xds.lambert_conformal_conic.attrs['GeoTransform'] = '-1708000.75 1000.0 0.0 -250000.5 0.0 -1000'
    xds = xds.assign_coords(x=(xds.x * 1000))
    xds = xds.assign_coords(y=(xds.y * 1000))
    xds.rio.write_crs(input_crs=input_crs, inplace=True)
    xds_UTM = xds.rio.reproject(output_crs)
    return xds_UTM



# directory
dir_in = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly'
dir_out = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly\tiff'
# file names:
f_ls = glob.glob(f'{dir_in}/*.nc')

# define the proper crs for the dataset:
input_crs = '+proj=lcc +datum=WGS84 +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +units=m +no_defs'

for f in f_ls:
    # extract var name info:
    base = os.path.basename(f)
    segs = base.split('monthly_')
    key = segs[-1][:-3]
    ds = xarray.open_dataset(f)

    if key == 'climate':

        # if it's climate data, extract tmin, tmax, and tmean:
        vars = ['tmin', 'tmax', 'tmean']
        # calculate tmean:
        tmean = (ds.tmin.values + ds.tmax.values)/2
        ds['tmean'] = (['time', 'y', 'x'], tmean)
        xds_UTM = manage_crs(ds, input_crs)
        for var in vars:
            print(f'Now exporting {var}')
            if 'grid_mapping' in xds_UTM[var].attrs.keys():
                del xds_UTM[var].attrs['grid_mapping']
            xds_UTM[var].rio.to_raster(f'{dir_out}/{var}.tif')
    else:
        if 'spei' in key:
            # if the file is calculated spei, rename the coords to y and x
            ds = ds.rename({'lat': 'y', 'lon': 'x'})

        temp = ds[key].values
        temp = np.moveaxis(temp, -1, 0)
        ds[key] = (['time', 'y', 'x'], temp)
        xds_UTM = manage_crs(ds, input_crs)
        print(f'Now exporting {key}')
        if 'grid_mapping' in xds_UTM[key].attrs.keys():
            del xds_UTM[var].attrs['grid_mapping']
        xds_UTM[key].rio.to_raster(f'{dir_out}/{key}.tif')


#%%--------------------- Get the date meta for the data, convert to .csv------------------
# directory
dir_in = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly'
dir_out = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly\tiff'
fn = 'monthly_prcp.nc'
xds = xarray.open_dataset(f'{dir_in}/{fn}')
dates = xds.time.data
months = dates.astype('datetime64[M]')
df = pd.DataFrame(data=months, columns=['Time'])
df.to_csv(f'{dir_out}/daymet_meta.csv', index=False)