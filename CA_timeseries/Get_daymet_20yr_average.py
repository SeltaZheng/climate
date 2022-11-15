"""
This script is used to retrieve monthly daymet data from 1980 to 1999
and calculate the 20yr average condition for the selected area.
"""
import pydaymet as daymet
import geopandas as gpd
import xarray
import numpy as np
import os, glob
import matplotlib.pyplot as plt
# directories:
dir_shp = r'D:\CA_TimeSeries\shp'
dir_out = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly\1980-1999'
shapefile = gpd.read_file(f'{dir_shp}/30m_extent_shp.shp')
geometry = shapefile.geometry[0]

# parameters:
# dates = ("1960-01-01", "2018-12-31")
# dates = ("1980-01-01", "2018-12-31")
dates = ("1980-01-01", "1999-12-31")

#%%------------------Step1: Get monthly data during the period -----------------------------------------


prcp = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables='prcp', time_scale="monthly")
# write to disk
prcp.to_netcdf(f'{dir_out}/monthly_prcp.nc', format='NETCDF4')

# get pet year by year to avoid memmory limits

# get pet, (only works for daily data)
i = 0
for y in range(1980, 2000): # (1980, 2019)
    print(y)
    dates = (f"{y}-01-01", f"{y}-12-31")
    pet = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables='dayl', pet="penman_monteith", time_scale="daily")
    # calculate monthly average of pet:
    l_time = pet.time.data
    months = l_time.astype('datetime64[M]')
    data = pet.pet.values
    for m in np.unique(months):
        idx = months == m
        # calculate mean for slice m:
        m_sum = np.sum(data[idx, :, :], axis=0)
        m_sum = m_sum[np.newaxis, :, :]
        if i == 0:
            pet_m = m_sum
        else:
            pet_m = np.concatenate((pet_m, m_sum), axis=0)
        i = i + 1


# construct new pet nc:
new_pet = prcp.copy()
new_pet["pet"] = (['time', 'y', 'x'],  pet_m)

new_pet = new_pet.drop_vars('prcp')
# new_pet = new_pet.drop_vars('lambert_conformal_conic')
new_pet['pet'].attrs = {'units': 'millimeters'}
# write to disk:
new_pet.to_netcdf(f'{dir_out}/monthly_pet.nc', format='NETCDF4')

del new_pet, prcp
# Download monthly data for other climate variables:
# parameters:
var = ["tmin", 'tmax', 'vp']
dates = ("1980-01-01", "1999-12-31")
# get data:
monthly = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables=var, time_scale="monthly")


i = 0
for y in range(1980, 2000):
    print(y)
    dates = (f"{y}-01-01", f"{y}-12-31")
    srad = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables=['srad', 'dayl'], time_scale="daily")
    l_time = srad.time.data
    months = l_time.astype('datetime64[M]')
    # convert srad from daylight average to daily total: MJ/m2/day
    data = srad.srad.values*srad.dayl.values/1000000
    for m in np.unique(months):
        idx = months == m
        # calculate mean for slice m:
        m_mean = np.mean(data[idx, :, :], axis=0)
        m_mean = m_mean[np.newaxis, :, :]
        if i == 0:
            srad_m = m_mean
        else:
            srad_m = np.concatenate((srad_m, m_mean), axis=0)
        i = i + 1

# add solar radiation to the climate object.
monthly["srad"] = (['time', 'y', 'x'],  srad_m)

monthly['srad'].attrs = srad.srad.attrs
# write to disk:
monthly.to_netcdf(f'{dir_out}/monthly_climate.nc', format='NETCDF4')


#%%----------------------- Step 2: calculate 20yr mean and convert to .tif----------------------------------------
def multi_year_avg(arr, interval=12):
    result = np.zeros((interval, arr.shape[1], arr.shape[2]))
    for i in range(interval):
        result[i, :, :] = np.nanmean(arr[i::interval, :, :], axis=0)

    return result



def manage_crs(xds, input_crs, output_crs="EPSG:32611"):
    # modify GeoTransform string to reflect 1000m spatial resolution
    # xds.lambert_conformal_conic.attrs['GeoTransform'] = '-1708000.75 1000.0 0.0 -250000.5 0.0 -1000'
    xds = xds.assign_coords(x=(xds.x * 1000))
    xds = xds.assign_coords(y=(xds.y * 1000))
    xds.rio.write_crs(input_crs=input_crs, inplace=True)
    xds_UTM = xds.rio.reproject(output_crs)
    return xds_UTM

f_ls = glob.glob(f'{dir_out}/*.nc')

# define the proper crs for the dataset:
input_crs = '+proj=lcc +datum=WGS84 +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +units=m +no_defs'

for f in f_ls:
    # extract var name info:
    base = os.path.basename(f)
    segs = base.split('monthly_')
    key = segs[-1][:-3]
    ds = xarray.open_dataset(f)

    if key == 'climate':
        vars = ['tmin', 'tmax', 'srad', 'vp']
    else:
        vars = [key]

    ds_new = ds.copy()
    ds_new = ds_new.drop_vars(vars)
    ds_new['time'] = ds_new.time[0:12]
    for var in vars:
        # calculate 20yr average for the monthly tmin and tmax:
        temp = multi_year_avg(ds[var].values)
        ds_new[var] = (['time', 'y', 'x'], temp)

    xds_UTM = manage_crs(ds_new, input_crs)
    for var in vars:
        print(f'Now exporting {var}')
        if 'grid_mapping' in xds_UTM[var].attrs.keys():
            del xds_UTM[var].attrs['grid_mapping']
        xds_UTM[var].rio.to_raster(f'{dir_out}/{var}.tif')

