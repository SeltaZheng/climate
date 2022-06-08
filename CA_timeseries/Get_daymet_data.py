"""
Retrieving daymet data by shape file
"""

import pydaymet as daymet
import geopandas as gpd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# directories:
dir_shp = r'D:\CA_TimeSeries\shp'
dir_out = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly'
shapefile = gpd.read_file(f'{dir_shp}/30m_extent_shp.shp')
geometry = shapefile.geometry[0]

# parameters:
# dates = ("1960-01-01", "2018-12-31")
# dates = ("1980-01-01", "2018-12-31")
dates = ("2010-01-01", "2018-12-31")

#%%------------------ GET prcp and pet data -----------------------------------------

prcp = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables='prcp', time_scale="monthly")

# modify prcp for climate indices:
prcp = prcp.drop_vars(['lat', 'lon'])
# rename dimensions:
prcp = prcp.rename({'y': 'lat', 'x': 'lon'})
# convert the dimension:
data = prcp.prcp.values
data = np.moveaxis(data, 0, -1)
prcp['prcp'] = (['lat', 'lon', 'time'], data)
prcp = prcp.drop_vars('lambert_conformal_conic')
prcp['prcp'].attrs= {'units': 'mm'}
# write to disk
prcp.to_netcdf(f'{dir_out}/monthly_prcp.nc', format='NETCDF4')

# get pet year by year to avoid memmory limits

# get pet, (only works for daily data)
i = 0
for y in range(2010, 2019): # (1980, 2019)
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
pet_m = np.moveaxis(pet_m, 0, -1)
new_pet["pet"]=(['lat', 'lon', 'time'],  pet_m)

new_pet = new_pet.drop_vars('prcp')
# new_pet = new_pet.drop_vars('lambert_conformal_conic')
new_pet['pet'].attrs = {'units': 'millimeters'}
# write to disk:
new_pet.to_netcdf(f'{dir_out}/monthly_pet.nc', format='NETCDF4')


# Download monthly data for other climate variables:
# parameters:
var = ["tmin", 'tmax', 'vp']
# dates = ("1960-01-01", "2018-12-31")
dates = ("2010-01-01", "2018-12-31")

# get data:
monthly = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables=var, time_scale="monthly")

monthly.to_netcdf(f'{dir_out}/monthly_climate.nc', format='NETCDF4')



#%%----------------------- get monthly srad data: -------------------------

i = 0
for y in range(2010, 2019):
    print(y)
    dates = (f"{y}-01-01", f"{y}-12-31")
    srad = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables='srad', time_scale="daily")
    # calculate monthly average of pet:
    l_time = srad.time.data
    months = l_time.astype('datetime64[M]')
    data = srad.srad.values
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

# construct monthly srad nc:

dates = ("2010-01-01", "2018-12-31")
monthly = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables='tmin', time_scale="monthly")

# get data:

new_srad = monthly.drop_vars('tmin')
srad_m = np.moveaxis(srad_m, 0, -1)
new_srad["srad"] = (['y', 'x', 'time'],  srad_m)
# new_pet = new_pet.drop_vars('lambert_conformal_conic')

new_srad['srad'].attrs = srad.srad.attrs
# write to disk:
new_srad.to_netcdf(f'{dir_out}/monthly_srad.nc', format='NETCDF4')


