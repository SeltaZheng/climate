"""
Retrieving daymet data by shape file
"""

import pydaymet as daymet
import geopandas as gpd


# directories:
dir_shp = r'D:\CA_TimeSeries\shp'
dir_out = r'D:\CA_TimeSeries\OtherData\Climate\daymet_monthly'
shapefile = gpd.read_file(f'{dir_shp}/30m_extent_shp.shp')
geometry = shapefile.geometry[0]

# parameters:
var = ["prcp", "tmin"]
dates = ("2000-01-01", "2000-06-30")
# get data:
monthly = daymet.get_bygeom(geometry, dates, crs='epsg:32611', variables=var, time_scale="monthly")
# ax = monthly.where(monthly.prcp > 0).prcp.plot(x="x", y="y", row="time", col_wrap=3)