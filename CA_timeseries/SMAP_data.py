"""
Example for reading the SMAP L4 soil moisture data 9km resolution
An easier SM data (NASA_USDA_SM) was downloaded from GEE later (see colab notebook for details), the SMAP L4 data was not used.
"""

import glob, os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# directories;

dir_in = r'D:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP_L4_3hour\227208495'

test = xr.open_mfdataset(f'{dir_in}/SMAP_L4_SM_gph_20150331T193000_Vv6032_001_HEGOUT.nc', group='Geophysical_Data')


