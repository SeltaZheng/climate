import os

import matplotlib.pyplot as plt
import rioxarray
import xarray
import numpy as np
import glob
import pandas as pd
# def manage_crs(xds, input_crs, output_crs="EPSG:32611"):
#     # modify GeoTransform string to reflect 1000m spatial resolution
#     # xds.lambert_conformal_conic.attrs['GeoTransform'] = '-1708000.75 1000.0 0.0 -250000.5 0.0 -1000'
#     xds = xds.assign_coords(x=(xds.x * 1000))
#     xds = xds.assign_coords(y=(xds.y * 1000))
#     xds.rio.write_crs(input_crs=input_crs, inplace=True)
#     xds_UTM = xds.rio.reproject(output_crs)
#     return xds_UTM
# # input crs: EPSG 4326
# input_crs = '+proj=longlat +datum=WGS84 +no_defs +type=crs'



#%%------------------------------ Quality check on the soil moisture data---------------------

dir_in = r'F:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP-HB\raw'
dir_out = r'F:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP-HB\QC'
fn_ls = glob.glob(f'{dir_in}/*.nc')
dates = []
for i, f in enumerate(fn_ls):
    basename = os.path.basename(f)
    print(f)
    ds = xarray.open_dataset(f)
    # count non nan pixels for each day
    days = ds['time'].data.astype('datetime64[D]')
    dates.append(days)
    mask = ~np.isnan(ds['SMAPHB_SM'])
    # save the daily pixel count into one file, and monthly count into one npz file
    c = np.count_nonzero(mask, axis=(1, 2))
    monthly_c = np.count_nonzero(mask, axis=0)
    y = days[0].astype('datetime64[Y]')
    # flip to get proper layout
    monthly_c = np.flip(monthly_c, axis=0)
    if i == 0:
        daily_c = c
    else:
        daily_c = np.concatenate((daily_c, c))

    fn_out = basename.replace('.nc', '.npz')
    np.savez(f'{dir_out}/{fn_out}', monthly_c)

    # plotting------------------
    if i == 0:
        fig, axes = plt.subplots(ncols=4, nrows=3, figsize=(12, 12))
        # flatten the axes:
        axs = axes.flat
        j = 0
    if i > 0:
        if y != y_old:
            # a different year, save the previous figure and create new fig
            # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            # fig.colorbar(im, cax=cbar_ax)
            fig.colorbar(im, ax=axes.ravel().tolist())
            fig.savefig(f'{dir_out}/{y_old}_valid_days.png', dpi=300)
            plt.close(fig)
            fig, axes = plt.subplots(ncols=4, nrows=3, figsize=(12, 12))
            # flatten the axes:
            axs = axes.flat
            j = 0

    im = axs[j].imshow(monthly_c, vmin=0, vmax=31)
    j = j+1
    y_old = y
    del ds
fig.colorbar(im, ax=axes.ravel().tolist())
fig.savefig(f'{dir_out}/{y}_valid_days.png', dpi=300)
plt.close(fig)
# concat dates lists into one list:
dates_com = dates[0]
for i in range(1, len(dates)):
    dates_com = np.concatenate((dates_com, dates[i]))
daily_valid = pd.DataFrame(data=dates_com, columns=['dates'])
daily_valid['counts'] = daily_c
daily_valid.to_csv(f'{dir_out}/daily_counts.csv', index=False)



#%%------------------------------ days per month with manximum valid pixels ----------------------------

dir_in = r'F:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP-HB\QC'
df_qc = pd.read_csv(f'{dir_in}/daily_counts.csv')
df_qc['mark'] = df_qc['dates'].str[:7]
valid = []
for m in df_qc['mark'].unique():
    sub = df_qc.loc[df_qc['mark'] == m, ['counts']]
    n = sum(sub['counts'] == max(sub['counts']))
    valid.append(n)

df_valid = pd.DataFrame(data=valid, columns=['valid'])
df_valid['month'] = df_qc['mark'].unique()
df_valid.to_csv(f'{dir_in}/valid_days.csv', index=False)

#%%------------------------------ calculate monthly mean for each year----------------------
dir_in = r'F:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP-HB\raw'
dir_qc = r'F:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP-HB\QC'
dir_out = r'F:\CA_TimeSeries\OtherData\Climate\SMAP\SMAP-HB\tiff'
df_qc = pd.read_csv(f'{dir_qc}/daily_counts.csv')

# extract year-month
df_qc['datetime'] = pd.to_datetime(df_qc['dates'])

df_qc['mark'] = df_qc['datetime'].dt.to_period('M')
df_qc['year'] = df_qc['datetime'].dt.to_period('Y')

base = 'SMAP-HB_surface-soil-moisture_30m_daily'

# loop through year
for y in df_qc['year'].unique()[1:]:
    df_y = df_qc.loc[df_qc['year'] == y, :]

    for i, m in enumerate(df_y['mark'].unique()):
        print(f'now processing {m}')
        sub = df_y.loc[df_y['mark'] == m, ['counts']].reset_index()
        idx = sub.index[sub['counts'] == max(sub['counts'])]
        f = f'{dir_in}/{base}_{m}.nc'
        ds = xarray.open_dataset(f)
        # only extract the days with maximum coverage:
        sm = ds['SMAPHB_SM'][idx, :, :].values

        # average non-nan pixels:
        mask = np.isnan(sm).any(axis=0)
        avg = np.nanmean(sm, axis=0)
        avg[mask] = -9999


        if i == 0:
            SMAPHB_SM_monthly = avg[np.newaxis, :, :]
        else:
            SMAPHB_SM_monthly = np.concatenate((SMAPHB_SM_monthly, avg[np.newaxis, :, :]), axis=0)

    # create a new nc and convert to tif:
    ds = ds.drop_vars('SMAPHB_SM')
    ds = ds.drop_vars('time')

    # add new time variable:
    ds['time'] = df_y['mark'].unique()
    ds = ds.rename({'lat': 'y', 'lon': 'x'})

    ds['SMAPHB_SM_monthly'] = (['time', 'y', 'x'], SMAPHB_SM_monthly)

    # to tiff:

    ds['SMAPHB_SM_monthly'].rio.to_raster(f'{dir_out}/SM_monthly_{y}.tif')


    del ds, SMAPHB_SM_monthly





