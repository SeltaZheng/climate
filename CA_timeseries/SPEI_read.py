"""
Read and visualize SPEI
"""
import xarray
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

f_input = 'SPEI_monthly_spei_gamma_03.nc'
SPEI = xarray.open_dataset(f_input)
data1 = SPEI.spei_gamma_03.values

plt.show()


f_spi = 'SPI_monthly_spi_gamma_03.nc'
SPI = xarray.open_dataset(f_spi)
data2 = SPI.spi_gamma_03.values
plt.plot(data1[100, 50, 370:], label='SPEI')
plt.plot(data2[100, 50, 370:], label='SPI')
plt.legend()
plt.show()

plt.imshow(data1[:, :, 100])
plt.show()

plt.imshow(data2[:, :, 100])
plt.show()