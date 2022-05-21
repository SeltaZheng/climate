"""
This script is a testing script to acquire daymet data using pydaymet
more intro can be found at: https://docs.hyriver.io/en/latest/notebooks/daymet.html
"""
import matplotlib.pyplot as plt
import pydaymet as daymet
from pynhd import NLDI

# vector file:
geometry = NLDI().get_basins("01031500").geometry[0]
dates = ("2000-01-01", "2000-01-06")
daily = daymet.get_bygeom(geometry, dates, variables="prcp", pet="hargreaves_samani")
ax = daily.where(daily.pet > 0).pet.plot(x="x", y="y", row="time", col_wrap=3)
# ax.fig.savefig("_static/daymet_grid.png", facecolor="w", bbox_inches="tight")

var = ["prcp", "tmin"]
dates = ("2000-01-01", "2000-06-30")
monthly = daymet.get_bygeom(geometry, dates, variables=var, time_scale="monthly")
ax = monthly.where(monthly.prcp > 0).prcp.plot(x="x", y="y", row="time", col_wrap=3)

