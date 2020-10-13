''' Script to further process the CP4A data saved by process_cp4.py
    into timeseries for individual gridboxes.
    Laura Burgin, August 2020 '''

import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas
from iris.exceptions import CoordinateNotFoundError

########## set up names of directories ################################

# directory to read data from and save time series to
data_dir = '/data/users/sburgan/CP4A/processed_data/Southern_Africa/'

data_save = '/data/users/sburgan/CP4A/processed_data/Southern_Africa/'

if __name__ == '__main__':

    #################   read in data  #####################
    
    cp4 = iris.load_cube(data_dir+'Southern Africa_CP4A_daily_precip_*')
    
    # for CP4 to convert kg/ms/s to mm/day
    #cp4.data = cp4.data*86400.0
    cp4 = cp4 * 86400.0

    # test plot to check the right area
    crs_latlon = ccrs.PlateCarree()
    ax = plt.axes(projection=crs_latlon)
    iplt.pcolormesh(cp4[0,:,:])   
    ax.set_extent((60, -50, 50, -50), crs=crs_latlon)
    ax.gridlines(crs=crs_latlon, linestyle=':')
    ax.add_feature(cfeature.BORDERS, linestyle='-')
    ax.add_feature(cfeature.COASTLINE, linestyle ='-')
    plt.plot(39.27, -6.80, marker='o', markersize=4.0, markerfacecolor='red',
                 transform=crs_latlon)
    iplt.show()

    #################   generate timeseries #####################
    
    # this section could probably be done better as a loop!

    # these were randomly selected gridboxes from an area equivalent to 
    # the corresponding GCM gridbox for the Lusaka region
