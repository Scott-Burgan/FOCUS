''' Plots the timeseries files as a 3x2 subplot
    Laura Burgin, August 2020 '''

import iris
import matplotlib
import iris
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import iris.plot as iplt
import pandas as pd


########## set up names of directories ################################

# directory to read the data from
data_dir = '/data/users/sburgan/CP4A/processed_data/Mongu/'

# directory to save the plots to
plot_dir = '/home/h05/sburgan/Documents/CP4A/plots/'

if __name__ == '__main__':

    #################   read in data  #####################

    cp4 = iris.load_cube(data_dir + 'Mongu_CP4A_daily_precip_*')

        # to convert kg/ms/s to mm/day
    cp4= cp4 * 86400.0


    cp4_1 = cp4[:, 1, 1]
    cp4_2 = cp4[:, 0, 0]
    cp4_3 = cp4[:, 3, 0]
    cp4_4 = cp4[:, 3, 3]
    cp4_5 = cp4[:, 2, 4]
    cp4_6 = cp4[:, 0, 2]

    #################   generate timeseries plots #####################
    
    # I'm not sure why I read the data in again like this and not the timeseries files

    fig = plt.figure(figsize=(8,10))
    plt.suptitle("Daily rainfall in CP4A gridboxes ", fontsize=16)
    plt.subplot(321)
    iplt.plot(cp4_1, 'bx', markersize=6)
    plt.ylim(0,350)
    plt.ylabel("mm/day", fontsize=12)
    plt.subplot(322)
    iplt.plot(cp4_2, 'gx', markersize=6)
    plt.ylim(0,350)
    plt.ylabel("mm/day", fontsize=12)
    plt.subplot(323)
    iplt.plot(cp4_3, 'rx', markersize=6)
    plt.ylim(0,350)
    plt.ylabel("mm/day", fontsize=12)
    plt.subplot(324)
    iplt.plot(cp4_4, 'yx', markersize=6)
    plt.ylim(0,350)
    plt.ylabel("mm/day", fontsize=12)
    plt.subplot(325)
    iplt.plot(cp4_5, 'mx', markersize=6)
    plt.ylim(0,350)
    plt.ylabel("mm/day", fontsize=12)
    plt.subplot(326)
    iplt.plot(cp4_6, 'cx', markersize=6)
    plt.ylim(0,350)
    plt.ylabel("mm/day", fontsize=12)
    plt.savefig(plot_dir+'CP4_daily_rainfall_ts_Mongu_TEST.png')

