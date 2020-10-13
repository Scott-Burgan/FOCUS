#! /bin/bash
''' Script to extract a smaller area of the CP4 dataset to make it small enough to work with
 and converts the data to daily means.
 Also creates a copy locally to avoid losing the data from /scratch
 Remember to delete the local version at an appropriate time
 Laura Burgin, August 2020 '''

########################################################################################
import iris
import iris.coord_categorisation
import numpy as np
from iris.experimental.equalise_cubes import equalise_attributes
from iris.analysis._interpolation import get_xy_dim_coords
from iris.util import unify_time_units
from os import path
import argparse

########## set up names of directories and other variables ################################

# local directory to save data to
data_dir ='/data/users/sburgan/CP4A/processed_data/Southern_Africa/'

# select a small region e.g. southern Zambia region
lat1 = -45
lat2 = 12
lon1 = 343
lon2 = 423

#Allow for runid to be used as a command line / bash argument

parser = argparse.ArgumentParser()
parser.add_argument('stashid')
args = parser.parse_args()
stash = args.stashid

# name of the model run
#stash = 'ac144'
region = 'Southern Africa'

#################   set-up functions    #####################
def load_CP4_cubelist(data_file, x_coord_name, y_coord_name):
    """ load a cp4 precip dataset using lat/lon constraints
        makes use of the time callback function below too """
    xy_constraint = iris.Constraint(
        coord_values={y_coord_name: lambda cell: lat1 < cell < lat2,
                      x_coord_name: lambda cell: lon1 < cell < lon2})
    data = iris.load(data_file, xy_constraint, callback=time_callback)
    return data


def time_callback(cube, field, filename):
    ''' adds year, month and day of month to the cube on loading '''
    iris.coord_categorisation.add_year(cube, 'time')
    iris.coord_categorisation.add_month(cube, 'time')
    iris.coord_categorisation.add_day_of_month(cube, 'time')


def convec_2_precipflux(cube_list):
    '''
    Function to add together large scale rain and snow to
    get the equivalent of 05216 (total precip flux) in a
    non convection resolving model. Input:
    cube_list = a list of cubes which contains 4203 and 4204 (large
    scale rain and snow)
    Output:
    cube_list = the same cube list but with the new precip field appended
    '''
    print('Downscaled model is convection permitting . . .')
    print('Correcting for stash difference between downscaled and driving model . . . ')
    print('large scale rainfall (m01s04i203) + large scale snowfall (m01s04i204) = total precipitation (m01s05i216)')
    ls_rain_constraint = iris.AttributeConstraint(STASH="m01s04i203")
    ls_rain_cube, = cube_list.extract(ls_rain_constraint)
    ls_snow_constraint = iris.AttributeConstraint(STASH="m01s04i204")
    ls_snow_cube, = cube_list.extract(ls_snow_constraint)
    var_cube = cube_list[0].copy()
    var_cube.data = ls_snow_cube.data + ls_rain_cube.data  # add snow and rain to get total precip
    precip_stash = iris.fileformats.pp.STASH(1, 5, 216)  # create instance of stash
    var_cube.attributes['STASH'] = precip_stash  # set stash that will be used for driving model
    var_cube.rename('precipitation_flux')
    cube_list.append(var_cube)  # add calculated precip flux to cube list
    cube_list.remove(ls_rain_cube)  # remove large scale rain from cube list
    cube_list.remove(ls_snow_cube)  # remove large scale snow from cube list
    return cube_list


if __name__ == '__main__':

    #################   read in data and extract domain #####################

    # This loops over the years and months of interest.
    # Wet season in Zambia is Oct to Mar so only months 1,2,3,10,11,12 included

    for yr in (1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007):
        # for yr in (2003,2004):
        # for months 1 to 3
        for mn in ('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'):
            cp4_data = '/data/users/sburgan/CP4A/precip/u-' + stash + '/' + stash + 'a.pa' + str(yr) + str(mn) + '*.pp'
            cp4_data_check = '/data/users/sburgan/CP4A/precip/u-' + stash + '/' + stash + 'a.pa' + str(yr) + str(
                mn) + '01_00.pp'
            if (path.exists(cp4_data_check) == True):
                cp4_cubelist = load_CP4_cubelist(cp4_data, 'longitude', 'latitude')
                # mask cube with shapefile
                # add together rain and snow (probably not necessary for Zambia but just in case)
                cp4_cube = convec_2_precipflux(cp4_cubelist)[0]
                # this narrows down even further to the Lusaka area
                # some experimentation needed to find the right area for other analyses
                # cp4_small = cp4_cube[:,34:64,17:64]                                               #Scott - Blanked out for now as already narrowed down

                # calculate daily means
                # print statement can be removed after running once to check cube metadata looks right
                cp4_cube_daily = cp4_cube.aggregated_by('day_of_month', iris.analysis.MEAN)
                print(cp4_cube_daily)

                # save out the daily files
                # for months 1-3
                iris.save(cp4_cube_daily, data_dir + '' + region + '_CP4A_daily_precip_' + str(yr) + str(mn) + '.pp')



