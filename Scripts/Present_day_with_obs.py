# create scatter plots of two variables from GCM models against each other and observations
# a new script using the new shared module
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import iris
import iris.coord_categorisation as iccat
import iris.analysis.cartography
import iris.experimental.equalise_cubes
import shared
import numpy as np

import ascend
from ascend import shape
import iris.quickplot as qplt


def get_constraints(years):
    # return the time and area constraints needed for this script
    time_con = iris.Constraint(time=lambda t: years[0] <= t.point.year <= years[-1])

    # This constraint accepts a gridbox if it's midpoint is within the stated coords
    area_con = iris.Constraint(
        latitude=lambda x: DOM_LAT[0] <= x <= DOM_LAT[1],
        longitude=lambda x: DOM_LON[0] <= x <= DOM_LON[1]
    )

    # this constraint only accepts a grid box if it all lies within the target co-ordinates stated
    #area_con = iris.Constraint(
    #    latitude=lambda x: (DOM_LAT[0] <= min(x.bound)) & (max(x.bound) <= DOM_LAT[1]),
    #    longitude=lambda x: (DOM_LON[0] <= min(x.bound)) & (max(x.bound) <= DOM_LON[1])
    #)

    # This constraint will accept a grid box if any part of it lies inside the target box specified.
    #area_con = iris.Constraint(
    #    latitude=lambda x: (DOM_LAT[0] <= min(x.bound) < DOM_LAT[1]) or (DOM_LAT[0] < max(x.bound) <= DOM_LAT[1]),
    #    longitude=lambda x: (DOM_LON[0] <= min(x.bound) < DOM_LON[1]) or (DOM_LON[0] < max(x.bound) <= DOM_LON[1])
    #)

    return (time_con, area_con)

def collapse_data(data, season):
    '''
    Collapse cube of data to a single point for that season
    :param data:
    :param season:
    :return: scalar cube of data with time and area co-ords collapsed
    '''
    iccat.add_season_membership(data, 'time', season=season)
    c_data = data.extract(iris.Constraint(season_membership=True))
    data.remove_coord('season_membership')

    # collapse time
    c_data = c_data.collapsed('time', iris.analysis.MEAN)

    # collapse area
    grid_weights = iris.analysis.cartography.area_weights(c_data)
    c_data = c_data.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_weights)

    return c_data

def compute_mean_tas_pr(models, years, seasons, shp=None, regrid_grid=None):
    '''

    :param models: CMIP5 models to use (list of strings)
    :param years: start and end year to work with (list of length 2)
    :param seasons: seasons to calculate statistics for (list of seasons)
    :param shp (optional): ascend shape object to use as mask for averaging
    :param regid_grid (optional): Grid to regrid model data to
    :return: Dict containing results
    '''
    vars = ['pr', 'tas']
    #vars = ['pr', 'ts']

    # initialise data structure to store results
    s_dict = dict.fromkeys(seasons)
    for s in seasons:
        s_dict[s] = dict.fromkeys(vars)
        for v in vars:
            s_dict[s][v] = dict.fromkeys(models)

    # loop over models doing calculations etc.
    for m in models:
        # get list of CMIP5 files for model
        for v in vars:
            if years[-1] > 2005:
                experiment = FUT_SCENARIO
            else:
                experiment = 'historical'
            print('looking for {}'.format(m))

            # load data
            ensemble = 'r1i1p1'
            cons = get_constraints(years)
            m_gridded_data = shared.ManageCMIP_Gridded_Data(m, v, experiment, ensemble)
            m_gridded_data.load(m_gridded_data.load(cons[0] & cons[1], regrid_grid))

            m_data = m_gridded_data.data

            # process data
            # mask with shapefile
            # TODO Should masking with shapefiles become part of the gridded_data class?
            if not shp == None:
                shp.mask_cube_inplace(m_data)

            # compute and store mean
            for s in seasons:
                mn = collapse_data(m_data, s)

                # store..
                s_dict[s][v][m] = np.asscalar(mn.data)

    return s_dict

def tas_pr_scatter(precip, temp, title, chosen_models=None, discarded_models=None):
    '''
    Make a scatter plot of precip vs temp
    :param precip: dict of precip values or tuple of two dicts (function will compute difference)
    :param temp: dict of temp values or tuple of two dicts (function will compute difference)
    :param title: title for plot
    :param save_path: filename to save as
    :param chosen_models: List of models to highlight in bold on the plot
    :return ax: axes object of the plot produced
    '''
    if type(precip) == dict:
        models = precip.keys()
    else:
        # precip should be a tuple or list
        models = precip[0].keys()

    plt.figure(figsize=(20, 10))
    for m in models:
        if type(precip) == dict:
            x = precip[m]
        else:
            x = precip[0][m] - precip[1][m]

        if type(temp) == dict:
            y = temp[m]
        else:
            y = temp[0][m] - temp[1][m]
        plt.scatter(x,y)
        if m in chosen_models:
            plt.annotate(m, (x, y), fontweight='bold')
        elif m in discarded_models:
            plt.annotate(m, (x, y), fontstyle='italic', fontweight='light', color='0.5')
        else:
            plt.annotate(m, (x, y))

    ax = plt.gca()
    ax.set_title(title)
    ax.set_xlabel('Precipitation (mm day-1)')
    ax.set_ylabel('Temperature ($^o$C)')

    return ax

## Main script.. User variables go here
# TODO ability to use two different domains. One to cover variable 1, and one for variable 2
#DOM_LAT = [3, 45]
#DOM_LON = [20, 70]
#DOM_LAT = [16.5, 20]
#DOM_LON = [-64, -59]
# The co-ords of our grid should line up with the 2.5x2.5 size grid being used for regridding.
DOM_LAT = [-45, 12]
DOM_LON = [-17, 63]

SEASONS = ['DJF', 'MAM', 'JJA', 'SON']
#SEASONS = ['MJ', 'JA', 'SON', 'DJFMA']
#SEASONS = ['JJAS', 'NDJFMA']
#SEASONS = ['JFMAMJJASOND']
FUT_SCENARIO = 'rcp85'

## End of user variables-------------------

# TODO change name of this variable and allow user to pick whatever country they want... Or in fact ANY area defined by a shapefile
# TODO also need to add flexibility to NOT use a shapefile and just stick to chosen area..
# get the KSA shape object
#country_shapes = shape.load_shp(ascend.EXAMPLE_NATURAL_EARTH)
#country = country_shapes.filter(BRK_NAME = 'Tanzania')[0]


models = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC',
          'CMCC-CM', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',
          'GISS-E2-H', 'GISS-E2-R', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR',
          'IPSL-CM5B-LR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC5', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'NorESM1-M']


chosen_models = []
discarded_models = ['ACCESS1-0']


# TODO - Consider regridding all models (and obs) to a common grid for analysis. (2.5 x 2.5 = CMAP grid)
# TODO - Check this works. and chosen domain aligns with 2.5x2.5 grid
# would need to extract a larger (extra 2.5 in each direction) than required domain (from that specified) first
# for extraction.. then do proper extraction considering grid bounds.

# load obs


cons = get_constraints([1979, 2005])

pr_gd = shared.cmap_pr_gd
pr_gd.load(cons[0] & cons[1])
#pr_gd.load(constraint=cons,shp_object = country)
pr_obs = pr_gd.data.copy()

temp_gd = shared.cru_tmp_gd
temp_gd.load(cons[0] & cons[1])
#temp_gd.load(constraint=cons, shp_object= country)
temp_obs = temp_gd.data.copy()

f_dict = compute_mean_tas_pr(models, [2070, 2099], SEASONS, shp=None, regrid_grid=pr_obs)
p_dict = compute_mean_tas_pr(models, [1979, 2005], SEASONS, shp=None, regrid_grid=pr_obs)

save_path = '/net/home/h05/sburgan/Documents/HORIZON2020_FOCUS/Plots/Scatter_plots'

for s in SEASONS:
	ax = tas_pr_scatter((f_dict[s]['pr'], p_dict[s]['pr']), (f_dict[s]['tas'], p_dict[s]['tas']),
	                        '{} {}: change 1979-2005 to 2070-2099'.format(s, FUT_SCENARIO),
	                        chosen_models=chosen_models, discarded_models=discarded_models)
	plt.savefig('{}/{}_{}_change_{}_{}.svg'.format(save_path, FUT_SCENARIO, s, DOM_LAT, DOM_LON))
	plt.close()

for s in SEASONS:
    # average over area and time
    mn_pr = collapse_data(pr_obs, s)
    mn_tmp = collapse_data(temp_obs, s)

    mn_pr = np.asscalar(mn_pr.data)
    mn_tmp = np.asscalar(mn_tmp.data)
    # make scatter plot
    ax = tas_pr_scatter(p_dict[s]['pr'], p_dict[s]['tas'], '{} - average 1979-2005'.format(s),
                        chosen_models=chosen_models, discarded_models=discarded_models)
    ax.scatter(mn_pr, mn_tmp, marker='X')
    ax.annotate('GCPC / CRU', (mn_pr, mn_tmp), fontweight='bold')
    
    if s == 'JJAS':
        # plot vertical line which is cut off for excluded models..
        ax.axvline(12.5, linestyle='--', color='k')
    plt.savefig('{}/{}_future_{}_{}.svg'.format(save_path, s, DOM_LAT, DOM_LON))
    plt.close()
