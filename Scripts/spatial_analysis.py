import matplotlib as mpl

# for working on SPICE or other displayless systems
mpl.use('Agg')

import iris
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units
import iris.coord_categorisation as iccat

# Documentation on installing and using the ucc libraries at
# http://www-hc/~ucc/ucc_python/html/installing.html
import ucc.contrib.fcfa.io_metadata as cmipio

import sys
import numpy as np
import re

# This script performs plotting of the CMIP5 model data over the domain specified and differences to observations

def extreme_x(data,min,max):

    # Finds cells with values outside of the colorbar range and mark them

    data.update_scalarmappable()
    ax = data.axes
    for p, color, value in zip(data.get_paths(), data.get_facecolors(), data.get_array()):
        x, y = p.vertices[:-2, :].mean(0)
        if np.all(color[:3] > 0.5):
            color = (0.0, 0.0, 0.0)
        else:
            color = (1.0, 1.0, 1.0)
        if(value < min):
            ax.text(x, y, r"$\circ$", ha="center", va="center", color=color)
        if(value > max):
            ax.text(x, y, r"$\circ$", ha="center", va="center", color=color)

def do_plot(data, title, nplots, subplot_n, plot_fn, gmin, gmax, dmin, dmax, type):
    '''
    Take cube and produce plot. This plots 1 plot in the subplot of the active figure
    :param data: iris cube to plot. 2 dimensional. However if wind. Supply a cube list. U and V cubes.
    :param title: title of plot (the function will use this to work out if it's a diff plot or not)
    :param nplots: The total number of rows in the figure that this plot is in..
    :param subplot_n: matplotlib subplot number for this plot
    :param plot_fn: function that will be called to do plotting (a function)
    :param gmin: model/ref scale minimum
    :param gmax: model/ref scale maximum
    :param dmin: diff scale minimum
    :param dmax: diff scale maximum
    :param type: model / ref / diff
    '''
    # TODO could nplots be calculated somehow from the plt.gcf() object?
    # select appropriate axes
    qplt.plt.subplot(nplots, 3, subplot_n)

    # Initialise cmap
    cmap = mpl.cm.get_cmap("bwr")

    if (type == "model"):
        if DO_MASK == True:
            cmap.set_bad('black',1.)
        else:
            cmap.set_bad('white',1.)

    if 'wind' in VAR:
        # calculate wind speed from vectors and plot
        ws = calc_wind_speed(data[0], data[1])
        cmap = mpl.cm.get_cmap("RdBu")
        plot_obj = plot_fn(ws,vmin=gmin,vmax=gmax,cmap=cmap)
        if title[:3] in ('ERA', 'Dif'):
            # don't plot every arrow because plot looks too messy with higher res ERA data
            # TODO: Celverer choice of density of arrows to plot based on number of grid boxes in figure or something
            iplt.quiver(data[0][::2, ::2], data[1][::2, ::2])
        else:
            iplt.quiver(data[0], data[1])
    # different color palette choices depending on variable chosen
    elif VAR == 'rainfall':
        cmap = mpl.cm.get_cmap("RdYlBu")
        if (type == "model"):
            if DO_MASK == True:
                cmap.set_bad('black',1.)
        plot_obj = plot_fn(data,vmin=gmin,vmax=gmax,cmap=cmap)
    else:
        plot_obj = plot_fn(data,vmin=gmin,vmax=gmax,cmap=cmap)

    # set scales
    # TODO (future...) cleverer scales based on data / common sesnse
    if type == 'diff':
        plot_obj.set_clim(vmin=dmin, vmax=dmax)

    # set appropriate contour label for MSLP plots
    if plot_fn == qplt.contour:
        qplt.plt.clabel(plot_obj, fmt='%0d')

    if VAR != 'MSLP':
        if(type == 'model' or type == 'ref'):
            extreme_x(plot_obj,gmin,gmax)
        if(type == 'diff'):
            extreme_x(plot_obj,dmin,dmax)

    # add coastlines and title
    #qplt.plt.gca().coastlines(alpha=0.5)
    qplt.plt.gca().coastlines('50m')
    qplt.plt.title(title)

def plot_all_seasons(m_data, r_data, diff_data, model, plot_type, gmin, gmax, dmin, dmax):
    '''
    This calls the do_plot function to plot all the subplots in the figure
    :param m_data: model data (iris cube)
    :param r_data: reference data, reanalysis or obs (iris cube)
    :param diff_data: difference data diff between model and reference (iris cube)
    :param model: model name (string)
    :param plot_type: plotting function to make plot (function)
    :param gmin: model/ref scale minimum
    :param gmax: model/ref scale maximum
    :param dmin: diff scale minimum
    :param dmax: diff scale maximum
    :return: Nothing. Side effect of running is plot saved to disk.
    '''
    plot_n = 1
    mx = MODELS.index(m)
    # setup figure size variable depending on number of seasons (taller for more seasons)
    qplt.plt.figure(figsize=(5 * len(SEASONS), 15))
    # loop over all the seasons
    # for each season, plot 3 plots, the model, the reference data and the differences
    for s in SEASONS:
        sx = SEASONS.index(s)
        # plot the model data
        type = "model"
        ptitle = 'Model: {}     [Season: {}]'.format(model, s)
        do_plot(m_data[s], ptitle, len(SEASONS), plot_n, plot_type, gmin[mx][sx], gmax[mx][sx], dmin, dmax, type)
        plot_n = plot_n + 1
        if VAR in ('rainfall', 'temperature', 'sst'):
            # rainfall and temperature come from obs datasets
            ptitle = 'Observations     [Season: {}]'.format(s)
        else:
            # other variables from ERA
            ptitle = 'Re-analysis     [Season: {}]'.format(s)
        # plot the reference data (obs or ERA)
        type = "ref"
        do_plot(r_data[s], ptitle, len(SEASONS), plot_n, plot_type, gmin[mx][sx], gmax[mx][sx], dmin, dmax, type)
        plot_n = plot_n + 1
        # plot the difference data
        type = "diff"
        if VAR in ('rainfall', 'temperature', 'sst'):
            ptitle = '(Model - Observations)     [Season: {}]'.format(s)
        else:
            ptitle = '(Model - Re-analysis)     [Season: {}]'.format(s)
        do_plot(diff_data[s], ptitle, len(SEASONS), plot_n, qplt.pcolormesh, gmin[mx][sx], gmax[mx][sx], dmin, dmax, type)
        plot_n = plot_n + 1

    # fix layout
    qplt.plt.tight_layout()
    # save
    #suptitle = "Model {} vs reference {}".format(something,something)
    #qplt.plt.suptitle(suptitle)
    qplt.plt.savefig('{}/{}.svg'.format(
        SAVE_PATH, '{}_{}_{}_{}'.format(
            model, VAR, years[0], years[1])))
    qplt.plt.close()

def calc_wind_speed(u, v):
    '''
    Calculate wind speed from u and v vectors
    :param u: u values (iris cube)
    :param v: v values (iris cube)
    :return: windspeed (iris cube)
    '''
    windspeed = (u ** 2 + v ** 2) ** 0.5
    windspeed.rename('windspeed')
    return windspeed

def apply_landmask(m, m_data):

    # Find and apply CMIP5 landmask to model data

    if (m != 0):
        print("Applying landmask data to model {}...".format(m))

        # Load landmask data: csm1-1-m model is found in piControl only
        if(m == 'bcc-csm1-1-m'):
            cm5mask = cmipio.RecentCmipData(model=m, var='sftlf', exp='piControl', freq='fx')
        else:
            cm5mask = cmipio.RecentCmipData(model=m, var='sftlf', exp='historical', freq='fx')

        mask_file = cm5mask.files()
        print(cm5mask.files())
        masks = iris.load(mask_file,'land_area_fraction')
        mask = masks[0] # Canadian model returns two landmasks...

        nlat = mask.shape[0]
        nlon = mask.shape[1]

        # How many time points does the model data have?
        nt = m_data.shape[0]

        bool_mask_3d = np.ones((nt, nlat, nlon), dtype=bool)

        for i in range(0,nlat):
            for j in range(0,nlon):
                val = mask.data[i][j]
                # Mask Sensitivity
                if (val == 0.0):
                    for t in range(0,nt):
                        bool_mask_3d[t][i][j] = 0
                else:
                    for t in range(0,nt):
                        # Mask these!
                        bool_mask_3d[t][i][j] = 1

        # Mask the model cube with landmask array
        m_data = iris.util.mask_cube(m_data,bool_mask_3d)

    return m_data

def get_range_diff(SEASONS, MODELS, diff_dict, senst):

    # Calculate sensible scale for diff plots

    xmin = np.zeros((len(MODELS),len(SEASONS)))
    xmax = np.zeros((len(MODELS),len(SEASONS)))

    # Find diff min and max

    for s in SEASONS:
        sx = SEASONS.index(s)
        for m in MODELS:
            mx = MODELS.index(m)
            if 'wind' in VAR:
                ws_d = calc_wind_speed(diff_dict[m][s][0], diff_dict[m][s][1])
                xmin[mx][sx] = np.percentile(ws_d.data.compressed(),5)
                xmax[mx][sx] = np.percentile(ws_d.data.compressed(),senst)
            else:
                xmin[mx][sx] = np.percentile(diff_dict[m][s].data.compressed(),5)
                xmax[mx][sx] = np.percentile(diff_dict[m][s].data,senst)

    dmin = np.amin(xmin)
    dmax = np.amax(xmax)

    # Adjust so scale symmetrical about zero point (for white at zero)

    if (dmin < 0 and dmax > 0):
        if (abs(dmin) > dmax):
            dmax = dmin * -1.0
        else:
            dmin = dmax * -1.0

    return dmin, dmax

def get_range(gmin, gmax, SEASONS, MODELS, all_dict, ref_dict, senst):

    # Find sensible scaling for model and ref for a given season and model

    obs_min = [0] * len(SEASONS)
    obs_max = [0] * len(SEASONS)

    # Find the min and max of the observations / re-analysis for each season

    for s in SEASONS:
        sx = SEASONS.index(s)
        if 'wind' in VAR:
            ws_r = calc_wind_speed(ref_dict[s][0], ref_dict[s][1])
            obs_min[sx] = np.percentile(ws_r.data.compressed(),5)
            obs_max[sx] = np.percentile(ws_r.data.compressed(),senst)
        else:
            obs_min[sx] = np.percentile(ref_dict[s].data.compressed(),5)
            obs_max[sx] = np.percentile(ref_dict[s].data.compressed(),senst)

        for m in MODELS:
            mx = MODELS.index(m)

            # For each season, compare model min/max to obs min/max
            # Set global min and max as appropriate

            for s in SEASONS:
                sx = SEASONS.index(s)

                if 'wind' in VAR:
                    ws_m = calc_wind_speed(all_dict[m][s][0],all_dict[m][s][1])
                    gmin[mx][sx] = np.percentile(ws_m.data.compressed(),5)
                    gmax[mx][sx] = np.percentile(ws_m.data.compressed(),senst)

                else:
                    gmin[mx][sx] = np.percentile(all_dict[m][s].data.compressed(),5)
                    gmax[mx][sx] = np.percentile(all_dict[m][s].data.compressed(),senst)

                    #print(np.amin(all_dict[m][s].data),gmin[mx][sx],np.amax(all_dict[m][s].data),gmax[mx][sx])

                    if (obs_min[sx] < gmin[mx][sx]):
                        gmin[mx][sx] = obs_min[sx]

                    if (obs_max[sx] > gmax[mx][sx]):
                        gmax[mx][sx] = obs_max[sx]

    return gmin, gmax

def calc_plot_data(cube_data, season):
    '''
    Calculate data to plot for given season.
    Essentially, extract just what we need and collapse time dimension to mean
    :param cube_data: iris cube of source data
    :param season: season as string
    :return: seasonal mean of data as iris cube
    '''
    # extract domain
    cube_data = cube_data.extract(area_con)
    assert not cube_data is None
    # extract season
    # the season membership approach is used because iris doens't like adding a season co-ordinate that
    # includes single month seasons, or doesn't span the whole year
    if type(cube_data) == iris.cube.CubeList:
        for c in cube_data:
            iccat.add_season_membership(c, 'time', season=season)
    else:
        iccat.add_season_membership(cube_data, 'time', season=season)
    plot_data = cube_data.extract(iris.Constraint(season_membership=True))

    # Units conversion if necessary
    if VAR == 'MSLP':
        plot_data.convert_units('hPa')
    elif VAR == 'rainfall':
        if(plot_data.units != 'mm month -1'):
            # convert_units can't handle mm/month.. but we at least check we are converting what we think we are here
            assert plot_data.units == cf_units.Unit('kg m-2 s-1'), 'Expect units kg m-2 s-1, got {}'.format(plot_data.units)
            plot_data = plot_data * 86400 * 30
            plot_data.units = 'mm month -1'
    elif VAR == 'temperature':
        # ensure units are kelvin
        if plot_data.units == 'degrees Celsius':
            plot_data.units = 'celsius'
            plot_data.convert_units('K')
    # Collapse time dimension to mean
    if type(cube_data) == iris.cube.CubeList:
        for i in range(len(cube_data)):
            plot_data[i] = plot_data[i].collapsed('time', iris.analysis.MEAN)
            cube_data[i].remove_coord('season_membership')
    else:
        plot_data = plot_data.collapsed('time', iris.analysis.MEAN)
        cube_data.remove_coord('season_membership')

    # return
    return plot_data

# TODO potentially refactor this function so it can be shared with the future_projections script
def load_data(data_path, cons, m):
    '''
    Load data for use in the scipt. Mainly contains fixes to deal with CMIP5 archive inconsistencies
    :param data_path: Path to gridded data files to load
    :param cons: iris constraints to apply on loading
    :return: iris cube(s) of data
    '''
    #print('Loading: {}'.format(data_path))
    # Callback used from FCFA library
    m_data = iris.load(data_path, callback=cmipio.callback_champ, constraints=cons)

    if 'wind' in VAR:
        for c in m_data:
            # there was a problem with some cubes having cell_methods whereas others didn't, so we remove them.
            c.cell_methods = None
        m_data = m_data.concatenate()
        if len(m_data) != 2:
            print(m_data)
            sys.exit('Error, expected 2 wind cubes. Got {}'.format(len(m_data)))
    elif VAR == 'temperature':
        # some of the models appear to have a bunch of extra unnecessary cubes.. just get air_temperature ones
        # 3 different possible cube names for temperature depending on if source is CMIP5, obs or ERA.
        m_data = m_data.extract(['air_temperature', '2 metre temperature', 'near-surface temperature'])
        m_data = m_data.concatenate_cube()
    elif VAR == 'sst':
        if(m == 0):
            m_data = m_data.extract(['surface_temperature'])
            m_data = m_data.concatenate_cube()
            m_data.coord('latitude').coord_system = None
            m_data.coord('longitude').coord_system = None
        else:
            m_data = m_data.extract(['surface_temperature'])
            m_data = m_data.concatenate_cube()
    elif VAR == 'rainfall':
        # similarly just extract the rainfall cubes from whatever has been loaded in this case
        m_data = m_data.extract(['precipitation_flux', 'GPCC Monthly total of precipitation', 'Average Monthly Rate of Precipitation'])
        # some models had the time coord as an aux coord, so promote
        for c in m_data:
            iris.util.promote_aux_coord_to_dim_coord(c, 'time')
        m_data = m_data.concatenate_cube()
        if(m == 0):
            m_data.data *= 30
            m_data.units = 'mm month -1'
    elif VAR == 'MSLP':
        # again, we only want the MSLP cubes in this case
        m_data = m_data.extract('air_pressure_at_sea_level')
        m_data = m_data.concatenate_cube()
    else:
        m_data = m_data.concatenate_cube()

    # Apply Landmask
    if (DO_MASK == True):
        m_data = apply_landmask(m, m_data)

    # remove longitude bounds to ensure rolling coords will work

    if type(m_data) == iris.cube.CubeList:
        for i,c in enumerate(m_data):
            for coord in m_data[i].dim_coords:
                coord.bounds = None
            # make sure co-ordinates are rolled to -180 to 180
            m_data[i] = m_data[i].intersection(longitude=(-180, 180))
    else:
        for coord in m_data.dim_coords:
            coord.bounds = None
        # make sure co-ordinates are rolled to -180 to 180
        m_data = m_data.intersection(longitude=(-180, 180))

    # make sure co-ordinates are rolled to -180 to 180
    # m_data = m_data.intersection(longitude=(-180, 180))

    # Add bounds to dim coords if they don't exist
    if type(m_data) == iris.cube.CubeList:
        for i,c in enumerate(m_data):
            for j,d in enumerate(m_data[i].dim_coords):
                if not m_data[i].dim_coords[j].has_bounds():
                    m_data[i].dim_coords[j].guess_bounds()
    else:
        for coord in m_data.dim_coords:
           if not coord.has_bounds():
               coord.guess_bounds()

    return m_data

#---------------------------------------------------------------------------------------
# Main script

# Dict to map user chosen variable to that expected in CMIP5 data, or ERA / Obs
# 1st element in value is CMIP5 var name, second is ERA-I or obs name

VARS_DICT = {
    'MSLP': [['psl'], ['mslp']],
    '850_wind': [['ua', 'va'], ['u850', 'v850']],
    '925_wind': [['ua', 'va'], ['u925', 'v925']],
    'temperature': [['tas'], ['t2m']],
    'rainfall': [['pr'], ['precip']],
    'sst': [['ts'],['surface_temperature']]
}

## User variables here
# -----------------------------
DOM_LAT = [-45, 12]
DOM_LON = [-17, 63]

PAD_DOM_LAT = [DOM_LAT[0]-7, DOM_LAT[1]+7]
PAD_DOM_LON = [DOM_LON[0]-7, DOM_LON[1]+7]

# most CMIP5 historical experiments end in 2005
years = [1981, 2005]

# SEASONS = ['JFMAMJJASOND']
#SEASONS = ['NDJFMA', 'JJAS']
SEASONS = ['DJF', 'MAM', 'JJA', 'SON']
#SEASONS = ['DJFMA', 'MJ', 'JA', 'SON']

# TODO use argparse to handle passing some arguments to script
#VAR = sys.argv[1]
VAR='temperature'

if (VAR == 'sst'):
    DO_MASK = True
else:
    DO_MASK = False

#SAVE_PATH = '/net/home/h01/slines/antigua/r637_antigua_gcm_selection/plots'
SAVE_PATH = '/net/home/h05/sburgan/Documents/HORIZON2020_FOCUS/Plots/Spatial_plots'

MODELS = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', \
          'CMCC-CM', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'EC-EARTH', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',\
          'GISS-E2-H', 'GISS-E2-R', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR',\
          'IPSL-CM5B-LR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC5', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'NorESM1-M']

#MODELS = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM']

# Initialise plot range variables
gmin = np.zeros((len(MODELS),len(SEASONS)))
gmax = np.zeros((len(MODELS),len(SEASONS)))

## End of user variables
# ----------------------------------

# to constrain to chosen years
time_con = iris.Constraint(time=lambda t: years[0] <= t.point.year <= years[-1])
# for extracting 850 winds from wind levels
if VAR == '850_wind':
    lev_con = iris.Constraint(air_pressure=85000)
elif VAR == '925_wind':
    ref_lev_con = iris.Constraint(pressure_level=925)
    lev_con = iris.Constraint(air_pressure=92500)
# to constrain to chosen domain
area_con = iris.Constraint(latitude=lambda x: PAD_DOM_LAT[0] <= x <= PAD_DOM_LAT[1],
                           longitude=lambda x: PAD_DOM_LON[0] <= x <= PAD_DOM_LON[1])

# create list of file(s) to load that contain reference data
DATA_FILES = []
print('Computing reference data\n')
if VAR == 'rainfall':
    # Data from CMAP
    DATA_FILES = '/project/ciid/obs_datasets/global/CMAP/monthly/pr/precip.mon.mean.1979-2019.nc'
elif VAR == 'temperature':
    # Data from CRU
    DATA_FILES = '/project/ciid/obs_datasets/global/CRU/ts4.02/tmp/cru_ts4.02.1901.2017.tmp.dat.nc'
elif VAR == 'sst':
    # Data from Reynolds SST
    DATA_FILES = '/project/ciid/obs_datasets/global/reynoldsSST/reynolds_sst_????_monthly.nc'
elif VAR == '925_wind':
    # Data from ERA-5
    for v in VARS_DICT[VAR][1]:
        DATA_FILES.append('/project/ciid/obs_datasets/global/ERA5/process/ERA5_{}_monthly.nc'.format(v))
else:
    # data comes from ERA-Interim
    for v in VARS_DICT[VAR][1]:
        DATA_FILES.append('/project/ciid/projects/ARAMCO/GCM_selection/ERA-Interim/ERAInterim.{}.monthly.nc'.format(v))

# load reference data, constraining by time
if VAR == '925_wind':
    ref_dat = load_data(DATA_FILES, time_con & ref_lev_con, 0)
else:
    ref_dat = load_data(DATA_FILES, time_con, 0)

# store ref data in a dictionary with seasons as keys
ref_dict = {}
for s in SEASONS:
    # the calc_plot_data function performs the area extraction..
    ref_dict[s] = calc_plot_data(ref_dat, s)
    if 'wind' in VAR:
        ref_res = ref_dict[s][0].shape[0] # Set reference resolution
    else:
        ref_res = ref_dict[s].shape[0] # Set reference resolution
# store model data in dictionary, models as keys
all_dict = {}
print('\nLoading model data\n')
for m in MODELS:
    print('Loading model {}'.format(m))
    DATA_FILES = []
    model_dict = {}
    all_dict[m] = {}

    # get list of CMIP5 files for MODEL, using cmipio library to get file list
    for v in VARS_DICT[VAR][0]:
        cmip5data = cmipio.RecentCmipData(model=m, var=v, ens='r1i1p1', exp='historical', freq='mon',
                                          domain='atmos')
        DATA_FILES.extend(cmip5data.files())

    # There are some rogue files in the champ archive that should not be there (ending 'YYYYMM_1x1.nc')
    # Their deletion has been requested, but...
    # this Regex guards against any files that do not end 'YYYYMM.nc' as expected
    DATA_FILES = [i for i in DATA_FILES if not re.search(r'\d{6}\.nc$',i) == None]

    # load data
    if 'wind' in VAR:
        # extract pressure level for winds
        m_data = load_data(DATA_FILES, time_con & lev_con, m)
        #print(m_data[0])
    else:
        m_data = load_data(DATA_FILES, time_con, m)

    # compute data to plot
    if len(SEASONS) > 0:
        for s in SEASONS:
            # calc plot data
            model_dict[s] = calc_plot_data(m_data, s)
            # store data
            all_dict[m][s] = model_dict[s]

# this block of code computes differences for all models ahead of plotting in a loop.
# TODO work out a way of using the all_abs_max data to compute a sensible scale
print('\nComputing model diffs\n')

diff_dict = {}
for m in MODELS:
    mc = MODELS.index(m)
    # store differences data in a dictionary
    diff_dict[m] = {}
    print("Regridding {}".format(m))
    for s in SEASONS:
        SWITCH_REGRID = False # Default is to regrid model to ref
        if 'wind' in VAR:
            diff = iris.cube.CubeList()
            # Smart regridding: if ref grid is coarser resolution than model
            # Then have diff on model resolution (by regridding ref to model)
            SWITCH_REGRID = True # Forced: wind looks better on coarse resolution
            for i in range(len(all_dict[m][s])):
                if (ref_res < all_dict[m][s][i].shape[0]):
                    SWITCH_REGRID = True
                if (SWITCH_REGRID):
                    obs_rg = ref_dict[s][i].regrid(all_dict[m][s][i], iris.analysis.AreaWeighted())
                    ref_dict[s][i] = obs_rg
                    diff.append(all_dict[m][s][i] - obs_rg)
                else:
                    model_rg = all_dict[m][s][i].regrid(ref_dict[s][i], iris.analysis.AreaWeighted())
                    diff.append(model_rg - ref_dict[s][i])
        else:
            if (ref_res < all_dict[m][s].shape[0]): # If the model is finer resolution than the ref data
                SWITCH_REGRID = True
            if (SWITCH_REGRID):
                obs_rg = ref_dict[s].regrid(all_dict[m][s], iris.analysis.AreaWeighted())
                diff = all_dict[m][s] - obs_rg
            else:
                model_rg = all_dict[m][s].regrid(ref_dict[s], iris.analysis.AreaWeighted())
                diff = model_rg - ref_dict[s]

        # add difference data to the dictionary
        diff_dict[m][s] = diff

area_con = iris.Constraint(latitude=lambda x: DOM_LAT[0] <= x <= DOM_LAT[1],
                           longitude=lambda x: DOM_LON[0] <= x <= DOM_LON[1])

for m in MODELS:
    mc = MODELS.index(m)
    for s in SEASONS:
        # Do a second area extraction to narrow down the wider padded domain
        # This solves the erroneous boundary cell issue.

        if 'wind' in VAR:
            for i in range(len(all_dict[m][s])):
                all_dict[m][s][i] = all_dict[m][s][i].extract(area_con)
                if (mc == 0):
                    ref_dict[s][i] = ref_dict[s][i].extract(area_con)
        all_dict[m][s] = all_dict[m][s].extract(area_con)
        if (mc == 0):
            ref_dict[s] = ref_dict[s].extract(area_con)
        diff_dict[m][s] = diff_dict[m][s].extract(area_con)


# Upper percentile probably needs adapting between data type

if VAR == 'sst':
    senst = 80
elif VAR == 'rainfall':
    senst = 70
else:
    senst = 90

# Get model and ref scale min and max
gmin, gmax = get_range(gmin, gmax, SEASONS, MODELS, all_dict, ref_dict, senst)
# Get diff scale min and max
dmin, dmax = get_range_diff(SEASONS, MODELS, diff_dict, senst)

# now plot
print('Plotting...')
if VAR == 'MSLP':
    plot_type = qplt.contour
else:
    plot_type = qplt.pcolormesh

# loop over all models and plot and save
for m in MODELS:
    print ('Processing... ', m)
    plot_all_seasons(all_dict[m], ref_dict, diff_dict[m], m, plot_type, gmin, gmax, dmin, dmax)
