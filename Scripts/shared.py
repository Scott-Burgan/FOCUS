# A module containing shared objects and functions for GCM selection work
import iris
import iris.coord_categorisation as iccat
import iris.experimental.equalise_cubes

import cf_units

# UCC libraries for loading CMIP models
# Documentation on installing and using the ucc libraries at
# http://www-hc/~ucc/ucc_python/html/installing.html
import ucc.contrib.fcfa.io_metadata as cmipio

import numpy as np

import re
import calendar
from collections import namedtuple
import ascend
from ascend import shape

variable_map = namedtuple('variable_map', ['ref_gd', 'cmip_name'])

class Domain:
    """
    Class to contain information about a Domain to be used for analysis of gridded data
    """
    def __init__(self, lat, lon, description=None):
        """Create an instance of the Domain class

        Arguments:
            lat {tuple} -- 2 element tuple of latitude in decimal
            lon {tuple} -- 2 element tuple of longitude in decimal
            description {string} -- A plain text description of the Domain
        """
        self.lat = lat
        self.lon = lon
        self.description = description

    @property
    def constraint(self):
        """Return an iris Constraint object for this Domain

        Returns:
            iris.Constraint -- iris.Constraint object for this Domain.
        """
        # return an iris constraint for this Domain
        # the logic requires bounds of gridboxes to be inside the Domain, not just centre points.
        area_con = iris.Constraint(
            latitude=lambda x: (self.lat[0] <= min(x.bound) < self.lat[1]) or (self.lat[0] < max(x.bound) <= self.lat[1]),
            longitude=lambda x: (self.lon[0] <= min(x.bound) < self.lon[1]) or (self.lon[0] < max(x.bound) <= self.lon[1])
        )
        return area_con


class Gridded_Data:
    """
    Class containing gridded data but also functionality around loading / post processing it properly
    """
    def __init__(self):
        self._load_fn = lambda: None
        self.description = None
        self.constraint = None

    def _extract_and_regrid(self, target, mdtol=0.5):
        '''
        Regrid cube onto target. Performing clever extraction first saves time / memory
        :param cube: cube to regrid (iris cube)
        :param target: Target grid to regrid to. (Iris cube with latitude and longitude co-ords)
        :return: Data loaded from path, regridded to target grid (iris cube)
        '''

        # create area constraint for loading, we use this lambda technique to accept if ANY part of the grid box overlaps
        # with the target grid
        in_bounds = lambda coord, value: np.min(coord.bounds) <= value <= np.max(coord.bounds)

        area_con = iris.Constraint(
            latitude=lambda x: in_bounds(target.coord('latitude'), np.min(x.bound)) or
                               in_bounds(target.coord('latitude'), np.max(x.bound)),
            longitude=lambda x: in_bounds(target.coord('longitude'), np.min(x.bound)) or
                                in_bounds(target.coord('longitude'), np.max(x.bound))
        )

        # crop target grid
        target = target.extract(area_con)
        # ideally would not need to extract from source cube here
        # (as regridding would shrink down) - but extraction helps to
        # save memory
        cube = self.data.extract(area_con)

        # now regrid to target
        if isinstance(cube, iris.cube.CubeList):
            for i in range(len(cube)):
                cube[i] = cube[i].regrid(target, iris.analysis.AreaWeighted(mdtol))
        else:
            cube = cube.regrid(target, iris.analysis.AreaWeighted(mdtol))

        return cube


    @staticmethod
    def _process_coords(data):
        '''
        perform standard adjustments to coords to aid with
        use in scripts later.
        '''
        # sort out longitude coord if necessary
        if data.coord('longitude').points[-1] > 180:
            data.coord('longitude').bounds = None
            data = data.intersection(longitude=(-180, 180))
            data.coord('longitude').guess_bounds()

        for c in ('latitude', 'longitude'):
            if not data.coord(c).has_bounds():
                data.coord(c).guess_bounds()

            # remove coord_system data
            data.coord(c).coord_system = None

        # add month_number co-ordinate if needed
        if not data.coords('month_number'):
            iccat.add_month_number(data, 'time')

        return data

    def get_area_weights(self):
        '''Calculate the area weights to calculate area averaged means'''

        if isinstance(self.data, iris.cube.CubeList):
            wts = iris.analysis.cartography.area_weights(self.data[0])
        else:
            wts = iris.analysis.cartography.area_weights(self.data)
        return wts

    def load(self, constraint=None, regrid=None, shp_object=None):
        '''
        Load the gridded data into memory
        :param constraint: any Iris constraints to apply to the cube
        :param regrid: If supplied, regrid the loaded data onto this grid (iris cube)
        :param shp_object: If supplied mask data according to this shape object
        :return:
        '''
        if self.constraint == None:
            self.constraint = constraint
        else:
            self.constraint = constraint & self.constraint

        self.regrid = regrid
        self.shp_object = shp_object

        # this sets the self.data attribute
        self._load_fn()

        if regrid:
            self.data = self._extract_and_regrid(regrid)
            if isinstance(self.data, iris.cube.CubeList):
                assert bounds_str(self.data[0]) == bounds_str(regrid)
            else:
                assert bounds_str(self.data) == bounds_str(regrid)

        # Execute this block if a shape file is given
        if shp_object:
            shp_object.mask_cube_inplace (self.data)

        if isinstance(self.data, iris.cube.CubeList):
            for c in self.data:
                c.attributes['gd_description'] = self.description
        else:
            self.data.attributes['gd_description'] = self.description


class Obs_Gridded_Data(Gridded_Data):
    '''
    Class for gridded data from observations
    '''
    def __init__(self, path, description, load_fn):
        self.path = path
        self.description = description
        self._load_fn = eval('self._load_{}'.format(load_fn))
        self.data = None
        self.constraint = None


    def _load_reynolds(self):
        '''
        Load reynolds SST data
        '''

        data = iris.load(self.path)
        data = data.concatenate_cube()

        data.convert_units('celsius')

        # make adjustments to co-ords
        data = self._process_coords(data)

        self.data = data

        if self.constraint:
            self.data = self.data.extract(self.constraint)


    def _load_cru_tmp(self):
        '''
        Load CRU temperature data
        '''

        data = iris.load_cube(self.path, 'near-surface temperature')

        # sort out units
        if data.units == 'degrees Celsius':
            data.units = 'Celsius'

        data.convert_units('celsius')

        data = self._process_coords(data)

        self.data = data

        # add year metadata
        iccat.add_year(self.data, 'time')

        if self.constraint:
            self.data = self.data.extract(self.constraint)


    def _load_cmap(self):
        self.data = iris.load_cube(self.path)

        # do standard co-ord processing
        self.data = self._process_coords(self.data)

        if self.constraint:
            self.data = self.data.extract(self.constraint)

        assert self.data.units == 'mm/day'


    def _load_gpcp(self):
        self.data = iris.load_cube(self.path, 'Average Monthly Rate of Precipitation')

        # do standard co-ord processing
        self.data = self._process_coords(self.data)

        if self.constraint:
            self.data = self.data.extract(self.constraint)

        assert self.data.units == 'mm/day'



class ERA_Gridded_Data(Gridded_Data):
    '''
    ERA based gridded data (from either ERA5 or ERA-Int)
    '''
    def __init__(self, path, vars, description, constraint=None):
        self.path = path
        self.vars = vars
        self.description = description
        self._load_fn = self._load_era
        self.data = None
        self.constraint = constraint


    def _load_era(self):
        # load the ERA data

        # check how many variables we are loading
        if len(self.vars) == 1:
            cube = iris.load(self.path.format(self.vars[0]), constraints=self.constraint)
            assert len(cube) == 1
            self.data = self._process_coords(cube[0])
        else:
            self.data = iris.cube.CubeList()
            for v in self.vars:
                cube = iris.load(self.path.format(v), constraints=self.constraint)
                assert len(cube) == 1

                self.data.append(self._process_coords(cube[0]))

        if self.vars[0] == 'mslp':
            self.data.convert_units('hPa')

    def convert_to_speed(self):
        if type(self.data) is iris.cube.CubeList:
            assert len(self.data) == 2
            # convert the u and v vectors to a single cube of windspeed
            self.data = calc_wind_speed(self.data[0], self.data[1])
        else:
            raise TypeError('Expecting a cubelist to calculate wind speed from')


class ManageCMIP_Gridded_Data(Gridded_Data):
    '''
    gridded_data class with special functionality for managecmip data
    '''
    def __init__(self, model, var, experiment, ensemble='r1i1p1'):
        # use the ucc library to get paths of cmip5 files to load
        cmip5data = cmipio.RecentCmipData(model=model, var=var, ens=ensemble, exp=experiment, freq='mon', domain='atmos')

        # get names of files, but apply safety regex to guard against rogue files that don't end YYMMDD.nc
        files = [i for i in cmip5data.files() if not re.search(r'\d{6}\.nc$', i) == None]

        self.path = files
        self.description = '{} {} {} data'.format(model, experiment, var)
        self.var = var
        self._load_fn = self._load_cmip
        self.data = None
        self.constraint = None

    def _load_cmip(self):
        '''
        Load champ cmip data
        '''
        var_name_dict = {
            'pr': 'precipitation_flux',
            'tas': 'air_temperature',
            'ts': 'surface_temperature',
            'ua': 'eastward_wind',
            'va': 'northward_wind',
            'psl': 'air_pressure_at_sea_level'
        }

        print('Loading:\n{}'.format(self.path))
        data = iris.load(self.path, callback=cmipio.callback_champ, constraints=var_name_dict[self.var])

        # adjust longitude co-ordinate if necessary
        for i in range(len(data)):
            data[i] = self._process_coords(data[i])

        # now, various other bodges to get cubes to concatenate into 1
        # equalise attributes and time
        iris.experimental.equalise_cubes.equalise_attributes(data)
        iris.util.unify_time_units(data)

        if self.constraint:
            data = data.extract(self.constraint)

        # some cubes (mainly precip) have time as an aux coord rather than dim, so promote
        for i in range(len(data)):
            try:
                iris.util.promote_aux_coord_to_dim_coord(data[i], 'time')
            except ValueError:
                # this occurs with HADGEM-2 (duplicate 12/2099)... try removing the last cube..
                print('WARNING - POSSIBLE PROBLEM')
                if data[-1].coord('time') == data[-2][-1].coord('time'):
                    del data[-1]
                elif len(data[-1].shape) != len(data[0].shape):
                    data[-1] = iris.util.new_axis(data[-1], 'time')
                else:
                    raise ValueError('Unknown problem with time coord for cube:\n{}'.format(data[i]))

        # now, attempt to concatenate
        try:
            data = data.concatenate_cube()
        except iris.exceptions.ConcatenateError:
            print('WARNING - POSSIBLE PROBLEM')
            # this checks if the time coord of the last cube is identical to that of the penultimate one
            # if so, it's a duplicate and can be removed
            if data[-1].coord('time') == data[-2].coord('time'):
                del data[-1]
            # another possible error is the final cube being only one timestep and therefore only 2 dimensions,
            # so it won't concantenate with the others..
            elif len(data[-1].shape) != len(data[0].shape):
                # check if this last time point is actually a duplicate of the last time point in the other cubes
                if data[-1].coord('time') == data[-2].coord('time')[-1]:
                    del data[-1]
                else:
                    data[-1] = iris.util.new_axis(data[-1], 'time')
            else:
                raise iris.exceptions.ConcatenateError('Unknown concatenate problem')
            # try again
            data = data.concatenate_cube()

        # at this point, data should now be a single cube.
        self.data = data

        if self.var == 'pr':
            # adjust precipitation units to mm day-1
            assert self.data.units == cf_units.Unit('kg m-2 s-1')
            self.data = self.data * 86400
            self.data.units = 'mm day-1'
        elif self.var == 'psl':
            self.data.convert_units('hPa')
        elif self.var in ('ts', 'tas'):
            self.data.convert_units('celsius')


class ManageCMIP_Wind_Gridded_Data(Gridded_Data):
    def __init__(self, u_gd, v_gd, plevel):
        """Contains two ManageCMIP_Gridded_Data objects. Representing u and v values.

        Arguments:
            u_gd {managecmip_gridded_data} -- managecmip u data
            v_gd {managecmip_gridded_data} -- managecmip v data
            plev {int} -- pressure level in Pa
        """
        self.u_data = u_gd
        self.v_data = v_gd

        assert u_gd.constraint == v_gd.constraint
        plev_con = iris.Constraint(air_pressure=plevel)
        self.constraint = u_gd.constraint & plev_con

        self.description = u_gd.description + ' + ' + v_gd.description

        self.path = u_gd.path + v_gd.path

        self._load_fn = self._load_wind


    def _load_wind(self):
        # load the data from both objects and put into a cubelist
        self.u_data.load(constraint=self.constraint, regrid=self.regrid)
        self.v_data.load(constraint=self.constraint, regrid=self.regrid)

        self.data = iris.cube.CubeList([self.u_data.data, self.v_data.data])


    def convert_to_speed(self):
        # convert the u and v vectors to a single cube of windspeed
        self.data = calc_wind_speed(self.data[0], self.data[1])

class HadGEM3RA_Gridded_Data(Gridded_Data):
    def __init__(self, model_id, var, datadir):
        """ Class for loading model data that has been downloaded from MASS

        """
        self.path = '{}/{}/*.pp'.format(datadir, model_id)
        self.description = '{} {} data'.format(model_id, var)
        self.var = var
        self._load_fn = self._load_HadGEM3RA
        self.data = None
        self.constraint = None

    def _load_HadGEM3RA(self):
        '''
        Load data
        '''

        var_name_dict = {
            'rainfall': 'precipitation_flux',
            'u_925': 'x_wind',
            'mslp': 'air_pressure_at_sea_level'
        }

        print('Loading:\n{}'.format(self.path))
        data = iris.load(
            self.path, callback=HadGEM3RA_load_callback, constraints=var_name_dict[self.var]
            )

        assert len(data) == 1
        self.data = data[0]

        if self.constraint:
            self.data = self.data.extract(self.constraint)

        if self.var == 'rainfall':
            # adjust precipitation units to mm day-1
            assert self.data.units == cf_units.Unit('kg m-2 s-1')
            self.data = self.data * 86400
            self.data.units = 'mm day-1'
        elif self.var == 'mslp':
            self.data.convert_units('hPa')


def mask_both(cube1, cube2):
    '''
    Compute a mask which masks both cube1 and cube2 and apply to both cubes
    Operate in place. Both cubes must have identical shapes
    '''
    mask_both = np.logical_or(cube1.data.mask, cube2.data.mask)
    cube1.data.mask = mask_both
    cube2.data.mask = mask_both


def bounds_str(cube):
    '''
    :return: String representation of the bounds of the gridded data
    '''
    return "[{}N {}N] [{}E {}E]".format(
        np.min(cube.coord('latitude').bounds),
        np.max(cube.coord('latitude').bounds),
        np.min(cube.coord('longitude').bounds),
        np.max(cube.coord('longitude').bounds)
    )


def calc_wind_speed(u, v, norm=False):
    '''
    Calculate wind speed from u and v vectors
    :param u: u values (iris cube)
    :param v: v values (iris cube)
    :return: windspeed (iris cube)
    '''
    windspeed = (u ** 2 + v ** 2) ** 0.5
    windspeed.rename('windspeed')

    if norm:
        windspeed = norm_vector(windspeed)

    return windspeed

def norm_vector(d):
    # normalize the data in vector d
    # d should be an Iris cube with 1d data
    # compute the mean
    res = d.copy()
    d_mean = np.mean(res.data)
    # subtract the mean from the data (the data is now an anomaly from the mean)
    res.data = d.data - d_mean
    # compute the standard deviation of the anomaly
    anom_sd = np.std(d.data)
    # divide by the standard deviation
    res.data = res.data / anom_sd

    return res

def HadGEM3RA_load_callback(cube, field, filename):
    if cube.name() == 'precipitation_flux':
        if cube.cell_methods[0].intervals[0] == '24 hour':
            raise iris.exceptions.IgnoreCubeException

    for c_name in ('latitude', 'longitude'):
        cube.coord(c_name).guess_bounds()
        cube.coord(c_name).coord_system = None

    iccat.add_month_number(cube, 'time')


# gridded_data objects... This data would probably be best eventually being loaded from a config file
reynolds_sst_gd = Obs_Gridded_Data(
    path='/project/ciid/obs_datasets/global/reynoldsSST/reynolds_sst_????_monthly.nc',
    description='Reynolds SST re-analysis',
    load_fn='reynolds'
)

cru_tmp_gd = Obs_Gridded_Data(
    path='/project/ciid/obs_datasets/global/CRU/ts4.02/tmp/cru_ts4.02.1901.2017.tmp.dat.nc',
    description='CRU mean temperature data. v4.02',
    load_fn='cru_tmp'
)

cmap_pr_gd = Obs_Gridded_Data(
    path='/project/ciid/obs_datasets/global/CMAP/monthly/pr/precip.mon.mean.1979-2019.nc',
    description='CMAP Precipitation data',
    load_fn='cmap'
)

gpcp_pr_gd = Obs_Gridded_Data(
    path='/project/ciid/obs_datasets/global/GPCP/monthly/pr/precip.mon.mean.1979-2019.nc',
    description='GPCP Precipitation data',
    load_fn='gpcp'
)

era5_925w = ERA_Gridded_Data(
    path='/project/ciid/obs_datasets/global/ERA5/process/ERA5_{}_monthly.nc',
    description='925 hPa winds from ERA5',
    vars=['u925', 'v925'],
    constraint=iris.Constraint(pressure_level=925)
)

era5_925_u = ERA_Gridded_Data(
    path='/project/ciid/obs_datasets/global/ERA5/process/ERA5_{}_monthly.nc',
    description='925 hPa winds from ERA5',
    vars=['u925'],
    constraint=iris.Constraint(pressure_level=925)
)

era_i_850w = ERA_Gridded_Data(
    path='/project/ciid/projects/ARAMCO/GCM_selection/ERA-Interim/ERAInterim.{}.monthly.nc',
    description='850 hPa winds from ERA-Int',
    vars=['u850', 'v850'],
    constraint=iris.Constraint(air_pressure=85000)
)

era_i_mslp = ERA_Gridded_Data(
    path='/project/ciid/projects/ARAMCO/GCM_selection/ERA-Interim/ERAInterim.{}.monthly.nc',
    description='MSLP from ERA-Int',
    vars=['mslp']
)

# dictionary of variable maps, to gd objects or CMIP5 variable names
variables = {
    'rainfall': variable_map(gpcp_pr_gd, 'pr'),
    'temperature': variable_map(cru_tmp_gd, 'tas'),
    'wind_925': variable_map(era5_925w, ('ua', 'va')),
    'u_925': variable_map(era5_925_u, 'ua'),
    'sst': variable_map(reynolds_sst_gd, 'ts'),
    'mslp': variable_map(era_i_mslp, 'psl')
}

def get_gd_data(var, model=None, experiment=None, path=None):
    """Function to setup data given variable info
    this is different depending on if we are loading model data or reference data....

    Arguments:
        var {string} -- Variable e.g. pr, psl ...

    Keyword Arguments:
        model {string} -- Model to use, either CMIP5 name, MO run ID (e.g. ba642)
        or none (observations) (default: {None})
        experiment {string} -- For CMIP5 models, experiment type.
        (historical, rcp85 etc.) (default: {None})
        path {string} -- For MO models, parent path to data files (default: {None})

    Returns:
        [Gridded_Data] -- Gridded_data class object
    """
    if model == None:
        gd_data = variables[var].ref_gd
    elif re.search(r'^[a-z]{2}-[a-z]{2}\d{3}$', model):
        # supplied model string matches MO run-id
        gd_data = HadGEM3RA_Gridded_Data(model, var, path)
        if '_' in var:
            plevel = int(var.split('_')[-1])
            gd_data.constraint = gd_data.constraint & iris.Constraint(pressure=plevel)
    else:
        # a CMIP5 model has been supplied, so we load the CMIP5 model data
        if type(variables[var].cmip_name) is str:
            gd_data = ManageCMIP_Gridded_Data(
                model=model, var=variables[var].cmip_name, experiment=experiment
                )
            if '_' in var:
                plevel = int(var.split('_')[-1]) * 100
                gd_data.constraint = gd_data.constraint & iris.Constraint(air_pressure=plevel)

        else:
            assert len(variables[var].cmip_name) == 2

            gd_data = []
            for v in variables[var].cmip_name:
                gd_data.append(
                    ManageCMIP_Gridded_Data(
                        model=model, var=v, experiment=experiment
                        )
                )
            plevel = int(var.split('_')[-1]) * 100
            gd_data = ManageCMIP_Wind_Gridded_Data(gd_data[0], gd_data[1], plevel)

    return gd_data

def get_constraints(years):
    # return the time and area constraints needed for this script
    time_con = iris.Constraint(time=lambda t: years[0] <= t.point.year <= years[-1])

    # This constraint accepts a gridbox if it's midpoint is within the stated coords
    #area_con = iris.Constraint(
    #    latitude=lambda x: DOM_LAT[0] <= x <= DOM_LAT[1],
    #    longitude=lambda x: DOM_LON[0] <= x <= DOM_LON[1]
    #)

    # this constraint only accepts a grid box if it all lies within the target co-ordinates stated
    #area_con = iris.Constraint(
    #    latitude=lambda x: (DOM_LAT[0] <= min(x.bound)) & (max(x.bound) <= DOM_LAT[1]),
    #    longitude=lambda x: (DOM_LON[0] <= min(x.bound)) & (max(x.bound) <= DOM_LON[1])
    #)

    # This constraint will accept a grid box if any part of it lies inside the target box specified.
    area_con = iris.Constraint(
        latitude=lambda x: (DOM_LAT[0] <= min(x.bound) < DOM_LAT[1]) or (DOM_LAT[0] < max(x.bound) <= DOM_LAT[1]),
        longitude=lambda x: (DOM_LON[0] <= min(x.bound) < DOM_LON[1]) or (DOM_LON[0] < max(x.bound) <= DOM_LON[1])
    )

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
    #vars = ['pr', 'tas']
    vars = ['pr', 'ts']

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
            if m == 'EC-EARTH':
                ensemble = 'r7i1p1'
            else:
                ensemble = 'r1i1p1'

            # load data
            cons = get_constraints(years)
            m_gridded_data = shared.ManageCMIP_Gridded_Data(m, v, experiment, ensemble)
            m_gridded_data.load(cons[0] & cons[1], regrid_grid)

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


