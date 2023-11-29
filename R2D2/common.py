from . import read
from . import google
from . import resolution
from . import models
from . import sync

class R2D2_data():    
    '''
    R2D2_data is a class for reading R2D2 data

    Attributes:
        p (dic): basic parameters
        qs (dic): 2D data at selected height
        qq (dic): 3D full data 
        qt (dic): 2D data at constant optical depths
        t (float): time
        vc (dic): data of on the fly analysis

        models (dic): Model S based stratification
    '''
R2D2_data.__init__         = read.init
R2D2_data.read_qq_select   = read.read_qq_select
R2D2_data.read_qq_select_z = read.read_qq_select_z
R2D2_data.read_qq          = read.read_qq
R2D2_data.read_qq_variable = read.read_qq_variable
R2D2_data.read_qq_tau      = read.read_qq_tau
R2D2_data.read_time        = read.read_time
R2D2_data.read_vc          = read.read_vc
R2D2_data.read_qq_check    = read.read_qq_check
R2D2_data.read_qq_2d       = read.read_qq_2d
R2D2_data.read_qq_slice    = read.read_qq_slice
R2D2_data.YinYangSet       = read.YinYangSet

R2D2_data.set_cells_gspread  = google.set_cells_gspread
R2D2_data.upgrade_resolution = resolution.upgrade_resolution

R2D2_data.models_init       = models.init

R2D2_data.sync_tau = sync.sync_tau
R2D2_data.sync_select = sync.sync_select
R2D2_data.sync_vc = sync.sync_vc
R2D2_data.sync_check = sync.sync_check
R2D2_data.sync_slice = sync.sync_slice
R2D2_data.sync_all = sync.sync_all

# color blind color set
# Obtained from https://github.com/matplotlib/matplotlib/issues/9460
blue    = '#7b85d4'
orange  = '#f37738'
green   = '#83c995'
magenta = '#d7369e'
gray    = '#c4c9d8'
ash     = '#859795'
yellow  = '#e9d943'
brown   = '#ad5d50'

# Solar parameters
msun = 1.988e33 # g
rsun = 6.957e10 # cm
lsun = 3.828e33 # erg/s
asun = 4.570e9 # yr