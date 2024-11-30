from . import read
from . import google
from . import resolution
from . import models
from . import sync

class R2D2_data():    
    '''
    Class for managing R2D2_data

    Attributes
    ----------
        p : dict
            basic parameters.
            
            .. seealso ::
                :class:`R2D2.R2D2_data.__init__`
        qs : dict
            2D data at selected height.
            See :class:`R2D2.R2D2_data.read_qq_select`
        qq : dict
            3D full data.
            See :class:`R2D2.R2D2_data.read_qq`
        qt : dict
            2D data at constant optical depths.
            See :class:`R2D2.R2D2_data.read_qq_tau`
        t : float
            time
            See :class:`R2D2.R2D2_data.read_time`
        vc : dict
            data of on the fly analysis
            See :class:`R2D2.R2D2_data.read_vc`

        models : dict
            Model S based stratification
            See :class:`R2D2.R2D2_data.models_init`
    '''
R2D2_data.__init__           = read.init
R2D2_data.read_qq_select     = read.read_qq_select
R2D2_data.read_qq_select_z   = read.read_qq_select_z
R2D2_data.read_qq            = read.read_qq
R2D2_data.read_qq_restricted = read.read_qq_restricted
R2D2_data.read_qq_all        = read.read_qq_all
R2D2_data.read_qq_tau        = read.read_qq_tau
R2D2_data.read_time          = read.read_time
R2D2_data.read_vc            = read.read_vc
R2D2_data.read_qq_check      = read.read_qq_check
R2D2_data.read_qq_2d         = read.read_qq_2d
R2D2_data.read_qq_slice      = read.read_qq_slice
R2D2_data.read_qq_ixr        = read.read_qq_ixr
R2D2_data.YinYangSet         = read.YinYangSet

R2D2_data.set_cells_gspread  = google.set_cells_gspread
R2D2_data.upgrade_resolution = resolution.upgrade_resolution

R2D2_data.models_init       = models.init
R2D2_data.eos               = models.eos

R2D2_data.sync_remap_qq = sync.sync_remap_qq
R2D2_data.sync_tau = sync.sync_tau
R2D2_data.sync_select = sync.sync_select
R2D2_data.sync_vc = sync.sync_vc
R2D2_data.sync_check = sync.sync_check
R2D2_data.sync_slice = sync.sync_slice
R2D2_data.sync_all = sync.sync_all