from .r2d2_read import R2D2_read
from .r2d2_sync import R2D2_sync
# from . import google
# from . import resolution
# from . import models

class R2D2_data:    
    '''
    Class for managing R2D2_data
    
    Notes
    -----
    self.read is an instance of R2D2.R2D2_read

    Attributes
    ----------
        read : R2D2.R2D2_read
            Instance of R2D2.R2D2_read
        sync : R2D2.R2D2_sync
            Instance of R2D2.R2D2_sync
    
        models : dict
            Model S based stratification.
            See :class:`R2D2.R2D2_read.models_init`
    '''

    def __init__(self,datadir,verbose=False,self_old=None):
        '''
        Initialize R2D2_data
        '''
        self.read = R2D2_read(datadir,verbose=verbose,self_old=self_old)
        self.sync = R2D2_sync(self.read)
        
    def summary(self):
        """
        Rapper of :meth:`R2D2.R2D2_read.summary`
        
        """
        R2D2_read.summary(self.read)

    def __getattr__(self, name):
        '''
        When an attribute is not found in R2D2_data, it is searched in R2D2_data.read        
        '''
        if hasattr(self.read, name):
            attr = getattr(self.read, name)
            return attr

        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

# R2D2_data.read_time = read.read_time

# R2D2_data.set_cells_gspread  = google.set_cells_gspread
# R2D2_data.upgrade_resolution = resolution.upgrade_resolution
# R2D2_data.models_init       = models.init
# R2D2_data.eos               = models.eos