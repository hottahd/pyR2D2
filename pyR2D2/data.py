import pyR2D2

class Data:
    '''
    Class for managing R2D2 data
    
    Attributes
    ----------
    read : R2D2.Read
        Instance of R2D2.Read
    sync : R2D2.Sync
        Instance of R2D2.Sync

    models : dict
        Model S based stratification.
        See :class:`R2D2.R2D2_read.models_init`
            
    '''
    
    def __init__(self,datadir,verbose=False,self_old=None):
        '''
        Initialize R2D2.Data
        '''
        self.read = pyR2D2.Read(datadir,verbose=verbose,self_old=self_old)
        self.sync = pyR2D2.Sync(self.read)
        
    def summary(self):
        """
        Rapper of :meth:`R2D2.Read.summary`
        
        """
        self.read.summary()

    def __getattr__(self, name):
        '''
        When an attribute is not found in R2D2.Data, it is searched in R2D2.Data.Read
        '''
        if hasattr(self.read, name):
            attr = getattr(self.read, name)
            if not callable(attr):  # 属性が関数でない場合のみ返す
                return attr

        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

# R2D2_data.read_time = read.read_time

# R2D2_data.set_cells_gspread  = google.set_cells_gspread
# R2D2_data.upgrade_resolution = resolution.upgrade_resolution
# R2D2_data.models_init       = models.init
# R2D2_data.eos               = models.eos