import pyR2D2

class Data:
    '''
    Class for managing R2D2 data
    
    Attributes
    ----------
    read : pyR2D2.Read
        Instance of pyR2D2.Read
    sync : pyR2D2.Sync
        Instance of pyR2D2.Sync
            
    '''
    
    def __init__(self,datadir,verbose=False,self_old=None):
        '''
        Initialize pyR2D2.Data
        '''
        self.read = pyR2D2.Read(datadir,verbose=verbose,self_old=self_old)
        self.sync = pyR2D2.Sync(self.read)
        
    def summary(self):
        """
        Rapper of :meth:`pyR2D2.Read.summary`
        
        """
        self.read.summary()

    def __getattr__(self, name):
        '''
        When an attribute is not found in pyR2D2.Data, it is searched in pyR2D2.Data.Read
        '''
        if hasattr(self.read, name):
            attr = getattr(self.read, name)
            if not callable(attr):  # 属性が関数でない場合のみ返す
                return attr

        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")