def init(self,):
    '''
    This method read Model S based stratification.
    '''
    import numpy as np

    self.models = {}
    
    f = open(self.p['datadir']+'../input_data/params.txt','r')
    self.models['ix'] = int(f.read())
    f.close()

    ix = self.models['ix']

    endian = '>'
    dtype = np.dtype
    dtyp=np.dtype([ \
                    ("head",endian+"i"),\
                    ("x",endian+str(ix)+"d"),\
                    ("pr0",endian+str(ix)+"d"),\
                    ("ro0",endian+str(ix)+"d"),\
                    ("te0",endian+str(ix)+"d"),\
                    ("se0",endian+str(ix)+"d"),\
                    ("en0",endian+str(ix)+"d"),\
                    ("op0",endian+str(ix)+"d"),\
                    ("tu0",endian+str(ix)+"d"),\
                    ("dsedr0",endian+str(ix)+"d"),\
                    ("dtedr0",endian+str(ix)+"d"),\
                    ("dprdro",endian+str(ix)+"d"),\
                    ("dprdse",endian+str(ix)+"d"),\
                    ("dtedro",endian+str(ix)+"d"),\
                    ("dtedse",endian+str(ix)+"d"),\
                    ("dendro",endian+str(ix)+"d"),\
                    ("dendse",endian+str(ix)+"d"),\
                    ("gx",endian+str(ix)+"d"),\
                    ("kp",endian+str(ix)+"d"),\
                    ("cp",endian+str(ix)+"d"),\
                    ("tail",endian+"i")\
    ])
    f = open(self.p['datadir']+"../input_data/value_cart.dac",'rb')
    back = np.fromfile(f,dtype=dtyp,count=1)
    f.close()
    
    for key in back.dtype.names:
        if back[key].size == ix:
            self.models[key] = back[key].reshape((ix),order='F')

def eos(self,ro,se,var):
    '''
    This method returns the Table equation of state.
    
    Parameters:
        ro (float): density
        se (float): entropy
        var (str): variable name, pr, te, en, op
    
    Returns:
       qq (float): corresponding variable
       
    Warning:
       This method is very slow for large numpy array. Please consider using fortraun code.
    '''
    
    import numpy as np
    
    iro = int((np.log(ro) - self.p['log_ro_e'][0])//self.p['dlogro_e'])
    ise = int((       se  - self.p['se_e'][0]    )//self.p['dse_e'])
    
    dlogro = np.log(ro) - self.p['log_ro_e'][iro]
    dse = se - self.p['se_e'][ise]
    
    qq = np.exp((\
        + self.p['log_'+var+'_e'][iro  ,ise  ]*(self.p['dlogro_e'] - dlogro)*(self.p['dse_e'] - dse) \
        + self.p['log_'+var+'_e'][iro+1,ise  ]*(                     dlogro)*(self.p['dse_e'] - dse) \
        + self.p['log_'+var+'_e'][iro  ,ise+1]*(self.p['dlogro_e'] - dlogro)*(                  dse) \
        + self.p['log_'+var+'_e'][iro+1,ise+1]*(                     dlogro)*(                  dse) \
        )/self.p['dlogro_e']/self.p['dse_e'])    
    return qq