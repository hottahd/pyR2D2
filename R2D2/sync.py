import os

def set(server,caseid,project=os.getcwd().split('/')[-2],dist='../run/'):
    '''
    This method downloads setting data from remote server

    Parameters:
        server (str): name of remote server
        caseid (str): caseid format of 'd001'
        project (str): name of project such as 'R2D2'
        dist (str): distination of data directory
    '''
    import os
        
    os.system('rsync -avP' \
              +' --exclude="a.out" ' \
              +' --exclude=".git" ' \
              +' --exclude="make/*.o" ' \
              +' --exclude="make/*.lst" ' \
              +' --exclude="make/*.mod" ' \
              +' --exclude="data/qq" ' \
              +' --exclude="data/remap/qq" ' \
              +' --exclude="data/remap/vl/vl*" ' \
              +' --exclude="data/slice/qq*" ' \
              +' --exclude="data/tau/qq*" ' \
              +' --exclude="output.*" ' \
              +server+':work/'+project+'/run/'+caseid+'/'+' '+dist+caseid+'/')
    
def sync_tau(self,server,project=os.getcwd().split('/')[-2]):
    '''
    This method downloads data at constant optical depth

    Parameters:
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
    '''

    import os

    caseid = self.p['datadir'].split('/')[-3]
    os.system('rsync -avP' \
              +' --exclude="param" ' \
              +' --exclude="qq" ' \
              +' --exclude="remap" ' \
              +' --exclude="slice" ' \
              +' --exclude="time/mhd" ' \
              +server+':work/'+project+'/run/'+caseid+'/data/ '+self.p['datadir'] )
    
def sync_remap_qq(self,n,server,project=os.getcwd().split('/')[-2]):
    '''
    This method downloads full 3D remap data

    Parameters:
        n (int): time step
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
    '''
    import os
    import numpy as np
    
    caseid = self.p['datadir'].split('/')[-3]
    
    # remapを行ったMPIランクの洗い出し
    nps = np.char.zfill(self.p['np_ijr'].flatten().astype(str),8)
    for ns in nps:
        par_dir = str(int(ns)//1000).zfill(5)+'/'
        chi_dir = str(int(ns)).zfill(8)+'/'
        
        os.makedirs(self.p['datadir']+'remap/qq/'+par_dir+chi_dir,exist_ok=True)
        os.system('rsync -avP ' \
            +server+':work/'+project+'/run/'+caseid+'/data/remap/qq/'+par_dir+chi_dir+'qq.dac.'+str(n).zfill(8)+'.'+ns \
                +' '+self.p['datadir']+'remap/qq/'+par_dir+chi_dir)
    
def sync_select(self,xs,server,project=os.getcwd().split('/')[-2]):
    '''
    This method downloads data at certain height

    Parameters:
        xs (float): height to be downloaded
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
    '''

    import os
    import numpy as np

    i0 = np.argmin(np.abs(self.p["x"]-xs))
    ir0 = self.p["i2ir"][i0]
    
    nps = np.char.zfill(self.p['np_ijr'][ir0-1,:].astype(str),8)

    files = ''
    caseid = self.p['datadir'].split('/')[-3]

    for ns in nps:
        par_dir = str(int(ns)//1000).zfill(5)+'/'
        chi_dir = str(int(ns)).zfill(8)+'/'
        
        os.makedirs(self.p['datadir']+'remap/qq/'+par_dir+chi_dir,exist_ok=True)        
        os.system('rsync -avP ' \
            +server+':work/'+project+'/run/'+caseid+'/data/remap/qq/'+par_dir+chi_dir+'qq.dac.*.'+ns \
                +' '+self.p['datadir']+'remap/qq/'+par_dir+chi_dir)
        
def sync_vc(self,server,project=os.getcwd().split('/')[-2]):
    '''
    This method downloads pre analyzed data

    Parameters:
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
    '''

    import os
    caseid = self.p['datadir'].split('/')[-3]

    set(server,caseid,project=project)
    os.system('rsync -avP' \
              +' --exclude="time/mhd" ' \
              +server+':work/'+project+'/run/'+caseid+'/data/remap/vl '
              +self.p['datadir']+'remap/' )
    
def sync_check(self,n,server,project=os.getcwd().split('/')[-2],end_step=False):
    '''
    This method downloads checkpoint data

    Parameters:
        n (int): step to be downloaded
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
        end_step (bool): If true, checkpoint of end step is read
    '''
    import numpy as np
    import os
    
    step = str(n).zfill(8)
    
    if end_step:
        if np.mod(self.p['nd'],2) == 0:
            step = 'e'
        if np.mod(self.p['nd'],2) == 1:
            step = 'o'
    
    caseid = self.p['datadir'].split('/')[-3]
    for ns in range(self.p['npe']):
        par_dir = str(int(ns)//1000).zfill(5)+'/'
        chi_dir = str(int(ns)).zfill(8)+'/'

        os.makedirs(self.p['datadir']+'qq/'+par_dir+chi_dir,exist_ok=True)
        os.system('rsync -avP ' \
              +server+':work/'+project+'/run/'+caseid+'/data/qq/'+par_dir+chi_dir+'qq.dac.'+step+'.'+str(int(ns)).zfill(8)+' ' \
              +self.p['datadir']+'qq/'+par_dir+chi_dir )

def sync_slice(self,n,server,project=os.getcwd().split('/')[-2]):
    '''
    This method downloads slice data

    Parameters:
        n (int): step to be downloaded
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
    '''
    import numpy as np
    import os

    step = str(n).zfill(8)
    
    caseid = self.p['datadir'].split('/')[-3]
    os.system('rsync -avP ' \
              +server+':work/'+project+'/run/'+caseid+'/data/slice/slice.dac ' \
              +self.p['datadir']+'/slice' )
    os.system('rsync -avP ' \
              +server+':work/'+project+'/run/'+caseid+'/data/slice/qq"*".dac.'+step+'."*" ' \
              +self.p['datadir']+'/slice' )
        
def sync_all(self,server,project=os.getcwd().split('/')[-2],dist='../run/'):
    '''
    This method downloads all the data

    Parameters:
        server (str): name of remote server
        project (str): name of project such as 'R2D2'
        dist (str): distination of data directory
    '''
    import os
    
    caseid = self.p['datadir'].split('/')[-3]
    os.system('rsync -avP ' \
              +server+':work/'+project+'/run/'+caseid+'/ ' \
              +dist+caseid)
