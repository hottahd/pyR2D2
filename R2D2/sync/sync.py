import os

class Sync:
    def __init__(self, read):
        """
        Initialize R2D2.Sync
        
        Parameters
        ----------
        read : R2D2.read
            Instance of R2D2.read
        """
        self.read = read
        
    def __getattr__(self, name):
        if hasattr(self.read, name):
            return getattr(self.read, name)

        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")
        
    @staticmethod
    def setup(server,caseid,ssh='ssh',project=os.getcwd().split('/')[-2],dist='../run/'):
        '''
        Downloads setting data from remote server

        Parameters
        ----------
        server : str
            Name of remote server
        caseid : str
            caseid format of 'd001'
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        dist :str
            Destination of data directory
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
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/'+' '+dist+caseid+'/')
            
    def tau(self,server,ssh='ssh',project=os.getcwd().split('/')[-2]):
        '''
        Downloads data at constant optical depth

        Parameters
        ----------
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        '''

        import os

        caseid = self.p['datadir'].split('/')[-3]
        os.system('rsync -avP' \
                +' --exclude="param" ' \
                +' --exclude="qq" ' \
                +' --exclude="remap" ' \
                +' --exclude="slice" ' \
                +' --exclude="time/mhd" ' \
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/ '+self.p['datadir'] )
        
    def remap_qq(self,n,server,ssh='ssh',project=os.getcwd().split('/')[-2]):
        '''
        Downloads full 3D remap data

        Parameters
        ----------
        n : int
            Target time step
        server : str
            Name of remote server
        project : str
            Name of project such as 'R2D2'
        ssh : str
            Type of ssh command
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
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/remap/qq/'+par_dir+chi_dir+'qq.dac.'+str(n).zfill(8)+'.'+ns \
                    +' '+self.p['datadir']+'remap/qq/'+par_dir+chi_dir)
        
    def select(self,xs,server,ssh='ssh',project=os.getcwd().split('/')[-2]):
        '''
        Downloads data at certain height

        Parameters
        ----------
            xs : float
                Height to be downloaded
            server : str
                Name of remote server
            ssh : str
                Type of ssh command
            project : str
                Name of project such as 'R2D2'
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
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/remap/qq/'+par_dir+chi_dir+'qq.dac.*.'+ns \
                    +' '+self.p['datadir']+'remap/qq/'+par_dir+chi_dir)
            
    def vc(self,server,ssh='ssh',project=os.getcwd().split('/')[-2]):
        '''
        Downloads pre analyzed data

        Parameters
        ----------
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        '''

        import os
        caseid = self.p['datadir'].split('/')[-3]

        R2D2_sync.setup(server,caseid,project=project)
        os.system('rsync -avP' \
                +' --exclude="time/mhd" ' \
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/remap/vl '
                +self.p['datadir']+'remap/' )
        
    def check(self,n,server,ssh='ssh',project=os.getcwd().split('/')[-2],end_step=False):
        '''
        Downloads checkpoint data

        Parameters
        ----------
        n : int
            Step to be downloaded
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        end_step : bool
            If true, checkpoint of end step is read
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
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/qq/'+par_dir+chi_dir+'qq.dac.'+step+'.'+str(int(ns)).zfill(8)+' ' \
                +self.p['datadir']+'qq/'+par_dir+chi_dir )

    def slice(self,n,server,ssh='ssh',project=os.getcwd().split('/')[-2]):
        '''
        Downloads slice data

        Parameters
        ----------
        n : int
            Step to be downloaded
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        '''
        import numpy as np
        import os

        step = str(n).zfill(8)
        
        caseid = self.p['datadir'].split('/')[-3]
        os.system('rsync -avP ' \
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/slice/slice.dac ' \
                +self.p['datadir']+'/slice' )
        os.system('rsync -avP ' \
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/data/slice/qq"*".dac.'+step+'."*" ' \
                +self.p['datadir']+'/slice' )
            
    def all(self,server,ssh='ssh',project=os.getcwd().split('/')[-2],dist='../run/'):
        '''
        This method downloads all the data

        Parameters
        ----------
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project :str
            Name of project such as 'R2D2'
        dist :str
            Destination of data directory
        '''
        import os
        
        caseid = self.p['datadir'].split('/')[-3]
        os.system('rsync -avP ' \
                +' -e "'+ssh+'" ' \
                +server+':work/'+project+'/run/'+caseid+'/ ' \
                +dist+caseid)
