import os
import numpy as np
import pyR2D2

class Parameters:
    '''
    '''
    def __init__(self, datadir: str):
        from scipy.io import FortranFile
        
        self.datadir = datadir
        
        # Time control parameters
        with open(self.datadir+"param/nd.dac","r") as f:
            nn = f.read().split()
            self.nd = int(nn[0]) # current maximum time step
            self.nd_tau = int(nn[1]) # current maximum time step for tau
                    
        # data/time/tau のファイルの数を数え、nd_tauと比較して大きい方をnd_tauとする
        if os.path.isdir(self.datadir+'time/tau'):
            self.nd_tau = max(len(os.listdir(datadir+'time/tau')) - 1,self.nd_tau)
            
        # Read basic parameters
        with open(self.datadir+"param/params.dac","r") as f:
            for line in f:
                parts = line.split()
                key = parts[1]
                value = parts[0]
                value_type = parts[2]
            
                if value_type == 'i': # integer
                    self.__dict__[key] = int(value)
                elif value_type == 'd': # double
                    self.__dict__[key] = float(value)
                elif value_type == 'c': # character
                    self.__dict__[key] = value
                elif value_type == 'l':
                    if value == 'F':
                        self.__dict__[key] = False
                    else:
                        self.__dict__[key] = True
        
        # self.rstarのない場合はself.rstar = pyR2D2.constant.RSUN
        if 'rstar' not in self.__dict__:
            self.rstar = pyR2D2.constant.RSUN

        # remapで出力している変数の種類
        self.remap_kind = ['ro','vx','vy','vz','bx','by','bz','se','ph']
        self.remap_kind_add = ['pr','te','op']

        # transform endiant for python
        if self.swap == 0: # little
            self.endian = "<"
        else:
            self.endian = ">"
            
        # MPI information
        f = FortranFile(self.datadir+'param/xyz.dac','r')
        shape = (self.ix0*self.jx0*self.kx0,3)
        self.xyz = f.read_reals(dtype=np.int32).reshape(shape,order='F')
        
        ## number of grid in each direction
        self.ix = self.ix0*self.nx
        self.jx = self.jx0*self.ny
        self.kx = self.kx0*self.nz
    
        self.ixg = self.ix + 2*self.margin*(self.xdcheck-1)
        self.jxg = self.jx + 2*self.margin*(self.ydcheck-1)
        self.kxg = self.kx + 2*self.margin*(self.zdcheck-1)
            
        # read background variables
        dtype=np.dtype([ \
                        ("head",self.endian+"i",),\
                        ("x",self.endian+"d",(self.ixg,)),\
                        ("y",self.endian+"d",(self.jxg,)),\
                        ("z",self.endian+"d",(self.kxg,)),\
                        ("pr0",self.endian+"d",(self.ixg,)),\
                        ("te0",self.endian+"d",(self.ixg,)),\
                        ("ro0",self.endian+"d",(self.ixg,)),\
                        ("se0",self.endian+"d",(self.ixg,)),\
                        ("en0",self.endian+"d",(self.ixg,)),\
                        ("op0",self.endian+"d",(self.ixg,)),\
                        ("tu0",self.endian+"d",(self.ixg,)),\
                        ("dsedr0",self.endian+"d",(self.ixg,)),\
                        ("dtedr0",self.endian+"d",(self.ixg,)),\
                        ("dprdro",self.endian+"d",(self.ixg,)),\
                        ("dprdse",self.endian+"d",(self.ixg,)),\
                        ("dtedro",self.endian+"d",(self.ixg,)),\
                        ("dtedse",self.endian+"d",(self.ixg,)),\
                        ("dendro",self.endian+"d",(self.ixg,)),\
                        ("dendse",self.endian+"d",(self.ixg,)),\
                        ("gx",self.endian+"d",(self.ixg,)),\
                        ("cp",self.endian+"d",(self.ixg,)),\
                        ("fa",self.endian+"d",(self.ixg,)),\
                        ("sa",self.endian+"d",(self.ixg,)),\
                        ("xi",self.endian+"d",(self.ixg,)),\
                        ("tail",self.endian+"i")\
        ])
        with open(self.datadir+"param/back.dac",'rb') as f:
            back = np.fromfile(f,dtype=dtype,count=1)

        # marginも含んだ座標
        self.xg = back['x'].reshape((self.ixg),order='F')
        self.yg = back['y'].reshape((self.jxg),order='F')
        self.zg = back['z'].reshape((self.kxg),order='F')

        # fluxはi+1/2で定義するので、そのための座標を定義
        self.x_flux = np.zeros(self.ix + 1)
        for i in range(0,self.ix + 1):
            self.x_flux[i] = 0.5*(self.xg[i+self.margin] + self.xg[i+self.margin-1])
        
        self.ro0g = back['ro0'].reshape((self.ixg),order='F')
        self.se0g = back['se0'].reshape((self.ixg),order='F')
        
        for key in back.dtype.names:
            if back[key].size == self.ixg:
                self.__dict__[key] = back[key].reshape((self.ixg),order="F")[self.margin:self.ixg-self.margin]
            elif back[key].size == self.jxg:
                self.__dict__[key] = back[key].reshape((self.jxg),order="F")[self.margin:self.jxg-self.margin]
            elif back[key].size == self.kxg:
                self.__dict__[key] = back[key].reshape((self.kxg),order="F")[self.margin:self.kxg-self.margin]
                
        self.xr = self.x/self.rstar

        if self.geometry == 'YinYang':            
            # original geometry in Yin-Yang grid
            self.jx_yy = self.jx
            self.kx_yy = self.kx

            self.jxg_yy = self.jx + 2*self.margin
            self.kxg_yy = self.kx + 2*self.margin
            
            self.y_yy = self.y - 0.5*(self.y[1] - self.y[0])
            self.z_yy = self.z - 0.5*(self.z[1] - self.z[0])

            self.yg_yy = self.yg
            self.zg_yy = self.zg
            
            # converted spherical geometry
            self.jx = self.jx*2
            self.kx = self.jx*2

            dy =   np.pi/self.jx
            dz = 2*np.pi/self.kx

            y = np.zeros(self.jx)
            z = np.zeros(self.kx)
            
            y[0] = 0.5*dy
            z[0] = 0.5*dz - np.pi

            for j in range(1,self.jx):
                y[j] = y[j-1] + dy

            for k in range(1,self.kx):
                z[k] = z[k-1] + dz

            self.y = y
            self.z = z
                
        if self.zdcheck == 2:
            dimension = '3d'
        else:
            dimension = '2d'
            
        ##############################
        # read value information
        if dimension == "3d":
            with open(self.datadir+"remap/vl/c.dac","r") as f:
                value = f.read().split()
                self.m2d_xy = int(value[0])
                self.m2d_xz = int(value[1])
                self.m2d_flux = int(value[2])
                if self.geometry == 'YinYang':
                    self.m2d_spex = int(value[3])
                del value[0:3]
                if self.geometry == 'YinYang':
                    del value[0]
                self.cl = list(map(str.strip,value)) ## strip space from character

            ##############################
            # read mpi information for remap
            dtyp=np.dtype([ \
                        ("iss", self.endian+str(self.npe)+"i4"),\
                        ("iee", self.endian+str(self.npe)+"i4"),\
                        ("jss", self.endian+str(self.npe)+"i4"),\
                        ("jee", self.endian+str(self.npe)+"i4"),\
                        ("iixl", self.endian+str(self.npe)+"i4"),\
                        ("jjxl", self.endian+str(self.npe)+"i4"),\
                        ("np_ijr", self.endian+str(self.ixr*self.jxr)+"i4"),\
                        ("ir", self.endian+str(self.npe)+"i4"),\
                        ("jr", self.endian+str(self.npe)+"i4"),\
                        ("i2ir", self.endian+str(self.ixg)+"i4"),\
                        ("j2jr", self.endian+str(self.jxg)+"i4"),\
            ])
            
            with open(datadir+"remap/remap_info.dac",'rb') as f:
                mpi = np.fromfile(f,dtype=dtyp,count=1)    
            
            for key in mpi.dtype.names:
                if key == "np_ijr":
                    self.__dict__[key] = mpi[key].reshape((self.ixr,self.jxr),order="F")
                else:
                    self.__dict__[key] = mpi[key].reshape((mpi[key].size),order="F")
        
            self.i2ir = self.i2ir[self.margin:self.ixg-self.margin]
            self.j2jr = self.j2jr[self.margin:self.jxg-self.margin]
            
            self.iss = self.iss - 1
            self.iee = self.iee - 1
            self.jss = self.jss - 1
            self.jee = self.jee - 1

            if os.path.isdir(self.datadir+'slice'):
                with open(self.datadir+"slice/params.dac","r") as f:
                    line = f.readline()
                    while line:
                        self.__dict__[line.split()[1]] = int(line.split()[0])
                        line = f.readline()

                dtype = np.dtype([ \
                    ('x_slice',self.endian + str(self.nx_slice)+'d'),\
                    ('y_slice',self.endian + str(self.ny_slice)+'d'),\
                    ('z_slice',self.endian + str(self.nz_slice)+'d'),\
                ])
                with open(self.datadir + 'slice/slice.dac','rb') as f:
                    slice = np.fromfile(f,dtype=dtype)

                    self.x_slice = slice['x_slice'].reshape(self.nx_slice,order='F')
                    self.y_slice = slice['y_slice'].reshape(self.ny_slice,order='F')
                    self.z_slice = slice['z_slice'].reshape(self.nz_slice,order='F')
            
        # read equation of state
        eosdir = self.datadir[:-5]+'input_data/'
        if os.path.exists(eosdir+'eos_table_sero.npz'):
            eos_d = np.load(eosdir+'eos_table_sero.npz')
            self.log_ro_e = eos_d['ro'] # density is defined in logarithmic scale
            self.se_e = eos_d['se']
            self.ix_e = len(self.log_ro_e)
            self.jx_e = len(self.se_e)
            
            self.log_pr_e = np.log(eos_d['pr']+ 1.e-200)
            self.log_en_e = np.log(eos_d['en']+ 1.e-200)
            self.log_te_e = np.log(eos_d['te']+ 1.e-200)
            self.log_op_e = np.log(eos_d['op']+ 1.e-200)
            
            self.dlogro_e = self.log_ro_e[1] - self.log_ro_e[0]
            self.dse_e = self.se_e[1] - self.se_e[0]
                    
        # read original data
        if os.path.exists(datadir+'cont_log.txt'):
            with open(datadir+'cont_log.txt') as f:
                self.origin = f.readlines()[6][-11:-7]
        else:
            self.origin = 'N/A'
            
class Read:
    '''
    Class for reading pyR2D2 data
    
    Attributes
    ----------
    p : pyR2D2.Prameters
        basic parameters.
        See :meth:`pyR2D2.Parameters`
    qs : dict
        2D data at selected height.
        See :meth:`pyR2D2.Read.qq_select`
    qq : dict
        3D full data.
        See :meth:`pyR2D2.Read.qq_3d`        
    qt : dict
        2D data at constant optical depths.
        See :meth:`pyR2D2.Read.qq_tau`
    qt_yin, qt_yan : dict
        2D data at constant optical depths in Yin-Yang grid.
        See :meth:`R2D2.Read.qq_tau`        
    qz : dict    
        2D data at constant z.
        See :meth:`R2D2.Read.qq_select_z`
    qr : dict
        3D data in a restricted region.
        See :meth:`R2D2.Read.qq_3d_restricted`
    qi : dict
        3D data in a MPI x-region.
        See :meth:`R2D2.Read.qq_ixr`
    ql : dict
        2D slice data.
        See :meth:`R2D2.Read.qq_slice`
    ql_yin, ql_yan : dict
        2D slice data in a Yin-Yang grid.
        See :meth:`pyR2D2.Read.qq_slice`
    q2 : dict
        Full data in 2D simulation.
        See :meth:`pyR2D2.Read.qq_2d`
    t : float
        time. See :meth:`pyR2D2.Read.time`
    vc : dict
        data of on the fly analysis.
        See :meth:`pyR2D2.Read.on_the_fly`
    
    '''
    
    def __init__(self, datadir: str, verbose: bool = False, self_old= None):
        '''
        This method reads basic data for calculation setting
        The data is stored in self.p dictionary

        Parameters
        ----------
        datadir : str
            directory of data
        verbose : bool
            True shows the information of data
        self_old : pyR2D2.Data
            if self_old is not None, datadir is compared with old one and if datadir is same as old one, self is updated with self_old
        '''
        from scipy.io import FortranFile
        
        # check if the data is already read2
        initialize_flag = True
        if self_old is not None:
            if self_old.datadir == datadir:
                initialize_flag = False
            
        if initialize_flag:  
            #self.p = {}
            self.p = Parameters(datadir)
            self.qs = {}
            self.qq = {}
            self.qt = {}
            self.qt_yin = {}
            self.qt_yan = {}            
            self.qz = {}
            self.qr = {}
            self.qi = {}
            self.ql = {}
            self.ql_yin = {}
            self.ql_yan = {}
            self.q2 = {}
            self.t = {}
            self.vc = {}
        else:
            self.__dict__.update(vars(self_old.read))
                                
        if verbose:
            self.summary()

    def __getattr__(self, name):
        '''
        When an attribute is not found in pyR2D2.Read, it is searched in pyR2D2.Read.p
        '''
        if hasattr(self.p, name):
            attr = getattr(self.p, name)
            if not callable(attr):  # 属性が関数でない場合のみ返す
                return attr

    def _get_filepath_remap_qq(self,n,np0):
        '''
        get file path of remap/qq/
        
        Parameters
        ----------
        n : int
            time step
        np0 : int
            MPI process number
            
        Returns
        -------
        filepath : str
            file path of remap/qq/
        '''            
        if os.path.isdir(self.datadir + "remap/qq/00000/"):
            cnou = '{0:05d}'.format(np0//1000)
            cno  = '{0:08d}'.format(np0)
            filepath = self.datadir+"remap/qq/"+cnou+"/"+cno+"/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0)
        else:
            # directoryを分けない古いバージョン対応
            filepath = self.datadir+"remap/qq/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0)
        return filepath            

    def _allocate_remap_qq(self,qq_type: str, ijk: tuple):
        """
        Paramters
        ---------
        qq_type : str
            type of data. 

        ijk : tuple
            size of data
        """
        memflag = True
        if 'ro' in self.__dict__[qq_type]:
            memflag = not self.__dict__[qq_type]['ro'].shape == ijk
        if 'ro' not in self.__dict__[qq_type] or memflag:
            for key in self.remap_kind + self.remap_kind_add:
                self.__dict__[qq_type][key] = np.zeros(ijk)

                
    def _dtype_remap_qq(self,np0):
        dtype=np.dtype([ \
                ("qq",self.endian + str(self.mtype*self.iixl[np0]*self.jjxl[np0]*self.kx)+"f"),\
                ("pr",self.endian + str(self.iixl[np0]*self.jjxl[np0]*self.kx)+"f"),\
                ("te",self.endian + str(self.iixl[np0]*self.jjxl[np0]*self.kx)+"f"),\
                ("op",self.endian + str(self.iixl[np0]*self.jjxl[np0]*self.kx)+"f"),\
        ])        

        return dtype

    # -- 
    def qq_select(self,xs: float,n: int, verbose: bool = True):
        '''
        Reads 2D slice data at a selected height.

        Parameters
        ----------
        xs : float 
            A selected height for data
        n : int 
            A setected time step for data
        verbose : bool
            False suppresses a message of store
            
        Notes
        -----
        The data is stored in self.qs dictionary
        '''

        i0 = np.argmin(np.abs(self.x - xs))
        ir0 = self.i2ir[i0]
                
        self._allocate_remap_qq(qq_type = 'qs', ijk = [self.jx, self.kx])
                
        for jr0 in range(1,self.jxr+1):
            np0 = self.np_ijr[ir0-1,jr0-1]
    
            if jr0 == self.jr[np0]:
                dtype = self._dtype_remap_qq(np0)
                filepath = self._get_filepath_remap_qq(n,np0)
                
                with open(filepath,'rb') as f:
                    qqq = np.fromfile(f,dtype=dtype,count=1)
                    for key, m in zip( self.remap_kind, range(self.mtype) ):
                        self.qs[key] = \
                                qqq["qq"].reshape((self.iixl[np0], self.jjxl[np0], self.kx, self.mtype), order="F")[i0-self.iss[np0],:,:,m]
        
                    for key in self.remap_kind_add:
                        self.qs[key] = \
                                qqq[key].reshape((self.iixl[np0], self.jjxl[np0], self.kx),order="F")[i0-self.iss[np0],:,:]

                info = {}
                info['xs'] = self.x[i0]
                self.qs['info'] = info

        if verbose:
            print('### variables are stored in self.qs ###')

    ##############################
    def qq_select_z(self,zs,n,verbose=True):
        '''
        Reads 2D slice data at constant z
        The data is stored in self.qz dictionary

        Parameters
        ----------
        zs : float
            A selected z for data
        n : int
            A selected time step for data
        verbose : bool
            False suppresses a message of store
            
        '''

        k0 = np.argmin(np.abs(self.z-zs))
        
        self._allocate_remap_qq(qq_type = 'qz', ijk = [self.ix, self.jx])
        
        for ir0 in range(1,self.ixr + 1):
            for jr0 in range(1,self.jxr + 1):
                np0 = self.np_ijr[ir0-1,jr0-1]

                dtype = self._dtype_remap_qq(np0)
                filepath = self._get_filepath_remap_qq(n,np0)
                
                with open(filepath,'rb') as f:
                    qqq = np.fromfile(f,dtype=dtype,count=1)
                    for key, m in zip( self.remap_kind, range(self.mtype) ):
                        self.qz[key][self.iss[np0]:self.iee[np0]+1,self.jss[np0]:self.jee[np0]+1] \
                                = qqq["qq"].reshape((self.iixl[np0], self.jjxl[np0], self.kx, self.mtype), order="F")[:,:,k0,m]

                    for key in self.remap_kind_add:
                        self.qz[key][self.iss[np0]:self.iee[np0]+1,self.jss[np0]:self.jee[np0]+1] \
                                = qqq[key].reshape((self.iixl[np0], self.jjxl[np0], self.kx), order="F")[:,:,k0]            
                    
            info = {}
            info['zs'] = self.z[k0]
            self.qz['info'] = info
            
        if verbose:
            print('### variables are stored in self.qz ###')
            
    ##############################
    def qq_3d(self,n,value,verbose=True):
        '''
        Reads 3D full data
        The data is stored in self.qq dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        value : str
            Kind of value. Options are:
                - "ro" : density
                - "vx", "vy", "vz": velocity
                - "bx", "by", "bz": magnetic field
                - "se": entropy
                - "ph": div B cleaning
                - "te": temperature
                - "op": Opacity
        verbose : bool
            False suppresses a message of store
        '''
        
        if type(value) == str:
            values_input = [value]
        if type(value) == list:
            values_input = value

        for value in values_input:
            if value not in self.p.remap_kind+self.p.remap_kind_add:
                print('######')
                print('value =',value)
                print('value should be one of ',self.p.remap_kind+self.p.remap_kind_add)
                print('return')
                return
        
        
        self._allocate_remap_qq(qq_type = 'qq', ijk = [self.ix, self.jx, self.kx])

        for ir0 in range(1,self.ixr + 1):
            for jr0 in range(1,self.jxr + 1):
                np0 = self.np_ijr[ir0-1,jr0-1]
                if ir0 == self.ir[np0] and jr0 == self.jr[np0]:

                    dtype = self._dtype_remap_qq(np0)                    
                    filepath = self._get_filepath_remap_qq(n,np0)
                    
                    with open(filepath,'rb') as f:
                        qqq = np.fromfile(f,dtype=dtype,count=1)

                        for value in values_input:
                            if value in self.p.remap_kind:
                                m = self.p.remap_kind.index(value)
                                self.qq[value][self.iss[np0]:self.iee[np0]+1,self.jss[np0]:self.jee[np0]+1,:] \
                                    = qqq["qq"].reshape((self.iixl[np0],self.jjxl[np0],self.kx,self.mtype),order="F")[:,:,:,m]
                            else:
                                self.qq[value][self.iss[np0]:self.iee[np0]+1,self.jss[np0]:self.jee[np0]+1,:] \
                                    = qqq[value].reshape((self.iixl[np0],self.jjxl[np0],self.kx),order="F")[:,:,:]

        if verbose:
            print('### variables are stored in self.qq ###')
            
    ##############################
    def qq_3d_restricted(self,n,value,x0,x1,y0,y1,z0,z1,verbose=True):
        '''
        Reads 3D restricted-area data
        The data is stored in self.qr dictionary

        Parameters
        ----------
        n : int 
            A selected time step for data
        value : str
            See :meth:`R2D2.Read.qq_3d` for options
    
        x0, y0, z0 : float
            Minimum x, y, z
        x1, y1, z1 : float
            Maximum x, y, z
        verbose : bool
            False suppresses a message of store
        '''
                
        i0, i1 = np.argmin(abs(self.x-x0)), np.argmin(abs(self.x-x1))
        j0, j1 = np.argmin(abs(self.y-y0)), np.argmin(abs(self.y-y1))
        k0, k1 = np.argmin(abs(self.z-z0)), np.argmin(abs(self.z-z1))
        ixr = i1 - i0 + 1
        jxr = j1 - j0 + 1
        kxr = k1 - k0 + 1
        

        if type(value) == str:
            values_input = [value]
        if type(value) == list:
            values_input = value
        
        for value in values_input:
            self.qr[value] = np.zeros((ixr,jxr,kxr),dtype=np.float32)
            
        for ir0 in range(1,self.ixr+1):
            for jr0 in range(1,self.jxr+1):
                np0 = self.np_ijr[ir0-1,jr0-1]
                
                if(not (self.iss[np0] > i1 or self.iee[np0] < i0 or self.jss[np0] > j1 or self.jee[np0] < j0) ):
                    
                    dtype = self._dtype_remap_qq(np0)
                    filepath = self._get_filepath_remap_qq(n,np0)

                    with open(filepath,'rb') as f:
                        qqq = np.fromfile(f,dtype=dtyp,count=1)

                        for value in values_input:
                            if value in self.p.remap_kind:
                                m = self.p.remap_kind.index(value)
                                isrt_rcv = max([0  ,self.iss[np0]-i0  ])
                                iend_rcv = min([ixr,self.iee[np0]-i0+1])
                                jsrt_rcv = max([0  ,self.jss[np0]-j0  ])
                                jend_rcv = min([jxr,self.jee[np0]-j0+1])
                                
                                isrt_snd = isrt_rcv - (self.iss[np0]-i0)
                                iend_snd = isrt_snd + (iend_rcv - isrt_rcv)
                                jsrt_snd = jsrt_rcv - (self.jss[np0]-j0)
                                jend_snd = jsrt_snd + (jend_rcv - jsrt_rcv)

                                self.qr[value][isrt_rcv:iend_rcv,jsrt_rcv:jend_rcv,:] = \
                                    qqq["qq"].reshape((self.iixl[np0],self.jjxl[np0],self.kx, self.mtype),order="F")[isrt_snd:iend_snd,jsrt_snd:jend_snd,k0:k1+1,m]
                            else:
                                self.qr[value][isrt_rcv:iend_rcv,jsrt_rcv:jend_rcv,:] = \
                                    qqq[value].reshape((self.iixl[np0],self.jjxl[np0],self.kx),order="F")[isrt_snd:iend_snd,jsrt_snd:jend_snd,k0:k1+1]                         
                    f.close()

        if verbose:
            print('### variales are stored in self.qr ###')

    ##############################
    def qq_3d_all(self,n,verbose=True):
        '''
        This method reads 3D full data all
        The data is stored in self.qq dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        verbose : bool
            False suppresses a message of store
        '''
        self.qq_3d(self,n,self.p.remap_kind+self.p.remap_kind_add,verbose=verbose)

    ##############################
    def qq_tau(self,n,verbose=True):
        '''
        Reads 2D data at certain optical depths.
        The data is stored in self.qt dictionary.
        In this version the selected optical depth is 1, 0.1, and 0.01

        Parameters
        ----------
        n : int
            A selected time step for data
        verbose : bool
            False suppresses a message of store
                
        '''

        if self.geometry == 'YinYang':
            files = ['_yin','_yan']
            qts = [self.qt_yin,self.qt_yan]
        else:
            files=['']
            qts = [self.qt]

        if(self.geometry == 'YinYang'):
            jx = self.jx_yy
            kx = self.kx_yy
        else:
            jx = self.jx
            kx = self.kx

        value_keys = ['in','ro','se','pr','te','vx','vy','vz','bx','by','bz','he','fr']
        for file,qt in zip(files,qts):
            with open(self.datadir+"tau/qq"+file+".dac."+'{0:08d}'.format(n),"rb") as f:
                qq_in0 = np.fromfile(f,self.endian + 'f',self.m_tu*self.m_in*jx*kx)
                
            for key, mk in zip(value_keys, range(13)):
                for tau, mt in zip(['','01','001'],range(3)):
                    qt[key+tau] = np.reshape(qq_in0,(self.m_tu,self.m_in,jx,kx),order="F")[mt,mk,:,:]
        
        if verbose:
            if self.geometry == 'YinYang':
                print('### variables are stored in self.qt_yin and self.qt_yan ###')
            else:
                print('### variables are stored in self.qt ###')
            
    ##############################
    def qq_ixr(self,ixrt,n,verbose=True):
        '''
        Reads 3D data at a selected a MPI process in x-direction.
        corrensponding to ixr-th region in remap coordinate
        
        Parameters
        ----------
        ixrt : int
            A selected x-region (0<=ixrt<=ixr-1)
            Note
            ----
            ixr is the number of MPI process in x-direction
            
        n : int
            A selected time step for data
        verbose : bool
            False suppresses a message of store
        '''
                
        # corresponding MPI process
        nps = np.where(self.ir - 1 == ixrt)[0]
        # correnponding i range
        i_ixrt = np.where(self.i2ir - 1 == ixrt)[0]
        
        self.qi['i_ixrt'] = i_ixrt
        self._allocate_remap_qq(qq_type = 'qi', ijk = [len(i_ixrt),self.jx,self.kx])
                
        for np0 in nps:
            if self.iixl[np0] != 0:
                dtype = self._dtype_remap_qq(np0)
                filepath = self._get_filepath_remap_qq(n,np0)
                                
                with open(filepath,'rb') as f:
                    qqq = np.fromfile(f,dtype=dtype,count=1)
                    for key, m in zip( self.remap_kind, range(self.mtype) ):
                        self.qi[key][:,self.jss[np0]:self.jee[np0]+1,:] = qqq["qq"].reshape((self.iixl[np0],self.jjxl[np0],self.kx,self.mtype),order="F")[:,:,:,m]

                    
                    for key in self.remap_kind_add:
                        self.qi[key][:,self.jss[np0]:self.jee[np0]+1,:] = qqq[key].reshape((self.iixl[np0],self.jjxl[np0],self.kx),order="F")[:,:,:]
                
        if verbose:
            print('### variables are stored in self.qi ###')
    ##############################
    def time(self,n,tau=False,verbose=False):
        '''
        Reads time at a selected time step
        The data is stored in self.t

        Parameters
        ----------
        n : int
            selected time step for data
        tau : bool
            if True time for optical depth (high cadence)
        verbose : bool
            False suppresses a message of store
                
        Returns
        -------
        time : float
            time at a selected time step
        '''

        if tau:
            with open(self.datadir+"time/tau/t.dac."+'{0:08d}'.format(n),"rb") as f:
                self.t = np.fromfile(f,self.endian+'d',1).reshape((1),order='F')[0]
        else:
            with open(self.datadir+"time/mhd/t.dac."+'{0:08d}'.format(n),"rb") as f:
                self.t = np.fromfile(f,self.endian+'d',1).reshape((1),order='F')[0]
        
        if verbose:
            print('### variales are stored in self.t ###')

        return self.t
    ##############################
    def on_the_fly(self,n,verbose=True):
        '''
        Reads on the fly analysis data from fortran.
        The data is stored in self.vc dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        verbose : bool
            False suppresses a message of store
        '''

        # read xy plane data
        with open(self.datadir + "remap/vl/vl_xy.dac."+'{0:08d}'.format(n),"rb") as f:
            vl = np.fromfile(f,self.endian + 'f', self.m2d_xy*self.ix*self.jx ) \
                .reshape((self.ix, self.jx, self.m2d_xy ),order="F")

        for m in range(self.m2d_xy):
            self.vc[self.cl[m]] = vl[:,:,m]

        # read xz plane data
        with open(self.datadir + "remap/vl/vl_xz.dac."+'{0:08d}'.format(n),"rb") as f:
            vl = np.fromfile(f,self.endian + 'f', self.m2d_xz*self.ix*self.kx) \
                .reshape((self.ix, self.kx, self.m2d_xz ),order="F")

        for m in range(self.m2d_xz ):
            self.vc[self.cl[m + self.m2d_xy]] = vl[:,:,m]
            
        # read flux related value
        with open(self.datadir + "remap/vl/vl_flux.dac."+'{0:08d}'.format(n),"rb") as f:
            vl = np.fromfile(f,self.endian + 'f',self.m2d_flux*(self.ix+1)*self.jx ) \
                .reshape((self.ix+1,self.jx,self.m2d_flux),order="F")

        for m in range(self.m2d_flux):
            self.vc[self.cl[m+self.m2d_xy + self.m2d_xz]] = vl[:,:,m]
        
        # read spectra
        if self.geometry == 'YinYang':
            with open(self.datadir + "remap/vl/vl_spex.dac."+'{0:08d}'.format(n),"rb") as f:
                vl = np.fromfile(f,self.endian + 'f', self.m2d_spex*self.ix*self.kx//4) \
                    .reshape((self.ix,self.kx//4,self.m2d_spex),order="F")

            for m in range(self.m2d_spex):
                self.vc[self.cl[m+self.m2d_xy + self.m2d_xz + self.m2d_flux]] = vl[:,:,m]
                
        if verbose:
            print('### variables are stored in self.vc ###')

    ##############################
    def qq_check(self,n,verbose=True,end_step=False):
        '''
        Reads 3D full data for checkpoint
        The data is stored in self.qc dictionary
    
        Parameters
        ----------
        n : int
            A selected time step for data
        verbose : bool
            False suppresses a message of store
        end_step : bool
            If true, checkpoint of end step is read.
        '''
        
        step = '{0:08d}'.format(n)
        if end_step:
            if np.mod(self.nd,2) == 0:
                step = 'e'
            if np.mod(self.nd,2) == 1:
                step = 'o'

        with open(self.datadir+"qq/qq.dac."+step,'rb') as f:
            self.qc = \
                    np.fromfile(f,self.endian+'d',self.mtype*self.ixg*self.jxg*self.kxg) \
                            .reshape((self.ixg,self.jxg,self.kxg,self.mtype),order="F")

        f.close()
        
        if verbose:
            print('### variables are stored in self.qc ###')

    ##############################
    def qq_slice(self,n_slice,direc,n,verbose=True):
        '''
        Reads 2D data of slice.
        The data is stored in self.ql dictionary

        Parameters
        ----------
        n_slice : int
            index of slice
        direc : str
            slice direction. 'x', 'y', or 'z'
        n : int
            a selected time step for data
        verbose : bool
            False suppresses a message of store
        '''

        if self.geometry == 'YinYang':
            files = ['_yin','_yan']
            qls = [self.ql_yin,self.ql_yan]
        else:
            files = ['']
            qls = [self.ql]
        
        for file,ql in zip(files,qls):
            with open(self.datadir + 'slice/qq'+direc+file
                    +'.dac.'+'{0:08d}'.format(n)+'.'
                    +'{0:08d}'.format(n_slice+1),'rb') as f:
                if direc == 'x':
                    if self.geometry == 'YinYang':
                        n1, n2 = self.jx_yy + 2*self.margin, self.kx_yy+2*self.margin
                    else:
                        n1, n2 = self.jx, self.kx
                if direc == 'y':
                    n1, n2 = self.ix, self.kx
                if direc == 'z':
                    n1, n2 = self.ix, self.jx
                qq_slice = np.fromfile(f,self.endian+'f',(self.mtype+2)*n1*n2)
        
            for key, m in zip( self.remap_kind + self.remap_kind_add[:-1] , range(self.mtype + 2) ):
                ql[key] = qq_slice.reshape((n1,n2,self.mtype+2),order='F')[:,:,m]
            info = {}
            info['direc'] = direc
            if direc == 'x': slice = self.x_slice
            if direc == 'y': slice = self.y_slice
            if direc == 'z': slice = self.z_slice
            info['slice'] = slice[n_slice]
            info['n_slice'] = n_slice
            ql['info'] = info
            
        if verbose:
            if self.geometry == 'YinYang':
                print('### variables are stored in self.ql_yin and self.ql_yan ###')
            else:
                print('### variables are stored in self.ql ###')
        
    ##############################
    def qq_2d(self,n,verbose=True):
        '''
        Reads full data of 2D calculation
        The data is stored in self.q2 dictionary

        Parameters
        ----------
        n : int
            A selected time step for data 
        verbose : bool
            False suppresses a message of store
        '''
        value_keys = ['ro','vx','vy','vz','bx','by','bz','se','pr','te','op','tu','fr']
        
        ### Only when memory is not allocated 
        ### and the size of array is different
        ### memory is allocated
        memflag = True        
        if 'ro' in self.q2:
            memflag = not self.q2['ro'].shape == (self.ix,self.jx)
        if 'ro' not in self.q2 or memflag:
            for key in value_keys:
                self.q2[key] = np.zeros((self.ix,self.jx))
        dtype = np.dtype([ ("qq",self.endian + str((self.mtype+5)*self.ix*self.jx)+"f") ])
        with open(self.datadir + "remap/qq/qq.dac."+'{0:08d}'.format(n),'rb') as f:
            qqq = np.fromfile(f,dtype=dtype,count=1)
            
        for key, m in zip(value_keys, range(self.mtype)):
            self.q2[key] = qqq["qq"].reshape((self.mtype+5,self.ix,self.jx),order="F")[m,:,:]

        if verbose:
            print('### variables are stored in self.q2 ###')

    def YinYangSet(self):
        '''
        YinYangSet function sets up the YinYang geometry for
        YinYang direct plot.
        '''
        if self.geometry == 'YinYang':
            if not 'Z_yy' in self.p.__dict__:
                print('Yes')
                self.p.Z_yy , self.p.Y_yy  = np.meshgrid(self.z_yy,self.y_yy)
                self.p.Zg_yy, self.p.Yg_yy = np.meshgrid(self.zg_yy,self.yg_yy)

                # Geometry in Yang grid
                self.p.Yo_yy = np.arccos(np.sin(self.p.Y_yy)*np.sin(self.p.Z_yy))
                self.p.Zo_yy = np.arcsin(np.cos(self.p.Y_yy)/np.sin(self.pYo_yy))
                
                self.p.Yog_yy = np.arccos(np.sin(self.p.Yg_yy)*np.sin(self.p.Zg_yy))
                self.p.Zog_yy = np.arcsin(np.cos(self.p.Yg_yy)/np.sin(self.p.Yog_yy))
                
                sct =  np.sin(self.p.Y_yy)*np.cos(self.p.Z_yy)
                sco = -np.sin(self.p.Yo_yy)*np.cos(self.p.Zo_yy)

                sctg =  np.sin(self.p.Yg_yy)*np.cos(self.p.Zg_yy)
                scog = -np.sin(self.p.Yog_yy)*np.cos(self.p.Zog_yy)
                
                tmp = sct*sco < 0
                self.p.Zo_yy[tmp] = np.sign(self.p.Zo_yy[tmp])*np.pi - self.p.Zo_yy[tmp]

                tmp = sctg*scog < 0
                self.p.Zog_yy[tmp] = np.sign(self.p.Zog_yy[tmp])*np.pi - self.p.Zog_yy[tmp]
    
    def models(self,verbose=True):
        '''
        This method read Model S based stratification.
        
        Parameters
        ----------
        verbose : bool
            False suppresses a message
        '''

        self.ms = {}
        
        with open(self.datadir+'../input_data/params.txt','r') as f:
            self.ms['ix'] = int(f.read())

        ix = self.ms['ix']

        endian = '>'
        dtype = np.dtype([ \
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
        with open(self.datadir+"../input_data/value_cart.dac",'rb') as f:
            back = np.fromfile(f,dtype=dtype,count=1)
        
        for key in back.dtype.names:
            if back[key].size == ix:
                self.ms[key] = back[key].reshape((ix),order='F')
                
        if verbose:
            print('### Model S based stratification is stored in self.ms ###')
    
    def summary(self):
        '''
        Show pyR2D2.Read summary.
        '''
            
        RED = '\033[31m'
        END = '\033[0m'

        print(RED + '### Star information ###' + END)
        if 'mstar' in self.p.__dict__:
            print('mstar = ','{:.2f}'.format(self.mstar/pyR2D2.constant.MSUN)+' msun')
        if 'astar' in self.p.__dict__:
            print('astar = ','{:.2e}'.format(self.astar)+' yr')

        print('')
        if self.p.geometry == 'Cartesian':
            print(RED + '### calculation domain ###' + END)
            print('xmax - rstar = ', '{:6.2f}'.format((self.xmax - self.rstar)*1.e-8),'[Mm], xmin - rstar = ', '{:6.2f}'.format((self.xmin - self.rstar)*1.e-8),'[Mm]')
            print('ymax         = ', '{:6.2f}'.format(self.ymax*1.e-8)       ,'[Mm], ymin         = ', '{:6.2f}'.format(self.ymin*1.e-8),'[Mm]' )
            print('zmax         = ', '{:6.2f}'.format(self.zmax*1.e-8)       ,'[Mm], zmin         = ', '{:6.2f}'.format(self.zmin*1.e-8),'[Mm]' )

        if self.geometry == 'Spherical':
            pi2rad = 180/np.pi
            print('### calculation domain ###')
            print('xmax - rstar = ', '{:6.2f}'.format((self.xmax - self.rstar)*1.e-8),'[Mm],  xmin - rstar = ', '{:6.2f}'.format((self.xmin - self.rstar)*1.e-8),'[Mm]')
            print('ymax        = ', '{:6.2f}'.format(self.ymax*pi2rad)        ,'[rad], ymin        = ', '{:6.2f}'.format(self.ymin*pi2rad), '[rad]' )
            print('zmax        = ', '{:6.2f}'.format(self.zmax*pi2rad)        ,'[rad], zmin        = ', '{:6.2f}'.format(self.zmin*pi2rad), '[rad]' )

        if self.geometry == 'YinYang':
            pi2rad = 180/np.pi
            print('### calculation domain ###')
            print('xmax - rstar = ', '{:6.2f}'.format((self.xmax - self.rstar)*1.e-8),'[Mm],  xmin - rstar = ', '{:6.2f}'.format((self.xmin - self.rstar)*1.e-8),'[Mm]')
            print('Yin-Yang grid is used to cover the whole sphere')

        print('')
        print(RED + '### number of grid ###' + END)
        print('(ix,jx,kx)=(',self.ix,',',self.jx,',',self.kx,')')

        print('')
        print(RED + '### calculation time ###' + END)
        print('time step (nd) =',self.nd)
        print('time step (nd_tau) =',self.nd_tau)
        t = self.time(self.nd)
        print('time =','{:.2f}'.format(t/3600),' [hour]')