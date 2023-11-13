def initialize_instance(main_locals,instance):
    if not instance in main_locals:
        main_locals[instance] = None

def caseid_select(main_locals,force_read=False):
    
    RED = '\033[31m'
    END = '\033[0m'

    if 'caseid' in main_locals and not force_read:
        caseid = main_locals['caseid']
    else:
       caseid = 'd'+str(input(RED + "input caseid id (3 digit): " + END)).zfill(3)
       
    return caseid

def locals_define(self,main_locals):
    for key, value in self.p.items():
        main_locals[key] = value    
    return

def define_n0(self,main_locals,nd_type='nd'):
    if 'n0' not in main_locals:
        main_locals['n0'] = 0
    if main_locals['n0'] > self.p[nd_type]:
        main_locals['n0'] = self.p[nd_type]
    

def show_information(self):
    import numpy as np
    import R2D2
    
    RED = '\033[31m'
    END = '\033[0m'

    print(RED + '### Star information ###' + END)
    print('mstar = ','{:.2f}'.format(self.p['mstar']/R2D2.msun)+' msun')    
    print('astar = ','{:.2e}'.format(self.p['astar'])+' yr')

    print('')
    if self.p['geometry'] == 'Cartesian':
        print(RED + '### calculation domain ###' + END)
        print('xmax - rstar = ', '{:6.2f}'.format((self.p['xmax'] - self.p['rstar'])*1.e-8),'[Mm], xmin - rstar = ', '{:6.2f}'.format((self.p['xmin'] - self.p['rstar'])*1.e-8),'[Mm]')
        print('ymax         = ', '{:6.2f}'.format(self.p['ymax']*1.e-8)       ,'[Mm], ymin         = ', '{:6.2f}'.format(self.p['ymin']*1.e-8),'[Mm]' )
        print('zmax         = ', '{:6.2f}'.format(self.p['zmax']*1.e-8)       ,'[Mm], zmin         = ', '{:6.2f}'.format(self.p['zmin']*1.e-8),'[Mm]' )

    if self.p['geometry'] == 'Spherical':
        pi2rad = 180/np.pi
        print('### calculation domain ###')
        print('xmax - rstar = ', '{:6.2f}'.format((self.p['xmax'] - self.p['rstar'])*1.e-8),'[Mm],  xmin - rstar = ', '{:6.2f}'.format((self.p['xmin'] - self.p['rstar'])*1.e-8),'[Mm]')
        print('ymax        = ', '{:6.2f}'.format(self.p['ymax']*pi2rad)        ,'[rad], ymin        = ', '{:6.2f}'.format(self.p['ymin']*pi2rad),'[rad]' )
        print('zmax        = ', '{:6.2f}'.format(self.p['zmax']*pi2rad)        ,'[rad], zmin        = ', '{:6.2f}'.format(self.p['zmin']*pi2rad),'[rad]' )

    if self.p['geometry'] == 'YinYang':
        pi2rad = 180/np.pi
        print('### calculation domain ###')
        print('xmax - rstar = ', '{:6.2f}'.format((self.p['xmax'] - self.p['rstar'])*1.e-8),'[Mm],  xmin - rstar = ', '{:6.2f}'.format((self.p['xmin'] - self.p['rstar'])*1.e-8),'[Mm]')
        print('Yin-Yang grid is used to cover the whole sphere')

    print('')
    print(RED + '### number of grid ###' + END)
    print('(ix,jx,kx)=(',self.p['ix'],',',self.p['jx'],',',self.p['kx'],')')

    print('')
    print(RED + '### calculation time ###' + END)
    print('time step (nd) =',self.p['nd'])
    print('time step (nd_tau) =',self.p['nd_tau'])
    t = self.read_time(self.p['nd'])
    print('time =','{:.2f}'.format(t/3600),' [hour]')