def initialize_instance(main_locals,instance):
    '''
    Initialize arbitrary instance of R2D2_data object in main program.
        
    Parameters:
       main_locals (dictionary): locals() in main program, which include local variables
       instance (str): name of instance of R2D2_data object
    '''
    if not instance in main_locals:
        main_locals[instance] = None

def caseid_select(main_locals,force_read=False):
    '''
    Choose caseid from input or from user input.
    
    Parameters:
        main_locals (dictionary): locals() in main program, which include local variables
        force_read (bool): if True, force to read caseid from user input, 
                           if False, read caseid from main_locals if exists
    '''
    
    RED = '\033[31m'
    END = '\033[0m'

    if 'caseid' in main_locals and not force_read:
        caseid = main_locals['caseid']
    else:
       caseid = 'd'+str(input(RED + "input caseid id (3 digit): " + END)).zfill(3)
       
    return caseid

def locals_define(self,main_locals):
    '''
    Substitute selp.p to main_locals in main program.
    
    Parameters:
       self (R2D2_data): instance of R2D2_data object
       main_locals (dictionary): locals() in main program, which include local variables
    '''
    for key, value in self.p.items():
        main_locals[key] = value    
    return

def define_n0(self,main_locals,nd_type='nd'):
    '''
    Define n0 in main_locals if not exists.
    
    Parameters:
       self (R2D2_data): instance of R2D2_data object
       main_locals (dictionary): locals() in main program, which include local variables
       nd_type (str): type of nd. 'nd' for MHD output, 'nd_tau' for high cadence output.
    '''
    if 'n0' not in main_locals:
        main_locals['n0'] = 0
    if main_locals['n0'] > self.p[nd_type]:
        main_locals['n0'] = self.p[nd_type]
    
    return main_locals['n0']

def show_information(self):
    '''
    Show data information
    
    Parameters:
       self (R2D2_data): instance of R2D2_data object
    '''
    
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

    return

def get_best_unit(size, unit_multipliers):
    for unit in reversed(['B', 'kB', 'MB', 'GB', 'TB', 'PB']):
        if size >= unit_multipliers[unit]:
            return unit
    return 'B'

def get_total_file_size(directory, unit=None):
    '''
    Evaluate total size of files in directory.
    
    Parameters:
       directory (str): directory path
       unit (str): unit of file size. Choose from 'B', 'kB', 'MB', 'GB', 'TB', 'PB'.
    Returns:
       total_size (int): total size of files in directory in bytes
    '''
    import os
    
    unit_multipliers = {
        'B': 1,
        'kB': 1024,
        'MB': 1024**2,
        'GB': 1024**3,
        'TB': 1024**4,
        'PB': 1024**5
    }
    
    if unit is not None and unit not in unit_multipliers:
        raise ValueError(f"Invalid unit: {unit}. Choose from 'B', 'kB', 'MB', 'GB', 'TB', 'PB'.")
    
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            try:
                total_size += os.path.getsize(filepath)
                
                if unit is None:
                    display_unit = get_best_unit(total_size, unit_multipliers)
                else:
                    display_unit = unit
                converted_size = total_size / unit_multipliers[display_unit]
                print(f"\r{' '*50}\rCurrent total size: {converted_size:.2f} {display_unit}", end='', flush=True)
                #print(f"\rCurrent total size: {converted_size:.2f} {display_unit} ", end='', flush=True)
            except FileNotFoundError:
                # ファイルが見つからない場合は無視
                continue
    
    if unit is None:
        unit = get_best_unit(total_size, unit_multipliers)
        
    final_size = total_size / unit_multipliers[unit]
    print(f"\nFinal total size: {final_size:.2f} {unit}")    
            
    return final_size, unit


def update_results_file(file_path, total_size, unit, caseid, dir_path):
    '''
    Update the results file with the size of a directory.
    
    Parameters:
       file_path (str): path to the results file
       total_size (float): total size of the directory
       unit (str): unit of the file size
       caseid (str): the case ID of the directory
       dir_path (str): the directory path
    '''
    import os
    from datetime import datetime
    # 現在の日時を取得
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

     # シンボリックリンクの場合の実体パスを取得
    if os.path.islink(dir_path):
        real_path = os.path.abspath(os.readlink(dir_path))
    else:
        real_path = os.path.abspath(dir_path)
    
    # 既存のデータを読み込む
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            existing_data = f.readlines()
    else:
        existing_data = []
    
    # データを辞書に変換
    existing_dict = {}
    for line in existing_data:
        parts = line.strip().rsplit(maxsplit=3)
        if len(parts) == 4:
            existing_caseid = parts[0]
            size = parts[1]
            timestamp = parts[2]
            path = parts[3]
            existing_dict[existing_caseid] = f"{size} {timestamp} {path}"
    
    # ディレクトリサイズの情報を更新
    directory_size_str = f"{total_size:6.2f} {unit}"
    existing_dict[caseid] = f"{directory_size_str} {now} {real_path}"
    
    # 更新されたデータを書き込む
    with open(file_path, 'w') as f:
        for caseid in sorted(existing_dict):
            size_timestamp_realpath = existing_dict[caseid]
            # フォーマット: ディレクトリ名、右揃えで6桁、小数点2桁、単位
            f.write(f"{caseid:<10} {size_timestamp_realpath}\n")