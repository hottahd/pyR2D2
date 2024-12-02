import R2D2

def init(main_locals,instance_name='d'):
    """
    Initialize R2D2.Data instance in main program.
    
    Parameters
    ----------
    main_locals : dict
        locals() in main program, which include local variables
    instance_name : str
        name of instance of R2D2.Data object

    
    Notes
    -----
    1. Select caseid from existing caseid in main_locals or from user input.
    2. Initialize instance of R2D2.Data object in main program.
    3. Assign the R2D2.Data.p dictionary keys to main_locals.
        
    Examples
    --------
    .. code-block:: python

        import R2D2
        R2D2.util.init(locals())
        
    """
    
    main_locals['caseid'] = caseid_select(main_locals)
    datadir = "../run/"+main_locals['caseid']+"/data/"
    
    initialize_instance(main_locals,instance_name)
    main_locals[instance_name] = R2D2.Data(datadir)
    
def initialize_instance(main_locals,instance_name):
    '''
    Initializes arbitrary instance of R2D2.Data object in main program.
        
    Parameters
    ----------
    main_locals : dict
        locals() in main program, which include local variables
    instance_name : str
        name of instance of R2D2.Data object
    '''
    if not instance_name in main_locals:
        main_locals[instance_name] = None

def caseid_select(main_locals,force_read=False):
    '''
    Choose caseid from input or from user input.
    
    Parameters
    ----------
    main_locals : dict 
        locals() in main program, which include local variables
    force_read : bool
        If True, force to read caseid from user input.
        If False, read caseid from main_locals if exists
    '''
    
    RED = '\033[31m'
    END = '\033[0m'

    if 'caseid' in main_locals and not force_read:
        caseid = main_locals['caseid']
    else:
       caseid = 'd'+str(input(RED + "input caseid id (3 digit): " + END)).zfill(3)
       
    return caseid

def locals_define(data,main_locals):
    '''
    Substitute selp.p to main_locals in main program.
    
    Parameters
    ----------
    data : R2D2.Data, or, R2D2.Read
        Instance of R2D2.Data or R2D2.Read classes
    main_locals : dict
        locals() in main program, which include local variables
    '''
    for key, value in data.p.items():
        main_locals[key] = value    
    return

def define_n0(data,main_locals,nd_type='nd'):
    '''
    Define n0 in main_locals if not exists.
    
    Parameters
    ----------
    data : R2D2.Data, or, R2D2.Read
        instance of R2D2.Data or R2D2.Read classes
    main_locals : dict
        locals() in main program, which include local variables
    nd_type : str
        type of nd. 'nd' for MHD output, 'nd_tau' for high cadence output.
    '''
    if 'n0' not in main_locals:
        main_locals['n0'] = 0
    if main_locals['n0'] > data.p[nd_type]:
        main_locals['n0'] = data.p[nd_type]
    
    return main_locals['n0']

def get_best_unit(size, unit_multipliers):
    """
    Get the best unit for displaying the size.

    Parameters
    ----------
    size : int
        size of the file
    unit_multipliers : dict
        dictionary of unit multipliers
    """
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
    Updates the results file with the size of a directory.
    
    Parameters
    ----------
    file_path : str
        Path to the results file
    total_size : float
        Total size of the directory
    unit : str
        Unit of the file size
    caseid : str
        The case ID of the directory
    dir_path : str
        The directory path
        
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