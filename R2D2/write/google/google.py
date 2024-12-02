'''
    Functions for push information to Google Spreadsheet
'''
import glob, os

def init_gspread(json_key, project):
    '''
    This function initialize the utility of google spread 

    Parameters
    ----------
    json_key : str
        file of json key to access Google API
    project : str
        Project name, typically name of upper directory

    Returns
    -------
    gc : gspread.client.Client
        Instance of Google API
    '''

    import gspread
    from oauth2client.service_account import ServiceAccountCredentials

    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(json_key, scope)
    gc = gspread.authorize(credentials)
    
    return gc
################################################################################
def fetch_URL_gspread(
    json_key=None,
    project=__file__.split('/')[-4],    
    ):
    '''
    Fetchs corresponding URL for google spreadsheet

    Parameters
    ----------
    json_key : str 
        File of json key to access Google API
    projet : str
        Project name, typically name of upper directory

    Returns
    -------
    URL : str
        Google spreadsheet URL
    '''
    
    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]
    gid = init_gspread(json_key,project).open(project).id
    
    return 'https://docs.google.com/spreadsheets/d/'+gid
################################################################################

def set_top_line(
    json_key=None,
    project=__file__.split('/')[-4],
    ):
    '''
    This function set top line of google spreadsheet

    Parameters
    ----------
    json_key : str
        File of json key to access Google API
    projet : str
        Project name, typically name of upper directory

    Returns
    -------
        None
    '''
    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]  
    gc = init_gspread(json_key,project)
    wks = gc.open(project).sheet1

    cells = wks.range('A1:T1')
    keys = [ 'Case ID' \
             ,'Mstar' \
             ,'(ix,jx,kx)' \
             ,'xmin [Mm]' \
             ,'xmax [Mm]' \
             ,'ymin' \
             ,'ymax' \
             ,'zmin' \
             ,'zmax' \
             ,'uni' \
             ,'dx [km]' \
             ,'m ray' \
             ,'dtout [s]' \
             ,'dtout_tau [s]' \
             ,'al' \
             ,'RSST' \
             ,'Om' \
             ,'Gemetry' \
             ,'origin'
             ,'update time' \
             ,'Server'
             ]

    for cell, key in zip(cells,keys):
        cell.value = key

    wks.update_cells(cells)
        
################################################################################
def set_cells_gspread(r2d2_data,
                json_key=None,
                project=__file__.split('/')[-4],
                caseid=None):
    '''
    Outputs parameters to Google spreadsheet

    Parameters
    ----------
    r2d2_data : R2D2.R2D2_data, or, R2D2.R2D2_read
        instance of R2D2_data or R2D2_read classes
    json_key : str
        File of json key to access Google API
    project : str
        Project name, typically name of upper directory
    caseid : str): caseid
    
    '''
    import datetime
    import numpy as np
    import R2D2

    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]  

    if caseid is None:
        caseid = r2d2_data.p['datadir'].split('/')[-3]
    
    gc = init_gspread(json_key,project)
    wks = gc.open(project).sheet1
    str_id = str(int(caseid[1:])+1)
    cells = wks.range('A'+str_id+':'+'T'+str_id)

    keys = [caseid]
    keys.append('{:.2f}'.format(r2d2_data.p['mstar']/R2D2.Constant.msun))
    keys.append(str(r2d2_data.p['ix'])+' '+str(r2d2_data.p['jx'])+' '+str(r2d2_data.p['kx']))
    keys.append( '{:6.2f}'.format((r2d2_data.p['xmin']-r2d2_data.p['rstar'])*1.e-8))
    keys.append( '{:6.2f}'.format((r2d2_data.p['xmax']-r2d2_data.p['rstar'])*1.e-8))
    

    if r2d2_data.p['geometry'] == 'Cartesian':
        keys.append( '{:6.2f}'.format(r2d2_data.p['ymin']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(r2d2_data.p['ymax']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(r2d2_data.p['zmin']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(r2d2_data.p['zmax']*1.e-8)+' [Mm]')

    if r2d2_data.p['geometry'] == 'Spherical':
        pi2rad = 180/np.pi
        keys.append( '{:6.2f}'.format(r2d2_data.p['ymin']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(r2d2_data.p['ymax']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(r2d2_data.p['zmin']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(r2d2_data.p['zmax']*pi2rad)+' [deg]')

    if r2d2_data.p['geometry'] == 'YinYang':
        pi2rad = 180/np.pi
        keys.append( '0 [deg]')
        keys.append( '180 [deg]')
        keys.append( '-180 [deg]')
        keys.append( '180 [deg]')
    
    if r2d2_data.p['ununiform_flag']:
        keys.append('F')
    else:
        keys.append('T')
    
    dx0 = (r2d2_data.p['x'][1] - r2d2_data.p['x'][0])*1.e-5
    dx1 = (r2d2_data.p['x'][r2d2_data.p['ix']-1] - r2d2_data.p['x'][r2d2_data.p['ix']-2])*1.e-5
    keys.append( '{:6.2f}'.format(dx0)+' '+'{:6.2f}'.format(dx1))
    keys.append( r2d2_data.p['rte'])
    keys.append( '{:6.2f}'.format(r2d2_data.p['dtout']))
    keys.append( '{:6.2f}'.format(r2d2_data.p['dtout_tau']))
    keys.append( '{:5.2f}'.format(r2d2_data.p['potential_alpha']))
    if r2d2_data.p['xi'].max() == 1.0:
        keys.append('F')
    else:
        keys.append('T')

    keys.append( '{:5.1f}'.format(r2d2_data.p['omfac']))
    keys.append(r2d2_data.p['geometry'])
    keys.append(r2d2_data.p['origin'])
    keys.append(str(datetime.datetime.now()).split('.')[0])
    keys.append(r2d2_data.p['server'])
    
    for cell, key in zip(cells,keys):
        cell.value = key
    
    wks.update_cells(cells)

