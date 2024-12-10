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
    from google.oauth2.service_account import Credentials

    scopes = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']

    credentials = Credentials.from_service_account_file(json_key, scopes=scopes)
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
def set_cells_gspread(data,
                json_key=None,
                project=__file__.split('/')[-4],
                caseid=None):
    '''
    Outputs parameters to Google spreadsheet

    Parameters
    ----------
    data : R2D2.Data, or, R2D2.Read
        instance of R2D2.Data or R2D2.Read classes
    json_key : str
        File of json key to access Google API
    project : str
        Project name, typically name of upper directory
    caseid : str
        Case ID
    
    '''
    import datetime
    import numpy as np
    import pyR2D2

    if json_key == None:
        json_key = glob.glob(os.environ['HOME']+'/json/*')[0]  

    if caseid is None:
        caseid = data.p['datadir'].split('/')[-3]
    
    gc = init_gspread(json_key,project)
    wks = gc.open(project).sheet1
    str_id = str(int(caseid[1:])+1)
    cells = wks.range('A'+str_id+':'+'T'+str_id)

    keys = [caseid]
    if 'mstar' in data.p:
        keys.append('{:.2f}'.format(data.p['mstar']/pyR2D2.constant.msun))
    else:
        keys.append('1.00') # solar mass
    keys.append(str(data.p['ix'])+' '+str(data.p['jx'])+' '+str(data.p['kx']))
    keys.append( '{:6.2f}'.format((data.p['xmin']-data.p['rstar'])*1.e-8))
    keys.append( '{:6.2f}'.format((data.p['xmax']-data.p['rstar'])*1.e-8))
    

    if data.p['geometry'] == 'Cartesian':
        keys.append( '{:6.2f}'.format(data.p['ymin']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(data.p['ymax']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(data.p['zmin']*1.e-8)+' [Mm]')
        keys.append( '{:6.2f}'.format(data.p['zmax']*1.e-8)+' [Mm]')

    if data.p['geometry'] == 'Spherical':
        pi2rad = 180/np.pi
        keys.append( '{:6.2f}'.format(data.p['ymin']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(data.p['ymax']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(data.p['zmin']*pi2rad)+' [deg]')
        keys.append( '{:6.2f}'.format(data.p['zmax']*pi2rad)+' [deg]')

    if data.p['geometry'] == 'YinYang':
        pi2rad = 180/np.pi
        keys.append( '0 [deg]')
        keys.append( '180 [deg]')
        keys.append( '-180 [deg]')
        keys.append( '180 [deg]')
    
    if data.p['ununiform_flag']:
        keys.append('F')
    else:
        keys.append('T')
    
    dx0 = (data.p['x'][1] - data.p['x'][0])*1.e-5
    dx1 = (data.p['x'][data.p['ix']-1] - data.p['x'][data.p['ix']-2])*1.e-5
    keys.append( '{:6.2f}'.format(dx0)+' '+'{:6.2f}'.format(dx1))
    keys.append( data.p['rte'])
    keys.append( '{:6.2f}'.format(data.p['dtout']))
    keys.append( '{:6.2f}'.format(data.p['dtout_tau']))
    keys.append( '{:5.2f}'.format(data.p['potential_alpha']))
    if data.p['xi'].max() == 1.0:
        keys.append('F')
    else:
        keys.append('T')

    keys.append( '{:5.1f}'.format(data.p['omfac']))
    keys.append(data.p['geometry'])
    keys.append(data.p['origin'])
    keys.append(str(datetime.datetime.now()).split('.')[0])
    keys.append(data.p['server'])
    
    for cell, key in zip(cells,keys):
        cell.value = key
    
    wks.update_cells(cells)

