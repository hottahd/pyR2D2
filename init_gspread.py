import R2D2
import glob
import os

caseid = R2D2.util.caseid_select(locals(),force_read=True)
datadir="../run/"+caseid+"/data/"

d = R2D2.R2D2_data(datadir,verbose=True)

json_key = glob.glob(os.environ['HOME']+'/json/*')[0]
project = os.getcwd().split('/')[-2]

R2D2.google.init_gspread(json_key,project)
R2D2.google.out_gspread(d,caseid,json_key,project)
print('\n### Google SpreadSheet write finished ###')
