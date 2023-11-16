import R2D2
import glob
import os

caseid = R2D2.util.caseid_select(locals(),force_read=True)
datadir="../run/"+caseid+"/data/"

d = R2D2.R2D2_data(datadir,verbose=True)

R2D2.google.set_top_line()
d.set_cells_gspread()
print('\n### Google SpreadSheet write finished ###')