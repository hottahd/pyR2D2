import R2D2
import glob
import os

caseid = R2D2.util.caseid_select(locals(),force_read=True)
datadir="../run/"+caseid+"/data/"

d = R2D2.Data(datadir,verbose=True)

R2D2.write.google.set_top_line()
R2D2.write.google.set_cells_gspread(d)
print('\n### Google SpreadSheet write finished ###')